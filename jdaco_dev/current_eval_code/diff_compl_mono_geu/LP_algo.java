package diff_compl_mono_geu;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.QuantDACOResultSet;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LPWizard;
import scpsolver.problems.LPWizardConstraint;
import scpsolver.problems.LinearProgram;
import mixed.abundance_estimation_algorithm_test;

public class LP_algo {
	private Map<HashSet<String>, Double> cached_abundance_of_complexes;
	private Map<String, Double> cached_remaining_abundance_of_proteins;
	private QuantDACOResultSet qdr;
	private LinearProgram lp;
	
	public LP_algo(QuantDACOResultSet qdr) {
		LPWizard lpw = new LPWizard(); 
		this.qdr = qdr;
		
		Map<String, List<HashSet<String>>> prot_compl_map = new HashMap<>();
		// define opt function without constant
		for (HashSet<String> compl:qdr.getResult()) {
			String c_j = String.join(",", compl);
			double factors = 0;
			for (String prot:compl) {
				double p_i = qdr.getProteinAbundance(prot);
				if (p_i == 0)
					continue;
				if (!prot_compl_map.containsKey(prot))
					prot_compl_map.put(prot, new LinkedList<>());
				prot_compl_map.get(prot).add(compl);
				double factor = -1.0 / p_i;
				factors += factor;
			}
			lpw.plus(c_j, factors);
		}
		
		// define constraints
		double N = 0;
		for (String prot:prot_compl_map.keySet()) {
			double p_i = qdr.getProteinAbundance(prot);
			if (p_i == 0)
				continue;
			N++;
			LPWizardConstraint c = lpw.addConstraint(prot, p_i,">=");
			
			for (HashSet<String> compl:prot_compl_map.get(prot)) {
				String c_j = String.join(",", compl);
				c.plus(c_j, 1.0);
			}
		}
		
		// add artificial constraint to mimick adding a constant
		lpw.plus("const", N);
		lpw.addConstraint("const1", 1.0,">=").plus("const", 1.0);
		lpw.addConstraint("const2", 1.0,"<=").plus("const", 1.0);
		
		// constraints to enforce C_j > 0
		for (HashSet<String> compl:qdr.getResult()) {
			String c_j = String.join(",", compl);
			lpw.addConstraint(c_j + "+", 0.0,"<=").plus(c_j, 1.0);
		}
		// setup solver
		this.lp = lpw.getLP();
		this.lp.setMinProblem(true);
	}
	
	public void solve(int time_constraint) {
		LinearProgramSolver solver  = SolverFactory.newDefault();
		solver.setTimeconstraint(time_constraint);
		
		// solve
		double[] sol = solver.solve(this.lp);
		this.cached_abundance_of_complexes = new HashMap<>();
		this.cached_remaining_abundance_of_proteins = new HashMap<>();
		
		if (sol == null || sol.length < this.qdr.getResult().size())
			return;
		
		for (HashSet<String> compl:this.qdr.getResult()) {
			this.cached_abundance_of_complexes.put(compl, sol[this.lp.getIndexmap().get(String.join(",", compl))]);
		}
	}
	
	public Map<HashSet<String>, Double> getAbundanceOfComplexes() {
		return cached_abundance_of_complexes;
	}


	public Map<String, Double> getRemainingAbundanceOfProteins() {
		return cached_remaining_abundance_of_proteins;
	}

	public LinearProgram getLP() {
		return lp;
	}

	public static void main(String[] args) {
	
		QuantDACOResultSet qdr = abundance_estimation_algorithm_test.alterating_example();
		
		LP_algo lp = new LP_algo(qdr);
		System.out.println(lp.getLP().getIndexmap());
		System.out.println(lp.getLP().convertToCPLEX());
		lp.solve(100);
		System.out.println(lp.getAbundanceOfComplexes());
		System.out.println(qdr.getAbundanceOfComplexes());
		
		QuantDACOResultSet dr = new QuantDACOResultSet("/Users/tho/Desktop/HG00171.csv.gz", "mixed_data/hocomoco_human_TFs_v10.txt.gz", "/Users/tho/Desktop/HG00171_major-transcripts.txt.gz");
		lp = new LP_algo(dr);
		lp.solve(1200);
		System.out.println(lp.getAbundanceOfComplexes());
		System.out.println(dr.getAbundanceOfComplexes());
	}

}
