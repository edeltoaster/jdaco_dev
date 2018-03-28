package diff_compl_mono_geu;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import framework.QuantDACOResultSet;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LPWizard;
import scpsolver.problems.LPWizardConstraint;
import scpsolver.problems.LinearProgram;

public class LP_algo {
	private Map<HashSet<String>, Double> cached_abundance_of_complexes;
	private Map<String, Double> cached_remaining_abundance_of_proteins;
	private QuantDACOResultSet qdr;
	private LinearProgram lp;
	
	public LP_algo(QuantDACOResultSet qdr) {
		LPWizard lpw = new LPWizard(); 
		this.qdr = qdr;
		
		// TODO: revise construction of linear program
		// define opt function
		for (HashSet<String> compl:qdr.getResult()) {
			String c_j = String.join(",", compl);
			lpw.plus(c_j, 1.0);
		}
		
		// define constraints
		for (String prot:qdr.getProteinToAssumedTranscript().keySet()) {
			double p_i = qdr.getProteinAbundance(prot);
			if (p_i == 0)
				continue;
			LPWizardConstraint c = lpw.addConstraint(prot, p_i,">=");
			
			for (HashSet<String> compl:qdr.getResult()) {
				if (!compl.contains(prot))
					continue;
				String c_j = String.join(",", compl);
				c.plus(c_j, 1.0);
			}
		}
		
		// setup solver
		this.lp = lpw.getLP();
		this.lp.setMinProblem(false);
	}
	
	public void solve(int time_constraint) {
		LinearProgramSolver solver  = SolverFactory.newDefault();
		solver.setTimeconstraint(time_constraint);
		
		// solve
		double[] sol = solver.solve(this.lp);
		int i = 0;
		this.cached_abundance_of_complexes = new HashMap<>();
		this.cached_remaining_abundance_of_proteins = new HashMap<>();
		
		if (sol == null || sol.length < this.qdr.getResult().size())
			return;
		
		for (HashSet<String> compl:this.qdr.getResult()) {
			this.cached_abundance_of_complexes.put(compl, sol[i]);
			i++;
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
		QuantDACOResultSet dr = new QuantDACOResultSet("/Users/tho/Desktop/HG00171.csv.gz", "mixed_data/hocomoco_human_TFs_v10.txt.gz", "/Users/tho/Desktop/HG00171_major-transcripts.txt.gz");
		LP_algo lp = new LP_algo(dr);
		lp.solve(1200);
		System.out.println(lp.getAbundanceOfComplexes());
	}

}
