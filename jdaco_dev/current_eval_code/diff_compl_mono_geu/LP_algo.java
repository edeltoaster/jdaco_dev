package diff_compl_mono_geu;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import com.joptimizer.exception.JOptimizerException;
import com.joptimizer.optimizers.LPOptimizationRequest;
import com.joptimizer.optimizers.LPPrimalDualMethod;

import framework.QuantDACOResultSet;
import mixed.abundance_estimation_algorithm_test;

public class LP_algo {
	private Map<HashSet<String>, Double> cached_abundance_of_complexes;
	private Map<String, Double> cached_remaining_abundance_of_proteins;
	private Map<HashSet<String>, Integer> compl_sol_map = new HashMap<>(); // stores relation of  complexes->j in c_j; c_0 is used to model the constant factor
	private double tolerance = 1.0E-7;
	LPOptimizationRequest or = new LPOptimizationRequest();
	private QuantDACOResultSet qdr;
	
	public LP_algo(QuantDACOResultSet qdr) {
		this.qdr = qdr;
		Logger.getRootLogger().setLevel(Level.OFF); // necesary to use JOptimizer
		
		Map<String, List<HashSet<String>>> prot_compl_map = new HashMap<>();
		Map<HashSet<String>, Double> compl_min_map = new HashMap<>();
		double[] c = new double[qdr.getResult().size()+1];
		
		// define opt function without constant
		int sol_i = 1; // 0: constant
		for (HashSet<String> compl:qdr.getResult()) {
			double factors = 0;
			for (String prot:compl) {
				double p_i = qdr.getProteinAbundance(prot);
				
				if (!prot_compl_map.containsKey(prot))
					prot_compl_map.put(prot, new LinkedList<>());
				prot_compl_map.get(prot).add(compl);
				
				// store min protein per complex to set tighter boundaries
				if (p_i < compl_min_map.getOrDefault(compl, Double.MAX_VALUE)) 
					compl_min_map.put(compl, p_i);
				
				double factor = -1.0 / p_i;
				factors += factor;
			}
			c[sol_i] = factors;
			this.compl_sol_map.put(compl, sol_i);
			sol_i++;
		}
		
		// define constraints
		double[] protein_levels = new double[prot_compl_map.keySet().size()];
		double[][] protein_constraints = new double[prot_compl_map.keySet().size()][qdr.getResult().size()+1];
		
		int i = 0;
		for (String prot:prot_compl_map.keySet()) {
			double p_i = qdr.getProteinAbundance(prot);
			double[] c_constr = new double[qdr.getResult().size() + 1];
			c_constr[0] = 0;
			protein_levels[i] = p_i;
			
			for (HashSet<String> compl:prot_compl_map.get(prot)) {
				c_constr[this.compl_sol_map.get(compl)] = 1.0;
			}
			
			protein_constraints[i] = c_constr;
			i++;
		}
		
		// add artificial constraint to mimick adding a constant
		c[0] = prot_compl_map.keySet().size();
		
		// upper/lower bounds on c_j
		double[] lower_bounds = new double[qdr.getResult().size() + 1];
		double[] upper_bounds = new double[qdr.getResult().size() + 1];
		
		// fix c_0 = constant in objective function
		lower_bounds[0] = 1.0;
		upper_bounds[0] = 1.0;
		
		// remaining lower bounds automatically 0.0, upper bounds to most abundant protein
		for (HashSet<String> compl:qdr.getResult()) {
			upper_bounds[compl_sol_map.get(compl)] = compl_min_map.get(compl);
		}
		
		//optimization problem
		this.or.setC(c);
		this.or.setG(protein_constraints);
		this.or.setH(protein_levels);
		this.or.setLb(lower_bounds);
		this.or.setUb(upper_bounds);
		this.or.setTolerance(this.tolerance);
		
		// solve
		this.solve();
	}
	
	private void solve() {
		
		LPPrimalDualMethod opt = new LPPrimalDualMethod();
		opt.setLPOptimizationRequest(or);
		
		try {
			opt.optimize();
		} catch (JOptimizerException e) {
			e.printStackTrace();
		}
		
		double[] sol = opt.getOptimizationResponse().getSolution();
		
		this.cached_abundance_of_complexes = new HashMap<>();
		this.cached_remaining_abundance_of_proteins = new HashMap<>();
		
		// set complex abundance from solution
		for (HashSet<String> compl:this.qdr.getResult()) {
			double res = sol[this.compl_sol_map.get(compl)];
			this.cached_abundance_of_complexes.put(compl, res);
		}
		
		// set remaining abundance
		Map<String, Double> used_abundance = new HashMap<>();
		for (Entry<HashSet<String>, Double> e:this.cached_abundance_of_complexes.entrySet()) {
			for (String p:e.getKey())
				used_abundance.put(p, used_abundance.getOrDefault(p, 0.0) + e.getValue());
		}
		for (String p:this.qdr.getProteinToAssumedTranscript().keySet()) {
			this.cached_remaining_abundance_of_proteins.put(p, this.qdr.getProteinAbundance(p) - used_abundance.getOrDefault(p, 0.0));
		}
		this.cached_remaining_abundance_of_proteins.keySet().removeIf(k -> this.cached_remaining_abundance_of_proteins.get(k) == 0.0);
		
	}
	
	public Map<HashSet<String>, Double> getAbundanceOfComplexes() {
		return cached_abundance_of_complexes;
	}


	public Map<String, Double> getRemainingAbundanceOfProteins() {
		return cached_remaining_abundance_of_proteins;
	}

	public Map<HashSet<String>, Integer> getComplSolMap() {
		return compl_sol_map;
	}

	public double getTolerance() {
		return tolerance;
	}

	public LPOptimizationRequest getOptimizationRequest() {
		return or;
	}

	public static void main(String[] args) {
		QuantDACOResultSet qdr = abundance_estimation_algorithm_test.simple_example();
		
		LP_algo lp = new LP_algo(qdr);
		
		System.out.println();
		System.out.println(lp.getAbundanceOfComplexes());
		System.out.println(qdr.getAbundanceOfComplexes());
		System.out.println();
		System.out.println(lp.getRemainingAbundanceOfProteins());
		System.out.println(qdr.getRemainingAbundanceOfProteins());
	}

}
