package PPIComp_hemato;

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.Utilities;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;

public class hemato_build_rand_diff_subnetworks {
	
	static int no_threads = 40;
	static double fdr = 0.05;
	static String t1 = "CLP";
	static String t2 = "CD4";
	
	// needs to be run on server
	static String network_folder = "BLUEPRINT_networks/";
	static String results_root = "hemato_eval/";
	
	public static Map<String, RewiringDetectorSample> getSubset(Map<String, RewiringDetectorSample> data, Set<String> subset) {
		Map<String, RewiringDetectorSample> subset_data = new HashMap<>();
		
		for (String s:subset)
			subset_data.put(s, data.get(s));
		
		return subset_data;
	}
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("read data");
		Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + t1 + "/");
		Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + t2 + "/");
		
		List<String> num_facts = new LinkedList<>();
		num_facts.add("run_id g1_size g2_size P_Rew P_rew_std diff_IAs");
		
		Set<Set<String>> g1_powerset = Utilities.getPowerSet(g1.keySet(), 2, 999);
		Set<Set<String>> g2_powerset = Utilities.getPowerSet(g2.keySet(), 2, 999);
		int overall_comp = g1_powerset.size()*g2_powerset.size();
		
		System.out.println("Compute " + overall_comp + " combinations");
		
		int i = 0;
		int overall = 1;
		for (Set<String> s1:g1_powerset) {
			Map<String, RewiringDetectorSample> g1s = getSubset(g1, s1);
			i++;
			int j = 0;
			for (Set<String> s2:g2_powerset) {
				Map<String, RewiringDetectorSample> g2s = getSubset(g2, s2);
				j++;
				String run_id = i+"_"+j;
				Utilities.writeEntries(g1s.keySet(), results_folder + run_id + "_g1s.txt.gz");
				Utilities.writeEntries(g2s.keySet(), results_folder + run_id + "_g2s.txt.gz");
				
				System.out.println("start RD calculation #" + overall + "/" + overall_comp +  " for " + run_id + " (" + g1s.size()*g2s.size() + " comparisons)");
				RewiringDetector rd = new RewiringDetector(g1s, g2s, fdr, no_threads, null, false, true);
				overall++;
				
				/*
				 * collect some facts and build output
				 */
				
				double P_rew = rd.getP_rew();
				double P_rew_std = rd.getP_rew_std();
				String out_temp = String.join(" ", new String[]{run_id, Integer.toString(g1s.size()), Integer.toString(g2s.size()), 
						Double.toString(P_rew), Double.toString(P_rew_std), Integer.toString(rd.getSignificantlyRewiredInteractions().size())});
				System.out.println(out_temp);
				System.out.flush();
				num_facts.add(out_temp);
				
				rd.writeDiffnet(results_folder + run_id + "_" + fdr + "_diffnet.txt.gz");
			}
			Utilities.writeEntries(num_facts, results_folder + "facts.txt");
			System.out.flush();
			System.gc();
		}
	}
	
	public static void main(String[] args) {
		
		System.out.println("hemato_rand_diff_subnetworks on " + new Date());
		System.out.println("no_threads:" + no_threads + ", FDR(s): 0.05");
		
		new File(results_root).mkdir();
		
		
		process(network_folder + "0.3/", results_root + "0.3/");
		
		System.out.println();
		System.out.println();
		
		process(network_folder + "0.35/", results_root + "0.35/");
	}
}
