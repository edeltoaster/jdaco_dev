package PPIComp_eval;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import framework.Utilities;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.StrPair;

public class build_rand_diff_subnetworks {
	
	static int no_threads = 48;
	static double FDR = 0.05;
	static int min_size = 3;
	static Map<Double, Integer> fraction_iteration = new TreeMap<>();
	
	static HashMap<Set<String>, Set<Set<String>>> sample_map;
	
	// needs to be run on server
	static String network_folder = "BRCA_networks/1.0/";
	static String results_root = "BRCA_rand_diffnets/";
	
	public static void defineParameters() {
		fraction_iteration.put(0.0025, 100); // 3 / 3
		fraction_iteration.put(0.005, 100); // 3 / 5
		fraction_iteration.put(0.01, 100); // 3 / 10
		fraction_iteration.put(0.05, 100); // 5 / 54
		fraction_iteration.put(0.1, 100); // 11 / 109
		fraction_iteration.put(0.25, 100); // 28 / 273
		fraction_iteration.put(0.5, 100); // 56 / 547
	}
	
	public static Map<String, RewiringDetectorSample> getRandomSubset(Map<String, RewiringDetectorSample> data, int min_size, double fraction) {
		Map<String, RewiringDetectorSample> subset = new HashMap<>();
		int sample_size = Math.max(min_size, (int) (fraction * data.size()));
		
		List<String> samples = new ArrayList<>(data.keySet());
		Collections.shuffle(samples);
		
		for (int i=0;i<sample_size;i++) {
			String sample = samples.get(i);
			subset.put(sample, data.get(sample));
		}
		
		return subset;
	}
	
	public static void process(String network_folder, String results_folder) {
		
		System.out.println("read data");
		Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + "normal/");
		Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + "tumor/");
		
		List<String> num_facts = new LinkedList<>();
		num_facts.add("fraction iteration g1_size g2_size P_Rew P_rew_std diff_IAs");
		
		for (Entry<Double, Integer> entry:fraction_iteration.entrySet()) {
			double fraction = entry.getKey();
			int iterations = entry.getValue();
			System.out.println();
			System.out.println("Analysis for " + fraction + " fraction:");
			sample_map = new HashMap<>(); // reset
			long start = System.currentTimeMillis();
			for (int i=1;i<=iterations;i++) {
				String run_id = fraction + "_" + i;
				
				// check that every pair never occurs twice
				Map<String, RewiringDetectorSample> g1s = null;
				Map<String, RewiringDetectorSample> g2s = null;
				
				while (g1s == null || g2s == null) {
					g1s = getRandomSubset(g1, min_size, fraction);
					g2s = getRandomSubset(g2, min_size, fraction);
					
					if (sample_map.containsKey(g1s.keySet()))
						if (sample_map.get(g1s.keySet()).contains(g2s.keySet())) {
							g1s = null;
							g2s = null;
						}
				}
				
				sample_map.getOrDefault(g1s.keySet(), new HashSet<>()).add(g2s.keySet());
				
				System.out.println("start RD calculations for " + i + " (" + g1s.size()*g2s.size() + " comparisons).");
				RewiringDetector rd = new RewiringDetector(g1s, g2s, FDR, no_threads, null, false, true);
				
				/*
				 * collect some facts and build output
				 */
				
				double P_rew = rd.getP_rew();
				double P_rew_std = rd.getP_rew_std();
				String out_temp = String.join(" ", new String[]{Double.toString(fraction), Integer.toString(i), Integer.toString(g1s.size()), Integer.toString(g2s.size()), 
						Double.toString(P_rew), Double.toString(P_rew_std), Integer.toString(rd.getSignificantlyRewiredInteractions().size())});
				System.out.println(out_temp);
				System.out.flush();
				num_facts.add(out_temp);
				
				List<String> temp = new LinkedList<>();
				for (StrPair pair:rd.getSignificantlyRewiredInteractions())
					temp.add(pair.getL() + " " + pair.getR());
				Utilities.writeEntries(temp, results_folder + run_id + "_rew_IAs.txt.gz");
				Utilities.writeEntries(g1s.keySet(), results_folder + run_id + "_g1s.txt.gz");
				Utilities.writeEntries(g2s.keySet(), results_folder + run_id + "_g2s.txt.gz");
			}
			
			long duration_minutes = (System.currentTimeMillis() - start) / 1000 / 60;
			
			Utilities.writeEntries(num_facts, results_folder + "facts.txt");
			System.out.println(iterations + " iterations for " + fraction + " took " + duration_minutes + " min.");
			System.out.flush();
			System.gc();
		}
	}
	
	public static void main(String[] args) {
		
		defineParameters();
		
		System.out.println("rand_diff_subnetworks on " + new Date());
		System.out.println("no_threads:" + no_threads + ", FDR:" + FDR + ", min_size:" + min_size);
		System.out.println("fractions/iterations: " + fraction_iteration);
		
		new File(results_root).mkdir();

		process(network_folder, results_root);
	}
}
