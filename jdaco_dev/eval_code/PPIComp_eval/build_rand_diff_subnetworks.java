package PPIComp_eval;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.Utilities;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.StrPair;

public class build_rand_diff_subnetworks {
	
	static int no_threads = 48;
	static double FDR = 0.05;
	static List<Double> fractions = Arrays.asList(0.01, 0.05, 0.1, 0.3, 0.5);
	static int min_size = 3;
	static int iterations = 3;
	
	// needs to be run on server
	static String network_folder = "BRCA_networks/1.0/";
	static String results_root = "BRCA_rand_diffnets/";
	
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
		
		new File(results_folder).mkdir();
		
		System.out.println("read data");
		Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + "normal/");
		Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + "tumor/");
		
		for (double fraction:fractions) {
			
			System.out.println("Analysis for " + fraction + " fraction:");
			List<String> num_facts = new ArrayList<>(iterations);
			num_facts.add("run_id g1_size g2_size P_Rew P_rew_std diff_IAs");
			
			for (int i=1;i<=iterations;i++) {
				String run_id = fraction + "_" + i;
				
				Map<String, RewiringDetectorSample> g1s = getRandomSubset(g1, min_size, fraction);
				Map<String, RewiringDetectorSample> g2s = getRandomSubset(g2, min_size, fraction);
				
				System.out.println("start RD calculations for " + i + " (" + g1s.size()*g2s.size() + " comparisons).");
				RewiringDetector rd = new RewiringDetector(g1s, g2s, FDR, no_threads, System.out, false, true);
				
				/*
				 * some output
				 */
				
				double P_rew_rounded = (double) Math.round( rd.getP_rew() * 1000d) / 1000d;
				System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getSignificantlyRewiredInteractions().size() + " dIAs" );
				
				
				/*
				 * collect some facts
				 */
				
				double P_rew = rd.getP_rew();
				double P_rew_std = rd.getP_rew_std();
				num_facts.add(String.join(" ", new String[]{run_id, Integer.toString(g1.size()), Integer.toString(g2.size()), 
						Double.toString(P_rew), Double.toString(P_rew_std), Integer.toString(rd.getSignificantlyRewiredInteractions().size())}));
				
				List<String> temp = new LinkedList<>();
				for (StrPair pair:rd.getSignificantlyRewiredInteractions())
					temp.add(pair.getL() + pair.getR());
				Utilities.writeEntries(temp, results_folder + run_id + "_rew_IAs.txt.gz");
				
			}
			
			Utilities.writeEntries(num_facts, results_folder + "facts_" + fraction + ".txt");
			System.out.println();
		}
	}
	
	public static void main(String[] args) {
		
		System.out.println("rand_diff_subnetworks on " + new Date());
		System.out.println("no_threads:" + no_threads + ", FDR:" + FDR + ", min_size:" + min_size);
		System.out.println("fractions:" + fractions + ", iterations:" + iterations);
		
		new File(results_root).mkdir();

		process(network_folder, results_root);
	}
}
