package PPIComp_eval;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.Utilities;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;

public class build_rand_diff_subnetworks {
	
	static int no_threads = 48;
	static double FDR = 0.05;
	static double fraction = 0.5;
	static int min_size = 3;
	static int iterations = 5;
	
	// needs to be run on server
	static String network_folder = "BRCA_networks/";
	static String results_root = "BRCA_diffnets/";
	
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
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder + ", " + fraction + " fraction used.");
		System.out.println();
		
		// define relations
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"normal", "tumor"});
		
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			System.out.println("read data");
			Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + state1 + "/");
			Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + state2 + "/");
			
			for (int i=1;i<=iterations;i++) {
				//String run_id = state1 + "_" + state2 + "_" + i;
				
				Map<String, RewiringDetectorSample> g1s = getRandomSubset(g1, min_size, fraction);
				Map<String, RewiringDetectorSample> g2s = getRandomSubset(g2, min_size, fraction);
				
				System.out.println("start RD calculations for " + i + "(" + g1s.size()*g2s.size() + "comparisons).");
				RewiringDetector rd = new RewiringDetector(g1s, g2s, FDR, no_threads, System.out, false, true);
				
				double P_rew_rounded = (double) Math.round( rd.getP_rew() * 1000d) / 1000d;
				System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getSignificantlyRewiredInteractions().size() + " dIAs" );
			}

		}
		
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		new File(results_root).mkdir();
		
		for (File f:Utilities.listDirectoriesAndFilesWithinFolder(new File(network_folder))) {
			
			if (!f.isDirectory())
				continue;
			
			String threshold_results = f.getName();
			process(f.getAbsolutePath() + "/", results_root + threshold_results + "_rand/");
			
		}
		
	}
}
