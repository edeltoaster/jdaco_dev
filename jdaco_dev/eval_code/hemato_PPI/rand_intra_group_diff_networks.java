package hemato_PPI;

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

public class rand_intra_group_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/";
	static String results_root = "/Users/tho/Desktop/BLUEPRINT_diffnets/";
	static int iterations = 5;
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder);
		
		for (File f:new File(network_folder).listFiles()) {
			
			if (!f.isDirectory())
				continue;
			
			String state1 = f.getName();
			String path = f.getAbsolutePath() + "/";
			
			Map<String, RewiringDetectorSample> all_networks = RewiringDetectorSample.readNetworks(path);
			List<String> network_list = new ArrayList<>(all_networks.keySet());
			
			System.out.println("Processing " + state1 + ", " + network_list.size() + " samples");
			
			// randomize
			List<Double> P_rews = new LinkedList<>();
			List<Double> diff_IAs = new LinkedList<>();
			for (int i=0; i<iterations; i++) {
				Map<String, RewiringDetectorSample> g1 = new HashMap<>();
				Map<String, RewiringDetectorSample> g2 = new HashMap<>();
				Collections.shuffle(network_list);
				
				for (int j=0;j<network_list.size();j++) {
					String s = network_list.get(j);
					if (j%2==0)
						g1.put(s, all_networks.get(s));
					else
						g2.put(s, all_networks.get(s));
				}
				
				RewiringDetector rd = new RewiringDetector(g1, g2, FDR);
				P_rews.add( rd.getP_rew() );
				diff_IAs.add( (double) rd.getInteractionReasonsCountMap().keySet().size());
			}
			
			double P_rew_rounded = (double) Math.round(Utilities.getMean(P_rews) * 1000d) / 1000d;
			double P_rew_std_rounded = (double) Math.round(Utilities.getStd(P_rews) * 1000d) / 1000d;
			double IAs_rounded = (double) Math.round(Utilities.getMean(diff_IAs) * 10d) / 10d;
			double IAs_std_rounded = (double) Math.round(Utilities.getStd(diff_IAs) * 10d) / 10d;
			System.out.println("P_rew: " + P_rew_rounded + "+-" + P_rew_std_rounded + " , " + IAs_rounded + "+-" + IAs_std_rounded);
			System.out.println();
		}
	}
	
	public static void main(String[] args) {
		
		new File(results_root).mkdir();
		
		System.out.println(iterations + " iterations");
		for (File f:new File(network_folder).listFiles()) {
			
			if (!f.isDirectory())
				continue;
			
			String threshold_results = f.getName();
			process(f.getAbsolutePath() + "/", results_root + threshold_results + "/");
			
		}
		
	}
}
