package hemato_PPI;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.Utilities;
import framework.RewiringDetector;

public class rand_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/";
	static String results_root = "/Users/tho/Desktop/BLUEPRINT_diffnets/";
	
	public static Map<String, ConstructedNetworks> readNetworks(String folder) {
		Map<String, ConstructedNetworks> data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(folder, "-ppin.txt.gz")) {
			String pre = f.getAbsolutePath().split("-ppin")[0];
			String donor = f.getName().split("-")[2];
			ConstructedNetworks cn = new ConstructedNetworks(pre + "-ppin.txt.gz", pre + "-ddin.txt.gz", pre + "-map.txt.gz", "homo_sapiens_core_81_38", true);
			
			data.put(donor, cn);
		}
		
		return data;
	}
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder);
		
		for (File f:new File(network_folder).listFiles()) {
			
			if (!f.isDirectory())
				continue;
			
			System.out.println(f.getAbsolutePath());
			String state1 = f.getName();
			String path = f.getAbsolutePath() + "/";
			System.out.println("Processing " + state1);
			Map<String, ConstructedNetworks> all_networks = readNetworks(path);
			List<String> network_list = new ArrayList<>(all_networks.keySet());
			
			// randomize
			for (int i=0; i<3;i++) {
				Map<String, ConstructedNetworks> g1 = new HashMap<>();
				Map<String, ConstructedNetworks> g2 = new HashMap<>();
				Collections.shuffle(network_list);
				
				for (int j=0;j<network_list.size();j++) {
					String s = network_list.get(j);
					if (j%2==0)
						g1.put(s, all_networks.get(s));
					else
						g2.put(s, all_networks.get(s));
				}
				RewiringDetector rd = new RewiringDetector(g1, g2, FDR);
				double P_rew_rounded = (double)Math.round(Utilities.getMean(rd.getP_rews().values()) * 1000d) / 1000d;
				System.out.println(rd.getP_rews().size()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getInteractionReasonsMap().size() + " dIAs" );
			}
			
			System.out.println();
		}
	}
	
	public static void main(String[] args) {
		
		new File(results_root).mkdir();
		
		for (File f:new File(network_folder).listFiles()) {
			
			if (!f.isDirectory())
				continue;
			
			String threshold_results = f.getName();
			process(f.getAbsolutePath() + "/", results_root + threshold_results + "/");
			
		}
		
	}
}
