package hemato_PPI;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.Utilities;
import framework.RewiringDetector;

public class build_diff_networks_stability_eval {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/";
	static int max_iter = 5;
	
	public static List<List<String>> permutateSamples(Map<String, ConstructedNetworks> data) {
		
		List<String> samples = new ArrayList<>(data.keySet());
		List<List<String>> sample_permutations = new LinkedList<>();
		
		for (List<Integer> perm:Utilities.getAllIntPermutations(samples.size())) {
			List<String> temp = new LinkedList<>();
			for (int i:perm) {
				temp.add(samples.get(i));
			}
			sample_permutations.add(temp);
		}
		
		return sample_permutations;
	}
	
	public static List<Map<String, ConstructedNetworks>> buildPermutatedData(String folder) {
		Map<String, ConstructedNetworks> all_data = ConstructedNetworks.readNetworks(folder);
		
		List<Map<String, ConstructedNetworks>> output = new LinkedList<>();
		
		// get sublists
		permutateSamples(all_data);
		
		return output;
	}
	
	public static Map<String, ConstructedNetworks> filterVenous(Map<String, ConstructedNetworks> input) {
		Map<String, ConstructedNetworks> filtered = new HashMap<String, ConstructedNetworks>();
		for (String sample:input.keySet()) {
			if (sample.startsWith("Venous"))
				continue;
			filtered.put(sample, input.get(sample));
		}
		return filtered;
	}
	
	public static void process(String network_folder) {
		
		System.out.println("Stability analysis for " + network_folder);
		System.out.println();
		
		// define relations
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"HSC", "MPP"});
		
		relations.add(new String[]{"MPP", "CMP"});
		relations.add(new String[]{"MPP", "CLP"});
		
		relations.add(new String[]{"CMP", "MEP"});
		relations.add(new String[]{"CMP", "GMP"});
		
		relations.add(new String[]{"MEP", "MK"});
		relations.add(new String[]{"MEP", "EB"});
		
		relations.add(new String[]{"GMP", "N"});
		relations.add(new String[]{"GMP", "M"});
		
		relations.add(new String[]{"CLP", "NK"});
		relations.add(new String[]{"CLP", "CD4"});
		
		// additional relations
		relations.add(new String[]{"CMP", "CLP"});
		relations.add(new String[]{"MEP", "GMP"});
		relations.add(new String[]{"MK", "EB"});
		relations.add(new String[]{"N", "M"});
		relations.add(new String[]{"NK", "CD4"});
		
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			Map<String, ConstructedNetworks> g1 = ConstructedNetworks.readNetworks(network_folder + state1 + "/");
			Map<String, ConstructedNetworks> g2 = ConstructedNetworks.readNetworks(network_folder + state2 + "/");
			
			System.out.print("Processing " + state1 + " (" + g1.keySet().size() + ") vs " + state2 + " (" + g2.keySet().size() + ") : ");
			RewiringDetector rd = new RewiringDetector(g1, g2, FDR, 4);
			System.out.println(rd.getSignificantlyRewiredInteractions().size());
			

			for (Map<String, ConstructedNetworks> s1 :buildPermutatedData(network_folder + state1 + "/"))
				for (Map<String, ConstructedNetworks> s2 :buildPermutatedData(network_folder + state2 + "/")) {
					rd = new RewiringDetector(s1, s2, FDR, 4);
					
					System.out.println(g1.keySet().size() + " " + g2.keySet().size() + " " + rd.getSignificantlyRewiredInteractions().size());
				}	
		}
		
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		for (File f:new File(network_folder).listFiles()) {
			
			if (!f.isDirectory() || f.getName().endsWith("0.0"))
				continue;
			
			process(f.getAbsolutePath() + "/");
			
		}
		
	}
}
