package hemato_PPI;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.RewiringDetector;
import framework.Utilities;

public class build_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Desktop/BLUEPRINT_networks_0.03125/";
	static String results_root = "/Users/tho/Desktop/test/";
	
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
	
	public static void main(String[] args) {
		
		// define relations
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"HSC", "MPP"});
		relations.add(new String[]{"MPP", "CMP"});
		relations.add(new String[]{"MPP", "CLP"});
		// TODO: finish
		
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			Map<String, ConstructedNetworks> g1 = readNetworks(network_folder + state1 + "/");
			Map<String, ConstructedNetworks> g2 = readNetworks(network_folder + state2 + "/");
			
			RewiringDetector rd = new RewiringDetector(g1, g2, 0.05);
			rd.writeDiffnet("/Users/tho/Desktop/" + state1 + "_" + state2 + ".txt");
		}
		
	}
}
