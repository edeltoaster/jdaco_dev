package hemato_PPI;

import java.io.File;
import java.util.HashMap;
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
		
		// read all data
		Map<String, ConstructedNetworks> HSC = readNetworks(network_folder + "MEP/");
		Map<String, ConstructedNetworks> MPP = readNetworks(network_folder + "MK/");
		
		RewiringDetector rd = new RewiringDetector(HSC, MPP, 0.05, results_root);
		rd.writeDiffnetAndReasons("/Users/tho/Desktop/diffnet.txt", "/Users/tho/Desktop/reasons.txt");
	}
}
