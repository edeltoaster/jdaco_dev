package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DACOResultSet;
import framework.Utilities;


public class compare_PrePPI_results_server {
	
	static String hc_results_folder = "DACO_PrePPIhc_TPMgene/res5/";
	static String all_results_folder = "DACO_PrePPI_TPMgene/res5/";
	static Set<String> seed = Utilities.readEntryFile("hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		
		System.out.println("reading hc data");
		Map<String, DACOResultSet> hc_data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(hc_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), seed);
			hc_data.put(sample, dr);
		}
		
		System.out.println("reading complete data");
		Map<String, DACOResultSet> all_data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(all_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), seed);
			all_data.put(sample, dr);
		}
		
		System.out.flush();
		
		System.out.println();
		System.out.println("Similarities:");
		for (String sample:hc_data.keySet()) {
			System.out.println(sample + " : " + hc_data.get(sample).getComplexSetsSimilarity(all_data.get(sample)));
		}
		
	}
}
