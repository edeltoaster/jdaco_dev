package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DACOResultSet;
import framework.Utilities;


public class compare_samples_server {
	
	static String daco_results_folder = "DACO_PrePPI_TPMgene/res5/";
	static Set<String> seed = Utilities.readEntryFile("hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "32");
		Map<String, DACOResultSet> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), seed);
			data.put(sample, dr);
		}
		
		for (String s1:data.keySet())
			for (String s2:data.keySet()) {
				if (s1.compareTo(s2) >= 0) 
					continue;
				System.out.println(s1 + " vs " + s2 + " : " + data.get(s1).getComplexSetsSimilarity(data.get(s2)));
			}
	}
}
