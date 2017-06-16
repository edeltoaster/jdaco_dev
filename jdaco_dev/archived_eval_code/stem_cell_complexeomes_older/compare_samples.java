package stem_cell_complexeomes_older;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DACOResultSet;
import framework.Utilities;


public class compare_samples {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPIhc_TPMgene/res5/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		
		// qualitative abundance
		Map<String, DACOResultSet> data = new HashMap<>();
		Set<HashSet<String>> reference_seed_variants = new HashSet<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), seed);
			data.put(sample, dr);
			reference_seed_variants.addAll(dr.getSeedToComplexMap().keySet());
		}
		
		for (String s1:data.keySet())
			for (String s2:data.keySet()) {
				if (s1.compareTo(s2) >= 0) 
					continue;
				System.out.println(s1 + " vs " + s2 + " : " + data.get(s1).getSeedVariantSetsDistance(reference_seed_variants, data.get(s2)));
			}
		
		// quantitative abundance
		// TODO: implement analysis?
	}
}
