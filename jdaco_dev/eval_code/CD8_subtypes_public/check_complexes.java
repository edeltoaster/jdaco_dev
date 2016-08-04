package CD8_subtypes_public;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DACOResultSet;
import framework.GOAnnotator;
import framework.Utilities;


public class check_complexes {
	
	static String results_folder = "/Users/tho/Dropbox/Work/projects/CD8_subsets_public/CD8_DACO/res7/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		
		Map<String, DACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			//String cell_type = sample.split("_")[1];
			
			//if (cell_type.equals("N") || cell_type.equals("TMNP"))
			results.put(sample, new DACOResultSet(f.getAbsolutePath(), seed));
		}
		
//		for (String sample1:results.keySet())
//			for (String sample2:results.keySet()) {
//				System.out.println(sample1 + " / " + sample2 + " : " + results.get(sample1).getComplexSimilarity(results.get(sample2)) + " " + results.get(sample1).getSeedVariantSimilarity(results.get(sample2)));
//			}
		
		GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/simple_tags_retrieved.txt.gz");
		for (String sample1:results.keySet()) {
			DACOResultSet daco_complexes = results.get(sample1);
			for (HashSet<String> tfs:daco_complexes.getSeedToComplexMap().keySet()) {
				for (HashSet<String> complex:daco_complexes.getSeedToComplexMap().get(tfs)) {
					
					System.out.println(sample1 + " : " + tfs + " -> " + complex + " : " + goa.rateProteins(complex));
				}
			}
		}
	}
}
