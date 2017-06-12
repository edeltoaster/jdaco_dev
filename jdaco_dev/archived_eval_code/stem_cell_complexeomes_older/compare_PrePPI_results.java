package stem_cell_complexeomes_older;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DACOResultSet;
import framework.Utilities;


public class compare_PrePPI_results {
	
	static String hc_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPIhc_TPMgene/res5/";
	static String all_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPI_TPMgene/res5/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		
		Set<HashSet<String>> reference_seed_variants = new HashSet<>();
		
		System.out.println("reading hc data");
		Map<String, DACOResultSet> hc_data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(hc_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), seed);
			hc_data.put(sample, dr);
			reference_seed_variants.addAll(dr.getSeedToComplexMap().keySet());
		}
		
		System.out.println("reading complete data");
		Map<String, DACOResultSet> all_data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(all_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), seed);
			all_data.put(sample, dr);
			reference_seed_variants.addAll(dr.getSeedToComplexMap().keySet());
		}
		
		System.out.flush();
		
		System.out.println();
		System.out.println("Similarities:");
		for (String sample:hc_data.keySet()) {
			System.out.println(sample + " : " + hc_data.get(sample).getResult().size() + " / " + all_data.get(sample).getResult().size());
			System.out.println(sample + " : " + hc_data.get(sample).getSeedToComplexMap().keySet().size() + " / " + all_data.get(sample).getSeedToComplexMap().keySet().size());
			System.out.println(sample + " : " + hc_data.get(sample).getSeedVariantSetsDistance(reference_seed_variants, all_data.get(sample)));
		}
		
	}
}
