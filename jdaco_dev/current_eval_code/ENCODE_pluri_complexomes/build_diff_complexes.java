package ENCODE_pluri_complexomes;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class build_diff_complexes { // intended to be run on a server

	static String seed_file = "mixed_data/hocomoco_human_core_TFs_v11.txt.gz";
	static String out_folder = "diff_compl/";
	static double FDR = 0.05;
	
	static String preppi_network_folder = "preppi_networks/";
	static String preppi_compl_folder = "res_preppi/";
	
	static String mentha_network_folder = "mentha_networks/";
	static String mentha_compl_folder = "res_mentha/";
	
	public static void process_data(String network) {
		
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("res_" + network + "/", ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed_file, network + "_networks/" + sample + "_major-transcripts.txt.gz");
		
			if (!sample.contains("H1-hESC") && !sample.contains("H7-hESC") && !sample.contains("induced-pluripotent-stem-cell"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		
		System.out.println();
		System.out.println("non-pluri samples : " + group1.size());
		System.out.println("pluri samples : " + group2.size());
		System.out.println();
		
		System.out.println("Determine differential complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, FDR, false, false, false, 0.75, 64, true);
		dcd.writeSignSortedComplexes(out_folder + network + "_sign.txt", false);
		dcd.writeSignSortedComplexes(out_folder + network + "_signh.txt", true);
		dcd.writeSignSortedVariants(out_folder + network + "_vsign.txt", false);
		dcd.writeSignSortedVariants(out_folder + network + "_vsignh.txt", true);
		
		System.out.println("Determine enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(FDR, 10000, 10);
		tfc_enrich.writeSignificantSeedProteinCombinations(out_folder + network + "_TFCs.txt");
		tfc_enrich.writeSignificantSeedProteinCombinations(out_folder + network + "_pos_TFCs.txt", out_folder + network + "_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(FDR, 10000, 10);
		tf_enrich2.writeSignificantSeedProteins(out_folder + network + "_TFs.txt");
		tf_enrich2.writeSignificantSeedProteins(out_folder + network + "_pos_TFs.txt", out_folder + network + "_neg_TFs.txt");
		
		System.out.println();
	}
	
	public static void main(String[] args) {
		File out_dir = new File(out_folder);
		if (!out_dir.exists())
			out_dir.mkdir();
		
		process_data("mentha");
		
		process_data("preppi");
	}
}
