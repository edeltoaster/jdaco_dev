package ENCODE_pluri_complexomes;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class build_diff_complexes { // intended to be run on a server

	static String seed_file = "mixed_data/hocomoco_human_core_TFs_v11.txt.gz";
	static String out_folder_pre = "diff_compl_";
	static double FDR = 0.05;
	static double[] min_variant_fractions = new double[] {0.75, 0.9, 0.8, 0.5};
	
	static String preppi_network_folder = "preppi_networks/";
	static String preppi_compl_folder = "res_preppi/";
	
	static String mentha_network_folder = "mentha_networks/";
	static String mentha_compl_folder = "res_mentha/";
	
	static String wmentha_network_folder = "wmentha_networks/";
	static String wmentha_compl_folder = "res_wmentha/";
	
	public static void process_data(String spec_out_folder, String network, double min_variant_fraction) {
		System.out.println("Processing " + network);
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
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, FDR, false, false, false, min_variant_fraction, 64, true);
		dcd.writeSignSortedComplexes(spec_out_folder + network + "_sign.txt", false);
		dcd.writeSignSortedComplexes(spec_out_folder + network + "_signh.txt", true);
		System.out.println(dcd.getSignificanceSortedComplexes().size() + " sign complexes");
		
		dcd.writeSignSortedVariants(spec_out_folder + network + "_vsign.txt", false);
		dcd.writeSignSortedVariants(spec_out_folder + network + "_vsignh.txt", true);
		System.out.println(dcd.getSignSortedVariants(false, false).size() + " sign variants");
		
		System.out.println("Determine enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(FDR, 10000, 10);
		tfc_enrich.writeSignificantSeedProteinCombinations(spec_out_folder + network + "_TFCs.txt");
		tfc_enrich.writeSignificantSeedProteinCombinations(spec_out_folder + network + "_pos_TFCs.txt", spec_out_folder + network + "_neg_TFCs.txt");
		System.out.println(tfc_enrich.getSignificanceSortedSeedProteinCombinations().size() + " sign seed protein combinations");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(FDR, 10000, 10);
		tf_enrich2.writeSignificantSeedProteins(spec_out_folder + network + "_TFs.txt");
		tf_enrich2.writeSignificantSeedProteins(spec_out_folder + network + "_pos_TFs.txt", spec_out_folder + network + "_neg_TFs.txt");
		System.out.println(tf_enrich2.getSignificanceSortedSeedProteins() + " sign seed proteins");
		
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		for (double min_var:min_variant_fractions) {
			System.out.println("Calculations for min_variant_fraction " + min_var);
			
			File out_dir = new File(out_folder_pre + min_var + "/");
			
			if (!out_dir.exists())
				out_dir.mkdir();
			
			process_data(out_dir.getAbsolutePath(), "wmentha", min_var);
			process_data(out_dir.getAbsolutePath(), "mentha", min_var);
			process_data(out_dir.getAbsolutePath(), "preppi", min_var);
			
			System.out.println();
			System.out.println();
		}

	}
}
