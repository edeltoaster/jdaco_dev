package diff_compl_BRCA;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DiffComplexDetector;
import framework.DiffSeedVarDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_BRCA_sub_complexes {

	static Map<String, String> subtype_map = BRCA_definitions.readSubtypes();
	
	public static void process(String subtype, Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2) {
		System.out.println("Processing BRCA subtype " + subtype);
		Map<String, QuantDACOResultSet> group1sub = new HashMap<>(group1);
		Map<String, QuantDACOResultSet> group2sub = new HashMap<>(group2);
		
		group1sub.keySet().removeIf(p -> !subtype_map.getOrDefault(p, "NA").equals(subtype));
		group2sub.keySet().removeIf(p -> !subtype_map.getOrDefault(p, "NA").equals(subtype));
		
		System.out.println(subtype + " normal samples : " + group1sub.size());
		System.out.println(subtype + " tumor samples : " + group2sub.size());
		
		String tfc_out = BRCA_definitions.diff_tfc_output_folder.substring(0, BRCA_definitions.diff_tfc_output_folder.length()-1) + "_" + subtype + "/";
		String compl_out = BRCA_definitions.diff_complex_output_folder.substring(0, BRCA_definitions.diff_complex_output_folder.length()-1) + "_" + subtype + "/";
		
		System.out.println();
		
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedVarDetector dsvd = new DiffSeedVarDetector(group1sub, group2sub, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.check_supersets, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		dsvd.diffTFComplAnalysis(tfc_out, BRCA_definitions.goa, BRCA_definitions.binding_data, 0.0001, BRCA_definitions.d_min, BRCA_definitions.d_max, true, null, null);
		
		System.out.println("Determining enriched TFs ...");
		DiffSeedVarDetector.SPEnrichment tf_enrich = dsvd.calculateSPEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(tfc_out + "enriched_pos_TFs.txt", tfc_out + "enriched_neg_TFs.txt");
		
		System.out.println();
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1sub, group2sub, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.check_supersets, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		dcd.diffTFComplAnalysis(compl_out, BRCA_definitions.goa, BRCA_definitions.binding_data, 0.0001, BRCA_definitions.d_min, BRCA_definitions.d_max, true, null, null);
		
		System.out.println("Determining enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(compl_out + "enriched_pos_TFCs.txt", compl_out + "enriched_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(compl_out + "enriched_pos_TFs.txt", compl_out + "enriched_neg_TFs.txt");
	
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		BRCA_definitions.printParameters();

		BRCA_definitions.goa.printTagInformation();
		
		System.out.println();
		
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(BRCA_definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String sample_prefix = sample.split("_")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), BRCA_definitions.seed, BRCA_definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("normal"))
				group1.put(sample_prefix, qdr);
			else {
				group2.put(sample_prefix, qdr);
			}
		}
		
		Set<String> subtypes = new HashSet<>(subtype_map.values());
		
		for (String subtype:subtypes)
			process(subtype, group1, group2);
	}
}
