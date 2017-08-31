package diff_compl_BRCA;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_BRCA_complexes {

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
		
		
		System.out.println("normal samples : " + group1.size());
		System.out.println("tumor samples : " + group2.size());
		
		System.out.println();
		
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.check_supersets, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		dsvd.diffTFComplAnalysis(BRCA_definitions.diff_tfc_output_folder, BRCA_definitions.goa, BRCA_definitions.binding_data, 0.0001, BRCA_definitions.d_min, BRCA_definitions.d_max, true, null, null);
		
		System.out.println("Determining enriched TFs ...");
		DiffSeedCombDetector.SPEnrichment tf_enrich = dsvd.calculateSPEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(BRCA_definitions.diff_tfc_output_folder + "enriched_pos_TFs.txt", BRCA_definitions.diff_tfc_output_folder + "enriched_neg_TFs.txt");
		
		System.out.println();
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.check_supersets, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		dcd.diffTFComplAnalysis(BRCA_definitions.diff_complex_output_folder, BRCA_definitions.goa, BRCA_definitions.binding_data, 0.0001, BRCA_definitions.d_min, BRCA_definitions.d_max, true, null, null);
		
		System.out.println("Determining enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(BRCA_definitions.diff_complex_output_folder + "enriched_pos_TFCs.txt", BRCA_definitions.diff_complex_output_folder + "enriched_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(BRCA_definitions.diff_complex_output_folder + "enriched_pos_TFs.txt", BRCA_definitions.diff_complex_output_folder + "enriched_neg_TFs.txt");
	}
}
