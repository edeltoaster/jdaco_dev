package diff_compl_ENC;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffSeedCombVarDetector;
import framework.GOAnnotator;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class check_MCF7_diff_complexes {
	
	public static void main(String[] args) {
		String compl_folder = "BRCA_complex_results_99_5_-25-25/";
		String tfc_folder = "BRCA_tfc_results_99_5_-25-25/";
		pluri_definitions.qvalue = 0.05;
		pluri_definitions.goa = new GOAnnotator("mixed_data/stem_tags_retrieved.txt.gz");
		
		pluri_definitions.printParameters();
		
		System.out.println();

		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(pluri_definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), pluri_definitions.seed, pluri_definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("MCF-7") )
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		
		
		System.out.println("all other samples : " + group1.size());
		System.out.println("MCF-7 samples : " + group2.size());
		
		System.out.println();
		
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedCombVarDetector dsvd = new DiffSeedCombVarDetector(group1, group2, pluri_definitions.qvalue, pluri_definitions.parametric, pluri_definitions.paired, pluri_definitions.check_supersets, pluri_definitions.min_variant_fraction, pluri_definitions.no_threads);
		dsvd.diffTFComplAnalysis(tfc_folder, pluri_definitions.goa, pluri_definitions.binding_data, 0.0001, pluri_definitions.d_min, pluri_definitions.d_max, true, null, null);
		
		System.out.println("Determining enriched TFs ...");
		DiffSeedCombVarDetector.SPEnrichment tf_enrich = dsvd.calculateSPEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(tfc_folder + "enriched_pos_TFs.txt", tfc_folder + "enriched_neg_TFs.txt");
		
		System.out.println();
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, pluri_definitions.qvalue, pluri_definitions.parametric, pluri_definitions.paired, pluri_definitions.check_supersets, pluri_definitions.min_variant_fraction, pluri_definitions.no_threads);
		dcd.diffTFComplAnalysis(compl_folder, pluri_definitions.goa, pluri_definitions.binding_data, 0.0001, pluri_definitions.d_min, pluri_definitions.d_max, true, null, null);
		
		System.out.println("Determine enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(compl_folder + "enriched_pos_TFCs.txt", compl_folder + "enriched_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(compl_folder + "enriched_pos_TFs.txt", compl_folder + "enriched_neg_TFs.txt");
	}
}
