package diff_compl_ENC;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class check_pluri_diff_complexes {
	
	public static void main(String[] args) {
		Set<String> allosome_proteins = DataQuery.getAllosomeProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		pluri_definitions.printParameters();
		
		pluri_definitions.goa.printTagInformation();
		
		System.out.println();

		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(pluri_definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), pluri_definitions.seed, pluri_definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("H1-hESC") && !sample.contains("H7-hESC") && !sample.contains("induced-pluripotent-stem-cell"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		
		
		System.out.println("non-pluri samples : " + group1.size());
		System.out.println("pluri samples : " + group2.size());
		
		System.out.println();
		
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, pluri_definitions.qvalue, pluri_definitions.parametric, pluri_definitions.paired, pluri_definitions.check_supersets, pluri_definitions.min_variant_fraction, pluri_definitions.no_threads);
		dsvd.diffTFComplAnalysis(pluri_definitions.output_folder_pre + "tfc/", pluri_definitions.goa, pluri_definitions.binding_data, 0.0001, pluri_definitions.d_min, pluri_definitions.d_max, true, allosome_proteins, pluri_definitions.pluri_factors);
		dsvd.writeSignSortedVariants(pluri_definitions.output_folder_pre + "tfc/pluri_sign.txt", true);
		
		System.out.println("Determining enriched TFs ...");
		DiffSeedCombDetector.SPEnrichment tf_enrich = dsvd.calculateSPEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(pluri_definitions.output_folder_pre + "tfc/" + "enriched_pos_TFs.txt", pluri_definitions.output_folder_pre + "tfc/" + "enriched_neg_TFs.txt");
		
		System.out.println();
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, pluri_definitions.qvalue, pluri_definitions.parametric, pluri_definitions.paired, pluri_definitions.check_supersets, pluri_definitions.min_variant_fraction, pluri_definitions.no_threads);
		dcd.diffTFComplAnalysis(pluri_definitions.output_folder_pre + "compl/", pluri_definitions.goa, pluri_definitions.binding_data, 0.0001, pluri_definitions.d_min, pluri_definitions.d_max, true, allosome_proteins, pluri_definitions.pluri_factors);
		dcd.writeSignSortedComplexes(pluri_definitions.output_folder_pre + "compl/pluri_sign.txt", true);
		
		System.out.println("Determine enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(pluri_definitions.output_folder_pre + "compl/" + "enriched_pos_TFCs.txt", pluri_definitions.output_folder_pre + "compl/" + "enriched_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(pluri_definitions.output_folder_pre + "compl/" + "enriched_pos_TFs.txt", pluri_definitions.output_folder_pre + "compl/" + "enriched_neg_TFs.txt");
		
		System.out.println();
		
		System.out.println("Determining diff. complexes (sub) ...");
		dcd = new DiffComplexDetector(group1, group2, pluri_definitions.qvalue, pluri_definitions.parametric, pluri_definitions.paired, true, pluri_definitions.min_variant_fraction, pluri_definitions.no_threads);
		dcd.diffTFComplAnalysis(pluri_definitions.output_folder_pre + "compl_sub/", pluri_definitions.goa, pluri_definitions.binding_data, 0.0001, pluri_definitions.d_min, pluri_definitions.d_max, true, allosome_proteins, pluri_definitions.pluri_factors);
		dcd.writeSignSortedComplexes(pluri_definitions.output_folder_pre + "compl_sub/pluri_sign.txt", true);
		
		System.out.println("Determine enriched TF combinations ...");
		tfc_enrich = dcd.calculateSPCEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(pluri_definitions.output_folder_pre + "compl_sub/" + "enriched_pos_TFCs.txt", pluri_definitions.output_folder_pre + "compl_sub/" + "enriched_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		tf_enrich2 = dcd.calculateSPEnrichment(pluri_definitions.qvalue, pluri_definitions.SPEnrich_iterations, pluri_definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(pluri_definitions.output_folder_pre + "compl_sub/" + "enriched_pos_TFs.txt", pluri_definitions.output_folder_pre + "compl_sub/" + "enriched_neg_TFs.txt");
	}
}
