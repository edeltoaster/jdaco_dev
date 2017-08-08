package diff_compl_ENC_pluri;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;
import framework.DiffComplexDetector.SPCEnrichment;
import framework.DiffComplexDetector.SPEnrichment;


public class check_pluri_diff_complexes {
	
	public static void main(String[] args) {
		Set<String> allosome_proteins = DataQuery.getAllosomeProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		definitions.printParameters();
		
		definitions.goa.printTagInformation();
		
		System.out.println();

		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("H1-hESC") && !sample.contains("H7-hESC") && !sample.contains("induced-pluripotent-stem-cell"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		
		System.out.println("non-pluri samples : " + group1.size());
		System.out.println("pluri samples : " + group2.size());
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dcd.diffTFComplAnalysis(definitions.diff_complex_output_folder, definitions.goa, definitions.binding_data, 0.0001, definitions.d_min, definitions.d_max, true, allosome_proteins, definitions.pluri_factors);
		
		System.out.println("Determine enriched TF combinations ...");
		SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(definitions.diff_complex_output_folder + "enriched_pos_TFCs.txt", definitions.diff_complex_output_folder + "enriched_neg_TFCs");
		
		System.out.println("Determining enriched TFs ...");
		SPEnrichment tf_enrich = dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(definitions.diff_complex_output_folder + "enriched_pos_TFs.txt", definitions.diff_complex_output_folder + "enriched_neg_TFs");
	}
}
