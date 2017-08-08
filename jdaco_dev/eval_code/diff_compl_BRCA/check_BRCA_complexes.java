package diff_compl_BRCA;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPEnrichment;
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
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.check_supersets, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		dcd.diffTFComplAnalysis(BRCA_definitions.diff_compl_output_folder, BRCA_definitions.goa, BRCA_definitions.binding_data, 0.0001, BRCA_definitions.d_min, BRCA_definitions.d_max, true, null, null);
		
		System.out.println("Determining enriched TFs ...");
		SPEnrichment tf_enrich = dcd.calculateSPEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(BRCA_definitions.diff_compl_output_folder + "enriched_pos_TFs.txt", BRCA_definitions.diff_compl_output_folder + "enriched_neg_TFs");
	}
}
