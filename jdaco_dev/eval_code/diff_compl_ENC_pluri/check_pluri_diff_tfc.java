package diff_compl_ENC_pluri;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.DiffSeedVarDetector;
import framework.DiffSeedVarDetector.SPEnrichment;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class check_pluri_diff_tfc {
	
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
		
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedVarDetector dsvd = new DiffSeedVarDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dsvd.diffTFComplAnalysis(definitions.diff_compl_output_folder, definitions.goa, definitions.binding_data, 0.0001, definitions.d_min, definitions.d_max, true, allosome_proteins, definitions.pluri_factors);
		
		System.out.println("Determining enriched TFs ...");
		SPEnrichment tf_enrich = dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(definitions.diff_compl_output_folder + "enriched_pos_TFs.txt", definitions.diff_compl_output_folder + "enriched_neg_TFs");
	}
}
