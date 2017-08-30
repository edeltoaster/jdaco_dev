package diff_compl_TCGA;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffSeedVarDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes {
	
	public static void process(String cancer_type) {
		System.out.println("Processing " +  cancer_type);
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String[] sample_spl = sample.split("_");
			String sample_pre = sample_spl[0] + "_" + sample_spl[1];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample_spl[0].equals(cancer_type))
				continue;
			
			if (sample.contains("normal"))
				group1.put(sample_pre, qdr);
			else {
				group2.put(sample_pre, qdr);
			}
		}
		
		System.out.println(cancer_type + " normal samples : " + group1.size());
		System.out.println(cancer_type + " tumor samples : " + group2.size());
		
		String tfc_out = definitions.diff_tfc_output_folder + cancer_type + "_";
		String compl_out = definitions.diff_complex_output_folder + cancer_type + "_";
		
		System.out.println();
		
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedVarDetector dsvd = new DiffSeedVarDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dsvd.diffTFComplAnalysis(tfc_out, definitions.goa, definitions.binding_data, 0.0001, definitions.d_min, definitions.d_max, true, null, null);
		dsvd.writeSignSortedVariants(tfc_out + "sign.txt", false);
		
		System.out.println("Determining enriched TFs ...");
		DiffSeedVarDetector.SPEnrichment tf_enrich = dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tf_enrich.writeSignificantSeedProteins(tfc_out + "enriched_pos_TFs.txt", tfc_out + "enriched_neg_TFs.txt");
		
		System.out.println();
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dcd.diffTFComplAnalysis(compl_out, definitions.goa, definitions.binding_data, 0.0001, definitions.d_min, definitions.d_max, true, null, null);
		dcd.writeSignSortedComplexes(compl_out + "sign.txt", false);
		
		System.out.println("Determining enriched TF combinations ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(compl_out + "enriched_pos_TFCs.txt", compl_out + "enriched_neg_TFCs.txt");
		
		System.out.println("Determining enriched TFs ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(compl_out + "enriched_pos_TFs.txt", compl_out + "enriched_neg_TFs.txt");
	
		// writing quantified results for later usage
		if (!new File(definitions.qr_output_folder).exists())
			new File(definitions.qr_output_folder).mkdir();
		
		for (String sample:group1.keySet())
			group1.get(sample).writeQuantifiedResult(definitions.qr_output_folder + sample + "_normal_qr.txt.gz");
		for (String sample:group2.keySet())
			group2.get(sample).writeQuantifiedResult(definitions.qr_output_folder + sample + "_tumor_qr.txt.gz");
		
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		definitions.printParameters();

		definitions.goa.printTagInformation();
		
		System.out.println();

		Map<String, Integer> cancer_types = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String cancer_type = f.getName().split("_")[0];
			cancer_types.put(cancer_type, cancer_types.getOrDefault(cancer_type, 0) + 1);
		}
		
		for (String cancer_type:cancer_types.keySet())
			process(cancer_type);
		
	}
}
