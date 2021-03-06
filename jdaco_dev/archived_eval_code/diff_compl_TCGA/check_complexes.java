package diff_compl_TCGA;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffSeedCombDetector;
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
		
		
		/*
		 * create output folder
		 */
		String compl_out = definitions.diff_complex_output_folder + cancer_type + "_";
		
		// first ensure that folder exists while allowing to put prefixes in like output_folder/bla_ ...
		File output_folder_obj = null;
		if (!compl_out.endsWith("/")) {
			String[] spl = compl_out.split("/");
			String actual_folder = String.join("/", Arrays.asList(spl).subList(0, spl.length-1));
			output_folder_obj = new File(actual_folder);
		} else 
			output_folder_obj = new File(compl_out);
			
		if (!output_folder_obj.exists())
			output_folder_obj.mkdir();
		
		System.out.println();
		
		
		/*
		 * conduct actual computations
		 */
		System.out.println("Determining diff. TF combinations ...");
		DiffSeedCombDetector dscd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dscd.writeSignSortedVariants(compl_out + "tfc_sign.txt.gz", false);
		dscd.writeSignSortedVariants(compl_out + "tfc_signh.txt.gz", true);
		
		System.out.println("Determining diff. complexes ...");
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dcd.writeSignSortedComplexes(compl_out + "compl_sign.txt.gz", false);
		dcd.writeSignSortedComplexes(compl_out + "compl_signh.txt.gz", true);
		dcd.writeSignSortedVariants(compl_out + "compl_vsign.txt.gz", false);
		dcd.writeSignSortedVariants(compl_out + "compl_vsignh.txt.gz", true);
		
		// in DiffSeedCompl
		System.out.println("Determining enriched TFs in diff. TF combinations ...");
		DiffSeedCombDetector.SPEnrichment tf_enrich1 = dscd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tf_enrich1.writeSignificantSeedProteins(compl_out + "tfc_pos_TFs.txt.gz", compl_out + "tfc_neg_TFs.txt.gz");
		
		// in DiffCompl
		System.out.println("Determining enriched TF combinations in diff. complexes ...");
		DiffComplexDetector.SPCEnrichment tfc_enrich = dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		tfc_enrich.writeSignificantSeedProteinCombinations(compl_out + "compl_pos_TFCs.txt.gz", compl_out + "compl_neg_TFCs.txt.gz");
		
		System.out.println("Determining enriched TFs in diff. complexes ...");
		DiffComplexDetector.SPEnrichment tf_enrich2 = dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		tf_enrich2.writeSignificantSeedProteins(compl_out + "compl_pos_TFs.txt.gz", compl_out + "compl_neg_TFs.txt.gz");
	
		
		/*
		 * write quantification results for later usage
		 */
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
		definitions.printInitParameters();

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
