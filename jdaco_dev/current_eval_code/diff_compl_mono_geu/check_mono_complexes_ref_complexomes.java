package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPCEnrichment;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_mono_complexes_ref_complexomes {

	static String diff_out_folder = "/Users/tho/Desktop/diff_results_corum/";
	static String res_folder = "/Users/tho/GDrive/Work/projects/CompleXChange/results/corum/";
	
	public static void main(String[] args) {
		definitions.printInitParameters();
		System.out.println();

		Map<String, QuantDACOResultSet> cm_data = new HashMap<>();
		Map<String, QuantDACOResultSet> ncm_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(res_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, "/Users/tho/GDrive/Work/projects/CompleXChange/results/mt_files/" + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				continue;
			else if (sample.startsWith("non"))
				ncm_data.put(sample, qdr);
			else
				cm_data.put(sample, qdr);
		}
		
		
		/**
		 * Diff. complexes (unpaired)
		 */
		
		new File(diff_out_folder).mkdir();
		System.out.println("Monocyte comparison cm->ncm (unpaired): " + cm_data.size() + " vs " + ncm_data.size());
		
		System.out.println("Diff. compl.");
		DiffComplexDetector mono_dcd = new DiffComplexDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcd.writeSignSortedComplexes(diff_out_folder + "mono_dcd_complh.txt", true);
		mono_dcd.writeSignSortedVariants(diff_out_folder + "mono_dcd_tfcsh.txt", true);
		mono_dcd.writeSignSortedComplexes(diff_out_folder + "mono_dcd_compl.txt", false);
		mono_dcd.writeSignSortedVariants(diff_out_folder + "mono_dcd_tfcs.txt", false);
		framework.DiffComplexDetector.SPEnrichment spe = mono_dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe.writeSignificantSeedProteins(diff_out_folder + "mono_dcd_SPenr.txt");
		SPCEnrichment spc = mono_dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spc.writeSignificantSeedProteinCombinations(diff_out_folder + "mono_dcd_SPCenr.txt");
		
		System.out.println("Diff. seed comb.");
		DiffSeedCombDetector mono_dsvd = new DiffSeedCombDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvd.writeSignSortedVariants(diff_out_folder + "mono_dsvcd_tfcsh.txt", true);
		mono_dsvd.writeSignSortedVariants(diff_out_folder + "mono_dsvcd_tfcs.txt", false);
		framework.DiffSeedCombDetector.SPEnrichment spe2 = mono_dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe2.writeSignificantSeedProteins(diff_out_folder + "mono_dsvcd_SPenr.txt");
		
		
		
		/**
		 * for hu.MAP
		 */
		diff_out_folder = "/Users/tho/Desktop/diff_results_humap/";
		res_folder = "/Users/tho/GDrive/Work/projects/CompleXChange/results/humap/";
		
		cm_data = new HashMap<>();
		ncm_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(res_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, "/Users/tho/GDrive/Work/projects/CompleXChange/results/mt_files/" + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				continue;
			else if (sample.startsWith("non"))
				ncm_data.put(sample, qdr);
			else
				cm_data.put(sample, qdr);
		}
		
		
		/**
		 * Diff. complexes (unpaired)
		 */
		
		new File(diff_out_folder).mkdir();
		System.out.println("Monocyte comparison cm->ncm (unpaired): " + cm_data.size() + " vs " + ncm_data.size());
		
		System.out.println("Diff. compl.");
		mono_dcd = new DiffComplexDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcd.writeSignSortedComplexes(diff_out_folder + "mono_dcd_complh.txt", true);
		mono_dcd.writeSignSortedVariants(diff_out_folder + "mono_dcd_tfcsh.txt", true);
		mono_dcd.writeSignSortedComplexes(diff_out_folder + "mono_dcd_compl.txt", false);
		mono_dcd.writeSignSortedVariants(diff_out_folder + "mono_dcd_tfcs.txt", false);
		spe = mono_dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe.writeSignificantSeedProteins(diff_out_folder + "mono_dcd_SPenr.txt");
		spc = mono_dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spc.writeSignificantSeedProteinCombinations(diff_out_folder + "mono_dcd_SPCenr.txt");
		
		System.out.println("Diff. seed comb.");
		mono_dsvd = new DiffSeedCombDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvd.writeSignSortedVariants(diff_out_folder + "mono_dsvcd_tfcsh.txt", true);
		mono_dsvd.writeSignSortedVariants(diff_out_folder + "mono_dsvcd_tfcs.txt", false);
		spe2 = mono_dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe2.writeSignificantSeedProteins(diff_out_folder + "mono_dsvcd_SPenr.txt");
	}
}
