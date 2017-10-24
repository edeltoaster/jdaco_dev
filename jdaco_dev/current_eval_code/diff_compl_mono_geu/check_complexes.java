package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPCEnrichment;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes {

	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> cm_data = new HashMap<>();
		Map<String, QuantDACOResultSet> ncm_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				continue;
			else if (sample.startsWith("non"))
				ncm_data.put(sample, qdr);
			else
				cm_data.put(sample, qdr);
		}
		
		new File(definitions.diff_out_folder).mkdir();
		
		System.out.println("Monocyte comparison cm->ncm (unpaired): " + cm_data.size() + " vs " + ncm_data.size());
		
		DiffComplexDetector mono_dcd = new DiffComplexDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcd_complh.txt", true);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcd_tfcsh.txt", true);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcd_compl.txt", false);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcd_tfcs.txt", false);
		framework.DiffComplexDetector.SPEnrichment spe = mono_dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dcd_SPenr.txt");
		SPCEnrichment spc = mono_dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spc.writeSignificantSeedProteinCombinations(definitions.diff_out_folder + "mono_dcd_SPCenr.txt");
		
		DiffSeedCombDetector mono_dsvd = new DiffSeedCombDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcd_tfcsh.txt", true);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcd_tfcs.txt", false);
		framework.DiffSeedCombDetector.SPEnrichment spe2 = mono_dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe2.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dsvcd_SPenr.txt");
		
		
		/**
		 * paired check
		 */
		
		Map<String, QuantDACOResultSet> pcm_data = new HashMap<>();
		Map<String, QuantDACOResultSet> pncm_data = new HashMap<>();
		
		System.out.println();
		
		// filter for paired analysis
		Set<String> not_in_both = new HashSet<>();
		
		// after this loop, not_in_both stores all numbered suffices of samples in cm_data
		for (String sample:cm_data.keySet()) {
			String suffix = sample.substring(sample.length()-2);
			not_in_both.add(suffix);
			pcm_data.put(suffix, cm_data.get(sample));
		}
		
		// after this loop, not_in_both stores all numbered suffices that are either in cm or in cnm
		for (String sample:ncm_data.keySet()) {
			String suffix = sample.substring(sample.length()-2);
			
			if (not_in_both.contains(suffix))
				not_in_both.remove(suffix);
			else
				not_in_both.add(suffix);
			pncm_data.put(suffix, ncm_data.get(sample));
		}
		
		System.out.println("Leaving out unmatched samples " + not_in_both);
		
		pcm_data.keySet().removeIf(s -> not_in_both.stream().anyMatch(b -> s.equals(b)));
		pncm_data.keySet().removeIf(s -> not_in_both.stream().anyMatch(b -> s.equals(b)));
		
		System.out.println("Monocyte comparison cm->ncm (paired): " + pcm_data.size() + " vs " + pncm_data.size());
		
		mono_dcd = new DiffComplexDetector(pcm_data, pncm_data, definitions.qvalue, definitions.parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcdp_complh.txt", true);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcdp_tfcsh.txt", true);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcdp_compl.txt", false);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcdp_tfcs.txt", false);
		spe = mono_dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dcdp_SPenr.txt");
		spc = mono_dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spc.writeSignificantSeedProteinCombinations(definitions.diff_out_folder + "mono_dcdp_SPCenr.txt");
		
		mono_dsvd = new DiffSeedCombDetector(pcm_data, pncm_data, definitions.qvalue, definitions.parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcdp_tfcsh.txt", true);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcdp_tfcs.txt", false);
		spe2 = mono_dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe2.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dsvcdp_SPenr.txt");
	}
}
