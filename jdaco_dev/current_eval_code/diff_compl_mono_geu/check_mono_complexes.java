package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPCEnrichment;
import framework.DiffComplexDetector.SPEnrichment;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_mono_complexes {

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
		
		
		/**
		 * Precompute quantified complexes
		 */
		
		System.out.println("Precompute quantified complexes ...");
		new File(definitions.qr_output_folder).mkdir();
		cm_data.keySet().parallelStream().forEach(s -> cm_data.get(s).writeQuantifiedResult(definitions.qr_output_folder + s + ".txt.gz"));
		ncm_data.keySet().parallelStream().forEach(s -> ncm_data.get(s).writeQuantifiedResult(definitions.qr_output_folder + s + ".txt.gz"));
		
		
		/**
		 * Diff. complexes (unpaired)
		 */
		
		new File(definitions.diff_out_folder).mkdir();
		System.out.println("Monocyte comparison cm->ncm (unpaired): " + cm_data.size() + " vs " + ncm_data.size());
		
		System.out.println("Diff. compl.");
		DiffComplexDetector mono_dcd = new DiffComplexDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcd_complh.txt", true);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcd_tfcsh.txt", true);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcd_compl.txt", false);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcd_tfcs.txt", false);
		framework.DiffComplexDetector.SPEnrichment spe = mono_dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dcd_SPenr.txt");
		SPCEnrichment spc = mono_dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spc.writeSignificantSeedProteinCombinations(definitions.diff_out_folder + "mono_dcd_SPCenr.txt");
		
		System.out.println("Diff. seed comb.");
		DiffSeedCombDetector mono_dsvd = new DiffSeedCombDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcd_tfcsh.txt", true);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcd_tfcs.txt", false);
		framework.DiffSeedCombDetector.SPEnrichment spe2 = mono_dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe2.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dsvcd_SPenr.txt");
		
		System.out.println("Monocyte comparison cm->ncm (unpaired, parametric): " + cm_data.size() + " vs " + ncm_data.size());
		
		System.out.println("Diff. compl.");
		mono_dcd = new DiffComplexDetector(cm_data, ncm_data, definitions.qvalue, true, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_pdcd_complh.txt", true);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdcd_tfcsh.txt", true);
		mono_dcd.writeSignSortedComplexes(definitions.diff_out_folder + "mono_pdcd_compl.txt", false);
		mono_dcd.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdcd_tfcs.txt", false);
		spe = mono_dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_pdcd_SPenr.txt");
		spc = mono_dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spc.writeSignificantSeedProteinCombinations(definitions.diff_out_folder + "mono_pdcd_SPCenr.txt");
		
		System.out.println("Diff. seed comb.");
		mono_dsvd = new DiffSeedCombDetector(cm_data, ncm_data, definitions.qvalue, true, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdsvcd_tfcsh.txt", true);
		mono_dsvd.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdsvcd_tfcs.txt", false);
		spe2 = mono_dsvd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spe2.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_pdsvcd_SPenr.txt");
		
		
		/**
		 * Diff. complexes (paired)
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
		
		System.out.println("Diff. compl.");
		DiffComplexDetector mono_dcdp = new DiffComplexDetector(pcm_data, pncm_data, definitions.qvalue, definitions.parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcdp.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcdp_complh.txt", true);
		mono_dcdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcdp_tfcsh.txt", true);
		mono_dcdp.writeSignSortedComplexes(definitions.diff_out_folder + "mono_dcdp_compl.txt", false);
		mono_dcdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_dcdp_tfcs.txt", false);
		SPEnrichment spep = mono_dcdp.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spep.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dcdp_SPenr.txt");
		SPCEnrichment spcp = mono_dcdp.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spcp.writeSignificantSeedProteinCombinations(definitions.diff_out_folder + "mono_dcdp_SPCenr.txt");
		
		System.out.println("Diff. seed comb.");
		DiffSeedCombDetector mono_dsvdp = new DiffSeedCombDetector(pcm_data, pncm_data, definitions.qvalue, definitions.parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcdp_tfcsh.txt", true);
		mono_dsvdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_dsvcdp_tfcs.txt", false);
		framework.DiffSeedCombDetector.SPEnrichment spep2 = mono_dsvdp.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spep2.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_dsvcdp_SPenr.txt");
		
		System.out.println("Monocyte comparison cm->ncm (paired, parametric): " + pcm_data.size() + " vs " + pncm_data.size());
		
		System.out.println("Diff. compl.");
		mono_dcdp = new DiffComplexDetector(pcm_data, pncm_data, definitions.qvalue, true, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dcdp.writeSignSortedComplexes(definitions.diff_out_folder + "mono_pdcdp_complh.txt", true);
		mono_dcdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdcdp_tfcsh.txt", true);
		mono_dcdp.writeSignSortedComplexes(definitions.diff_out_folder + "mono_pdcdp_compl.txt", false);
		mono_dcdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdcdp_tfcs.txt", false);
		spep = mono_dcdp.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spep.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_pdcdp_SPenr.txt");
		spcp = mono_dcdp.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold);
		spcp.writeSignificantSeedProteinCombinations(definitions.diff_out_folder + "mono_pdcdp_SPCenr.txt");
		
		System.out.println("Diff. seed comb.");
		mono_dsvdp = new DiffSeedCombDetector(pcm_data, pncm_data, definitions.qvalue, true, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		mono_dsvdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdsvcdp_tfcsh.txt", true);
		mono_dsvdp.writeSignSortedVariants(definitions.diff_out_folder + "mono_pdsvcdp_tfcs.txt", false);
		spep2 = mono_dsvdp.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		spep2.writeSignificantSeedProteins(definitions.diff_out_folder + "mono_pdsvcdp_SPenr.txt");
	}
}
