package diff_compl_mono_geu;


import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DiffComplexDetector;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class test_mono_complexes {
	
	public static int iterations = 5000;
	public static int[] samples_per_group = {15, 10, 5}; // 17 samples in both groups
	
	public static int paired_iterations = 500;
	public static int[] paired_samples_per_group = {13, 10, 5}; // 16 paired samples in both groups
	
	// current reference to compare with
	public static Set<HashSet<String>> ref_dcd_complexes;
	public static Set<HashSet<String>> ref_dcd_tfcs;
	public static Set<HashSet<String>> ref_dsvd_tfcs;

	
	/**
	 * Computes #detected events, Precision, Recall and F-Score for 
	 * dcd_complexes, dcd_tfcs and dsvd_tfcs given results compared to currently stored global reference:
	 * dcd_complexes_events dcd_complexes_prec dcd_complexes_rec dcd_complexes_F dcd_tfcs_events dcd_tfcs_prec dcd_tfcs_rec dcd_tfcs_F dsvd_tfcs_events dsvd_tfcs_prec dsvd_tfcs_rec dsvd_tfcs_F.
	 * If NaN values are found (0 results) -> return 0.0.
	 * @return
	 */
	public static String[] computeMetrics(Collection<HashSet<String>> test_dcd_complexes, Collection<HashSet<String>> test_dcd_tfcs, Collection<HashSet<String>> test_dsvd_tfcs) {
		
		// determine matching data
		Set<HashSet<String>> matching_dcd_complexes = new HashSet<>(test_dcd_complexes);
		matching_dcd_complexes.retainAll(ref_dcd_complexes);
		
		Set<HashSet<String>> matching_dcd_tfcs = new HashSet<>(test_dcd_tfcs);
		matching_dcd_tfcs.retainAll(ref_dcd_tfcs);
		
		Set<HashSet<String>> matching_dsvd_tfcs = new HashSet<>(test_dsvd_tfcs);
		matching_dsvd_tfcs.retainAll(ref_dsvd_tfcs);
		
		// compute all metrics
		int dcd_complexes_events = test_dcd_complexes.size();
		double dcd_complexes_prec = matching_dcd_complexes.size() / (double) dcd_complexes_events;
		double dcd_complexes_rec = matching_dcd_complexes.size() / (double) ref_dcd_complexes.size();
		double dcd_complexes_F = 2 * dcd_complexes_prec * dcd_complexes_rec / (dcd_complexes_prec + dcd_complexes_rec);
		
		int dcd_tfcs_events = test_dcd_tfcs.size();
		double dcd_tfcs_prec = matching_dcd_tfcs.size() / (double) dcd_tfcs_events;
		double dcd_tfcs_rec = matching_dcd_tfcs.size() / (double) ref_dcd_tfcs.size();
		double dcd_tfcs_F = 2 * dcd_tfcs_prec * dcd_tfcs_rec / (dcd_tfcs_prec + dcd_tfcs_rec);
		
		int dsvd_tfcs_events = test_dsvd_tfcs.size();
		double dsvd_tfcs_prec = matching_dsvd_tfcs.size() / (double) dsvd_tfcs_events;
		double dsvd_tfcs_rec = matching_dsvd_tfcs.size() / (double) ref_dsvd_tfcs.size();
		double dsvd_tfcs_F = 2 * dsvd_tfcs_prec * dsvd_tfcs_rec / (dsvd_tfcs_prec + dsvd_tfcs_rec);
		
		String[] temp = new String[]{Integer.toString(dcd_complexes_events), Double.toString(dcd_complexes_prec), Double.toString(dcd_complexes_rec), Double.toString(dcd_complexes_F), 
				Integer.toString(dcd_tfcs_events), Double.toString(dcd_tfcs_prec), Double.toString(dcd_tfcs_rec), Double.toString(dcd_tfcs_F), 
				Integer.toString(dsvd_tfcs_events), Double.toString(dsvd_tfcs_prec), Double.toString(dsvd_tfcs_rec), Double.toString(dsvd_tfcs_F)};
		
		// prevent problems
		for (int i = 0; i < temp.length; i++) {
			String s = temp[i];
			if (s.equals("NaN"))
				temp[i] = "0.0";
		}
		
		return temp;
	}
	
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
		cm_data.keySet().parallelStream().forEach(s -> cm_data.get(s).getAbundanceOfComplexes());
		ncm_data.keySet().parallelStream().forEach(s -> ncm_data.get(s).getAbundanceOfComplexes());
		
		List<String> cm_samples = new ArrayList<>(cm_data.keySet());
		List<String> ncm_samples = new ArrayList<>(ncm_data.keySet());
		
		List<String> to_write = new LinkedList<String>();
		
		to_write.add("paired parametric samples_total smaller_group samples_group1 samples_group2 iteration "
				+ "dcd_complexes_events dcd_complexes_prec dcd_complexes_rec dcd_complexes_F dcd_tfcs_events dcd_tfcs_prec dcd_tfcs_rec dcd_tfcs_F dsvd_tfcs_events dsvd_tfcs_prec dsvd_tfcs_rec dsvd_tfcs_F");
		
		boolean[] parametric = {false, true};
		
		System.out.println("Starting unpaired analysis ...");
		for (boolean assume_parametric:parametric) {
			System.out.println("Assume paramtric: " + assume_parametric);
			// compute reference
			Map<String, QuantDACOResultSet> group1 = new HashMap<>();
			Map<String, QuantDACOResultSet> group2 = new HashMap<>();
			
			cm_data.keySet().stream().forEach(s -> group1.put(s, cm_data.get(s)));
			ncm_data.keySet().stream().forEach(s -> group2.put(s, ncm_data.get(s)));
			
			DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, assume_parametric, false, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
			DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, assume_parametric, false, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
			
			ref_dcd_complexes = new HashSet<>(dcd.getSignificanceSortedComplexes());
			ref_dcd_tfcs = dcd.getSignificantSeedCombVariants();
			ref_dsvd_tfcs = new HashSet<>(dsvd.getSignificanceSortedVariants());
			
			String[] metrics = computeMetrics(ref_dcd_complexes, ref_dcd_tfcs, ref_dsvd_tfcs);
			
			to_write.add(false + " " + assume_parametric + " " + (group1.size() + group2.size()) + " 17 " + group1.size() + " " + group2.size() + " 0 " + String.join(" ", metrics));
			
			for (int no_samples_cm:samples_per_group) {
				for (int no_samples_ncm:samples_per_group) {
					
					System.out.println("Doing " + no_samples_cm + " vs " + no_samples_ncm);
					int smaller_group = Math.min(no_samples_cm, no_samples_ncm);
					Set<String> permutations_used = new HashSet<>();
					String permutation;
					boolean permutate = true;

					for (int i = 0; i < iterations; i++) {
						
						while (permutate) {
							group1.clear();
							group2.clear();
							Collections.shuffle(cm_samples);
							Collections.shuffle(ncm_samples);
							
							for (int j = 0; j < no_samples_cm; j++) {
								String sample = cm_samples.get(j);
								group1.put(sample, cm_data.get(sample));
							}
							
							for (int j = 0; j < no_samples_ncm; j++) {
								String sample = ncm_samples.get(j);
								group2.put(sample, ncm_data.get(sample));
							}
							
							permutation = String.join("/", group1.keySet()) + "," + String.join("/", group2.keySet());
							
							if (!permutations_used.contains(permutation)) {
								permutate = false;
								permutations_used.add(permutation);
							}
						}
						
						dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, assume_parametric, false, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
						dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, assume_parametric, false, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
						metrics = computeMetrics(dcd.getSignificanceSortedComplexes(), dcd.getSignificantSeedCombVariants(), dsvd.getSignificanceSortedVariants());
						
						to_write.add(false + " " + assume_parametric + " " + (group1.size() + group2.size()) + " " + smaller_group + " " + group1.size() + " " + group2.size() + " " + i + " " + String.join(" ", metrics));
						
						permutate = true;
					}
					
					Utilities.writeEntries(to_write, "mono_subsets_test.txt.gz");
					
				}
			}
			
		}
		
		
		/**
		 * same check for paired
		 */
		
		Map<String, QuantDACOResultSet> pcm_data = new HashMap<>();
		Map<String, QuantDACOResultSet> pncm_data = new HashMap<>();
		
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
		
		List<String> paired_samples = new ArrayList<>(pcm_data.keySet());
		
		System.out.println("Starting paired analysis ...");
		for (boolean assume_parametric:parametric) {
			System.out.println("Assume paramtric: " + assume_parametric);
			// compute reference
			Map<String, QuantDACOResultSet> group1 = new HashMap<>();
			Map<String, QuantDACOResultSet> group2 = new HashMap<>();
			
			pcm_data.keySet().stream().forEach(s -> group1.put(s, pcm_data.get(s)));
			pncm_data.keySet().stream().forEach(s -> group2.put(s, pncm_data.get(s)));
			
			DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, assume_parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
			DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, assume_parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
			
			ref_dcd_complexes = new HashSet<>(dcd.getSignificanceSortedComplexes());
			ref_dcd_tfcs = dcd.getSignificantSeedCombVariants();
			ref_dsvd_tfcs = new HashSet<>(dsvd.getSignificanceSortedVariants());
			
			String[] metrics = computeMetrics(ref_dcd_complexes, ref_dcd_tfcs, ref_dsvd_tfcs);
			
			to_write.add(true + " " + assume_parametric + " " + (group1.size() + group2.size()) + " 16 " + group1.size() + " " + group2.size() + " 0 " + String.join(" ", metrics));
			
			for (int no_samples:paired_samples_per_group) {
				System.out.println("Doing " + no_samples);
				Set<String> permutations_used = new HashSet<>();
				String permutation;
				boolean permutate = true;
				int smaller_group = no_samples; // for completeness, not really needed here
				for (int i = 0; i < paired_iterations; i++) {
					
					while (permutate) {
						group1.clear();
						group2.clear();
						Collections.shuffle(paired_samples);
						
						for (int j = 0; j < no_samples; j++) {
							String sample = paired_samples.get(j);
							group1.put(sample, pcm_data.get(sample));
							group2.put(sample, pncm_data.get(sample));
						}
						
						permutation = String.join("/", group1.keySet());
						
						if (!permutations_used.contains(permutation)) {
							permutate = false;
							permutations_used.add(permutation);
						}
					}
					
					dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, assume_parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
					dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, assume_parametric, true, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
					metrics = computeMetrics(dcd.getSignificanceSortedComplexes(), dcd.getSignificantSeedCombVariants(), dsvd.getSignificanceSortedVariants());
					
					to_write.add(true + " " + assume_parametric + " " + (group1.size() + group2.size()) + " " + smaller_group + " " + group1.size() + " " + group2.size() + " " + i + " " + String.join(" ", metrics));
					
					permutate = true;
				}
				
				Utilities.writeEntries(to_write, "mono_subsets_test.txt.gz");
			}
			
		}
	}
}
