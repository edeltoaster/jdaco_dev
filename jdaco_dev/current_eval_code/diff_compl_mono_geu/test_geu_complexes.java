package diff_compl_mono_geu;


import java.io.File;
import java.util.ArrayList;
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


public class test_geu_complexes {

	public static int iterations = 200;
	public static int[] equal_samples_per_group = {29, 25, 20, 15, 10, 5}; // 58 total samples
	public static int[] unequal_samples_per_group = {25, 20, 15, 10, 5};
	
	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> geu_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				geu_data.put(sample, qdr);
		}
		
		// equal sample size
		List<String> to_write = new LinkedList<String>();
		to_write.add("groupsize parametric samples_total samples_group1 samples_group2 iteration no_complexes no_compl_tfcs no_tfcs");
		List<String> samples = new ArrayList<>(geu_data.keySet());
		
		Set<String> permutations_used = new HashSet<>();
		String permutation;
		boolean permutate = true;
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (int samples_per_group:equal_samples_per_group) {
			for (int i = 0; i < iterations; i++) {
				
				if (permutate) {
					group1.clear();
					group2.clear();
					Collections.shuffle(samples);
					
					int j = 0;
					while (j < samples_per_group) {
						String sample = samples.get(j);
						group1.put(sample, geu_data.get(sample));
						j++;
					}
					
					while (j < samples_per_group*2) {
						String sample = samples.get(j);
						group2.put(sample, geu_data.get(sample));
						j++;
					}
					
					permutation = String.join("/", group1.keySet()) + "," + String.join("/", group2.keySet());
					
					if (!permutations_used.contains(permutation)) {
						permutate = false;
						permutations_used.add(permutation);
					}
				}
				
				DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_complexes = dcd.getSignificanceSortedComplexes().size();
				int no_compl_tfcs = dcd.getSignificantSeedCombVariants().size();
				DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_tfcs = dsvd.getSignificanceSortedVariants().size();
				
				int total = group1.size() + group2.size();
				to_write.add("equal nonparametric " + total + " " + group1.size() + " " + group2.size() + " " + i + " " + no_complexes + " " + no_compl_tfcs + " " + no_tfcs);
				permutate = true;
			}
			
			Utilities.writeEntries(to_write, "geu_rand_out.txt.gz");
		}
		
		// unequal sample size
		permutations_used.clear();
		for (int samples_first_group:unequal_samples_per_group) {
			for (int i = 0; i < iterations; i++) {
				
				if (permutate) {
					group1.clear();
					group2.clear();
					Collections.shuffle(samples);
					
					int j = 0;
					while (j < samples_first_group) {
						String sample = samples.get(j);
						group1.put(sample, geu_data.get(sample));
						j++;
					}
					
					while (j < geu_data.size()) {
						String sample = samples.get(j);
						group2.put(sample, geu_data.get(sample));
						j++;
					}
					
					permutation = String.join("/", group1.keySet()) + "," + String.join("/", group2.keySet());
					
					if (!permutations_used.contains(permutation)) {
						permutate = false;
						permutations_used.add(permutation);
					}
				}
				
				DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_complexes = dcd.getSignificanceSortedComplexes().size();
				int no_compl_tfcs = dcd.getSignificantSeedCombVariants().size();
				DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_tfcs = dsvd.getSignificanceSortedVariants().size();
				
				int total = group1.size() + group2.size();
				to_write.add("unequal nonparametric " + total + " " + group1.size() + " " + group2.size() + " " + i + " " + no_complexes + " " + no_compl_tfcs + " " + no_tfcs);
				permutate = true;
			}
			
			Utilities.writeEntries(to_write, "geu_rand_out.txt.gz");
		}
		
		/**
		 * parametric testing
		 */
		
		// equal sample size
		permutations_used.clear();
		for (int samples_per_group:equal_samples_per_group) {
			for (int i = 0; i < iterations; i++) {
				
				if (permutate) {
					group1.clear();
					group2.clear();
					Collections.shuffle(samples);
					
					int j = 0;
					while (j < samples_per_group) {
						String sample = samples.get(j);
						group1.put(sample, geu_data.get(sample));
						j++;
					}
					
					while (j < samples_per_group*2) {
						String sample = samples.get(j);
						group2.put(sample, geu_data.get(sample));
						j++;
					}
					
					permutation = String.join("/", group1.keySet()) + "," + String.join("/", group2.keySet());
					
					if (!permutations_used.contains(permutation)) {
						permutate = false;
						permutations_used.add(permutation);
					}
				}
				
				DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, true, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_complexes = dcd.getSignificanceSortedComplexes().size();
				int no_compl_tfcs = dcd.getSignificantSeedCombVariants().size();
				DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, true, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_tfcs = dsvd.getSignificanceSortedVariants().size();
				
				int total = group1.size() + group2.size();
				to_write.add("equal parametric " + total + " " + group1.size() + " " + group2.size() + " " + i + " " + no_complexes + " " + no_compl_tfcs + " " + no_tfcs);
				permutate = true;
			}
			
			Utilities.writeEntries(to_write, "geu_rand_out.txt.gz");
		}
		
		// unequal sample size
		permutations_used.clear();
		for (int samples_first_group:unequal_samples_per_group) {
			for (int i = 0; i < iterations; i++) {
				
				if (permutate) {
					group1.clear();
					group2.clear();
					Collections.shuffle(samples);
					
					int j = 0;
					while (j < samples_first_group) {
						String sample = samples.get(j);
						group1.put(sample, geu_data.get(sample));
						j++;
					}
					
					while (j < geu_data.size()) {
						String sample = samples.get(j);
						group2.put(sample, geu_data.get(sample));
						j++;
					}
					
					permutation = String.join("/", group1.keySet()) + "," + String.join("/", group2.keySet());
					
					if (!permutations_used.contains(permutation)) {
						permutate = false;
						permutations_used.add(permutation);
					}
				}
				
				DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, true, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_complexes = dcd.getSignificanceSortedComplexes().size();
				int no_compl_tfcs = dcd.getSignificantSeedCombVariants().size();
				DiffSeedCombDetector dsvd = new DiffSeedCombDetector(group1, group2, definitions.qvalue, true, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
				int no_tfcs = dsvd.getSignificanceSortedVariants().size();
				
				int total = group1.size() + group2.size();
				to_write.add("unequal parametric " + total + " " + group1.size() + " " + group2.size() + " " + i + " " + no_complexes + " " + no_compl_tfcs + " " + no_tfcs);
				permutate = true;
			}
			
			Utilities.writeEntries(to_write, "geu_rand_out.txt.gz");
		}
	}
}
