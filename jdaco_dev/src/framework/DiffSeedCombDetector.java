package framework;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

/**
 * Class for differential seed combination variants in DACO complexes
 * @author Thorsten Will
 */
public class DiffSeedCombDetector {
	
	// given data
	private final Map<String, QuantDACOResultSet> group1;
	private final Map<String, QuantDACOResultSet> group2;
	private final double FDR;
	private final boolean parametric;
	private final boolean paired;
	private final boolean incorporate_supersets;
	private int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	
	// processed / determined data
	private final Map<HashSet<String>, LinkedList<HashSet<String>>> seed_combination_variants;
	private final Map<HashSet<String>, LinkedList<Double>> group1_abundances;
	private final Map<HashSet<String>, LinkedList<Double>> group2_abundances;
	private Map<HashSet<String>, Double> group1_means;
	private Map<HashSet<String>, Double> group2_means;
	private Map<HashSet<String>, Double> fold_changes;
	
	// differential analysis results
	private Map<HashSet<String>, Double> variants_raw_pvalues;
	private final Map<HashSet<String>, Double> significant_variants_qvalues;
	private final List<HashSet<String>> significance_sorted_variants;
	private final Map<HashSet<String>, String> significance_variants_directions;
	
	// helper objects
	private final ForkJoinPool pool;
	private Map<String, String> up_to_gene_map; // computed on demand, use getters
	
	/**
	 * Constructor
	 * @param group1
	 * @param group2
	 * @param FDR
	 * @param parametric
	 * @param paired
	 * @param incorporate_supersets
	 * @param min_variant_fraction
	 * @param no_threads
	 */
	public DiffSeedCombDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR, boolean parametric, boolean paired, boolean incorporate_supersets, double min_variant_fraction, int no_threads, boolean filter_allosome) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.parametric = parametric;
		this.paired = paired;
		this.incorporate_supersets = incorporate_supersets;
		this.no_threads = no_threads;
		
		pool = new ForkJoinPool(this.no_threads);
		
		// ensure every sample has a matching counterpart in the other group
		if (this.paired)
			if (!group1.keySet().equals(group2.keySet())) {
				System.err.println("Non-equal sample groups supplied as matched data.");
				System.exit(1);
			}
		
		// determine potentially relevant seed variant combinations ...
		Map<HashSet<String>, Integer> count_map_g1 = new HashMap<>();
		group1.values().stream().forEach(sample -> sample.getSeedToComplexMap().keySet().stream().forEach(var -> count_map_g1.put(var, count_map_g1.getOrDefault(var, 0) + 1)));
		Map<HashSet<String>, Integer> count_map_g2 = new HashMap<>();
		group2.values().stream().forEach(sample -> sample.getSeedToComplexMap().keySet().stream().forEach(var -> count_map_g2.put(var, count_map_g2.getOrDefault(var, 0) + 1)));
		
		Set<HashSet<String>> unfiltered_seed_comb_var = new HashSet<>(count_map_g1.keySet());
		unfiltered_seed_comb_var.addAll(count_map_g2.keySet());
		
		// filter complexes with allosome proteins
		if (filter_allosome) {
			String db = DataQuery.getEnsemblOrganismDatabaseFromProteins(count_map_g1.keySet().iterator().next());
			Set<String> allosome = DataQuery.getAllosomeProteins(db);
			
			Set<HashSet<String>> to_remove = new HashSet<>();
			for (HashSet<String> tfc:unfiltered_seed_comb_var) {
				Set<String> overlap = new HashSet<>(tfc);
				overlap.retainAll(allosome);
				if (overlap.size() > 0)
					to_remove.add(tfc);
 			}
			unfiltered_seed_comb_var.removeAll(to_remove);
		}

		
		// check if above consideration threshold and determine subsets (if necessary)
		this.seed_combination_variants = new HashMap<>();
		double min_count_g1 = min_variant_fraction * group1.size();
		double min_count_g2 = min_variant_fraction * group2.size();
		for (HashSet<String> tfc:unfiltered_seed_comb_var) {
			
			this.seed_combination_variants.put(tfc, new LinkedList<HashSet<String>>());
			this.seed_combination_variants.get(tfc).add(tfc);
			int count_g1 = count_map_g1.getOrDefault(tfc, 0);
			int count_g2 = count_map_g2.getOrDefault(tfc, 0);
			
			if (this.incorporate_supersets)
				for (HashSet<String> current_tfc:unfiltered_seed_comb_var) 
					if (current_tfc.containsAll(tfc) && !current_tfc.equals(tfc)) {
						this.seed_combination_variants.get(tfc).add(current_tfc);
						count_g1 += count_map_g1.getOrDefault(current_tfc, 0);
						count_g2 += count_map_g2.getOrDefault(current_tfc, 0);
					}
			
			// filter
			if (count_g1 < min_count_g1 && count_g2 < min_count_g2)
				this.seed_combination_variants.remove(tfc);
		}
		
		
		
		// determine abundance values
		this.group1_abundances = this.determineAbundanceOfSeedVariantsComplexes(group1);
		this.group2_abundances = this.determineAbundanceOfSeedVariantsComplexes(group2);
		
		// determine medians
		ForkJoinTask<Map<HashSet<String>, Double>> group1_task = pool.submit(() -> this.group1_abundances.entrySet().parallelStream().collect(Collectors.toMap(e -> e.getKey(), e -> Utilities.getMean(e.getValue()))));
		ForkJoinTask<Map<HashSet<String>, Double>> group2_task = pool.submit(() -> this.group2_abundances.entrySet().parallelStream().collect(Collectors.toMap(e -> e.getKey(), e -> Utilities.getMean(e.getValue()))));
		try {
			this.group1_means = group1_task.get();
			this.group2_means = group2_task.get();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// determine differential abundance, apply multiple hypothesis correction and filter to significant seed variants
		if (!this.paired) {
			
			if (this.parametric) // parametric Welch test
				this.significant_variants_qvalues = this.determineUnpairedPValuesParametric();
			else // non-parametric MWU test
				this.significant_variants_qvalues = this.determineUnpairedPValuesNonParametric();
			
		} else {
			
			if (this.parametric) // paired t-test
				this.significant_variants_qvalues = this.determinePairedPValuesParametric();
			else // Wilcoxon signed-rank test
				this.significant_variants_qvalues = this.determinePairedPValuesNonParametric();
			
		}
		// precompute fold-changes
		this.fold_changes = new HashMap<HashSet<String>, Double>();
		this.variants_raw_pvalues.keySet().forEach(c -> this.fold_changes.put(c, Utilities.calcFoldChange(this.group1_means.get(c), this.group2_means.get(c))));
		
		// sort from most to least significant
		this.significance_sorted_variants = new ArrayList<>(this.significant_variants_qvalues.keySet());
		this.significance_sorted_variants.sort( (v1, v2) -> diffCompareTo(v1, v2));
		
		this.significance_variants_directions = this.determineDirections();
		
		pool.shutdownNow();
	}
	
	/**
	 * Compatibility constructor
	 * @param group1
	 * @param group2
	 * @param FDR
	 * @param parametric
	 * @param paired
	 * @param incorporate_supersets
	 * @param min_variant_fraction
	 * @param no_threads
	 */
	public DiffSeedCombDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR, boolean parametric, boolean paired, boolean incorporate_supersets, double min_variant_fraction, int no_threads) {
		this(group1, group2, FDR, parametric, paired, incorporate_supersets, min_variant_fraction, no_threads, false);
	}
	
	/**
	 * Minimal constructor
	 * @param group1
	 * @param group2
	 * @param FDR
	 * @param parametric
	 * @param paired
	 */
	public DiffSeedCombDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR, boolean parametric, boolean paired) {
		this(group1, group2, FDR, parametric, paired, false, 0.0, Runtime.getRuntime().availableProcessors(), false);
	}
	
	/**
	 * Custom compareTo function that first sorts by the adjusted p-value/q-value and breaks ties using the amount of fold-change and absolute differences of the mean between groups
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int diffCompareTo(HashSet<String> v1, HashSet<String> v2) {
		int sign_compareTo = this.significant_variants_qvalues.get(v1).compareTo(this.significant_variants_qvalues.get(v2));
		
		if (sign_compareTo == 0) {
			// determining fold-chance for v1
			double fold_change1 = this.fold_changes.get(v1);
			// ... and its extend
			fold_change1 = Utilities.amountFoldChange(fold_change1);
			
			// determining fold-chance for v2
			double fold_change2 = this.fold_changes.get(v2);
			// ... and its extend
			fold_change2 = Utilities.amountFoldChange(fold_change2);
			
			int sign_compareTo2 = Double.compare(fold_change2, fold_change1);
			
			if (sign_compareTo2 == 0) {
				double v1_mean_diff = Math.abs(this.group1_means.get(v1) - this.group2_means.get(v1));
				double v2_mean_diff = Math.abs(this.group1_means.get(v2) - this.group2_means.get(v2));
				return Double.compare(v2_mean_diff, v1_mean_diff); // note that we want to have the bigger difference as the smaller entry
			} else
				return sign_compareTo2;
			
		} else
			return sign_compareTo;
	}
	
	/**
	 * Determines the abundance values of each seed variant found with the overall abundance of complexes containing the exact variant in the samples of the groups. 
	 * Ordering of samples is ensured to order paired data correctly.
	 * @param group
	 * @return map of seed variant and abundance in each sample
	 */
	private Map<HashSet<String>, LinkedList<Double>> determineAbundanceOfSeedVariantsComplexes(Map<String, QuantDACOResultSet> group) {
		// init empty data structure for abundance data
		Map<HashSet<String>, LinkedList<Double>> group_abundances = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants.keySet()) {
			group_abundances.put(variant, new LinkedList<Double>());
		}
		
		// abundance of a seed_comb_variant is the sum of all complexes it is contained in
		ForkJoinTask<Map<String, Map<HashSet<String>, Double>>> task = this.pool.submit(() -> group.entrySet().parallelStream().collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue().getAbundanceOfSeedVariantsComplexes())));
		Map<String, Map<HashSet<String>, Double>> precomputed_sample_abundances = null;
		try {
			precomputed_sample_abundances = task.get();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// ensures ordering of samples for paired data even when different datastructures are used
		List<String> samples = new ArrayList<>(group.keySet());
		samples.sort(String::compareTo);
		
		for (String sample:samples) {
			Map<HashSet<String>, Double> sample_abundances = precomputed_sample_abundances.get(sample);
			for (HashSet<String> variant:this.seed_combination_variants.keySet()) {
				double abundance_values = 0.0;
				for (HashSet<String> current_variant:this.seed_combination_variants.get(variant))
					abundance_values += sample_abundances.getOrDefault(current_variant, 0.0);
				
				group_abundances.get(variant).add(abundance_values);
			}
		}
		
		return group_abundances;
	}
	
	/**
	 * Applies parametric heteroscedastic two-sample t-test/Welch test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determineUnpairedPValuesParametric() {
		TTest tt = new TTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants.keySet()) {
			// two-sample two-tailed Welch test, t-test with differing variances (heteroscedastic)
			double pm = tt.tTest(Utilities.getDoubleArray(this.group1_abundances.get(variant)), Utilities.getDoubleArray(this.group2_abundances.get(variant)));
			test_results.put(variant, pm);
		}
		
		this.variants_raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Applies non-parametric MWU test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determineUnpairedPValuesNonParametric() {
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants.keySet()) {
			// MWU test
			double pm = mwu.mannWhitneyUTest(Utilities.getDoubleArray(this.group1_abundances.get(variant)), Utilities.getDoubleArray(this.group2_abundances.get(variant)));
			test_results.put(variant, pm);
		}
		
		this.variants_raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Applies paired t-test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determinePairedPValuesParametric() {
		TTest tt = new TTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants.keySet()) {
			// paired t-test
			double pm = tt.pairedTTest(Utilities.getDoubleArray(this.group1_abundances.get(variant)), Utilities.getDoubleArray(this.group2_abundances.get(variant)));
			test_results.put(variant, pm);
		}
		
		this.variants_raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Applies non-parametric Wilcoxon signed-rank test test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determinePairedPValuesNonParametric() {
		WilcoxonSignedRankTest wsrt = new WilcoxonSignedRankTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		
		// above 20 samples, test statistic W is approx. normally distributed
		boolean compute_exact_p = true;
		if (this.group1.size() > 20)
			compute_exact_p = false;
		
		for (HashSet<String> variant:this.seed_combination_variants.keySet()) {
			// Wilcoxon signed-rank test
			double pm = wsrt.wilcoxonSignedRankTest(Utilities.getDoubleArray(this.group1_abundances.get(variant)), Utilities.getDoubleArray(this.group2_abundances.get(variant)), compute_exact_p);
			test_results.put(variant, pm);
		}
		
		this.variants_raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Precomputes directions of change
	 * @return
	 */
	private Map<HashSet<String>, String> determineDirections() {
		
		Map<HashSet<String>, String> significance_variants_directions = new HashMap<>();
		for (HashSet<String> variant:this.significance_sorted_variants) {
			double mean_g1 = this.group1_means.get(variant);
			double mean_g2 = this.group2_means.get(variant);
			
			String sign = "-";
			if (mean_g2 > mean_g1)
				sign = "+";
			
			if (mean_g1 == mean_g2) {
				// get medians instead
				double median_g1 = Utilities.getMedian(this.group1_abundances.get(variant));
				double median_g2 = Utilities.getMedian(this.group2_abundances.get(variant));
				
				sign = "-";
				if (median_g2 > median_g1)
					sign = "+";
			}
			
			significance_variants_directions.put(variant, sign);
		}
		
		return significance_variants_directions;
	}
	
	/**
	 * Does a standard analysis for the special case of transcription factors as seed proteins.
	 * All output is written to the specified folder.
	 * @param output_folder
	 */
	public void diffTFComplAnalysis(String output_folder, GOAnnotator goa, String binding_data_path, double binding_data_threshold, int binding_d_min, int binding_d_max, boolean also_compute_SCC, Set<String> proteins_to_remove, Set<String> proteins_of_interest) {

		// some first output
		System.out.println(this.getNumberOfTests() + " TF combinations tested.");
		System.out.println(this.getSignificanceSortedVariants().size() + " diff. TF combinations overall.");
		
		// no need to look further if there is nothing to tell about
		if (this.getSignificanceSortedVariants().size() == 0) {
			System.out.println();
			return;
		}
		
		// gather some helping data
		Map<String, String> up_name_map = this.getUniprotToGeneMap();
		
		/*
		 * build diff. complexes overview
		 */
		
		Set<String> involved_tfs = new HashSet<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		Map<HashSet<String>, List<HashSet<String>>> tfc_to_complexes = new HashMap<>();
		Map<String, String> directions_map = new HashMap<>();
		for (HashSet<String> tfc:this.getSignificanceSortedVariants()) {
			
			String sign = this.getSignificantVariantsDirections().get(tfc);
			involved_tfs.addAll(tfc);
			directions_map.put(tfc.toString(), sign);
			
			String compl_string = tfc.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			double pval = this.getSignificantVariantsQValues().get(tfc);
			String out_string = sign + " " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(out_string);
			} else {
				res_pos_all.add(out_string);
			}
			
			// determine actual complexes
			List<HashSet<String>> complexes = new LinkedList<>();
			for (HashSet<String> actual_tfc:this.seed_combination_variants.get(tfc)) {
				for (QuantDACOResultSet qdr:group1.values())
					if (qdr.getSeedToComplexMap().containsKey(actual_tfc))
						complexes.addAll(qdr.getSeedToComplexMap().get(actual_tfc));
				for (QuantDACOResultSet qdr:group2.values())
					if (qdr.getSeedToComplexMap().containsKey(actual_tfc))
						complexes.addAll(qdr.getSeedToComplexMap().get(actual_tfc));
			}
			tfc_to_complexes.put(tfc, complexes);
		}
		
		System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. TF combinations.");
		
		/*
		 * write results 
		 */
		
		// first ensure that folder exists while allowing to put prefixes in like output_folder/bla_ ...
		File output_folder_obj = null;
		if (!output_folder.endsWith("/")) {
			String[] spl = output_folder.split("/");
			String actual_folder = String.join("/", Arrays.asList(spl).subList(0, spl.length-1));
			output_folder_obj = new File(actual_folder);
		} else 
			output_folder_obj = new File(output_folder);
			
			output_folder_obj.mkdir();
		
		// actually writing something
		Utilities.writeEntries(res_pos_all, output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, output_folder + "res_neg_all.txt");
		
		// do the same with some given proteins not considered
		if (proteins_to_remove != null) {
			res_pos_all.clear();
			res_neg_all.clear();
			for (HashSet<String> tfc:this.getSignificanceSortedVariants()) {
				
				// skip TF comb that involve a protein that should be removed/not be considered
				if (tfc.stream().anyMatch(p -> proteins_to_remove.contains(p)))
					continue;
				
				String sign = this.getSignificantVariantsDirections().get(tfc);
				
				String compl_string = tfc.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				double pval = this.getSignificantVariantsQValues().get(tfc);
				String out_string = sign + " " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
				
				// distinguish between increased/positive abundance and diminishing/negative abundance
				if (sign.equals("-")) {
					res_neg_all.add(out_string);
				} else {
					res_pos_all.add(out_string);
				}
			}
			
			System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. TF combinations (pruned).");
			
			Utilities.writeEntries(res_pos_all, output_folder + "res_pos_pruned.txt");
			Utilities.writeEntries(res_neg_all, output_folder + "res_neg_pruned.txt");
		}
		
		// do the same with only some given proteins considered
		Set<HashSet<String>> POI_sign_tfcs = new HashSet<>();
		if (proteins_of_interest != null) {
			res_pos_all.clear();
			res_neg_all.clear();
			for (HashSet<String> tfc:this.getSignificanceSortedVariants()) {
				
				// skip TF comb that do not involve a protein that should be considered
				if (!tfc.stream().anyMatch(p -> proteins_of_interest.contains(p)))
					continue;
				
				POI_sign_tfcs.add(tfc);
				String sign = this.getSignificantVariantsDirections().get(tfc);
				
				String compl_string = tfc.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				double pval = this.getSignificantVariantsQValues().get(tfc);
				String out_string = sign + " " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
				
				// distinguish between increased/positive abundance and diminishing/negative abundance
				if (sign.equals("-")) {
					res_neg_all.add(out_string);
				} else {
					res_pos_all.add(out_string);
				}
			}
			
			System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. TF combinations (POI).");
			
			Utilities.writeEntries(res_pos_all, output_folder + "res_pos_POI.txt");
			Utilities.writeEntries(res_neg_all, output_folder + "res_neg_POI.txt");
		}
		
		
		/*
		 * build regulatory network
		 */
		
		// build information for regulatory network
		//Map<String, String> occ_across_samples = new HashMap<>();
		Map<String, String> mean_abun_g1 = new HashMap<>();
		Map<String, String> mean_abun_g2 = new HashMap<>();
		Map<String, String> GO_details = new HashMap<>();
		Map<String, String> GO_overall = new HashMap<>();
		for (HashSet<String> tf_comb:tfc_to_complexes.keySet()) {
			String[] annotation_output = this.getSPCAnnotations(tf_comb, tfc_to_complexes.get(tf_comb), goa);
			String tfc = tf_comb.toString();
			//occ_across_samples.put(tfc, annotation_output[0]);
			mean_abun_g1.put(tfc, annotation_output[1]);
			mean_abun_g2.put(tfc, annotation_output[2]);
			GO_details.put(tfc, annotation_output[3]);
			GO_overall.put(tfc, annotation_output[4]);
		}
		
		// read binding data
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(binding_data_path, involved_tfs, binding_data_threshold, involved_tfs);
		
		// build regulatory network
		System.out.println("Building regulatory network ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(tfc_to_complexes.keySet(), involved_tfs, bdh, binding_d_min, binding_d_max, no_threads, 1);
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(output_folder + "regnet.txt");
		
		// write annotation data
		Map<String, Map<String,String>> annotational_data = new TreeMap<>();
		annotational_data.put("Mean_abundance_G1", mean_abun_g1);
		annotational_data.put("Mean_abundance_G2", mean_abun_g2);
		annotational_data.put("Direction", directions_map);
		annotational_data.put("GO_details", GO_details);
		annotational_data.put("GO_overall", GO_overall);
		regnet.writeNodeTable(output_folder + "nodetable.txt", annotational_data);
		
		// prune if interesting
		if (also_compute_SCC) {
			regnet.pruneToLargestSCCs();
			System.out.println("SCC: " + regnet.getSizesStr());
			regnet.writeRegulatoryNetwork(output_folder + "regnet_SCC.txt");
			regnet.writeNodeTable(output_folder + "nodetable_SCC.txt", annotational_data);
		}
		
		if (proteins_to_remove != null) {
			regnet.removeProteinSet(proteins_to_remove);
			System.out.println("pruned: " + regnet.getSizesStr());
			regnet.writeRegulatoryNetwork(output_folder + "regnet_pruned.txt");
			regnet.writeNodeTable(output_folder + "nodetable_pruned.txt", annotational_data);
		}
		
		if (proteins_of_interest != null && POI_sign_tfcs.size() > 0) {
			regnet = new RegulatoryNetwork(POI_sign_tfcs, involved_tfs, bdh, binding_d_min, binding_d_max, no_threads, 1);
			System.out.println("POI: " + regnet.getSizesStr());
			regnet.writeRegulatoryNetwork(output_folder + "regnet_POI.txt");
			regnet.writeNodeTable(output_folder + "nodetable_POI.txt", annotational_data);
		}
	}
	
	/**
	 * Returns parsable output in the space separated format (optionally Uniprot Accs are converted to gene identifiers):
	 * (sub)complex direction q-value fold-change mean-change member_seed_comb member_complexes
	 * @param include_header
	 * @param human_readable
	 */
	public List<String> getSignSortedVariants(boolean include_header, boolean human_readable) {
		List<String> to_write = new LinkedList<>();
		
		if (include_header)
			to_write.add("(sub)seed_comb direction q-value fold-change mean-change member_seed_comb member_complexes");
		
		if (human_readable)
			this.getUniprotToGeneMap();
		
		for (HashSet<String> seed_comb:this.significance_sorted_variants) {
			
			// translate seed combination if necessary
			Set<String> seed_comb_temp = null;
			if (human_readable)
				seed_comb_temp = seed_comb.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet());
			else 
				seed_comb_temp = seed_comb;
			String seed_comb_string = String.join("/", seed_comb_temp);
			
			// determine fold-change
			double fold_change = this.fold_changes.get(seed_comb);
			
			// median change calculation
			double mean_change = this.group2_means.get(seed_comb) - this.group1_means.get(seed_comb);
			
			// determine member seed combinations if subset was used
			List<HashSet<String>> member_seed_comb = this.seed_combination_variants.get(seed_comb);
			
			if (human_readable) {
				List<HashSet<String>> member_seed_comb_temp = new ArrayList<>(member_seed_comb.size());
				member_seed_comb.stream().forEach(l -> member_seed_comb_temp.add(new HashSet<String>(l.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()))));
				member_seed_comb = member_seed_comb_temp;
			}
			
			String member_seed_comb_string = String.join(",", member_seed_comb.stream().map(l -> String.join("/", l)).collect(Collectors.toSet()));
			
			
			// determine actual complexes and translate if necessary
			List<HashSet<String>> complexes = new LinkedList<>();
			for (HashSet<String> actual_tfc:this.seed_combination_variants.get(seed_comb)) {
				for (QuantDACOResultSet qdr:group1.values())
					if (qdr.getSeedToComplexMap().containsKey(actual_tfc))
						complexes.addAll(qdr.getSeedToComplexMap().get(actual_tfc));
				for (QuantDACOResultSet qdr:group2.values())
					if (qdr.getSeedToComplexMap().containsKey(actual_tfc))
						complexes.addAll(qdr.getSeedToComplexMap().get(actual_tfc));
			}
			
			if (human_readable) {
				List<HashSet<String>> member_compl_temp = new ArrayList<>(complexes.size());
				complexes.stream().forEach(l -> member_compl_temp.add(new HashSet<String>(l.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()))));
				complexes = member_compl_temp;
			}
			
			String member_complexes_string = String.join(",", complexes.stream().map(l -> String.join("/", l)).collect(Collectors.toSet()));
			
			// shortening numbers depending on human/machine-usage
			String qval_string = null;
			String foldc_string = null;
			String mean_change_string = null;
			if (human_readable) {
				qval_string = String.format(Locale.US, "%.3g", this.significant_variants_qvalues.get(seed_comb));
				foldc_string = String.format(Locale.US, "%.3g", fold_change);
				mean_change_string = String.format(Locale.US, "%.3g", mean_change);
			} else {
				qval_string = Double.toString(this.significant_variants_qvalues.get(seed_comb));
				foldc_string = Double.toString(fold_change);
				mean_change_string = Double.toString(mean_change);
			}
			
			// write in format: (sub)seed_comb direction q-value fold-change median-change member_seed_comb member_complexes
			List<String> line = Arrays.asList(seed_comb_string, this.significance_variants_directions.get(seed_comb), qval_string, foldc_string, mean_change_string, member_seed_comb_string, member_complexes_string);
			to_write.add(String.join(" ", line));
		}

		return to_write;
	}
	
	/**
	 * Writes parsable output in the space separated format (optionally Uniprot Accs are converted to gene identifiers):
	 * (sub)complex direction q-value fold-change mean-change member_seed_comb member_complexes
	 * @param out_file
	 * @param human_readable
	 */
	public void writeSignSortedVariants(String out_file, boolean human_readable) {
		
		Utilities.writeEntries(this.getSignSortedVariants(true, human_readable), out_file);
	}
	
	
	/*
	 * getters
	 */
	
	/**
	 * Returns data of group1
	 * @return
	 */
	public Map<String, QuantDACOResultSet> getGroup1() {
		return group1;
	}

	/**
	 * Returns data of group2
	 * @return
	 */
	public Map<String, QuantDACOResultSet> getGroup2() {
		return group2;
	}

	/**
	 * Returns FDR used
	 * @return
	 */
	public double getFDR() {
		return FDR;
	}

	/**
	 * Returns the if parametric tests or non-parametric tests were used
	 * @return
	 */
	public boolean getParametricTestsUsed() {
		return this.parametric;
	}
	
	/**
	 * Returns the if a paired tests were used
	 * @return
	 */
	public boolean getPairedTestsUsed() {
		return this.paired;
	}
	
	/**
	 * Returns if supersets have been incorporated in the significance calculations
	 * @return
	 */
	public boolean getIncorporateSupersets() {
		return this.incorporate_supersets;
	}
	
	/**
	 * Returns all seed combinations that have been found in all samples
	 * @return
	 */
	public Set<HashSet<String>> getSeedCombinationVariants() {
		return seed_combination_variants.keySet();
	}

	/**
	 * Returns abundance data of group1
	 * @return
	 */
	public Map<HashSet<String>, LinkedList<Double>> getGroup1Abundances() {
		return group1_abundances;
	}

	/**
	 * Returns abundance data of group2
	 * @return
	 */
	public Map<HashSet<String>, LinkedList<Double>> getGroup2Abundances() {
		return group2_abundances;
	}

	/**
	 * Returns mean abundance data of group1
	 * @return
	 */
	public Map<HashSet<String>, Double> getGroup1MeanAbundances() {
		return group1_means;
	}

	/**
	 * Returns mean abundance data of group2
	 * @return
	 */
	public Map<HashSet<String>, Double> getGroup2MeanAbundances() {
		return group2_means;
	}
	
	/**
	 * Returns map of significant seed combination variants and associated adjusted p-values / q-values
	 * @return
	 */
	public Map<HashSet<String>, Double> getSignificantVariantsQValues() {
		return significant_variants_qvalues;
	}

	/**
	 * Returns map of seed combination variants and associated raw p-values
	 * @return
	 */
	public Map<HashSet<String>, Double> getVariantsRawPValues() {
		return variants_raw_pvalues;
	}
	
	/**
	 * Returns the number of statistical tests made.
	 * @return
	 */
	public int getNumberOfTests() {
		return variants_raw_pvalues.size();
	}
	
	/**
	 * Returns a list of significantly deregulated seed combination variants in their order of significance
	 * @return
	 */
	public List<HashSet<String>> getSignificanceSortedVariants() {
		return significance_sorted_variants;
	}
	
	/**
	 * Returns a map of sign. seed variants to their direction of change as + / -
	 * @return
	 */
	public Map<HashSet<String>, String> getSignificantVariantsDirections() {
		return this.significance_variants_directions;
	}
	
	
	/*
	 * General utility
	 */
	
	/**
	 * Returns a map of Uniprot accession to gene name map matching the organism.
	 * @return
	 */
	public Map<String, String> getUniprotToGeneMap() {
		if (this.up_to_gene_map != null)
			return this.up_to_gene_map;
		
		this.up_to_gene_map = DataQuery.getUniprotToGeneNameMap(this.group1.values().iterator().next().getAbundantSeedProteins());
		
		return this.up_to_gene_map;
	}
	
	
	/*
	 * generic helper functions
	 */
	
	/**
	 * Given the seed protein combinations occurring across samples, direction, sample data and a GOA definition, count and sort their appearances as well as their inferred GO annotations.
	 * Returns [occurrences across samples, mean abundances group1 (rounded to 2 positions), mean abundances group2 (rounded to 2 positions), GO annotations details string, set of all GO annotations found]
	 * @param complexes
	 * @param goa
	 * @return
	 */
	private String[] getSPCAnnotations(HashSet<String> seed_protein_comb, List<HashSet<String>> complexes, GOAnnotator goa) {
		
		// precompute data
		Map<HashSet<String>, Integer> count_map = new HashMap<>();
		complexes.stream().forEach(c -> count_map.put(c, count_map.getOrDefault(c, 0) + 1));
		List<HashSet<String>> occ_sorted_complexes = new ArrayList<>(count_map.keySet());
		occ_sorted_complexes.sort( (c1, c2) -> count_map.get(c2).compareTo(count_map.get(c1)));
		
		Map<String, String> up_name = this.getUniprotToGeneMap();
		Map<HashSet<String>, String> names_map = new HashMap<>();
		occ_sorted_complexes.stream().forEach(c -> names_map.put(c, String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList()))));
		String occ_sorted_complexes_string = String.join(",", occ_sorted_complexes.stream().map(c -> names_map.get(c) + ":" + count_map.get(c)).collect(Collectors.toList()));
		
		// determine mean abundance values
		Map<HashSet<String>, Double> group1_mean_abundances = new HashMap<>();
		occ_sorted_complexes.stream().forEach(c -> group1_mean_abundances.put(c, Utilities.getMean(this.group1.values().stream().map(qdr -> qdr.getAbundanceOfComplexes().getOrDefault(c, 0.0)).collect(Collectors.toList()))));
		String abun_g1_string = String.join(",", occ_sorted_complexes.stream().map(c -> names_map.get(c) + ":" + String.format(Locale.US, "%.2g", group1_mean_abundances.get(c)) ).collect(Collectors.toList()));
		
		Map<HashSet<String>, Double> group2_mean_abundances = new HashMap<>();
		occ_sorted_complexes.stream().forEach(c -> group2_mean_abundances.put(c, Utilities.getMean(this.group2.values().stream().map(qdr -> qdr.getAbundanceOfComplexes().getOrDefault(c, 0.0)).collect(Collectors.toList()))));
		String abun_g2_string = String.join(",", occ_sorted_complexes.stream().map(c -> names_map.get(c) + ":" + String.format(Locale.US, "%.2g", group2_mean_abundances.get(c)) ).collect(Collectors.toList()));
		
		// annotate with GO annotations
		Map<Set<String>, String> GOA_map = new HashMap<>();
		occ_sorted_complexes.stream().forEach(c -> GOA_map.put(c, goa.rateProteins(c)));
		GOA_map.entrySet().removeIf(e -> e.getValue().equals("/"));
		
		String GOAnnotation_string_details = "/";
		String GOAnnotation_string_overall = "/";
		List<HashSet<String>> GO_complexes = new ArrayList<>(occ_sorted_complexes);
		GO_complexes.retainAll(GOA_map.keySet());
		if (!GO_complexes.isEmpty()) {
			GOAnnotation_string_details = String.join(",", GO_complexes.stream().map(c -> names_map.get(c) + ":" + GOA_map.get(c)).collect(Collectors.toList()));
			GOAnnotation_string_overall = String.join(",", GO_complexes.stream().map(c -> GOA_map.get(c)).collect(Collectors.toSet()));
		}
		
		// build output datastructure
		String[] output = new String[5];
		output[0] = occ_sorted_complexes_string;
		output[1] = abun_g1_string;
		output[2] = abun_g2_string;
		output[3] = GOAnnotation_string_details;
		output[4] = GOAnnotation_string_overall;
		
		return output;
	}
	
	/**
	 * Given the seed variant occurring across samples, direction, sample data and a GOA definition, count and sort their appearances as well as their inferred GO annotations.
	 * Returns [sorted actual complexes string, sorted GO annotations string, set of all GO annotations found, sorted mean abundances (rounded to 2 positions)]
	 * @param complexes
	 * @param goa
	 * @return
	 */
	public static String[] getSortedComplexesAnnotations(HashSet<String> variant, String sign, GOAnnotator goa, Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2) {
		Map<String, QuantDACOResultSet> group_to_check = group1;
		if (sign.equals("+"))
			group_to_check = group2;
		
		// determine actual complexes
		List<Set<String>> complexes = new LinkedList<>();
		for (QuantDACOResultSet qdr:group_to_check.values())
			if (qdr.getSeedToComplexMap().containsKey(variant))
				complexes.addAll(qdr.getSeedToComplexMap().get(variant));
		
		// get naming data
		Map<String, String> up_name = DataQuery.getUniprotToGeneNameMap(complexes.get(0));
		
		// count occurrences across group samples
		Map<Set<String>, Integer> count_map = new HashMap<>();
		complexes.stream().forEach(c -> count_map.put(c, count_map.getOrDefault(c, 0) + 1));
		List<Set<String>> occ_sorted_complexes = new ArrayList<>(count_map.keySet());
		occ_sorted_complexes.sort( (c1, c2) -> count_map.get(c2).compareTo(count_map.get(c1)));
		String occ_sorted_complexes_string = String.join(",", occ_sorted_complexes.stream().map(c -> String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList())) + ":" + count_map.get(c)).collect(Collectors.toList()));
		
		// annotate with GO annotations
		Map<Set<String>, String> GOA_map = new HashMap<>();
		occ_sorted_complexes.stream().forEach(c -> GOA_map.put(c, goa.rateProteins(c)));
		GOA_map.entrySet().removeIf(e -> e.getValue().equals("/"));
		
		String sorted_GOAnnotation_string = "/";
		String GOAnnotation_string = "/";
		occ_sorted_complexes.retainAll(GOA_map.keySet());
		if (!occ_sorted_complexes.isEmpty()) {
			sorted_GOAnnotation_string = String.join(",", occ_sorted_complexes.stream().map(c -> String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList())) + ":" + GOA_map.get(c)).collect(Collectors.toList()));
			GOAnnotation_string = String.join(",", occ_sorted_complexes.stream().map(c->GOA_map.get(c)).collect(Collectors.toSet()));
		}
		
		// determine abundance values
		Map<Set<String>, Double> mean_abundances = new HashMap<>();
		List<Set<String>> abun_sorted_complexes = new ArrayList<>(new HashSet<>(complexes));
		final Map<String, QuantDACOResultSet> group_final = group_to_check;
		abun_sorted_complexes.stream().forEach(c -> mean_abundances.put(c, Utilities.getMean(group_final.values().stream().map(qdr -> qdr.getAbundanceOfComplexes().getOrDefault((HashSet<String>) c, 0.0)).collect(Collectors.toList()))));
		abun_sorted_complexes.sort( (c1, c2) -> mean_abundances.get(c2).compareTo(mean_abundances.get(c1)));
		String abun_sorted_complexes_string = String.join(",", abun_sorted_complexes.stream().map(c -> String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList())) + ":" + String.format(Locale.US, "%.2g", mean_abundances.get(c)) ).collect(Collectors.toList()));
		
		// build output datastructure
		String[] output = new String[4];
		output[0] = occ_sorted_complexes_string;
		output[1] = sorted_GOAnnotation_string;
		output[2] = GOAnnotation_string;
		output[3] = abun_sorted_complexes_string;
		
		return output;
	}
	
	/*
	 * Seed protein enrichment calculations
	 */
	
	/**
	 * Compute seed protein enrichment in the fashion of GSEA, 
	 * see Subramanian et al. (2005)
	 * @param FDR
	 * @param iterations
	 * @param complex_participation_cutoff
	 * @return
	 */
	public SPEnrichment calculateSPEnrichment(double FDR, int iterations, int complex_participation_cutoff) {
		return new SPEnrichment(FDR, iterations, complex_participation_cutoff);
	}
	
	/**
	 * Compute seed protein enrichment in the fashion of GSEA, 
	 * see Subramanian et al. (2005)
	 */
	public final class SPEnrichment {
		
		private final double[] sorted_scores;
		private final Map<String, Double> sp_in_complex_count;
		private final double no_complexes;
		private final List<String> sign_sorted_sps;
		private final Map<String, Double> sign_sp_qvalue;
		private final Map<String, String> sign_sp_direction;
		
		public SPEnrichment(double FDR, int iterations, int complex_participation_cutoff) {
			
			// -log p-values -> scores
			Map<HashSet<String>, Double> scores = new HashMap<>();
			for (HashSet<String> variant:variants_raw_pvalues.keySet()) {
				double mean_diff = group2_means.get(variant) - group1_means.get(variant);
				double score = -Math.log(variants_raw_pvalues.get(variant));
				
				if (mean_diff < 0)
					score *= -1;
				
				scores.put(variant, score);
			}
			
			// construct sorted matching arrays of variant/score
			List<HashSet<String>> sorted_complex_list = new ArrayList<>(scores.keySet());
			sorted_complex_list.sort((v1, v2) -> diffScoresCompareTo(v2, v1, scores));
			sorted_scores = new double[sorted_complex_list.size()];
			for (int i = 0; i< sorted_scores.length; i++)
				sorted_scores[i] = scores.get(sorted_complex_list.get(i));
			
			// determine complex participation of each seed protein (analogous to N_H in org. paper)
			sp_in_complex_count = new HashMap<>();
			variants_raw_pvalues.keySet().stream().forEach(sps -> sps.stream().forEach(sp -> sp_in_complex_count.put(sp, sp_in_complex_count.getOrDefault(sp, 0.0) + 1)));
			// remove SPs that are not in complexes very often
			sp_in_complex_count.keySet().removeIf(sp -> sp_in_complex_count.get(sp).doubleValue() < complex_participation_cutoff);
			
			// precompute #complexes (analogous to N in org. paper)
			no_complexes = sorted_scores.length;
			
			// for each seed protein: calc. enrichment score
			Map<String, Double> sp_enrichment_scores = new HashMap<>();
			for (String sp:sp_in_complex_count.keySet())
				sp_enrichment_scores.put(sp, calculateEnrichmentScore(sp, sorted_complex_list));

			// get distr of NULL distributions by permutation
			List<List<HashSet<String>>> randomized_complex_lists = new ArrayList<>(iterations);
			for (int i = 0; i< iterations; i++) {
				List<HashSet<String>> complex_list = new ArrayList<>(sorted_complex_list);
				Collections.shuffle(complex_list);
				randomized_complex_lists.add(complex_list);
			}
			
			// compute ES and distributions for each seed protein of interest
			ForkJoinPool pool = new ForkJoinPool(no_threads);
			Map<String, Double> normalized_pos_sp_ES = new HashMap<>();
			Map<String, Double> normalized_neg_sp_ES = new HashMap<>();
			List<Double> normalized_pos_rnd_data = new LinkedList<>();
			List<Double> normalized_neg_rnd_data = new LinkedList<>();
			for (String sp:sp_in_complex_count.keySet()) {
				double sp_ES = sp_enrichment_scores.get(sp);
				ForkJoinTask<List<Double>> rnd_counts = pool.submit(() -> randomized_complex_lists.parallelStream().map(l -> calculateEnrichmentScore(sp, l)).collect(Collectors.toList()));
				List<Double> rnd_data = null;
				try {
					rnd_data = rnd_counts.get();
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(1);
				}
				
				// towards multiple hypothesis correction
				double pos_mean = Utilities.getMean(rnd_data.stream().filter(d -> d.doubleValue() >= 0).collect(Collectors.toList()));
				double neg_mean = Utilities.getMean(rnd_data.stream().filter(d -> d.doubleValue() <= 0).map(d -> Math.abs(d)).collect(Collectors.toList()));
				
				for (double d:rnd_data) {
					if (d > 0)
						normalized_pos_rnd_data.add(d / pos_mean);
					else if (d < 0)
						normalized_neg_rnd_data.add(d / neg_mean); // neg_mean is absolute -> will still be negative
					else {// 0 in both, just for correctness
						normalized_pos_rnd_data.add(d);
						normalized_neg_rnd_data.add(d);
					}
				}
				
				if (sp_ES >= 0)
					normalized_pos_sp_ES.put(sp, sp_ES / pos_mean);
				else
					normalized_neg_sp_ES.put(sp, sp_ES / neg_mean);
				
			}
			
			// corrected FDR q-values
			sign_sp_qvalue = new HashMap<>();
			sign_sp_direction = new HashMap<>();
			for (String sp:normalized_pos_sp_ES.keySet()) {
				double norm_ES = normalized_pos_sp_ES.get(sp);
				// calculate q-value
				double FDR_q = Math.max(normalized_pos_rnd_data.stream().filter(d -> d.doubleValue() >= norm_ES).count(), 1) / (double) normalized_pos_rnd_data.size();
				FDR_q /= normalized_pos_sp_ES.keySet().stream().filter(sp2 -> normalized_pos_sp_ES.get(sp2).doubleValue() >= norm_ES).count() / (double) normalized_pos_sp_ES.keySet().size(); // will be >1 since norm_ES is in there
				sign_sp_qvalue.put(sp, FDR_q);
				sign_sp_direction.put(sp, "+");
			}
			for (String sp:normalized_neg_sp_ES.keySet()) {
				double norm_ES = normalized_neg_sp_ES.get(sp);
				// calculate q-value
				double FDR_q = Math.max(normalized_neg_rnd_data.stream().filter(d -> d.doubleValue() <= norm_ES).count(), 1) / (double) normalized_neg_rnd_data.size();
				FDR_q /= normalized_neg_sp_ES.keySet().stream().filter(sp2 -> normalized_neg_sp_ES.get(sp2).doubleValue() <= norm_ES).count() / (double) normalized_neg_sp_ES.keySet().size(); // will be >1 since norm_ES is in there
				sign_sp_qvalue.put(sp, FDR_q);
				sign_sp_direction.put(sp, "-");
			}
			
			// only retain those below threshold
			sign_sp_qvalue.keySet().removeIf(tf -> sign_sp_qvalue.get(tf).doubleValue() >= FDR);
			sign_sorted_sps = new ArrayList<>(sign_sp_qvalue.keySet());
			sign_sorted_sps.sort((v1, v2) -> sign_sp_qvalue.get(v1).compareTo(sign_sp_qvalue.get(v2)));
		}
		
		/**
		 * Custom compareTo function that first sorts by the log-scores and breaks ties using the differences of the median between groups
		 * @param v1
		 * @param v2
		 * @param scores
		 * @return
		 */
		private int diffScoresCompareTo(HashSet<String> v1, HashSet<String> v2, Map<HashSet<String>, Double> scores) {
			int sign_compareTo = scores.get(v1).compareTo(scores.get(v2));
			
			if (sign_compareTo == 0) {
				// determining fold-chance for v1
				double fold_change1 = fold_changes.get(v1);
				// determining fold-chance for v2
				double fold_change2 = fold_changes.get(v2);
				int sign_compareTo2 = Double.compare(fold_change1, fold_change2);
				
				if (sign_compareTo2 == 0) {
					double v1_mean_diff = group2_means.get(v1) - group1_means.get(v1);
					double v2_mean_diff = group2_means.get(v2) - group1_means.get(v2);
					return Double.compare(v1_mean_diff, v2_mean_diff);
				} else
					return sign_compareTo2;
				
			} else
				return sign_compareTo;
		}
		
		/**
		 * Calculates enrichment score in the fashion of GSEA, 
		 * handles enrichment and depletion independently, returns [max_ES, min_ES]
		 * @param query_SP
		 * @param scores
		 * @param complex_list
		 * @return
		 */
		private double calculateEnrichmentScore(String query_SP, List<HashSet<String>> complex_list) {
			// determine sum of scores (N_R in org. paper)
			int i = 0;
			double N_R = 0.0;
			for (HashSet<String> complex:complex_list) {
				if (complex.contains(query_SP))
					N_R += Math.abs(sorted_scores[i]);
				++i;
			}
			
			double N_H = sp_in_complex_count.get(query_SP);
			double running_sum = 0.0;
			double max_dev = 0.0;
			i = 0;
			for (HashSet<String> complex:complex_list) {
				if (complex.contains(query_SP))
					running_sum += (Math.abs(sorted_scores[i]) / N_R);
				else
					running_sum -= 1.0 / (no_complexes - N_H);
				
				if (Math.abs(running_sum) > Math.abs(max_dev))
					max_dev = running_sum;

				
				++i;
			}
			
			return max_dev;
		}
		
		/**
		 * Writes result into two specified text files.
		 * @param pos_out_path
		 * @param neg_out_path
		 */
		public void writeSignificantSeedProteins(String pos_out_path, String neg_out_path) {
			
			List<String> pos_sp_enrich_out = new LinkedList<>();
			List<String> neg_sp_enrich_out = new LinkedList<>();
			Map<String, String> up_to_gene_map = getUniprotToGeneMap();
			for (String sp:this.getSignificanceSortedSeedProteins()) {
				String dir = this.getSignificantSeedProteinDirections().get(sp);
				String out = sp + " " + up_to_gene_map.getOrDefault(sp, sp) + " " + this.getSignificantSeedProteinQvalues().get(sp);
				
				if (dir.equals("+"))
					pos_sp_enrich_out.add(out);
				else
					neg_sp_enrich_out.add(out);
			}
			Utilities.writeEntries(pos_sp_enrich_out, pos_out_path);
			Utilities.writeEntries(neg_sp_enrich_out, neg_out_path);
		}
		
		/**
		 * Writes result into a specified text file.
		 * @param out_path
		 */
		public void writeSignificantSeedProteins(String out_path) {
			
			List<String> sp_enrich_out = new LinkedList<>();
			Map<String, String> up_to_gene_map = getUniprotToGeneMap();
			
			for (String sp:this.getSignificanceSortedSeedProteins()) {
				String dir = this.getSignificantSeedProteinDirections().get(sp);
				String out = dir + " " + sp + " " + up_to_gene_map.getOrDefault(sp, sp) + " " + this.getSignificantSeedProteinQvalues().get(sp);
				sp_enrich_out.add(out);
			}
			
			Utilities.writeEntries(sp_enrich_out, out_path);
		}
		
		/**
		 * Returns sign. seed variants sorted by q-value
		 * @return
		 */
		public List<String> getSignificanceSortedSeedProteins() {
			return sign_sorted_sps;
		}
		
		/**
		 * Returns map of sign. seed proteins and their respective q-values
		 * @return
		 */
		public Map<String, Double> getSignificantSeedProteinQvalues() {
			return sign_sp_qvalue;
		}
		
		/**
		 * Returns map of sign. seed proteins and their direction of change as + / -
		 * @return
		 */
		public Map<String, String> getSignificantSeedProteinDirections() {
			return sign_sp_direction;
		}
	}
}
