package framework;

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
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;

/**
 * Class for differential DACO complexes
 * @author Thorsten Will
 */
public class DiffComplexDetector {
	
	// given data
	private final Map<String, QuantDACOResultSet> group1;
	private final Map<String, QuantDACOResultSet> group2;
	private final double FDR;
	private final boolean parametric;
	private final boolean incorporate_supersets;
	private int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	
	// processed / determined data
	private final Set<HashSet<String>> seed_combination_variants;
	private final Map<HashSet<String>, LinkedList<HashSet<String>>> seed_combination_variant_superset;
	private final Map<HashSet<String>, LinkedList<Double>> group1_abundances;
	private final Map<HashSet<String>, LinkedList<Double>> group2_abundances;
	private Map<HashSet<String>, Double> group1_medians;
	private Map<HashSet<String>, Double> group2_medians;
	
	// differential analysis results
	private Map<HashSet<String>, Double> variants_raw_pvalues;
	private final Map<HashSet<String>, Double> significant_variants_qvalues;
	private final List<HashSet<String>> significance_sorted_variants;
	private final Map<HashSet<String>, String> significance_variants_directions;
	
	// helper objects
	private final ForkJoinPool pool;
	
	public DiffComplexDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR, boolean parametric, boolean incorporate_supersets, int no_threads) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.parametric = parametric;
		this.incorporate_supersets = incorporate_supersets;
		this.no_threads = no_threads;
		
		pool = new ForkJoinPool(this.no_threads);
		
		// determine potentially relevant seed variant combinations ...
		this.seed_combination_variants = new HashSet<>();
		for (QuantDACOResultSet qdr:group1.values())
			this.seed_combination_variants.addAll(qdr.getSeedToComplexMap().keySet());
		for (QuantDACOResultSet qdr:group2.values())
			this.seed_combination_variants.addAll(qdr.getSeedToComplexMap().keySet());
		
		// ... and subsets (if necessary)
		if (this.incorporate_supersets) {
			this.seed_combination_variant_superset = new HashMap<>();
			for (HashSet<String> variant:this.seed_combination_variants)
				for (HashSet<String> current_variant:this.seed_combination_variants) 
					if (current_variant.containsAll(variant) && !current_variant.equals(variant)) {
						if (!this.seed_combination_variant_superset.containsKey(variant))
							this.seed_combination_variant_superset.put(variant, new LinkedList<HashSet<String>>());
						this.seed_combination_variant_superset.get(variant).add(current_variant);
					}
		} else {
			this.seed_combination_variant_superset = null;
		}

		// determine abundance values
		this.group1_abundances = this.determineAbundanceOfSeedVariantsComplexes(group1);
		this.group2_abundances = this.determineAbundanceOfSeedVariantsComplexes(group2);
		
		// determine medians
		ForkJoinTask<Map<HashSet<String>, Double>> group1_task = pool.submit(() -> this.group1_abundances.entrySet().parallelStream().collect(Collectors.toMap(e -> e.getKey(), e -> Utilities.getMedian(e.getValue()))));
		ForkJoinTask<Map<HashSet<String>, Double>> group2_task = pool.submit(() -> this.group2_abundances.entrySet().parallelStream().collect(Collectors.toMap(e -> e.getKey(), e -> Utilities.getMedian(e.getValue()))));
		try {
			this.group1_medians = group1_task.get();
			this.group2_medians = group2_task.get();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// determine differential abundance, apply multiple hypothesis correction and filter to significant seed variants
		if (this.parametric) // parametric Welch test
			this.significant_variants_qvalues = this.determinePValuesParametric();
		else // non-parametric MWU test
			this.significant_variants_qvalues = this.determinePValuesNonParametric();
		
		// sort from most to least significant
		this.significance_sorted_variants = new ArrayList<>(this.significant_variants_qvalues.keySet());
		this.significance_sorted_variants.sort( (v1, v2) -> diffCompareTo(v1, v2));
		
		this.significance_variants_directions = this.determineDirections();
		
		pool.shutdownNow();
	}
	
	/**
	 * Constructor only for testing
	 * @param raw_pvalues_medians_file
	 * @param FDR
	 * @param no_threads
	 */
	public DiffComplexDetector(String raw_pvalues_medians_file, double FDR, int no_threads) {
		// empty parameters
		this.group1 = null;
		this.group2 = null;
		this.FDR = FDR;
		this.parametric = false;
		this.incorporate_supersets = false;
		this.no_threads = no_threads;
		
		// empty results
		this.significant_variants_qvalues = null;
		this.seed_combination_variants = null;
		this.seed_combination_variant_superset = null;
		this.pool = null;
		this.group1_abundances = null;
		this.group2_abundances = null;
		
		// read data
		this.variants_raw_pvalues = new HashMap<>();
		this.group1_medians = new HashMap<>();
		this.group2_medians = new HashMap<>();
		
		for (String line:Utilities.readFile(raw_pvalues_medians_file)) {
			String[] spl = line.trim().split("\\s+");
			HashSet<String> variant = new HashSet<>(Arrays.asList(spl[0].split(",")));
			double raw_p = Double.parseDouble(spl[1]);
			double group1_med = Double.parseDouble(spl[2]);
			double group2_med = Double.parseDouble(spl[3]);
			this.variants_raw_pvalues.put(variant, raw_p);
			this.group1_medians.put(variant, group1_med);
			this.group2_medians.put(variant, group2_med);
		}
		
		this.significance_sorted_variants = new ArrayList<>(this.variants_raw_pvalues.keySet()); // only for testing
		this.significance_variants_directions = this.determineDirections();
	}
	
	/**
	 * Custom compareTo function that first sorts by the adjusted p-value/q-value and breaks ties using the absolute differences of the median between groups
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int diffCompareTo(HashSet<String> v1, HashSet<String> v2) {
		int sign_compareTo = this.significant_variants_qvalues.get(v1).compareTo(this.significant_variants_qvalues.get(v2));
		
		if (sign_compareTo == 0) {
			double v1_med_diff = Math.abs(this.group1_medians.get(v1) - this.group2_medians.get(v1));
			double v2_med_diff = Math.abs(this.group1_medians.get(v2) - this.group2_medians.get(v2));
			return Double.compare(v2_med_diff, v1_med_diff); // note that we want to have the bigger difference as the smaller entry
		} else
			return sign_compareTo;
	}
	
	/**
	 * Determines the abundance values of each seed variant found with the overall abundance of complexes containing the exact variant in the samples of the groups.
	 * @param group
	 * @return map of seed variant and abundance in each sample
	 */
	private Map<HashSet<String>, LinkedList<Double>> determineAbundanceOfSeedVariantsComplexes(Map<String, QuantDACOResultSet> group) {
		// init empty data structure for abundance data
		Map<HashSet<String>, LinkedList<Double>> group_abundances = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants) {
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
		
		for (String sample:group.keySet()) {
			Map<HashSet<String>, Double> sample_abundances = precomputed_sample_abundances.get(sample);
			for (HashSet<String> variant:this.seed_combination_variants) {
				double abundance_values = sample_abundances.getOrDefault(variant, 0.0);
				
				// if intended, also take supersets into account
				if (this.incorporate_supersets && this.seed_combination_variant_superset.containsKey(variant)) {
					for (HashSet<String> current_variant:this.seed_combination_variant_superset.get(variant)) { // sufficient to sum over sample as it is zero anyhow
						abundance_values += sample_abundances.getOrDefault(current_variant, 0.0);
					}
				}
				
				group_abundances.get(variant).add(abundance_values);
			}
		}
		
		return group_abundances;
	}
	
	/**
	 * Applies parametric heteroscedastic two-sample t-test/Welch test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determinePValuesParametric() {
		TTest tt = new TTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants) {
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
	private Map<HashSet<String>, Double> determinePValuesNonParametric() {
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variants) {
			// MWU test
			double pm = mwu.mannWhitneyUTest(Utilities.getDoubleArray(this.group1_abundances.get(variant)), Utilities.getDoubleArray(this.group2_abundances.get(variant)));
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
			double median_g1 = this.group1_medians.get(variant);
			double median_g2 = this.group2_medians.get(variant);
			
			String sign = "-";
			if (median_g2 > median_g1)
				sign = "+";
			
			significance_variants_directions.put(variant, sign);
		}
		
		return significance_variants_directions;
	}
	
	public void printResults() {
		for (HashSet<String> variant: this.significance_sorted_variants) {
			System.out.println(variant + " " + this.significant_variants_qvalues.get(variant));
		}
	}

	// TODO: think about what is helpful, probably inherit classes that are more specific, like TF complexes
	
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
	 * Returns the if a parametric test or a non-parametric test was used
	 * @return
	 */
	public boolean getParametricTestUsed() {
		return this.parametric;
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
		return seed_combination_variants;
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
	 * Returns median abundance data of group1
	 * @return
	 */
	public Map<HashSet<String>, Double> getGroup1MedianAbundances() {
		return group1_medians;
	}

	/**
	 * Returns median abundance data of group2
	 * @return
	 */
	public Map<HashSet<String>, Double> getGroup2MedianAbundances() {
		return group2_medians;
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
	 * generic helper functions
	 */
	
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
	public SPEnrichment calculateTFEnrichment(double FDR, int iterations, int complex_participation_cutoff) {
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
				double med_diff = group2_medians.get(variant) - group1_medians.get(variant);
				double score = -Math.log(variants_raw_pvalues.get(variant));
				
				if (med_diff < 0)
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
			
			// for each TF: calc. enrichment score
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
				FDR_q /= normalized_pos_sp_ES.keySet().stream().filter(tf2 -> normalized_pos_sp_ES.get(tf2).doubleValue() >= norm_ES).count() / (double) normalized_pos_sp_ES.keySet().size(); // will be >1 since norm_ES is in there
				sign_sp_qvalue.put(sp, FDR_q);
				sign_sp_direction.put(sp, "+");
			}
			for (String sp:normalized_neg_sp_ES.keySet()) {
				double norm_ES = normalized_neg_sp_ES.get(sp);
				// calculate q-value
				double FDR_q = Math.max(normalized_neg_rnd_data.stream().filter(d -> d.doubleValue() <= norm_ES).count(), 1) / (double) normalized_neg_rnd_data.size();
				FDR_q /= normalized_neg_sp_ES.keySet().stream().filter(tf2 -> normalized_neg_sp_ES.get(tf2).doubleValue() <= norm_ES).count() / (double) normalized_neg_sp_ES.keySet().size(); // will be >1 since norm_ES is in there
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
				double v1_med_diff = group2_medians.get(v1) - group1_medians.get(v1);
				double v2_med_diff = group2_medians.get(v2) - group1_medians.get(v2);
				return Double.compare(v1_med_diff, v2_med_diff);
			} else
				return sign_compareTo;
		}
		
		/**
		 * Calculates enrichment score in the fashion of GSEA, 
		 * handles enrichment and depletion independently, returns [max_ES, min_ES]
		 * @param query_TF
		 * @param scores
		 * @param complex_list
		 * @return
		 */
		private double calculateEnrichmentScore(String query_TF, List<HashSet<String>> complex_list) {
			// determine sum of scores (N_R in org. paper)
			int i = 0;
			double N_R = 0.0;
			for (HashSet<String> complex:complex_list) {
				if (complex.contains(query_TF))
					N_R += Math.abs(sorted_scores[i]);
				++i;
			}
			
			double N_H = sp_in_complex_count.get(query_TF);
			double running_sum = 0.0;
			double max_dev = 0.0;
			i = 0;
			for (HashSet<String> complex:complex_list) {
				if (complex.contains(query_TF))
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
