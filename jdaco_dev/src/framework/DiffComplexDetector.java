package framework;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
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
	private final Map<HashSet<String>, Double> significance_variants_pvalues;
	private final List<HashSet<String>> significance_sorted_variants;
	
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
			this.significance_variants_pvalues = this.determinePValuesParametric();
		else // non-parametric MWU test
			this.significance_variants_pvalues = this.determinePValuesNonParametric();
		
		// sort from most to least significant
		this.significance_sorted_variants = new ArrayList<>(this.significance_variants_pvalues.keySet());
		this.significance_sorted_variants.sort( (v1, v2) -> diffCompareTo(v1, v2));
		
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
		this.significance_variants_pvalues = null;
		this.significance_sorted_variants = null;
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
		
	}
	
	/**
	 * Custom compareTo function that first sorts by the adjusted p-value and breaks ties using the absolute differences of the median between groups
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int diffCompareTo(HashSet<String> v1, HashSet<String> v2) {
		int sign_compareTo = this.significance_variants_pvalues.get(v1).compareTo(this.significance_variants_pvalues.get(v2));
		
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
	
	public void printResults() {
		for (HashSet<String> variant: this.significance_sorted_variants) {
			System.out.println(variant + " " + this.significance_variants_pvalues.get(variant));
		}
	}

	// TODO: think about what is helpful, probably inherit classes that are more specific, like TF complexes
	
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
			double v1_med_diff = this.group2_medians.get(v1) - this.group1_medians.get(v1);
			double v2_med_diff = this.group2_medians.get(v2) - this.group1_medians.get(v2);
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
	private Double[] calculateEnrichmentScore(String query_TF, double[] scores, List<HashSet<String>> complex_list) {
		double running_sum = 0.0;
		double max_sum = 0.0;
		double min_sum = 0.0;
		
		int i = 0;
		for (HashSet<String> complex:complex_list) {
			if (complex.contains(query_TF))
				running_sum += scores[i];
			else
				running_sum -= scores[i];
			// TODO: do differently
			if (running_sum > max_sum)
				max_sum = running_sum;
			else if (running_sum < min_sum)
				min_sum = running_sum;
			
			++i;
		}
		
		return new Double[]{max_sum, min_sum};
	}
	
	public Map<String, Double> determineGSEAStyle(double FDR, int iterations) {
		// TODO: extra class
		Map<HashSet<String>, Double> scores = new HashMap<>();
		for (HashSet<String> variant:this.variants_raw_pvalues.keySet()) {
			double med_diff = this.group2_medians.get(variant) - this.group1_medians.get(variant);
			double score = -Math.log(this.variants_raw_pvalues.get(variant));
			
			if (med_diff < 0)
				score *= -1;
			scores.put(variant, score);
		}
		
		// construct sorted matching arrays of variant/score
		List<HashSet<String>> sorted_complex_list = new ArrayList<>(scores.keySet());
		sorted_complex_list.sort((v1, v2) -> diffScoresCompareTo(v2, v1, scores));
		double[] sorted_scores = new double[sorted_complex_list.size()];
		for (int i = 0; i< sorted_scores.length; i++)
			sorted_scores[i] = scores.get(sorted_complex_list.get(i));
		
		// determine all TFs
		Set<String> TFs = new HashSet<>();
		this.variants_raw_pvalues.keySet().stream().forEach(v -> TFs.addAll(v));
		
		// for each TF: calc. enrichment score
		Map<String, Double> tf_enrichment_scores = new HashMap<>();
		Map<String, Double> tf_depletion_scores = new HashMap<>();
		for (String tf:TFs) {
			Double[] tf_ES = calculateEnrichmentScore(tf, sorted_scores, sorted_complex_list);
			tf_enrichment_scores.put(tf, tf_ES[0].doubleValue());
			tf_depletion_scores.put(tf, tf_ES[1].doubleValue());
		}
		
		// get distr of NULL distributions by permutation
		List<List<HashSet<String>>> randomized_complex_lists = new ArrayList<>(iterations);
		for (int i = 0; i< iterations; i++) {
			List<HashSet<String>> complex_list = new ArrayList<>(sorted_complex_list);
			Collections.shuffle(complex_list);
			randomized_complex_lists.add(complex_list);
		}
		
		// check p-val for each
		Map<String, Double> tfdir_raw_pvalues = new HashMap<>();
		ForkJoinPool pool = new ForkJoinPool(this.no_threads);
		for (String tf:TFs) {
			double tf_ES = tf_enrichment_scores.get(tf);
			double tf_DS = tf_depletion_scores.get(tf);
			ForkJoinTask<List<Double[]>> rnd_counts = pool.submit(() -> randomized_complex_lists.parallelStream().map(l -> calculateEnrichmentScore(tf, sorted_scores, l)).collect(Collectors.toList()));
			List<Double[]> rnd_data = null;
			try {
				rnd_data = rnd_counts.get();
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			
			long enrich_counts = Math.max(rnd_data.stream().filter(d -> d[0].doubleValue() > tf_ES).count(), 1);
			long depl_counts = Math.max(rnd_data.stream().filter(d -> d[1].doubleValue() < tf_DS).count(), 1);
			
			double enrich_p = enrich_counts / (double) iterations;
			double depl_p = depl_counts / (double) iterations;
			tfdir_raw_pvalues.put("+" + tf, enrich_p);
			tfdir_raw_pvalues.put("-" + tf, depl_p);
		}
		
		// p-value correction and return
		return Utilities.convertRawPValuesToBHFDR(tfdir_raw_pvalues, FDR);
		
	}
	
	
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
	 * Returns map of significant seed combination variants and associated adjusted p-values
	 * @return
	 */
	public Map<HashSet<String>, Double> getSignificanceVariantsPValues() {
		return significance_variants_pvalues;
	}

	/**
	 * Returns map of seed combination variants and associated raw p-values
	 * @return
	 */
	public Map<HashSet<String>, Double> getVariantsRawPValues() {
		return variants_raw_pvalues;
	}
	
	/**
	 * Returns a list of significantly deregulated seed combination variants in their order of significance
	 * @return
	 */
	public List<HashSet<String>> getSignificanceSortedVariants() {
		return significance_sorted_variants;
	}
	
}
