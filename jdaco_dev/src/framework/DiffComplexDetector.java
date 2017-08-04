package framework;

import java.util.ArrayList;
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
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

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
	private final boolean paired;
	private int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	
	// processed / determined data
	private final Set<HashSet<String>> relevant_complexes;
	private final Map<HashSet<String>, LinkedList<Double>> group1_abundances;
	private final Map<HashSet<String>, LinkedList<Double>> group2_abundances;
	private Map<HashSet<String>, Double> group1_medians;
	private Map<HashSet<String>, Double> group2_medians;
	
	// differential analysis results
	private Map<HashSet<String>, Double> raw_pvalues;
	private final Map<HashSet<String>, Double> qvalues;
	private final List<HashSet<String>> significance_sorted_complexes;
	private final Map<HashSet<String>, String> directions;
	
	// helper objects
	private final ForkJoinPool pool;
	
	public DiffComplexDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR, boolean parametric, boolean paired, double min_variant_fraction, int no_threads) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.parametric = parametric;
		this.paired = paired;
		this.no_threads = no_threads;
		
		pool = new ForkJoinPool(this.no_threads);
		
		// ensure every sample has a matching counterpart in the other group
		if (this.paired)
			if (!group1.keySet().equals(group2.keySet())) {
				System.err.println("Non-equal sample groups supplied as matched data.");
				System.exit(1);
			}
		
		// determine potentially relevant seed variant combinations ...
		Map<HashSet<String>, Integer> count_map = new HashMap<>();
		group1.values().stream().forEach(sample -> sample.getSeedToComplexMap().values().stream().forEach(cl -> cl.stream().forEach(c -> count_map.put(c, count_map.getOrDefault(c, 0) + 1))));
		double min_count_g1 = min_variant_fraction * group1.size();
		count_map.entrySet().removeIf(e -> e.getValue().intValue() < min_count_g1);
		this.relevant_complexes = new HashSet<>(count_map.keySet());
		count_map.clear();
		group2.values().stream().forEach(sample -> sample.getSeedToComplexMap().values().stream().forEach(cl -> cl.stream().forEach(c -> count_map.put(c, count_map.getOrDefault(c, 0) + 1))));
		double min_count_g2 = min_variant_fraction * group2.size();
		count_map.entrySet().removeIf(e -> e.getValue().intValue() < min_count_g2);
		this.relevant_complexes.addAll(count_map.keySet());

		// determine abundance values
		this.group1_abundances = this.determineComplexesAbundances(group1);
		this.group2_abundances = this.determineComplexesAbundances(group2);
		
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
		if (!this.paired) {
			
			if (this.parametric) // parametric Welch test
				this.qvalues = this.determineUnpairedPValuesParametric();
			else // non-parametric MWU test
				this.qvalues = this.determineUnpairedPValuesNonParametric();
			
		} else {
			
			if (this.parametric) // paired t-test
				this.qvalues = this.determinePairedPValuesParametric();
			else // Wilcoxon signed-rank test
				this.qvalues = this.determinePairedPValuesNonParametric();
			
		}
		
		// sort from most to least significant
		this.significance_sorted_complexes = new ArrayList<>(this.qvalues.keySet());
		this.significance_sorted_complexes.sort( (v1, v2) -> diffCompareTo(v1, v2));
		
		this.directions = this.determineDirections();
		
		pool.shutdownNow();
	}
	
	/**
	 * Custom compareTo function that first sorts by the adjusted p-value/q-value and breaks ties using the absolute differences of the median between groups
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int diffCompareTo(HashSet<String> v1, HashSet<String> v2) {
		int sign_compareTo = this.qvalues.get(v1).compareTo(this.qvalues.get(v2));
		
		if (sign_compareTo == 0) {
			double v1_med_diff = Math.abs(this.group1_medians.get(v1) - this.group2_medians.get(v1));
			double v2_med_diff = Math.abs(this.group1_medians.get(v2) - this.group2_medians.get(v2));
			return Double.compare(v2_med_diff, v1_med_diff); // note that we want to have the bigger difference as the smaller entry
		} else
			return sign_compareTo;
	}
	
	/**
	 * Determines the abundance of complexes in the samples of the groups. 
	 * Ordering of samples is ensured to order paired data correctly.
	 * @param group
	 * @return map abundance in each sample
	 */
	private Map<HashSet<String>, LinkedList<Double>> determineComplexesAbundances(Map<String, QuantDACOResultSet> group) {
		// init empty data structure for abundance data
		Map<HashSet<String>, LinkedList<Double>> group_abundances = new HashMap<>();
		for (HashSet<String> complex:this.relevant_complexes) {
			group_abundances.put(complex, new LinkedList<Double>());
		}
		
		ForkJoinTask<Map<String, Map<HashSet<String>, Double>>> task = this.pool.submit(() -> group.entrySet().parallelStream().collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue().getAbundanceOfComplexes())));
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
			for (HashSet<String> complex:this.relevant_complexes) {
				double abundance_values = sample_abundances.getOrDefault(complex, 0.0);
				group_abundances.get(complex).add(abundance_values);
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
		for (HashSet<String> complex:this.relevant_complexes) {
			// two-sample two-tailed Welch test, t-test with differing variances (heteroscedastic)
			double pm = tt.tTest(Utilities.getDoubleArray(this.group1_abundances.get(complex)), Utilities.getDoubleArray(this.group2_abundances.get(complex)));
			test_results.put(complex, pm);
		}
		
		this.raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Applies non-parametric MWU test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determineUnpairedPValuesNonParametric() {
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> complex:this.relevant_complexes) {
			// MWU test
			double pm = mwu.mannWhitneyUTest(Utilities.getDoubleArray(this.group1_abundances.get(complex)), Utilities.getDoubleArray(this.group2_abundances.get(complex)));
			test_results.put(complex, pm);
		}
		
		this.raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Applies paired t-test and FDR-correction
	 * @return
	 */
	private Map<HashSet<String>, Double> determinePairedPValuesParametric() {
		TTest tt = new TTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> complex:this.relevant_complexes) {
			// paired t-test
			double pm = tt.pairedTTest(Utilities.getDoubleArray(this.group1_abundances.get(complex)), Utilities.getDoubleArray(this.group2_abundances.get(complex)));
			test_results.put(complex, pm);
		}
		
		this.raw_pvalues = test_results;
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
		
		for (HashSet<String> complex:this.relevant_complexes) {
			// Wilcoxon signed-rank test
			double pm = wsrt.wilcoxonSignedRankTest(Utilities.getDoubleArray(this.group1_abundances.get(complex)), Utilities.getDoubleArray(this.group2_abundances.get(complex)), compute_exact_p);
			test_results.put(complex, pm);
		}
		
		this.raw_pvalues = test_results;
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}
	
	/**
	 * Precomputes directions of change
	 * @return
	 */
	private Map<HashSet<String>, String> determineDirections() {
		
		Map<HashSet<String>, String> significance_variants_directions = new HashMap<>();
		for (HashSet<String> variant:this.significance_sorted_complexes) {
			double median_g1 = this.group1_medians.get(variant);
			double median_g2 = this.group2_medians.get(variant);
			
			String sign = "-";
			if (median_g2 > median_g1)
				sign = "+";
			
			if (median_g1 == median_g2) {
				// get means instead
				double mean_g1 = Utilities.getMean(this.group1_abundances.get(variant));
				double mean_g2 = Utilities.getMean(this.group2_abundances.get(variant));
				
				sign = "-";
				if (mean_g2 > mean_g1)
					sign = "+";
			}
			
			significance_variants_directions.put(variant, sign);
		}
		
		return significance_variants_directions;
	}
	
	public void printResults() {
		for (HashSet<String> variant: this.significance_sorted_complexes) {
			System.out.println(variant + " " + this.qvalues.get(variant));
		}
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
	 * Returns all complexes assessed
	 * @return
	 */
	public Set<HashSet<String>> getRelevantComplexes() {
		return relevant_complexes;
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
	 * Returns map of significant complexes and associated adjusted p-values / q-values
	 * @return
	 */
	public Map<HashSet<String>, Double> getSignificantVariantsQValues() {
		return qvalues;
	}

	/**
	 * Returns map of complexes and associated raw p-values
	 * @return
	 */
	public Map<HashSet<String>, Double> getRawPValues() {
		return raw_pvalues;
	}
	
	/**
	 * Returns the number of statistical tests made.
	 * @return
	 */
	public int getNumberOfTests() {
		return raw_pvalues.size();
	}
	
	/**
	 * Returns a list of significantly deregulated complexes in their order of significance
	 * @return
	 */
	public List<HashSet<String>> getSignificanceSortedVariants() {
		return significance_sorted_complexes;
	}
	
	/**
	 * Returns a map of sign. complexes to their direction of change as + / -
	 * @return
	 */
	public Map<HashSet<String>, String> getSignificantVariantsDirections() {
		return this.directions;
	}
}
