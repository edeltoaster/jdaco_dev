package framework;

import java.io.File;
import java.util.ArrayList;
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
	private final boolean incorporate_supersets;
	private int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	
	// processed / determined data
	private final Map<HashSet<String>, LinkedList<HashSet<String>>> relevant_complexes;
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
	private Map<String, String> up_to_gene_map; // computed on demand, use getters
	private Set<String> seed_proteins; // computed on demand, use getters
	
	public DiffComplexDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR, boolean parametric, boolean paired, boolean incorporate_supersets, double min_variant_fraction, int no_threads) {
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
		
		// bookmark and countall potentially relevant complexes ...
		Map<HashSet<String>, Integer> count_map_g1 = new HashMap<>();
		group1.values().stream().forEach(sample -> sample.getSeedToComplexMap().values().stream().forEach(cl -> cl.stream().forEach(c -> count_map_g1.put(c, count_map_g1.getOrDefault(c, 0) + 1))));
		Map<HashSet<String>, Integer> count_map_g2 = new HashMap<>();
		group2.values().stream().forEach(sample -> sample.getSeedToComplexMap().values().stream().forEach(cl -> cl.stream().forEach(c -> count_map_g2.put(c, count_map_g2.getOrDefault(c, 0) + 1))));
		Set<HashSet<String>> unfiltered_complexes = new HashSet<>(count_map_g1.keySet());
		unfiltered_complexes.addAll(count_map_g2.keySet());
		
		// check if above consideration threshold and determine subsets (if necessary)
		this.relevant_complexes = new HashMap<>();
		double min_count_g1 = min_variant_fraction * group1.size();
		double min_count_g2 = min_variant_fraction * group2.size();
		for (HashSet<String> complex:unfiltered_complexes) {
			this.relevant_complexes.put(complex, new LinkedList<HashSet<String>>());
			this.relevant_complexes.get(complex).add(complex);
			int count_g1 = count_map_g1.getOrDefault(complex, 0);
			int count_g2 = count_map_g2.getOrDefault(complex, 0);
			
			if (this.incorporate_supersets)
				for (HashSet<String> current_complex:unfiltered_complexes) 
					if (current_complex.containsAll(complex) && !current_complex.equals(complex)) {
						this.relevant_complexes.get(complex).add(current_complex);
						count_g1 += count_map_g1.getOrDefault(current_complex, 0);
						count_g2 += count_map_g2.getOrDefault(current_complex, 0);
					}
			// filter
			if (count_g1 < min_count_g1 && count_g2 < min_count_g2)
				this.relevant_complexes.remove(complex);
		}
		
		// remove subsets if an exact set is also above threshold
		Set<HashSet<String>> to_remove = this.relevant_complexes.keySet().stream().filter(rc -> this.relevant_complexes.get(rc).stream().anyMatch(c -> !rc.equals(c) && this.relevant_complexes.containsKey(c))).collect(Collectors.toSet());
		this.relevant_complexes.keySet().removeAll(to_remove);
		
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
		for (HashSet<String> complex:this.relevant_complexes.keySet()) {
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
		
		// ensures ordering of samples for paired data even when different data structures are used
		List<String> samples = new ArrayList<>(group.keySet());
		samples.sort(String::compareTo);
		
		for (String sample:samples) {
			Map<HashSet<String>, Double> sample_abundances = precomputed_sample_abundances.get(sample);
			for (HashSet<String> complex:this.relevant_complexes.keySet()) {
				double abundance_values = sample_abundances.getOrDefault(complex, 0.0);
				
				// if intended, also take supersets into account
				if (this.incorporate_supersets) {
					for (HashSet<String> current_complex:this.relevant_complexes.get(complex)) { // sufficient to sum over sample as it is zero anyhow
						abundance_values += sample_abundances.getOrDefault(current_complex, 0.0);
					}
				}
				
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
		for (HashSet<String> complex:this.relevant_complexes.keySet()) {
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
		for (HashSet<String> complex:this.relevant_complexes.keySet()) {
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
		for (HashSet<String> complex:this.relevant_complexes.keySet()) {
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
		
		for (HashSet<String> complex:this.relevant_complexes.keySet()) {
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
	
	
	/*
	 * standard analysis
	 */
	
	/**
	 * Does a standard analysis for the special case of transcription factors as seed proteins.
	 * All output is written to the specified folder.
	 * @param output_folder
	 */
	public void diffTFComplAnalysis(String output_folder, GOAnnotator goa, String binding_data_path, double binding_data_threshold, int binding_d_min, int binding_d_max, boolean also_compute_SCC, Set<String> proteins_to_remove) {

		// some first output
		System.out.println(this.getNumberOfTests() + " complexes tested.");
		System.out.println(this.getSignificanceSortedVariants().size() + " diff. complexes overall.");
		
		// no need to look further if there is nothing to tell about
		if (this.getSignificanceSortedVariants().size() == 0) {
			System.out.println();
			return;
		}
		
		// gather some helping data
		Set<String> seed_tfs = this.getSeedProteins();
		Map<String, String> up_name_map = this.getUniprotToGeneMap();
		
		
		/*
		 * build diff. complexes overview
		 */
		
		Set<String> involved_tfs = new HashSet<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		Map<HashSet<String>, List<HashSet<String>>> tfc_to_complexes = new HashMap<>();
		for (HashSet<String> complex:this.getSignificanceSortedVariants()) {
			
			String sign = this.getSignificantVariantsDirections().get(complex);
			
			// determine TFs that are involved and note in relevant data structures
			HashSet<String> tfs = new HashSet<>(complex);
			tfs.retainAll(seed_tfs);
			involved_tfs.addAll(tfs);
			if (!tfc_to_complexes.containsKey(tfs))
				tfc_to_complexes.put(tfs, new LinkedList<HashSet<String>>());
			tfc_to_complexes.get(tfs).add(complex);
			
			String compl_string = complex.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			String tfs_string = tfs.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			double pval = this.getSignificantVariantsQValues().get(complex);
			String out_string = sign + " " + tfs_string + " : " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(out_string);
			} else {
				res_pos_all.add(out_string);
			}
		}
		
		System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. complexes.");
		System.out.println(tfc_to_complexes.size() + " TF combinations in diff. complexes.");
		
		// write results
		if (!new File(output_folder).exists())
			new File(output_folder).mkdir();
		Utilities.writeEntries(res_pos_all, output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, output_folder + "res_neg_all.txt");
		
		// do the same with some given proteins not considered
		if (proteins_to_remove != null) {
			res_pos_all.clear();
			res_neg_all.clear();
			for (HashSet<String> complex:this.getSignificanceSortedVariants()) {
				
				// skip complexes that involve a protein that should be removed/not be considered
				if (complex.stream().anyMatch(p -> proteins_to_remove.contains(p)))
					continue;
				
				String sign = this.getSignificantVariantsDirections().get(complex);
				
				// determine TFs that are involved
				HashSet<String> tfs = new HashSet<>(complex);
				tfs.retainAll(seed_tfs);
				
				String compl_string = complex.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				String tfs_string = tfs.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				double pval = this.getSignificantVariantsQValues().get(complex);
				String out_string = sign + " " + tfs_string + " : " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
				
				// distinguish between increased/positive abundance and diminishing/negative abundance
				if (sign.equals("-")) {
					res_neg_all.add(out_string);
				} else {
					res_pos_all.add(out_string);
				}
			}
			
			System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. complexes (pruned).");
			
			// write results
			if (!new File(output_folder).exists())
				new File(output_folder).mkdir();
			Utilities.writeEntries(res_pos_all, output_folder + "res_pos_pruned.txt");
			Utilities.writeEntries(res_neg_all, output_folder + "res_neg_pruned.txt");
		}
		
		/*
		 * build regulatory network
		 */
		
		// build information for regulatory network
		Map<String, String> med_abun_g1 = new HashMap<>();
		Map<String, String> med_abun_g2 = new HashMap<>();
		Map<String, String> directions = new HashMap<>();
		Map<String, String> GO_details = new HashMap<>();
		Map<String, String> GO_overall = new HashMap<>();
		for (HashSet<String> tf_comb:tfc_to_complexes.keySet()) {
			String[] annotation_output = this.getSPCAnnotations(tf_comb, tfc_to_complexes.get(tf_comb), goa);
			String tfc = tf_comb.toString();
			med_abun_g1.put(tfc, annotation_output[0]);
			med_abun_g2.put(tfc, annotation_output[1]);
			directions.put(tfc, annotation_output[2]);
			GO_details.put(tfc, annotation_output[3]);
			GO_overall.put(tfc, annotation_output[4]);
		}
		
		// read binding data
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(binding_data_path, involved_tfs, binding_data_threshold, involved_tfs);
		
		// build regulatory network
		System.out.println("Building regulatory network ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(tfc_to_complexes.keySet(), bdh, binding_d_min, binding_d_max, no_threads, 1);
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(output_folder + "regnet.txt");
		
		// write annotation data
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Med_abundance_G1", med_abun_g1);
		annotational_data.put("Med_abundance_G2", med_abun_g2);
		annotational_data.put("Direction", directions);
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
	 * Returns all (sub)complexes assessed and which variants of them are found
	 * @return
	 */
	public Map<HashSet<String>, LinkedList<HashSet<String>>> getRelevantComplexes() {
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
	
	/**
	 * Returns the set of seed proteins used in the complex predictions.
	 * @return
	 */
	public Set<String> getSeedProteins() {
		
		if (this.seed_proteins != null)
			return this.seed_proteins;
		
		this.seed_proteins = new HashSet<>();
		group1.values().stream().forEach(qdr -> seed_proteins.addAll(qdr.getAbundantSeedProteins()));
		group2.values().stream().forEach(qdr -> seed_proteins.addAll(qdr.getAbundantSeedProteins()));
		
		return this.seed_proteins;
	}
	
	
	/*
	 * helper functions
	 */
	
	/**
	 * Given the seed protein combination occurring across samples, direction, sample data and a GOA definition, count and sort their appearances as well as their inferred GO annotations.
	 * Returns [mean abundances group1 (rounded to 2 positions), mean abundances group2 (rounded to 2 positions), directions of change per complex, GO annotations details string, set of all GO annotations found]
	 * @param complexes
	 * @param goa
	 * @return
	 */
	private String[] getSPCAnnotations(HashSet<String> seed_protein_comb, List<HashSet<String>> complexes, GOAnnotator goa) {
		// precompute naming data
		Map<String, String> up_name = this.getUniprotToGeneMap();
		Map<HashSet<String>, String> names_map = new HashMap<>();
		complexes.stream().forEach(c -> names_map.put(c, String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList()))));
		
		// determine median abundance values
		String abun_g1_string = String.join(",", complexes.stream().map(c -> names_map.get(c) + ":" + String.format(Locale.US, "%.2g", this.group1_medians.getOrDefault(c, 0.0)) ).collect(Collectors.toList()));
		String abun_g2_string = String.join(",", complexes.stream().map(c -> names_map.get(c) + ":" + String.format(Locale.US, "%.2g", this.group2_medians.getOrDefault(c, 0.0)) ).collect(Collectors.toList()));
		
		// get direction of change
		String directions_string = String.join(",", complexes.stream().map(c -> names_map.get(c) + ":" + this.getSignificantVariantsDirections().get(c)).collect(Collectors.toList()));
		
		// annotate with GO annotations
		Map<Set<String>, String> GOA_map = new HashMap<>();
		complexes.stream().forEach(c -> GOA_map.put(c, goa.rateProteins(c)));
		GOA_map.entrySet().removeIf(e -> e.getValue().equals("/"));
		
		String GOAnnotation_string_details = "/";
		String GOAnnotation_string_overall = "/";
		List<HashSet<String>> GO_complexes = new ArrayList<>(complexes);
		GO_complexes.retainAll(GOA_map.keySet());
		if (!GO_complexes.isEmpty()) {
			GOAnnotation_string_details = String.join(",", GO_complexes.stream().map(c -> names_map.get(c) + ":" + GOA_map.get(c)).collect(Collectors.toList()));
			GOAnnotation_string_overall = String.join(",", GO_complexes.stream().map(c -> GOA_map.get(c)).collect(Collectors.toSet()));
		}
		
		// build output datastructure
		String[] output = new String[5];
		output[0] = abun_g1_string;
		output[1] = abun_g2_string;
		output[2] = directions_string;
		output[3] = GOAnnotation_string_details;
		output[4] = GOAnnotation_string_overall;
		
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
			for (HashSet<String> complexes:raw_pvalues.keySet()) {
				double med_diff = group2_medians.get(complexes) - group1_medians.get(complexes);
				double score = -Math.log(raw_pvalues.get(complexes));
				
				if (med_diff < 0)
					score *= -1;
				scores.put(complexes, score);
			}
			
			// construct sorted matching arrays of variant/score
			List<HashSet<String>> sorted_complex_list = new ArrayList<>(scores.keySet());
			sorted_complex_list.sort((v1, v2) -> diffScoresCompareTo(v2, v1, scores));
			sorted_scores = new double[sorted_complex_list.size()];
			for (int i = 0; i< sorted_scores.length; i++)
				sorted_scores[i] = scores.get(sorted_complex_list.get(i));
			
			// determine seed proteins
			Set<String> seed_proteins = getSeedProteins();
			
			sp_in_complex_count = new HashMap<>();
			raw_pvalues.keySet().stream().forEach(proteins -> proteins.stream().filter(p -> seed_proteins.contains(p)).forEach(sp -> sp_in_complex_count.put(sp, sp_in_complex_count.getOrDefault(sp, 0.0) + 1)));
			
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
		
		public void writeSignificantSeedProteins(String pos_out_path, String neg_out_path) {
			
			List<String> pos_tf_enrich_out = new LinkedList<>();
			List<String> neg_tf_enrich_out = new LinkedList<>();
			Map<String, String> up_to_gene_map = getUniprotToGeneMap();
			for (String tf:this.getSignificanceSortedSeedProteins()) {
				String dir = this.getSignificantSeedProteinDirections().get(tf);
				String out = tf + " " + up_to_gene_map.getOrDefault(tf, tf) + " " + this.getSignificantSeedProteinQvalues().get(tf);
				
				if (dir.equals("+"))
					pos_tf_enrich_out.add(out);
				else
					neg_tf_enrich_out.add(out);
			}
			Utilities.writeEntries(pos_tf_enrich_out, pos_out_path);
			Utilities.writeEntries(neg_tf_enrich_out, neg_out_path);
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
