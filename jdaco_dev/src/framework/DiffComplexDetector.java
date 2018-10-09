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
		
		// bookmark and count all potentially relevant complexes ...
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
	 * Custom compareTo function that first sorts by the adjusted p-value/q-value and breaks ties using the amount of fold-change and absolute differences of the median between groups
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int diffCompareTo(HashSet<String> v1, HashSet<String> v2) {
		int sign_compareTo = this.qvalues.get(v1).compareTo(this.qvalues.get(v2));
		
		if (sign_compareTo == 0) {
			// determining fold-chance for v1
			double fold_change1 = Utilities.calcFoldChange(this.group1_medians.get(v1), this.group2_medians.get(v1));
			// ... and its extend
			fold_change1 = Utilities.amountFoldChange(fold_change1);
			
			// determining fold-chance for v2
			double fold_change2 = Utilities.calcFoldChange(this.group1_medians.get(v2), this.group2_medians.get(v2));
			// ... and its extend
			fold_change2 = Utilities.amountFoldChange(fold_change2);
			
			int sign_compareTo2 = Double.compare(fold_change2, fold_change1);
			
			if (sign_compareTo2 == 0) {
				double v1_med_diff = Math.abs(this.group1_medians.get(v1) - this.group2_medians.get(v1));
				double v2_med_diff = Math.abs(this.group1_medians.get(v2) - this.group2_medians.get(v2));
				return Double.compare(v2_med_diff, v1_med_diff); // note that we want to have the bigger difference as the smaller entry
			} else
				return sign_compareTo2;
			
		} else
			return sign_compareTo;
	}
	
	/**
	 * Custom compareTo function that first sorts by the adjusted p-value/q-value and breaks ties using the amount of fold-change and absolute differences of the median between groups
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int diffCompareTo2(HashSet<String> v1, HashSet<String> v2, Map<HashSet<String>, Double> all_qvalues) {
		int sign_compareTo = all_qvalues.get(v1).compareTo(all_qvalues.get(v2));
		
		if (sign_compareTo == 0) {
			// determining fold-chance for v1
			double fold_change1 = Utilities.calcFoldChange(this.group1_medians.get(v1), this.group2_medians.get(v1));
			// ... and its extend
			fold_change1 = Utilities.amountFoldChange(fold_change1);
			
			// determining fold-chance for v2
			double fold_change2 = Utilities.calcFoldChange(this.group1_medians.get(v2), this.group2_medians.get(v2));
			// ... and its extend
			fold_change2 = Utilities.amountFoldChange(fold_change2);
			
			int sign_compareTo2 = Double.compare(fold_change2, fold_change1);
			
			if (sign_compareTo2 == 0) {
				double v1_med_diff = Math.abs(this.group1_medians.get(v1) - this.group2_medians.get(v1));
				double v2_med_diff = Math.abs(this.group1_medians.get(v2) - this.group2_medians.get(v2));
				return Double.compare(v2_med_diff, v1_med_diff); // note that we want to have the bigger difference as the smaller entry
			} else
				return sign_compareTo2;
			
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
				double abundance_values = 0.0;
				for (HashSet<String> current_complex:this.relevant_complexes.get(complex)) // sufficient to sum over sample as it is zero anyhow
					abundance_values += sample_abundances.getOrDefault(current_complex, 0.0);
				
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
		
		Map<HashSet<String>, String> directions = new HashMap<>();
		for (HashSet<String> variant:this.raw_pvalues.keySet()) {
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
			
			directions.put(variant, sign);
		}
		
		return directions;
	}
	
	
	/*
	 * standard analysis
	 */
	
	/**
	 * Does a standard analysis for the special case of transcription factors as seed proteins.
	 * All output is written to the specified folder.
	 * @param output_folder
	 */
	public void diffTFComplAnalysis(String output_folder, GOAnnotator goa, String binding_data_path, double binding_data_threshold, int binding_d_min, int binding_d_max, boolean also_compute_SCC, Set<String> proteins_to_remove, Set<String> proteins_of_interest) {

		// some first output
		System.out.println(this.getNumberOfTests() + " complexes tested.");
		System.out.println(this.getSignificanceSortedComplexes().size() + " diff. complexes overall.");
		
		// no need to look further if there is nothing to tell about
		if (this.getSignificanceSortedComplexes().size() == 0) {
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
		for (HashSet<String> complex:this.getSignificanceSortedComplexes()) {
			
			String sign = this.getSignificantComplexesDirections().get(complex);
			
			// determine TFs that are involved and note in relevant data structures
			HashSet<String> tfs = new HashSet<>(complex);
			tfs.retainAll(seed_tfs);
			involved_tfs.addAll(tfs);
			if (!tfc_to_complexes.containsKey(tfs))
				tfc_to_complexes.put(tfs, new LinkedList<HashSet<String>>());
			tfc_to_complexes.get(tfs).add(complex);
			
			String compl_string = complex.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			String tfs_string = tfs.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			double pval = this.getSignificantComplexesQValues().get(complex);
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
			
		if (!output_folder_obj.exists())
			output_folder_obj.mkdir();
		
		// actually writing something
		Utilities.writeEntries(res_pos_all, output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, output_folder + "res_neg_all.txt");
		
		// do the same with some given proteins not considered
		if (proteins_to_remove != null) {
			res_pos_all.clear();
			res_neg_all.clear();
			for (HashSet<String> complex:this.getSignificanceSortedComplexes()) {
				
				// skip complexes that involve a protein that should be removed/not be considered
				if (complex.stream().anyMatch(p -> proteins_to_remove.contains(p)))
					continue;
				
				String sign = this.getSignificantComplexesDirections().get(complex);
				
				// determine TFs that are involved
				HashSet<String> tfs = new HashSet<>(complex);
				tfs.retainAll(seed_tfs);
				
				String compl_string = complex.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				String tfs_string = tfs.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				double pval = this.getSignificantComplexesQValues().get(complex);
				String out_string = sign + " " + tfs_string + " : " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
				
				// distinguish between increased/positive abundance and diminishing/negative abundance
				if (sign.equals("-")) {
					res_neg_all.add(out_string);
				} else {
					res_pos_all.add(out_string);
				}
			}
			
			System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. complexes (pruned).");
			
			Utilities.writeEntries(res_pos_all, output_folder + "res_pos_pruned.txt");
			Utilities.writeEntries(res_neg_all, output_folder + "res_neg_pruned.txt");
		}
		
		// do the same with only some given proteins considered
		Set<HashSet<String>> POI_sign_tfcs = new HashSet<>();
		if (proteins_of_interest != null) {
			res_pos_all.clear();
			res_neg_all.clear();
			for (HashSet<String> complex:this.getSignificanceSortedComplexes()) {
				
				// skip complexes that involve a protein that should be removed/not be considered
				if (!complex.stream().anyMatch(p -> proteins_of_interest.contains(p)))
					continue;
				
				String sign = this.getSignificantComplexesDirections().get(complex);
				
				// determine TFs that are involved
				HashSet<String> tfs = new HashSet<>(complex);
				tfs.retainAll(seed_tfs);
				POI_sign_tfcs.add(tfs);
				
				String compl_string = complex.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				String tfs_string = tfs.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
				double pval = this.getSignificantComplexesQValues().get(complex);
				String out_string = sign + " " + tfs_string + " : " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
				
				// distinguish between increased/positive abundance and diminishing/negative abundance
				if (sign.equals("-")) {
					res_neg_all.add(out_string);
				} else {
					res_pos_all.add(out_string);
				}
			}
			
			System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "- diff. complexes (POI).");
			
			Utilities.writeEntries(res_pos_all, output_folder + "res_pos_POI.txt");
			Utilities.writeEntries(res_neg_all, output_folder + "res_neg_POI.txt");
		}
		
		
		/*
		 * build regulatory network
		 */
		
		// build information for regulatory network
		Map<String, String> med_abun_g1 = new HashMap<>();
		Map<String, String> med_abun_g2 = new HashMap<>();
		Map<String, String> directions = new HashMap<>();
		Map<String, String> direction = new HashMap<>();
		Map<String, String> GO_details = new HashMap<>();
		Map<String, String> GO_overall = new HashMap<>();
		for (HashSet<String> tf_comb:tfc_to_complexes.keySet()) {
			String[] annotation_output = this.getSPCAnnotations(tf_comb, tfc_to_complexes.get(tf_comb), goa);
			String tfc = tf_comb.toString();
			med_abun_g1.put(tfc, annotation_output[0]);
			med_abun_g2.put(tfc, annotation_output[1]);
			directions.put(tfc, annotation_output[2]);
			direction.put(tfc, annotation_output[3]);
			GO_details.put(tfc, annotation_output[4]);
			GO_overall.put(tfc, annotation_output[5]);
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
		annotational_data.put("Med_abundance_G1", med_abun_g1);
		annotational_data.put("Med_abundance_G2", med_abun_g2);
		annotational_data.put("Direction_details", directions);
		annotational_data.put("Direction_overall", direction);
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
	 * Returns parsable output on significantly deregulated complexes in the space separated format (optionally Uniprot Accs are converted to gene identifiers and numbers are shortened):
	 * (sub)complex direction q-value fold-change median-change member_seed_comb member_complexes
	 * @param include_header
	 * @param human_readable
	 */
	public List<String> getSignSortedComplexes(boolean include_header, boolean human_readable) {
		List<String> to_write = new LinkedList<>();
		
		if (include_header)
			to_write.add("(sub)complex direction q-value fold-change median-change member_seed_comb member_complexes");
		
		if (human_readable)
			this.getUniprotToGeneMap();
		
		for (HashSet<String> sign_complex:this.significance_sorted_complexes) {
			
			// id conversion is necessary
			Set<String> sign_compl_temp = null;
			if (human_readable)
				sign_compl_temp = sign_complex.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet());
			else 
				sign_compl_temp = sign_complex;
			String compl_string = String.join("/", sign_compl_temp);
			
			// fold-change calculation
			double fold_change = Utilities.calcFoldChange(this.group1_medians.get(sign_complex), this.group2_medians.get(sign_complex));
			
			// median change calculation
			double median_change = this.group2_medians.get(sign_complex) - this.group1_medians.get(sign_complex);
			
			// determine member seed combination
			HashSet<String> seed_sub = new HashSet<>(sign_complex);
			seed_sub.retainAll(this.getSeedProteins());
			
			if (human_readable) {
				seed_sub = new HashSet<String>(seed_sub.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()));
			}
			String member_seed_comb_string = String.join("/", seed_sub);
			
			// determine member complexes
			List<HashSet<String>> members = this.relevant_complexes.get(sign_complex);
			if (human_readable) {
				List<HashSet<String>> members_temp = new ArrayList<>(members.size());
				members.stream().forEach(l -> members_temp.add(new HashSet<String>(l.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()))));
				members = members_temp;
			}
			String members_string = String.join(",", members.stream().map(l -> String.join("/", l)).collect(Collectors.toList()));
			
			// shortening numbers depending on human/machine-usage
			String qval_string = null;
			String foldc_string = null;
			String med_change_string = null;
			if (human_readable) {
				qval_string = String.format(Locale.US, "%.3g", this.qvalues.get(sign_complex));
				foldc_string = String.format(Locale.US, "%.3g", fold_change);
				med_change_string = String.format(Locale.US, "%.3g", median_change);
			} else {
				qval_string = Double.toString(this.qvalues.get(sign_complex));
				foldc_string = Double.toString(fold_change);
				med_change_string = Double.toString(median_change);
			}
			
			//format: (sub)complex direction q-value fold-change median-change member_seed_comb member_complexes
			List<String> line = Arrays.asList(compl_string, this.directions.get(sign_complex), qval_string, foldc_string, med_change_string, member_seed_comb_string, members_string);
			to_write.add(String.join(" ", line));
		}
		
		return to_write;
	}
	
	/**
	 * Returns parsable output on all complexes in the space separated format (optionally Uniprot Accs are converted to gene identifiers and numbers are shortened):
	 * (sub)complex direction q-value fold-change median-change member_seed_comb member_complexes
	 * @param include_header
	 * @param human_readable
	 */
	public List<String> getAllComplexes(boolean include_header, boolean human_readable) {
		List<String> to_write = new LinkedList<>();
		
		if (include_header)
			to_write.add("(sub)complex direction q-value fold-change median-change member_seed_comb member_complexes");
		
		if (human_readable)
			this.getUniprotToGeneMap();
		
		Map<HashSet<String>, Double> all_qvalues = Utilities.convertRawPValuesToBHFDRall(this.raw_pvalues);
		List<HashSet<String>> sorted_complexes = new ArrayList<>(all_qvalues.keySet());
		sorted_complexes.sort( (v1, v2) -> diffCompareTo2(v1, v2, all_qvalues));
		for (HashSet<String> complex:sorted_complexes) {
			
			// id conversion is necessary
			Set<String> sign_compl_temp = null;
			if (human_readable)
				sign_compl_temp = complex.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet());
			else 
				sign_compl_temp = complex;
			String compl_string = String.join("/", sign_compl_temp);
			
			// fold-change calculation
			double fold_change = Utilities.calcFoldChange(this.group1_medians.get(complex), this.group2_medians.get(complex));
			
			// median change calculation
			double median_change = this.group2_medians.get(complex) - this.group1_medians.get(complex);
			
			// determine member seed combination
			HashSet<String> seed_sub = new HashSet<>(complex);
			seed_sub.retainAll(this.getSeedProteins());
			
			if (human_readable) {
				seed_sub = new HashSet<String>(seed_sub.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()));
			}
			String member_seed_comb_string = String.join("/", seed_sub);
			
			// determine member complexes
			List<HashSet<String>> members = this.relevant_complexes.get(complex);
			if (human_readable) {
				List<HashSet<String>> members_temp = new ArrayList<>(members.size());
				members.stream().forEach(l -> members_temp.add(new HashSet<String>(l.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()))));
				members = members_temp;
			}
			String members_string = String.join(",", members.stream().map(l -> String.join("/", l)).collect(Collectors.toList()));
			
			// shortening numbers depending on human/machine-usage
			String qval_string = null;
			String foldc_string = null;
			String med_change_string = null;
			if (human_readable) {
				qval_string = String.format(Locale.US, "%.3g", all_qvalues.get(complex));
				foldc_string = String.format(Locale.US, "%.3g", fold_change);
				med_change_string = String.format(Locale.US, "%.3g", median_change);
			} else {
				qval_string = Double.toString(all_qvalues.get(complex));
				foldc_string = Double.toString(fold_change);
				med_change_string = Double.toString(median_change);
			}
			
			//format: (sub)complex direction q-value fold-change median-change member_seed_comb member_complexes
			List<String> line = Arrays.asList(compl_string, this.directions.get(complex), qval_string, foldc_string, med_change_string, member_seed_comb_string, members_string);
			to_write.add(String.join(" ", line));
		}
		
		return to_write;
	}
	
	/**
	 * Writes parsable output in the space separated format (optionally Uniprot Accs are converted to gene identifiers and numbers are shortened):
	 * (sub)complex direction q-value fold-change median-change member_seed_comb member_complexes
	 * @param out_file
	 * @param human_readable
	 */
	public void writeSignSortedComplexes(String out_file, boolean human_readable) {
		
		Utilities.writeEntries(this.getSignSortedComplexes(true, human_readable), out_file);
	}
	
	/**
	 * Writes parsable output in the space separated format (optionally Uniprot Accs are converted to gene identifiers and numbers are shortened):
	 * (sub)complex direction q-value fold-change median-change member_seed_comb member_complexes
	 * @param out_file
	 * @param human_readable
	 */
	public void writeAllComplexes(String out_file, boolean human_readable) {
		
		Utilities.writeEntries(this.getAllComplexes(true, human_readable), out_file);
	}
	
	/**
	 * Reads written output of writeSignSortedComplexes into a convenient data structure.
	 * @param output_file
	 * @return
	 */
	public static SignSortedComplexesResult readSignSortedComplexResult(String output_file) {
		return new SignSortedComplexesResult(output_file);
	}
	
	/**
	 * Reads written output of writeSignSortedComplexes into a convenient data structure (only until certain q-value is met).
	 * @param output_file
	 * @param qval_threshold
	 * @return
	 */
	public static SignSortedComplexesResult readSignSortedComplexResult(String output_file, double qval_threshold) {
		return new SignSortedComplexesResult(output_file, qval_threshold);
	}
	
	/**
	 * Returns parsable output in the space separated format (optionally Uniprot Accs are converted to gene identifiers and numbers are shortened):
	 * seed_comb direction direction_coarse avg_q-value sum_fold-change sum_median-change member_complexes direction_details q-value_details fold-change_details; in the order of avg_q-value, amount of fold-change and abs. sum. median change
	 * @param include_header
	 * @param human_readable
	 */
	public List<String> getSignSortedVariants(boolean include_header, boolean human_readable) {
		
		HashMap<HashSet<String>, List<Double>> sign_comb_qs_map = new HashMap<>();
		HashMap<HashSet<String>, List<Double>> sign_comb_medc_map = new HashMap<>();
		HashMap<HashSet<String>, List<HashSet<String>>> sign_comb_sign_compl_map = new HashMap<>();
		HashMap<HashSet<String>, Double> sign_compl_foldc = new HashMap<>();
		HashMap<HashSet<String>, Double> sign_compl_medc = new HashMap<>();
		HashMap<HashSet<String>, Double> sign_comb_g1_medsum = new HashMap<>();
		HashMap<HashSet<String>, Double> sign_comb_g2_medsum = new HashMap<>();
		
		// determine seed combination variants, collect p-values as well as note dereg. complexes, fold changes
		for (HashSet<String> sign_complex:this.significance_sorted_complexes) {
			double q_value = this.qvalues.get(sign_complex);
			
			// determine seed subset
			HashSet<String> seed_sub = new HashSet<>(sign_complex);
			seed_sub.retainAll(this.getSeedProteins());
			
			// store q-value(s)
			if (!sign_comb_qs_map.containsKey(seed_sub)) {
				sign_comb_qs_map.put(seed_sub, new LinkedList<Double>());
				sign_comb_sign_compl_map.put(seed_sub, new LinkedList<HashSet<String>>());
				sign_comb_medc_map.put(seed_sub, new LinkedList<Double>());
			}
			sign_comb_qs_map.get(seed_sub).add(q_value);
			sign_comb_sign_compl_map.get(seed_sub).add(sign_complex);
			
			// store fold change
			double fold_change = Utilities.calcFoldChange(this.group1_medians.get(sign_complex), this.group2_medians.get(sign_complex));
			
			sign_compl_foldc.put(sign_complex, fold_change);
			double medc = this.group2_medians.get(sign_complex) - this.group1_medians.get(sign_complex);
			sign_compl_medc.put(sign_complex, medc);
			
			// store median change
			sign_comb_medc_map.get(seed_sub).add(medc);
			
			// store medians of groups
			sign_comb_g1_medsum.put(seed_sub, this.group1_medians.get(sign_complex) + sign_comb_g1_medsum.getOrDefault(seed_sub, 0.0));
			sign_comb_g2_medsum.put(seed_sub, this.group2_medians.get(sign_complex) + sign_comb_g2_medsum.getOrDefault(seed_sub, 0.0));
		}
		
		// avg q-value
		HashMap<HashSet<String>, Double> sign_comb_q_map = new HashMap<>();
		sign_comb_qs_map.keySet().forEach(c -> sign_comb_q_map.put(c, Utilities.getMean(sign_comb_qs_map.get(c))));
		
		// get fold-change of summarized group medians
		HashMap<HashSet<String>, Double> sign_comb_fc_map = new HashMap<>();
		for (HashSet<String> seed_sub:sign_comb_g1_medsum.keySet()) {
			double sum_fold_change = Utilities.calcFoldChange(sign_comb_g1_medsum.get(seed_sub), sign_comb_g2_medsum.get(seed_sub));
			sign_comb_fc_map.put(seed_sub, sum_fold_change);
		}
		
		// summarized signed median-abundance change
		HashMap<HashSet<String>, Double> sign_comb_sumabun_map = new HashMap<>();
		sign_comb_medc_map.keySet().forEach(c -> sign_comb_sumabun_map.put(c, sign_comb_medc_map.get(c).stream().reduce(0.0, Double::sum)));
		
		// sort by min q-value, then by amount of fold change and median change
		List<HashSet<String>> sign_sorted_combs = new ArrayList<>(sign_comb_q_map.keySet());
		sign_sorted_combs.sort( (c1, c2) -> diffCompareToVar(c1, c2, sign_comb_q_map, sign_comb_fc_map, sign_comb_sumabun_map));
		
		// preface
		List<String> to_write = new LinkedList<>();
		
		if (include_header)
			to_write.add("seed_comb direction direction_coarse avg_q-value sum_fold-change sum_median-change member_complexes direction_details q-value_details fold-change_details");
		
		if (human_readable)
			this.getUniprotToGeneMap();
		
		for (HashSet<String> seed_comb:sign_sorted_combs) {
			
			// translate seed combination if necessary
			Set<String> seed_comb_temp = null;
			if (human_readable)
				seed_comb_temp = seed_comb.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet());
			else 
				seed_comb_temp = seed_comb;
			String seed_comb_string = String.join("/", seed_comb_temp);
			
			// determine direction and individual values for each complex
			List<String> member_complexes = new LinkedList<>();
			List<String> direction_details = new LinkedList<>();
			List<String> qval_details = new LinkedList<>();
			List<String> foldc_details = new LinkedList<>();
			List<String> medc_details = new LinkedList<>();
			
			Set<String> direction_set = new HashSet<>();
			for (HashSet<String> complex:sign_comb_sign_compl_map.get(seed_comb)) {
				
				// member complexes identifier conversion (if wanted)
				String compl_string = null;
				String qval_string = null;
				String foldc_string = null;
				String medc_string = null;
				if (human_readable) {
					compl_string = String.join("/", complex.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toSet()));
					qval_string = String.format(Locale.US, "%.3g", this.qvalues.get(complex));
					foldc_string = String.format(Locale.US, "%.3g", sign_compl_foldc.get(complex));
					medc_string = String.format(Locale.US, "%.3g", sign_compl_medc.get(complex));
				} else {
					compl_string = String.join("/", complex);
					qval_string = Double.toString(this.qvalues.get(complex));
					foldc_string = Double.toString(sign_compl_foldc.get(complex));
					medc_string = Double.toString(sign_compl_medc.get(complex));
				}
				member_complexes.add(compl_string);
				
				// direction stuff
				String direction = this.directions.get(complex);
				direction_set.add(direction);
				direction_details.add(compl_string + ":" + direction);
				
				// q-value details
				qval_details.add(compl_string + ":" + qval_string); 
				// fold-change details
				foldc_details.add(compl_string + ":" + foldc_string);
				// median-change details
				medc_details.add(compl_string + ":" + medc_string);
			}
			
			// determine overall direction
			String direction_string = "+";
			if (!direction_set.contains(direction_string))
				direction_string = "-";
			else if (direction_set.size() > 1)
				direction_string = "+-";
			
			String direction_coarse_string = null;
			if (direction_string.equals("+-")) {
				double foldc = sign_comb_fc_map.get(seed_comb);
				if (foldc > 2.0)
					direction_coarse_string = "+";
				else if (foldc < 0.5)
					direction_coarse_string = "-";
				else 
					direction_coarse_string = direction_string;
			} else {
				direction_coarse_string = direction_string;
			}
			// merge details strings
			String direction_details_string = String.join(",", direction_details);
			String member_complexes_string = String.join(",", member_complexes);
			String qval_details_string = String.join(",", qval_details);
			String foldc_details_string = String.join(",", foldc_details);
			String medc_details_string = String.join(",", medc_details);
			
			// shortening numbers depending on human/machine-usage
			String qval_string = null;
			String foldc_string = null;
			String sum_med_change_string = null;
			if (human_readable) {
				qval_string = String.format(Locale.US, "%.3g", sign_comb_q_map.get(seed_comb));
				foldc_string = String.format(Locale.US, "%.3g", sign_comb_fc_map.get(seed_comb));
				sum_med_change_string = String.format(Locale.US, "%.3g", sign_comb_sumabun_map.get(seed_comb));
			} else {
				qval_string = Double.toString(sign_comb_q_map.get(seed_comb));
				foldc_string = Double.toString(sign_comb_fc_map.get(seed_comb));
				sum_med_change_string = Double.toString(sign_comb_sumabun_map.get(seed_comb));
			}
			
			// write in format: seed_comb direction direction_coarse avg_q-value sum_fold-change sum_median-change member_complexes direction_details q-value_details fold-change_details median-change_details
			List<String> line = Arrays.asList(seed_comb_string, direction_string, direction_coarse_string, qval_string, foldc_string, sum_med_change_string, member_complexes_string, direction_details_string, qval_details_string, foldc_details_string, medc_details_string);
			to_write.add(String.join(" ", line));
		}
		
		return to_write;
	}
		
	/**
	 * Writes parsable output in the space separated format (optionally Uniprot Accs are converted to gene identifiers and numbers are shortened):
	 * seed_comb direction direction_coarse avg_q-value sum_fold-change sum_median-change member_complexes direction_details q-value_details fold-change_details; in the order of avg_q-value, amount of fold-change and abs. sum. median change
	 * @param out_file
	 * @param human_readable
	 */
	public void writeSignSortedVariants(String out_file, boolean human_readable) {
		
		Utilities.writeEntries(this.getSignSortedVariants(true, human_readable), out_file);
	}
	
	/**
	 * Custom compareTo function that first sorts by the avg. q-value, then by the extend of fold-change between groups and then by the absolute median change
	 * @param v1
	 * @param v2
	 * @param sign_comb_q_map
	 * @param sign_comb_foldc_map
	 * @param sign_comb_medc_map
	 * @return
	 */
	private int diffCompareToVar(HashSet<String> v1, HashSet<String> v2, HashMap<HashSet<String>, Double> sign_comb_q_map, HashMap<HashSet<String>, Double> sign_comb_foldc_map, HashMap<HashSet<String>, Double> sign_comb_medc_map) {
		int sign_compareTo = sign_comb_q_map.get(v1).compareTo(sign_comb_q_map.get(v2));
		
		if (sign_compareTo == 0) {
			double f1 = sign_comb_foldc_map.get(v1);
			f1 = Utilities.amountFoldChange(f1);
			
			double f2 = sign_comb_foldc_map.get(v2);
			f2 = Utilities.amountFoldChange(f2);
			
			int sign_compareTo2 = Double.compare(f2, f1); // note that we want to have the bigger difference as the smaller entry
			
			if (sign_compareTo2 == 0) {
				return Double.compare(Math.abs(sign_comb_medc_map.get(v2)), Math.abs(sign_comb_medc_map.get(v1))); // note that we want to have the bigger difference as the smaller entry
			} else
				return sign_compareTo2;
			
		} else
			return sign_compareTo;
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
	public Map<HashSet<String>, Double> getSignificantComplexesQValues() {
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
	public List<HashSet<String>> getSignificanceSortedComplexes() {
		return significance_sorted_complexes;
	}
	
	/**
	 * Returns a set of seed protein combination variants that are part of significantly deregulated complexes
	 * @return
	 */
	public Set<HashSet<String>> getSignificantSeedCombVariants() {
		
		Set<HashSet<String>> sign_seed_combs = new HashSet<>();
		for (HashSet<String> sign_complex:this.significance_sorted_complexes) {
			// determine seed subset
			HashSet<String> seed_sub = new HashSet<>(sign_complex);
			seed_sub.retainAll(this.getSeedProteins());
			
			sign_seed_combs.add(seed_sub);
		}
		
		return sign_seed_combs;
	}
	
	/**
	 * Returns a map of sign. complexes to their direction of change as + / -
	 * @return
	 */
	public Map<HashSet<String>, String> getSignificantComplexesDirections() {
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
	 * Returns [mean abundances group1 (rounded to 2 positions), mean abundances group2 (rounded to 2 positions), directions of change per complex, overall direction change, GO annotations details string, set of all GO annotations found]
	 * @param complexes
	 * @param goa
	 * @return
	 */
	private String[] getSPCAnnotations(HashSet<String> seed_protein_comb, List<HashSet<String>> complexes, GOAnnotator goa) {
		// precompute naming data
		Map<String, String> up_name = this.getUniprotToGeneMap();
		Map<HashSet<String>, String> names_map = new HashMap<>();
		// complexes sorted by p-val
		complexes.stream().forEach(c -> names_map.put(c, String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList()))));
		
		// determine median abundance values
		String abun_g1_string = String.join(",", complexes.stream().map(c -> names_map.get(c) + ":" + String.format(Locale.US, "%.2g", this.group1_medians.getOrDefault(c, 0.0)) ).collect(Collectors.toList()));
		String abun_g2_string = String.join(",", complexes.stream().map(c -> names_map.get(c) + ":" + String.format(Locale.US, "%.2g", this.group2_medians.getOrDefault(c, 0.0)) ).collect(Collectors.toList()));
		
		// get direction of change
		String directions_string = String.join(",", complexes.stream().map(c -> names_map.get(c) + ":" + this.getSignificantComplexesDirections().get(c)).collect(Collectors.toList()));
		Set<String> direction_set = new HashSet<>();
		for (HashSet<String> complex:complexes)
			direction_set.add(this.getSignificantComplexesDirections().get(complex));
		
		// determine overall direction
		String direction_string = "+";
		if (!direction_set.contains(direction_string))
			direction_string = "-";
		else if (direction_set.size() > 1)
			direction_string = "+-";
		
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
		String[] output = new String[6];
		output[0] = abun_g1_string;
		output[1] = abun_g2_string;
		output[2] = directions_string;
		output[3] = direction_string;
		output[4] = GOAnnotation_string_details;
		output[5] = GOAnnotation_string_overall;
		
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
	
	
	/*
	 * Seed protein combination enrichment
	 */
	
	/**
	 * Compute seed protein combination enrichment in the fashion of GSEA, 
	 * see Subramanian et al. (2005)
	 * @param FDR
	 * @param iterations
	 * @param complex_participation_cutoff
	 * @return
	 */
	public SPCEnrichment calculateSPCEnrichment(double FDR, int iterations, int complex_participation_cutoff) {
		return new SPCEnrichment(FDR, iterations, complex_participation_cutoff);
	}
	
	/**
	 * Compute seed protein enrichment in the fashion of GSEA, 
	 * see Subramanian et al. (2005)
	 */
	public final class SPCEnrichment {
		
		private final double[] sorted_scores;
		private final Map<HashSet<String>, Double> spc_in_complex_count;
		private final double no_complexes;
		private final Map<HashSet<String>, HashSet<String>> complex_to_spc;
		private final List<HashSet<String>> sign_sorted_sps;
		private final Map<HashSet<String>, Double> sign_spc_qvalue;
		private final Map<HashSet<String>, String> sign_spc_direction;
		
		public SPCEnrichment(double FDR, int iterations, int complex_participation_cutoff) {
			
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
			complex_to_spc = new HashMap<>();
			spc_in_complex_count = new HashMap<>();
			for (HashSet<String> complex:raw_pvalues.keySet()) {
				HashSet<String> spc_in_complex = new HashSet<>(complex);
				spc_in_complex.retainAll(seed_proteins);
				
				if (spc_in_complex.size() > 0) {
					spc_in_complex_count.put(spc_in_complex, spc_in_complex_count.getOrDefault(spc_in_complex, 0.0) + 1);
					complex_to_spc.put(complex, spc_in_complex);
				}
			}
			
			// remove SPCs that are not in complexes very often
			spc_in_complex_count.keySet().removeIf(spc -> spc_in_complex_count.get(spc).doubleValue() < complex_participation_cutoff);
			
			// precompute #complexes (analogous to N in org. paper)
			no_complexes = sorted_scores.length;
			
			// for each seed protein combination: calc. enrichment score
			Map<HashSet<String>, Double> spc_enrichment_scores = new HashMap<>();
			for (HashSet<String> spc:spc_in_complex_count.keySet())
				spc_enrichment_scores.put(spc, calculateEnrichmentScore(spc, sorted_complex_list));

			// get distr of NULL distributions by permutation
			List<List<HashSet<String>>> randomized_complex_lists = new ArrayList<>(iterations);
			for (int i = 0; i< iterations; i++) {
				List<HashSet<String>> complex_list = new ArrayList<>(sorted_complex_list);
				Collections.shuffle(complex_list);
				randomized_complex_lists.add(complex_list);
			}
			
			// compute ES and distributions for each seed protein combination of interest
			ForkJoinPool pool = new ForkJoinPool(no_threads);
			Map<HashSet<String>, Double> normalized_pos_sp_ES = new HashMap<>();
			Map<HashSet<String>, Double> normalized_neg_sp_ES = new HashMap<>();
			List<Double> normalized_pos_rnd_data = new LinkedList<>();
			List<Double> normalized_neg_rnd_data = new LinkedList<>();
			for (HashSet<String> spc:spc_in_complex_count.keySet()) {
				double spc_ES = spc_enrichment_scores.get(spc);
				ForkJoinTask<List<Double>> rnd_counts = pool.submit(() -> randomized_complex_lists.parallelStream().map(l -> calculateEnrichmentScore(spc, l)).collect(Collectors.toList()));
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
				
				if (spc_ES >= 0)
					normalized_pos_sp_ES.put(spc, spc_ES / pos_mean);
				else
					normalized_neg_sp_ES.put(spc, spc_ES / neg_mean);
				
			}
			
			// corrected FDR q-values
			sign_spc_qvalue = new HashMap<>();
			sign_spc_direction = new HashMap<>();
			for (HashSet<String> spc:normalized_pos_sp_ES.keySet()) {
				double norm_ES = normalized_pos_sp_ES.get(spc);
				// calculate q-value
				double FDR_q = Math.max(normalized_pos_rnd_data.stream().filter(d -> d.doubleValue() >= norm_ES).count(), 1) / (double) normalized_pos_rnd_data.size();
				FDR_q /= normalized_pos_sp_ES.keySet().stream().filter(sp2 -> normalized_pos_sp_ES.get(sp2).doubleValue() >= norm_ES).count() / (double) normalized_pos_sp_ES.keySet().size(); // will be >1 since norm_ES is in there
				sign_spc_qvalue.put(spc, FDR_q);
				sign_spc_direction.put(spc, "+");
			}
			for (HashSet<String> spc:normalized_neg_sp_ES.keySet()) {
				double norm_ES = normalized_neg_sp_ES.get(spc);
				// calculate q-value
				double FDR_q = Math.max(normalized_neg_rnd_data.stream().filter(d -> d.doubleValue() <= norm_ES).count(), 1) / (double) normalized_neg_rnd_data.size();
				FDR_q /= normalized_neg_sp_ES.keySet().stream().filter(sp2 -> normalized_neg_sp_ES.get(sp2).doubleValue() <= norm_ES).count() / (double) normalized_neg_sp_ES.keySet().size(); // will be >1 since norm_ES is in there
				sign_spc_qvalue.put(spc, FDR_q);
				sign_spc_direction.put(spc, "-");
			}
			
			// only retain those below threshold
			sign_spc_qvalue.keySet().removeIf(tf -> sign_spc_qvalue.get(tf).doubleValue() >= FDR);
			sign_sorted_sps = new ArrayList<>(sign_spc_qvalue.keySet());
			sign_sorted_sps.sort((v1, v2) -> sign_spc_qvalue.get(v1).compareTo(sign_spc_qvalue.get(v2)));
		}
		
		/**
		 * Custom compareTo function that first sorts by the log-scores and breaks ties using the differences of the median between groups.
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
		 * @param query_SPC
		 * @param scores
		 * @param complex_list
		 * @return
		 */
		private double calculateEnrichmentScore(HashSet<String> query_SPC, List<HashSet<String>> complex_list) {
			// determine sum of scores (N_R in org. paper)
			int i = 0;
			double N_R = 0.0;
			for (HashSet<String> complex:complex_list) {
				if (complex_to_spc.get(complex).equals(query_SPC))
					N_R += Math.abs(sorted_scores[i]);
				++i;
			}
			
			double N_H = spc_in_complex_count.get(query_SPC);
			double running_sum = 0.0;
			double max_dev = 0.0;
			i = 0;
			for (HashSet<String> complex:complex_list) {
				if (complex_to_spc.get(complex).equals(query_SPC))
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
		public void writeSignificantSeedProteinCombinations(String pos_out_path, String neg_out_path) {
			
			List<String> pos_spc_enrich_out = new LinkedList<>();
			List<String> neg_spc_enrich_out = new LinkedList<>();
			Map<String, String> up_to_gene_map = getUniprotToGeneMap();
			for (HashSet<String> spc:this.getSignificanceSortedSeedProteinCombinations()) {
				String dir = this.getSignificantSeedProteinCombDirections().get(spc);
				String out = String.join("/", spc) + " " + String.join("/", spc.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toList())) + " " + this.getSignificantSeedProteinCombQvalues().get(spc);
				
				if (dir.equals("+"))
					pos_spc_enrich_out.add(out);
				else
					neg_spc_enrich_out.add(out);
			}
			
			Utilities.writeEntries(pos_spc_enrich_out, pos_out_path);
			Utilities.writeEntries(neg_spc_enrich_out, neg_out_path);
		}
		
		/**
		 * Writes result into a specified text file.
		 * @param out_path
		 */
		public void writeSignificantSeedProteinCombinations(String out_path) {
			
			List<String> spc_enrich_out = new LinkedList<>();
			Map<String, String> up_to_gene_map = getUniprotToGeneMap();
			for (HashSet<String> spc:this.getSignificanceSortedSeedProteinCombinations()) {
				String dir = this.getSignificantSeedProteinCombDirections().get(spc);
				String out = dir + " " + String.join("/", spc) + " " + String.join("/", spc.stream().map(p -> up_to_gene_map.getOrDefault(p, p)).collect(Collectors.toList())) + " " + this.getSignificantSeedProteinCombQvalues().get(spc);
				
				spc_enrich_out.add(out);
			}
			
			Utilities.writeEntries(spc_enrich_out, out_path);
		}
		
		/**
		 * Returns sign. seed protein combinations sorted by q-value
		 * @return
		 */
		public List<HashSet<String>> getSignificanceSortedSeedProteinCombinations() {
			return sign_sorted_sps;
		}
		
		/**
		 * Returns map of sign. seed protein combinations and their respective q-values
		 * @return
		 */
		public Map<HashSet<String>, Double> getSignificantSeedProteinCombQvalues() {
			return sign_spc_qvalue;
		}
		
		/**
		 * Returns map of sign. seed proteins and their direction of change as + / -
		 * @return
		 */
		public Map<HashSet<String>, String> getSignificantSeedProteinCombDirections() {
			return sign_spc_direction;
		}
	}
	
	
	/*
	 * Convenient helper data structure to transfer results
	 */
	
	public static final class SignSortedComplexesResult {
		
		private final List<HashSet<String>> significance_sorted_complexes = new LinkedList<>();
		private final Map<HashSet<String>, String> directions = new HashMap<>();
		private final Map<HashSet<String>, Double> qvalues = new HashMap<>();
		private final Map<HashSet<String>, Double> fold_change = new HashMap<>();
		private final Map<HashSet<String>, Double> median_change = new HashMap<>();
		private final Map<HashSet<String>, HashSet<String>> member_seed_comb = new HashMap<>();
		private final Map<HashSet<String>, List<HashSet<String>>> member_seed_complexes = new HashMap<>();

		//format: 0 (sub)complex 1 direction 2 q-value 3 fold-change 4 median-change 5 member_seed_comb 6 member_complexes
		public SignSortedComplexesResult(String in_file) {
			for (String line:Utilities.readFile(in_file)) {
				if (line.startsWith("("))
					continue;
				String[] spl = line.split("\\s+");
				HashSet<String> complex = new HashSet<String>(Arrays.asList(spl[0].split("/")));
				this.significance_sorted_complexes.add(complex);
				this.directions.put(complex, spl[1]);
				this.qvalues.put(complex, Double.parseDouble(spl[2]));
				this.fold_change.put(complex, Double.parseDouble(spl[3]));
				this.median_change.put(complex, Double.parseDouble(spl[4]));
				this.member_seed_comb.put(complex, new HashSet<String>(Arrays.asList(spl[5].split("/"))));
				this.member_seed_complexes.put(complex, Arrays.stream(spl[6].split(",")).map(s -> new HashSet<String>(Arrays.asList(s.split("/")))).collect(Collectors.toList()));
			}
		}
		
		public SignSortedComplexesResult(String in_file, double qval_threshold) {
			for (String line:Utilities.readFile(in_file)) {
				if (line.startsWith("("))
					continue;
				String[] spl = line.split("\\s+");
				double qval = Double.parseDouble(spl[2]);
				
				// stops with the last entry below the threshold
				if (qval >= qval_threshold)
					break;
				
				HashSet<String> complex = new HashSet<String>(Arrays.asList(spl[0].split("/")));
				this.significance_sorted_complexes.add(complex);
				this.directions.put(complex, spl[1]);
				this.qvalues.put(complex, qval);
				this.fold_change.put(complex, Double.parseDouble(spl[3]));
				this.median_change.put(complex, Double.parseDouble(spl[4]));
				this.member_seed_comb.put(complex, new HashSet<String>(Arrays.asList(spl[5].split("/"))));
				this.member_seed_complexes.put(complex, Arrays.stream(spl[6].split(",")).map(s -> new HashSet<String>(Arrays.asList(s.split("/")))).collect(Collectors.toList()));
			}
		}
		
		public List<HashSet<String>> getSignificanceSortedComplexes() {
			return significance_sorted_complexes;
		}

		public Map<HashSet<String>, String> getDirections() {
			return directions;
		}

		public Map<HashSet<String>, Double> getQValues() {
			return qvalues;
		}

		public Map<HashSet<String>, Double> getFoldChange() {
			return fold_change;
		}

		public Map<HashSet<String>, Double> getMedianChange() {
			return median_change;
		}

		public Map<HashSet<String>, HashSet<String>> getMemberSeedComb() {
			return member_seed_comb;
		}

		public Map<HashSet<String>, List<HashSet<String>>> getMemberSeedComplexes() {
			return member_seed_complexes;
		}
	}
}
