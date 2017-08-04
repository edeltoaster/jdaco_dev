package diff_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_BRCA_complexes {
	
	public static Map<String, String> up_name_map;
	
	public static void check_paired_diff_compl(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2) {
		String output_folder = "cdiffnet_results_99_5_-25-25/";
	
		Set<String> interesting_targets = new HashSet<>();
		
		System.out.println("Determine differential complexomes ...");
		Set<String> involved_tfs = new HashSet<>();
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		
		Map<String, String> effect = new HashMap<>();
		Map<String, String> summarized_effect = new HashMap<>();
		Map<String, String> abundances = new HashMap<>();
		Map<String, String> actual_complexes = new HashMap<>();
		Map<String, String> directions = new HashMap<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		
		System.out.println(dcd.getNumberOfTests() + " TF complexes tested.");
		System.out.println(dcd.getSignificanceSortedVariants().size() + " diff. TF complexes.");
		
		if (dcd.getSignificanceSortedVariants().size() == 0) {
			System.out.println();
			return;
		}
		
		Map<HashSet<String>, List<HashSet<String>>> tfs_to_complexes = new HashMap<>();
		for (HashSet<String> complex:dcd.getSignificanceSortedVariants()) {
			
			String sign = dcd.getSignificantVariantsDirections().get(complex);
			HashSet<String> tfs = new HashSet<>(complex);
			tfs.retainAll(BRCA_definitions.seed);
			
			if (!tfs_to_complexes.containsKey(tfs))
				tfs_to_complexes.put(tfs, new LinkedList<HashSet<String>>());
			tfs_to_complexes.get(tfs).add(complex);
			
			String compl_string = complex.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			String tfs_string = tfs.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			double pval = dcd.getSignificantVariantsQValues().get(complex);
			involved_tfs.addAll(tfs);
			
			String out_string = sign + " " + tfs_string + " : " + compl_string + " -> " + String.format(Locale.US, "%.4g", pval);
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(out_string);
			} else {
				// everything in network
				res_pos_all.add(out_string);
			}
			
			directions.put(tfs.toString(), sign);
		}
		
		for (HashSet<String> tf_comb:tfs_to_complexes.keySet()) {
			actual_complexes.put(tf_comb.toString(), String.join(",", tfs_to_complexes.get(tf_comb).stream().map(c -> String.join("/", c.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()))).collect(Collectors.toList())));
		}
		
		System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "-");
		
		new File(output_folder).mkdir();
		
		// write results
		Utilities.writeEntries(res_pos_all, output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, output_folder + "res_neg_all.txt");
		
		interesting_targets.addAll(involved_tfs);
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(BRCA_definitions.binding_data, involved_tfs, 0.0001, interesting_targets);
		
		/**
		 *  writing network data
		 */
		
		System.out.println("Building ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(tfs_to_complexes.keySet(), bdh, BRCA_definitions.d_min, BRCA_definitions.d_max, BRCA_definitions.no_threads, 1);
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(output_folder + "regnet.txt");
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Epi_effect", summarized_effect);
		annotational_data.put("Epi_effect_details", effect);
		annotational_data.put("Actual_complexes", actual_complexes);
		annotational_data.put("Mean_abundances", abundances);
		annotational_data.put("Direction", directions);
		regnet.writeNodeTable(output_folder + "nodetable.txt", annotational_data);
		// pruning
		regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(output_folder + "regnet_pruned.txt");
		regnet.writeNodeTable(output_folder + "nodetable_pruned.txt", annotational_data);
	}

	public static void main(String[] args) {
		BRCA_definitions.printParameters();
		
		BRCA_definitions.goa.printTagInformation();
		
		System.out.println();
		new File(BRCA_definitions.diff_compl_output_folder).mkdir();
		
		up_name_map = DataQuery.getUniprotToGeneNameMap(BRCA_definitions.seed);
		
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(BRCA_definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String sample_prefix = sample.split("_")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), BRCA_definitions.seed, BRCA_definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("normal"))
				group1.put(sample_prefix, qdr);
			else {
				group2.put(sample_prefix, qdr);
			}
		}
		System.out.println("normal samples : " + group1.size());
		System.out.println("tumor samples : " + group2.size());
		
		check_paired_diff_compl(group1, group2);
	}
}
