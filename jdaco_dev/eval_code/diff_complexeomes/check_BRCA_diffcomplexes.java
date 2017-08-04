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
import framework.DiffSeedVarDetector;
import framework.DiffSeedVarDetector.SPEnrichment;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_BRCA_diffcomplexes {
	
	public static Map<String, String> up_name_map;
	
	public static void check_paired_diff_compl(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2) {
		String output_folder = BRCA_definitions.diff_compl_output_folder;
	
		Set<String> interesting_targets = new HashSet<>();
		
		System.out.println("Determine differential complexomes ...");
		Set<String> involved_tfs = new HashSet<>();
		DiffSeedVarDetector dcd = new DiffSeedVarDetector(group1, group2, BRCA_definitions.qvalue, BRCA_definitions.parametric, BRCA_definitions.paired, BRCA_definitions.check_supersets, BRCA_definitions.min_variant_fraction, BRCA_definitions.no_threads);
		
		List<HashSet<String>> tf_variants = new LinkedList<>();
		Map<String, String> effect = new HashMap<>();
		Map<String, String> summarized_effect = new HashMap<>();
		Map<String, String> abundances = new HashMap<>();
		Map<String, String> actual_complexes = new HashMap<>();
		Map<String, String> directions = new HashMap<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		
		System.out.println(dcd.getNumberOfTests() + " TF variants tested.");
		System.out.println(dcd.getSignificanceSortedVariants().size() + " diff. TF variants.");
		
		if (dcd.getSignificanceSortedVariants().size() == 0) {
			System.out.println();
			return;
		}
		
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			
			String sign = dcd.getSignificantVariantsDirections().get(variant);
			
			String tf_comb = variant.stream().map(p -> up_name_map.getOrDefault(p, p)).collect(Collectors.toList()).toString();
			double pval = dcd.getSignificantVariantsQValues().get(variant);
			involved_tfs.addAll(variant);
			String out_string = sign + " " + tf_comb + " -> " + String.format(Locale.US, "%.4g", pval);
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(out_string);
			} else {
				// everything in network
				res_pos_all.add(out_string);
			}
			
			tf_variants.add(variant);
			directions.put(variant.toString(), sign);
			
			String[] annotation_data = DiffSeedVarDetector.getSortedComplexesAnnotations(variant, sign, BRCA_definitions.goa, group1, group2);
			effect.put(variant.toString(), annotation_data[1]);
			summarized_effect.put(variant.toString(), annotation_data[2]);
			actual_complexes.put(variant.toString(), annotation_data[0]);
			abundances.put(variant.toString(), annotation_data[3]);
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
		RegulatoryNetwork pluri_regnet = new RegulatoryNetwork(tf_variants, bdh, BRCA_definitions.d_min, BRCA_definitions.d_max, BRCA_definitions.no_threads, 1);
		System.out.println(pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(output_folder + "regnet.txt");
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Epi_effect", summarized_effect);
		annotational_data.put("Epi_effect_details", effect);
		annotational_data.put("Actual_complexes", actual_complexes);
		annotational_data.put("Mean_abundances", abundances);
		annotational_data.put("Direction", directions);
		pluri_regnet.writeNodeTable(output_folder + "nodetable.txt", annotational_data);
		// pruning
		pluri_regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(output_folder + "regnet_pruned.txt");
		pluri_regnet.writeNodeTable(output_folder + "nodetable_pruned.txt", annotational_data);
		
		/**
		 * Playing around with seed protein enrichment
		 */
		
		System.out.println("Calculating TF enrichment ...");
		SPEnrichment tf_enrich = dcd.calculateTFEnrichment(BRCA_definitions.qvalue, BRCA_definitions.SPEnrich_iterations, BRCA_definitions.SPEnrich_compl_part_threshold);
		List<String> pos_tf_enrich_out = new LinkedList<>();
		List<String> neg_tf_enrich_out = new LinkedList<>();
		for (String tf:tf_enrich.getSignificanceSortedSeedProteins()) {
			String dir = tf_enrich.getSignificantSeedProteinDirections().get(tf);
			String out = tf + " " + up_name_map.getOrDefault(tf, tf) + " " + tf_enrich.getSignificantSeedProteinQvalues().get(tf);
			
			if (dir.equals("+"))
				pos_tf_enrich_out.add(out);
			else
				neg_tf_enrich_out.add(out);
		}
		Utilities.writeEntries(pos_tf_enrich_out, output_folder + "tf_enrich_pos.txt");
		Utilities.writeEntries(neg_tf_enrich_out, output_folder + "tf_enrich_neg.txt");
		System.out.println();
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
