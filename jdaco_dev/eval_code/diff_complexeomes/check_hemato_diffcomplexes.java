package diff_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import diff_complexeomes.definitions;
import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPEnrichment;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_hemato_diffcomplexes {
	
	public static String daco_results_folder = "res_99_5/";
	public static String networks_folder = "BP_networks/"; // intended to be run on the servers
	public static String diff_compl_output_folder = "BP_out/";
	
	public static void main(String[] args) {
		definitions.printParameters();
		System.out.println("folders overwritten for BLUEPRINT stuff");
		System.out.println("parametric");
		
		System.out.println("Reading data ...");
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("CLP")) // pairwise or more?
				group1.put(sample, qdr);
			else if (sample.contains("CD4")) {
				group2.put(sample, qdr);
			}
		}
		System.out.println("other hemato-samples : " + group1.size());
		System.out.println("CD4 samples : " + group2.size());
		
		Set<String> allosome_proteins = DataQuery.getAllosomeProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		
		System.out.println("Determine differential complexomes ...");
		Set<String> involved_tfs = new HashSet<>();
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, true, definitions.check_supersets, definitions.no_threads);
		
		List<HashSet<String>> cd4_tf_variants = new LinkedList<>();
		Map<String, String> cd4_effect = new HashMap<>();
		Map<String, String> cd4_summarized_effect = new HashMap<>();
		Map<String, String> cd4_abundances = new HashMap<>();
		Map<String, String> cd4_actual_complexes = new HashMap<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_pos_all_noallo = new LinkedList<>();
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			
			String sign = dcd.getSignificantVariantsDirections().get(variant);
			
			String hgncs = DataQuery.batchHGNCNamesFromProteins(variant).toString();
			double pval = dcd.getSignificantVariantsQValues().get(variant);
			involved_tfs.addAll(variant);
			String out_string = sign + " " + hgncs + ", " + pval;
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				continue;
			} else {
				// everything in network
				res_pos_all.add(out_string);
				cd4_tf_variants.add(variant);
				
				if (variant.stream().noneMatch(p -> allosome_proteins.contains(p)))
					res_pos_all_noallo.add(out_string);
				
				String[] annotation_data = DiffComplexDetector.getSortedComplexesAnnotations(variant, sign, definitions.goa, group1, group2);
				cd4_effect.put(variant.toString(), annotation_data[1]);
				cd4_summarized_effect.put(variant.toString(), annotation_data[2]);
				cd4_actual_complexes.put(variant.toString(), annotation_data[0]);
				cd4_abundances.put(variant.toString(), annotation_data[3]);
			}
		}
		
		new File(diff_compl_output_folder).mkdir();
		
		// write results
		Utilities.writeEntries(res_pos_all, diff_compl_output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_pos_all_noallo, diff_compl_output_folder + "res_pos_all_noallo.txt");
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, involved_tfs, 0.0001, involved_tfs);
		
		/**
		 *  writing pluri network data
		 */
		
		System.out.println("Building CD4 regnet ...");
		RegulatoryNetwork pluri_regnet = new RegulatoryNetwork(cd4_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		System.out.println(pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(diff_compl_output_folder + "cd4_regnet.txt");
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Epi_effect", cd4_summarized_effect);
		annotational_data.put("Epi_effect_details", cd4_effect);
		annotational_data.put("Actual_complexes", cd4_actual_complexes);
		annotational_data.put("Mean_abundances", cd4_abundances);
		pluri_regnet.writeNodeTable(diff_compl_output_folder + "cd4_nodetable.txt", annotational_data);
		// pruning
		pluri_regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + pluri_regnet.getSizesStr());
		pluri_regnet.removeProteinSet(allosome_proteins);
		System.out.println("allo: " + pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(diff_compl_output_folder + "cd4_regnet_pruned.txt");
		pluri_regnet.writeNodeTable(diff_compl_output_folder + "cd4_nodetable_pruned.txt", annotational_data);
		
		/**
		 * Playing around with seed protein enrichment
		 */
		
		System.out.println("Calculating TF enrichment ...");
		SPEnrichment tf_enrich = dcd.calculateTFEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		List<String> pos_tf_enrich_out = new LinkedList<>();
		List<String> neg_tf_enrich_out = new LinkedList<>();
		for (String tf:tf_enrich.getSignificanceSortedSeedProteins()) {
			String dir = tf_enrich.getSignificantSeedProteinDirections().get(tf);
			String out = tf + " " + DataQuery.getHGNCNameFromProtein(tf) + " " + tf_enrich.getSignificantSeedProteinQvalues().get(tf);
			
			if (dir.equals("+"))
				pos_tf_enrich_out.add(out);
			else
				neg_tf_enrich_out.add(out);
		}
		Utilities.writeEntries(pos_tf_enrich_out, diff_compl_output_folder + "tf_enrich_pos.txt");
		Utilities.writeEntries(neg_tf_enrich_out, diff_compl_output_folder + "tf_enrich_neg.txt");
	}
}
