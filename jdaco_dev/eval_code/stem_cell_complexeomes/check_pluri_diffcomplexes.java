package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPEnrichment;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;
import stem_cell_complexeomes.definitions;


public class check_pluri_diffcomplexes {
	
	public static void main(String[] args) {
		definitions.printParameters();
		
		System.out.println("Reading data ...");
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("H1-hESC") && !sample.contains("H7-hESC") && !sample.contains("induced-pluripotent-stem-cell"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("non-stem : " + group1.size());
		System.out.println("stem samples : " + group2.size());
		
		Set<String> allosome_proteins = DataQuery.getAllosomeProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		
		System.out.println("Determine differential complexomes ...");
		Set<String> involved_tfs = new HashSet<>();
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.qvalue, definitions.parametric, definitions.check_supersets, definitions.no_threads);
		
		List<HashSet<String>> pluri_tf_variants = new LinkedList<>();
		Map<String, String> pluri_effect = new HashMap<>();
		List<HashSet<String>> plurisub_tf_variants = new LinkedList<>();
		Map<String, String> plurisub_effect = new HashMap<>();
		Map<String, String> plurisub_actual_complexes = new HashMap<>();
		Map<String, String> pluri_actual_complexes = new HashMap<>();
		List<HashSet<String>> nonpluri_tf_variants = new LinkedList<>();
		Map<String, String> nonpluri_effect = new HashMap<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_pos_pluri = new LinkedList<>();
		List<String> res_pos_all_noallo = new LinkedList<>();
		List<String> res_pos_pluri_noallo = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			
			String sign = dcd.getSignificantVariantsDirections().get(variant);
			
			String hgncs = DataQuery.batchHGNCNamesFromProteins(variant).toString();
			double pval = dcd.getSignificantVariantsQValues().get(variant);
			involved_tfs.addAll(variant);
			String out_string = sign + " " + hgncs + ", " + pval;
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(out_string);
				nonpluri_tf_variants.add(variant);
				
				// determine actual complexes
				List<Set<String>> complexes = new LinkedList<>();
				for (QuantDACOResultSet qdr:group1.values()) // get examples from group1 since they are underrepresented in group2
					if (qdr.getSeedToComplexMap().containsKey(variant))
						complexes.addAll(qdr.getSeedToComplexMap().get(variant));
				nonpluri_effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
			} else {
				// everything in pluri-network
				res_pos_all.add(out_string);
				pluri_tf_variants.add(variant);
				
				if (variant.stream().noneMatch(p -> allosome_proteins.contains(p)))
					res_pos_all_noallo.add(out_string);
				
				// determine actual complexes
				List<Set<String>> complexes = new LinkedList<>();
				for (QuantDACOResultSet qdr:group2.values())
					if (qdr.getSeedToComplexMap().containsKey(variant))
						complexes.addAll(qdr.getSeedToComplexMap().get(variant));
				pluri_effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
				String complexes_string = String.join(",", complexes.stream().map(c -> String.join("/", DataQuery.batchHGNCNamesFromProteins(c))).collect(Collectors.toList()));
				pluri_actual_complexes.put(variant.toString(), complexes_string);
				
				// filter for those including pluri factors
				Set<String> overlap = new HashSet<>(definitions.pluri_factors);
				overlap.retainAll(variant);
				if (overlap.size() == 0)
					continue;
				
				res_pos_pluri.add(out_string);
				plurisub_tf_variants.add(variant);
				
				if (variant.stream().noneMatch(p -> allosome_proteins.contains(p)))
					res_pos_pluri_noallo.add(out_string);
				
				// determine actual complexes
				complexes = new LinkedList<>();
				for (QuantDACOResultSet qdr:group2.values())
					if (qdr.getSeedToComplexMap().containsKey(variant))
						complexes.addAll(qdr.getSeedToComplexMap().get(variant));
				plurisub_effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
				complexes_string = String.join(",", complexes.stream().map(c -> String.join("/", DataQuery.batchHGNCNamesFromProteins(c))).collect(Collectors.toList()));
				plurisub_actual_complexes.put(variant.toString(), complexes_string);
			}
		}
		
		new File(definitions.diff_compl_output_folder).mkdir();
		
		// write results
		Utilities.writeEntries(res_pos_all, definitions.diff_compl_output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, definitions.diff_compl_output_folder + "res_neg_all.txt");
		Utilities.writeEntries(res_pos_pluri, definitions.diff_compl_output_folder + "res_pos_pluri.txt");
		Utilities.writeEntries(res_pos_all_noallo, definitions.diff_compl_output_folder + "res_pos_all_noallo.txt");
		Utilities.writeEntries(res_pos_pluri_noallo, definitions.diff_compl_output_folder + "res_pos_pluri_noallo.txt");
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, involved_tfs, 0.0001, involved_tfs);
		
		/**
		 *  writing pluri network data
		 */
		
		System.out.println("Building pluri sub-regnet ...");
		RegulatoryNetwork plurisub_regnet = new RegulatoryNetwork(plurisub_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		System.out.println(plurisub_regnet.getSizesStr());
		plurisub_regnet.writeRegulatoryNetwork(definitions.diff_compl_output_folder + "plurisub_regnet.txt");
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", plurisub_effect);
		annotational_data.put("Actual_complexes", plurisub_actual_complexes);
		plurisub_regnet.writeNodeTable(definitions.diff_compl_output_folder + "plurisub_nodetable.txt", annotational_data);
		// pruning
		plurisub_regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + plurisub_regnet.getSizesStr());
		plurisub_regnet.removeProteinSet(allosome_proteins);
		System.out.println("allo: " + plurisub_regnet.getSizesStr());
		plurisub_regnet.writeRegulatoryNetwork(definitions.diff_compl_output_folder + "plurisub_regnet_pruned.txt");
		plurisub_regnet.writeNodeTable(definitions.diff_compl_output_folder + "plurisub_nodetable_pruned.txt", annotational_data);
		// more very concise output
		List<String> res_pos_pluri_network = new LinkedList<>();
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			if (!plurisub_regnet.getComplexToTargets().keySet().contains(variant))
				continue;
			String sign = dcd.getSignificantVariantsDirections().get(variant);
			String hgncs = DataQuery.batchHGNCNamesFromProteins(variant).toString();
			double pval = dcd.getSignificantVariantsQValues().get(variant);
			involved_tfs.addAll(variant);
			String out_string = sign + " " + hgncs + ", " + pval;
			res_pos_pluri_network.add(out_string);
		}
		Utilities.writeEntries(res_pos_pluri_network, definitions.diff_compl_output_folder + "res_pos_pluri_network.txt");
		
		System.out.println("Building pluri regnet ...");
		RegulatoryNetwork pluri_regnet = new RegulatoryNetwork(pluri_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		System.out.println(pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(definitions.diff_compl_output_folder + "pluri_regnet.txt");
		annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", pluri_effect);
		annotational_data.put("Actual_complexes", pluri_actual_complexes);
		pluri_regnet.writeNodeTable(definitions.diff_compl_output_folder + "pluri_nodetable.txt", annotational_data);
		// pruning
		pluri_regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + pluri_regnet.getSizesStr());
		pluri_regnet.removeProteinSet(allosome_proteins);
		System.out.println("allo: " + pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(definitions.diff_compl_output_folder + "pluri_regnet_pruned.txt");
		pluri_regnet.writeNodeTable(definitions.diff_compl_output_folder + "pluri_nodetable_pruned.txt", annotational_data);
		
		/**
		 *  writing non-pluri network data
		 */
		
		System.out.println("Building non-pluri regnet ...");
		RegulatoryNetwork nonpluri_regnet = new RegulatoryNetwork(nonpluri_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		System.out.println(nonpluri_regnet.getSizesStr());
		nonpluri_regnet.writeRegulatoryNetwork(definitions.diff_compl_output_folder + "nonpluri_regnet.txt");
		annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", nonpluri_effect);
		pluri_regnet.writeNodeTable(definitions.diff_compl_output_folder + "pluri_nodetable.txt", annotational_data);
		// pruning
		nonpluri_regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + nonpluri_regnet.getSizesStr());
		nonpluri_regnet.removeProteinSet(allosome_proteins);
		System.out.println("allo: " + nonpluri_regnet.getSizesStr());
		nonpluri_regnet.writeRegulatoryNetwork(definitions.diff_compl_output_folder + "nonpluri_regnet_pruned.txt");
		nonpluri_regnet.writeNodeTable(definitions.diff_compl_output_folder + "nonpluri_nodetable_pruned.txt", annotational_data);
		
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
		Utilities.writeEntries(pos_tf_enrich_out, definitions.diff_compl_output_folder + "tf_enrich_pos.txt");
		Utilities.writeEntries(neg_tf_enrich_out, definitions.diff_compl_output_folder + "tf_enrich_neg.txt");
	}
}
