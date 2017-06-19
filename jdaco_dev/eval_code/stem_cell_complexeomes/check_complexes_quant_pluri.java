package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;
import stem_cell_complexeomes.definitions;


public class check_complexes_quant_pluri {
	
	public static void main(String[] args) {
		definitions.printParameters();
		
		System.out.println("Reading data ...");
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			//qdr.removeOpposinglyAnnotatedComplexes(definitions.goa);
			
			if (!sample.contains("H1-hESC") && !sample.contains("H7-hESC") && !sample.contains("induced-pluripotent-stem-cell"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("non-stem : " + group1.size());
		System.out.println("stem samples : " + group2.size());
		
		System.out.println("Determine differential complexomes ...");
		Set<String> involved_tfs = new HashSet<>();
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, definitions.pvalue, definitions.no_threads, definitions.check_supersets);
		
		List<HashSet<String>> pluri_tf_variants = new LinkedList<>();
		Map<String, String> pluri_effect = new HashMap<>();
		List<HashSet<String>> plurisub_tf_variants = new LinkedList<>();
		Map<String, String> plurisub_effect = new HashMap<>();
		List<HashSet<String>> nonpluri_tf_variants = new LinkedList<>();
		Map<String, String> nonpluri_effect = new HashMap<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_pos_pluri = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			double median_tissues = Utilities.getMedian(dcd.getGroup1Abundances().get(variant));
			double median_ESCs = Utilities.getMedian(dcd.getGroup2Abundances().get(variant));
			
			String sign = "-";
			if (median_ESCs > median_tissues)
				sign = "+";
			
			String hgncs = DataQuery.batchHGNCProteinsGenes(variant).toString();
			double pval = dcd.getSignificanceVariantsPValues().get(variant);
			involved_tfs.addAll(variant);
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(sign + " " + hgncs + ", " + pval);
				
				nonpluri_tf_variants.add(variant);
				
				// determine actual complexes
				List<Set<String>> complexes = new LinkedList<>();
				for (QuantDACOResultSet qdr:group1.values()) // get examples from group1 since they are underrepresented in group2
					if (qdr.getSeedToComplexMap().containsKey(variant))
						complexes.addAll(qdr.getSeedToComplexMap().get(variant));
				nonpluri_effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
			} else {
				// everything in pluri-network
				res_pos_all.add(sign + " " + hgncs + ", " + pval);
				
				pluri_tf_variants.add(variant);
				
				// determine actual complexes
				List<Set<String>> complexes = new LinkedList<>();
				for (QuantDACOResultSet qdr:group2.values())
					if (qdr.getSeedToComplexMap().containsKey(variant))
						complexes.addAll(qdr.getSeedToComplexMap().get(variant));
				pluri_effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
				
				
				// filter for those including pluri factors
				Set<String> overlap = new HashSet<>(definitions.pluri_factors);
				overlap.retainAll(variant);
				if (overlap.size() == 0)
					continue;
				
				res_pos_pluri.add(sign + " " + hgncs + ", " + pval);
				
				plurisub_tf_variants.add(variant);
				
				// determine actual complexes
				complexes = new LinkedList<>();
				for (QuantDACOResultSet qdr:group2.values())
					if (qdr.getSeedToComplexMap().containsKey(variant))
						complexes.addAll(qdr.getSeedToComplexMap().get(variant));
				plurisub_effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
			}
		}
		
		// write results
		Utilities.writeEntries(res_pos_all, "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, "res_neg_all.txt");
		Utilities.writeEntries(res_pos_pluri, "res_pos_pluri.txt");
		
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, involved_tfs, 0.0001, involved_tfs);
		
		/**
		 *  writing pluri network data
		 */
		
		System.out.println("Building pluri sub-regnet ...");
		RegulatoryNetwork plurisub_regnet = new RegulatoryNetwork(plurisub_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		plurisub_regnet.writeRegulatoryNetwork("plurisub_regnet_only.txt");
		plurisub_regnet.writeRegulatoryNetwork("plurisub_regnet_only_min2.txt", 2);
		
		// adding annotational data
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", plurisub_effect);
		plurisub_regnet.writeNodeTable("plurisub_node_table_only.txt", annotational_data);
		
		
		System.out.println("Building pluri regnet ...");
		RegulatoryNetwork pluri_regnet = new RegulatoryNetwork(pluri_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		pluri_regnet.writeRegulatoryNetwork("pluri_regnet_only.txt");
		pluri_regnet.writeRegulatoryNetwork("pluri_regnet_only_min2.txt", 2);
		
		// adding annotational data
		annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", plurisub_effect);
		pluri_regnet.writeNodeTable("pluri_node_table_only.txt", annotational_data);
		
		/**
		 *  writing non-pluri network data
		 */
		
		System.out.println("Building non-pluri regnet ...");
		RegulatoryNetwork nonpluri_regnet = new RegulatoryNetwork(nonpluri_tf_variants, bdh, definitions.d_min, definitions.d_max, definitions.no_threads, 1);
		nonpluri_regnet.writeRegulatoryNetwork("nonpluri_regnet_only.txt");
		nonpluri_regnet.writeRegulatoryNetwork("nonpluri_regnet_only_min2.txt", 2);
		
		// also adding annotational data here
		annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", nonpluri_effect);
		nonpluri_regnet.writeNodeTable("nonpluri_node_table_only.txt", annotational_data);
	}
}
