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


public class check_complexes_quant_pluri {
	
	public static void main(String[] args) {
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
//			if (sample.startsWith("BM_"))
//				group1.put(sample, qdr);
//			else if (sample.contains("hESC")) {
//				group2.put(sample, qdr);
//			}
			
			if (!sample.contains("hESC"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		
		Set<String> involved_tfs = new HashSet<>();
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, 0.05);
		List<HashSet<String>> pluri_tf_variants = new LinkedList<>();
		Map<String, String> effect = new HashMap<>();
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			double median_tissues = Utilities.getMedian(dcd.getGroup1Abundances().get(variant));
			double median_ESCs = Utilities.getMedian(dcd.getGroup2Abundances().get(variant));
			
			String sign = "-";
			if (median_ESCs > median_tissues)
				sign = "+";
			
			String hgncs = DataQuery.batchHGNCProteinsGenes(variant).toString();
			double pval = dcd.getSignificanceVariantsPValues().get(variant);
			
			Set<String> overlap = new HashSet<>(definitions.pluri_factors);
			overlap.retainAll(variant);
			
			// filter for increases abundance and those including pluri factors
			if (sign.equals("-") || overlap.size() == 0)
				continue;
			
			involved_tfs.addAll(variant);
			pluri_tf_variants.add(variant);
			
			// determine actual complexes
			List<Set<String>> complexes = new LinkedList<>();
			for (QuantDACOResultSet qdr:group2.values())
				if (qdr.getSeedToComplexMap().containsKey(variant))
					complexes.addAll(qdr.getSeedToComplexMap().get(variant));
			effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
			
			System.out.println(sign + " " + hgncs + ", " + pval);
		}
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, involved_tfs, 0.0001, involved_tfs);
		
		System.out.println("Building regnet ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(pluri_tf_variants, bdh, -50, 50, 4, 1);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only.txt");
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only_min2.txt", 2);
		
		// adding annotational data
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", effect);
		
		regnet.writeNodeTable("/Users/tho/Desktop/node_table_only.txt", annotational_data);
	}
}
