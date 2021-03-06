package stem_cell_complexomes_older;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;
import stem_cell_complexomes_older.definitions;


public class check_complexes_quant {
	
	public static void main(String[] args) {
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		boolean check_supersets = false;
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("hESC"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		
		Set<String> involved_tfs = new HashSet<>();
		DiffSeedCombDetector dcd = new DiffSeedCombDetector(group1, group2, 0.01, false, false, check_supersets, 0.0, 4);
		Map<String, String> effect = new HashMap<>();
		List<String> res_pos = new LinkedList<>();
		List<HashSet<String>> tf_variants = new LinkedList<>();
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			double median_tissues = Utilities.getMedian(dcd.getGroup1Abundances().get(variant));
			double median_ESCs = Utilities.getMedian(dcd.getGroup2Abundances().get(variant));
			
			String sign = "-";
			if (median_ESCs > median_tissues)
				sign = "+";
			
			String hgncs = DataQuery.batchHGNCNamesFromProteins(variant).toString();
			double pval = dcd.getSignificantVariantsQValues().get(variant);
			
			// filter for increases abundance
			if (sign.equals("-"))
				continue;
			
			res_pos.add(sign + " " + hgncs + ", " + pval);

			involved_tfs.addAll(variant);
			tf_variants.add(variant);
			
			// determine actual complexes
			List<Set<String>> complexes = new LinkedList<>();
			for (QuantDACOResultSet qdr:group2.values())
				if (qdr.getSeedToComplexMap().containsKey(variant))
					complexes.addAll(qdr.getSeedToComplexMap().get(variant));
			effect.put(variant.toString(), definitions.goa.rateCollectionOfProteins(complexes));
		}
		
		// write results
		Utilities.writeEntries(res_pos, "/Users/tho/Desktop/res_all.txt");
		
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, involved_tfs, 0.0001, involved_tfs);
		
		System.out.println("Building regnet ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(tf_variants, involved_tfs, bdh, -30, 30, 4, 1);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only.txt");
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only_min2.txt", 2);
		
		// adding annotational data
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", effect);
		
		regnet.writeNodeTable("/Users/tho/Desktop/node_table_only.txt", annotational_data);
	}
}
