package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPIhc_TPMgene/res5/";
	static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/PrePPIhc_TPMgene_networks/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
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
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			double median_tissues = Utilities.getMedian(dcd.getGroup1Abundances().get(variant));
			double median_ESCs = Utilities.getMedian(dcd.getGroup2Abundances().get(variant));
			involved_tfs.addAll(variant);
			
			String sign = "-";
			if (median_ESCs > median_tissues)
				sign = "+";
			
			String hgncs = DataQuery.batchHGNCProteinsGenes(variant).toString();
			double pval = dcd.getSignificanceVariantsPValues().get(variant);
			
			if (sign.equals("-"))
				continue;
			
			System.out.println(sign + " " + hgncs + ", " + pval);
		}
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10/hocomoco_v10_EPD_2.5k.txt.gz", involved_tfs, 0.0001, involved_tfs);
		System.out.println("Building regnet ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(dcd.getSignificanceSortedVariants(), bdh, -50, 50, 4, 1);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only.txt");
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only_min2.txt", 2);
		regnet.writeHumanNodeTable("/Users/tho/Desktop/node_table_only.txt");
	}
}
