package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPI_95_wTPM/res5/";
	static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/PrePPI_TPMgene_networks/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("BM_"))
				group1.put(sample, qdr);
			else if (sample.contains("hESC")) {
				group2.put(sample, qdr);
			}
		}
		
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, 0.05);
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			double median_tissues = Utilities.getMedian(dcd.getGroup1Abundances().get(variant));
			double median_ESCs = Utilities.getMedian(dcd.getGroup2Abundances().get(variant));
			
			String sign = "-";
			if (median_ESCs > median_tissues)
				sign = "+";
			
			String hgncs = DataQuery.batchHGNCProteinsGenes(variant).toString();
			double pval = dcd.getSignificanceVariantsPValues().get(variant);
			
			if (sign.equals("-"))
				continue;
			
			System.out.println(sign + " " + hgncs + ", " + pval);
		}
	}
}
