package stem_cell_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPI_95_95_TPM/res5/";
	static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/PrePPI_TPM_networks/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void main(String[] args) {
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("BM"))
				group1.put(sample, qdr);
			else if (sample.contains("hESC")) {
				group2.put(sample, qdr);
			}
		}
		
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, 0.05);
		dcd.printResults();
		String sox2 = "P48431";
		String oct4 = "Q01860";
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		System.out.println(dcd.getGroup2().keySet());
		for (HashSet<String> combs:dcd.getGroup2Abundances().keySet()) {
			if (combs.contains(sox2) && combs.contains(oct4)) {
				double pm = mwu.mannWhitneyUTest(Utilities.getDoubleArray(dcd.getGroup2Abundances().get(combs)), Utilities.getDoubleArray(dcd.getGroup1Abundances().get(combs)));
				System.out.println(combs + " " + pm + "/" + dcd.getSignificanceVariantsPValues().get(combs) + ":"+ dcd.getGroup2Abundances().get(combs) + " vs " + dcd.getGroup1Abundances().get(combs));
			}
		}
	}
}
