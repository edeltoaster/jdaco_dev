package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes {

	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> geu_data = new HashMap<>();
		Map<String, QuantDACOResultSet> cm_data = new HashMap<>();
		Map<String, QuantDACOResultSet> ncm_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				geu_data.put(sample, qdr);
			else if (sample.startsWith("non"))
				ncm_data.put(sample, qdr);
			else
				cm_data.put(sample, qdr);
		}
		
		System.out.println();
		
		DiffComplexDetector mono_dcd = new DiffComplexDetector(cm_data, ncm_data, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		System.out.println(mono_dcd.getSignSortedComplexes(false, false).size());
		System.out.println(mono_dcd.getSignSortedVariants(false, false).size());
		mono_dcd.writeSignSortedComplexes("mono_compl.txt", true);
		mono_dcd.writeSignSortedVariants("mono_tfcs.txt", true);
	}
}
