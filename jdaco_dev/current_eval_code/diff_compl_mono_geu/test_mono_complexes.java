package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.QuantDACOResultSet;
import framework.Utilities;


public class test_mono_complexes {

	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> cm_data = new HashMap<>();
		Map<String, QuantDACOResultSet> ncm_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				continue;
			else if (sample.startsWith("non"))
				ncm_data.put(sample, qdr);
			else
				cm_data.put(sample, qdr);
		}
		
		
		/**
		 * Precompute quantified complexes
		 */
		
		System.out.println("Precompute quantified complexes ...");
		new File(definitions.qr_output_folder).mkdir();
		cm_data.keySet().parallelStream().forEach(s -> cm_data.get(s).writeQuantifiedResult(definitions.qr_output_folder + s + ".txt.gz"));
		ncm_data.keySet().parallelStream().forEach(s -> ncm_data.get(s).writeQuantifiedResult(definitions.qr_output_folder + s + ".txt.gz"));
	
		// TODO: implement test
		
	}
}
