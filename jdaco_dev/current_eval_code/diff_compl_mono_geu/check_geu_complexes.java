package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_geu_complexes {

	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> geu_data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG"))
				geu_data.put(sample, qdr);
		}
		
		// TODO: code eval
	}
}
