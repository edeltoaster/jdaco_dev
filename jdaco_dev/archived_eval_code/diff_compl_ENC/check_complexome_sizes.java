package diff_compl_ENC;


import java.io.File;

import framework.DACOResultSet;
import framework.Utilities;


public class check_complexome_sizes {
	
	public static void main(String[] args) {
		pluri_definitions.printParameters();
		
		double n = 0;
		double total = 0;
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(pluri_definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			DACOResultSet dr = new DACOResultSet(f.getAbsolutePath(), pluri_definitions.seed_file);
			System.out.println(sample + " : " + dr.getResult().size() + " complexes.");
			n += dr.getResult().size();
			total += 1;
		}
		
		System.out.println();
		System.out.println("average: " + n/total + " per sample.");
	}
}
