package breastcancer_proteome_complexome;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.QuantDACOResultSet;
import framework.Utilities;

public class diff_complexes_rel_mentha {
	static String net_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/mentha_networks/";
	static String compl_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/res/";
	static String seed_file = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/JDACO_run_stuff/hocomoco_human_core_TFs_v11.txt.gz";
	static String qr_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/relquant_compl/";
	
	static Map<String, QuantDACOResultSet> data = new HashMap<>();
	
	public static void load_data(boolean write) {
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(compl_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			System.out.println("Processing " + sample);
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed_file, net_folder + sample + "_major-transcripts_rel.txt.gz");
			
			if (write)
				qdr.writeQuantifiedResult(qr_folder + sample + ".txt.gz");
			
			data.put(sample, qdr);
		}
	}
	
	public static void main(String[] args) {
		load_data(true);
	}
}
