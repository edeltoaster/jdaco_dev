package eval_biofilm;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.ConstructedNetworks;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.Utilities;

public class builtSalivaNetworks {
	static String project_folder = "/Users/tho/Dropbox/Work/biofilm/";
	
	public static void main(String[] args) {
		
		Map<String, Set<String>> data_prot_map = new HashMap<String, Set<String>>();
		
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders(project_folder + "proteins/", "speichel.txt")) {
			String filename = file.getName();
			String sample = filename.split("_")[0];
			System.out.println(file.getAbsolutePath());
		}
		
		// initialize stuff
		PPIN ppi = new PPIN("mixed_data/human_ppi.tsv.gz");
		
		NetworkBuilder builder = new NetworkBuilder(ppi);
		
		ConstructedNetworks constr = builder.constructAssociatedNetworksFromProteinSet(null);
	}

}
