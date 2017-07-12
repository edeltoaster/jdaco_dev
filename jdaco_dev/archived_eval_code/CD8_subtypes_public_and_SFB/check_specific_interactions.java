package CD8_subtypes_public_and_SFB;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.PPIN;
import framework.StrPair;
import framework.Utilities;

public class check_specific_interactions {
	
	static String CD8_networks_folder = "/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/CD8_networks_0.13/";//"/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/CD8_networks_0.0/"; // 0.13
	static String hemato_networks_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/0.31/";//"/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/quant_hemo_networks_0.0/";
	
	public static boolean containsIA(String ppin_path, StrPair IA) {
		PPIN ppin = new PPIN(ppin_path);
		return ppin.getInteractionsFast().contains(IA);
	}
	
	public static void main(String[] args) {
		Set<StrPair> IA_of_interest = new HashSet<>();
		IA_of_interest.add(new StrPair("P17947", "P15976")); // PU.1 / GATA-1
		IA_of_interest.add(new StrPair("Q13422", "P50750")); // Ik / Cdk9 
		IA_of_interest.add(new StrPair("Q13422", "P15976")); // Ik / GATA-1
		IA_of_interest.add(new StrPair("Q13422", "P23769")); // Ik / GATA-2
		IA_of_interest.add(new StrPair("Q13422", "P23771")); // Ik / GATA-3
		IA_of_interest.add(new StrPair("P29590", "P06400")); // PML / retinoblastoma pRb 105
		IA_of_interest.add(new StrPair("P17542", "P25791")); // SCL / LMO2
		
		for (StrPair pair:IA_of_interest) {
			List<String> list = new ArrayList<>(2);
			list.add(pair.getL());
			list.add(pair.getR());
			System.out.println("Checking " + DataQuery.batchHGNCNamesFromProteins(list));
			
			// check reference
			if (!containsIA("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/human_mentha_8_jul.txt.gz", pair)) { // strictly speaking more than the old hematopoiesis data
				System.out.println("not in reference");
				continue;
			}
			
			Map<String, Integer> cell_type_count = new HashMap<>();
			Map<String, Integer> cell_type_overall = new HashMap<>();
			
			for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(hemato_networks_folder, "_ppin.txt.gz")) {
				String cell_type = f.getName().split("_")[0];
				cell_type = f.getParentFile().getName(); // for old data
				
				if (!cell_type_count.containsKey(cell_type)) {
					cell_type_count.put(cell_type, 0);
					cell_type_overall.put(cell_type, 0);
				}
				
				if (containsIA(f.getAbsolutePath(), pair))
					cell_type_count.put(cell_type, cell_type_count.get(cell_type) + 1);
				cell_type_overall.put(cell_type, cell_type_overall.get(cell_type) + 1);
				
			}
			
			
			for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(CD8_networks_folder, "_ppin.txt.gz")) {
				String cell_type = f.getName().split("_")[1];
				
				if (!cell_type_count.containsKey(cell_type)) {
					cell_type_count.put(cell_type, 0);
					cell_type_overall.put(cell_type, 0);
				}
				
				if (containsIA(f.getAbsolutePath(), pair))
					cell_type_count.put(cell_type, cell_type_count.get(cell_type) + 1);
				cell_type_overall.put(cell_type, cell_type_overall.get(cell_type) + 1);
				
			}
			
			// output for each cell_type
			for (String ct:cell_type_count.keySet()) {
				System.out.println(ct + ":" + cell_type_count.get(ct) + " / " + cell_type_overall.get(ct));
			}
			
			System.out.println();
		}
		
	}
}
