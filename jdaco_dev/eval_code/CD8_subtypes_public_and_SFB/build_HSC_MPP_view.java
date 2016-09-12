package CD8_subtypes_public_and_SFB;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import framework.PPIN;
import framework.StrPair;
import framework.Utilities;

public class build_HSC_MPP_view {
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/0.31/";
	
	public static List<PPIN> loadNetworks(String path) {
		List<PPIN> ppins = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(path, "_ppin.txt.gz")) {
			ppins.add( new PPIN(f.getAbsolutePath()) );
		}
		return ppins;
	}
	
	public static void main(String[] args) {
		PPIN original_ppin = new PPIN("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/human_mentha_17_jan.txt.gz");
		List<PPIN> HSCs = loadNetworks(network_folder + "HSC/");
		List<PPIN> MPPs = loadNetworks(network_folder + "MPP/");
		
		System.out.println(HSCs.size() + " -> " + MPPs.size());
		
//		List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
//		Map<String, String> up_hgnc = new HashMap<>();
//		for (String[] s:HGNC_map) {
//			up_hgnc.put(s[1], s[0]);
//		}
		
		List<String> to_write = new LinkedList<>();
		to_write.add("Protein1 Protein2 HSC_percentage MPP_percentage Membership");
		for (StrPair pair:original_ppin.getInteractionsFast()) {
			String p1 = pair.getL();
			String p2 = pair.getR();
			
			double HSCc = 0;
			for (PPIN ppin:HSCs)
				if (ppin.getInteractionsFast().contains(pair))
					HSCc += 1;
			HSCc /= HSCs.size();
			
			double MPPc = 0;
			for (PPIN ppin:MPPs)
				if (ppin.getInteractionsFast().contains(pair))
					MPPc += 1;
			
			MPPc /= MPPs.size();
			
			String color = "none";
			if (HSCc > 0.5 && MPPc > 0.5) 
				color = "HSCs/MPPs";
			else if (HSCc > 0.5)
				color = "HSCs";
			else if (MPPc > 0.5)
				color = "MPPs";
			
			to_write.add(p1 + " " + p2 + " " + HSCc + " " + MPPc + " " + color);
		}
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/ref_HSC_MPP.txt");
		
	}

}
