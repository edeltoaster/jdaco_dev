package old_mixed;

import java.io.File;
import framework.PPIN;
import framework.Utilities;

public class compTissueSTRING {

	public static void computeForNetwork(String input_folder, String results_out_folder) {
		
		// make output folder
		new File(results_out_folder).mkdir();
		new File(results_out_folder+"networks").mkdir();
		
		System.out.println("Reading and processing samples ...");
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders(input_folder, ".tsv")) {
			
			String sample_path = file.getAbsolutePath();
			String file_name = file.getName();
			
			System.out.println("Processing "+ file_name);
			PPIN ppi = new PPIN(sample_path);
			PPIN STRING_network = ppi.getAsSTRINGWeighted(false);
			STRING_network.writePPIN(results_out_folder + "networks/" + file_name);
		}
	}
	
	public static void main(String[] args) {
		new File("/Users/tho/Desktop/STRING/").mkdir();
		
		System.out.println("IntAct");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/IntAct/networks/", "/Users/tho/Desktop/STRING/IntAct/");
		PPIN ppi = new PPIN("mixed_data/human_intact.tsv");
		PPIN STRING_network = ppi.getAsSTRINGWeighted(false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/IntAct/unpruned.tsv");
		
		System.out.println("PrePPI");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/PrePPI/networks/", "/Users/tho/Desktop/STRING/PrePPI/");
		ppi = new PPIN("mixed_data/human_ppi.tsv");
		STRING_network = ppi.getAsSTRINGWeighted(false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/PrePPI/unpruned.tsv");
		
		System.out.println("PrePPI95");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/PrePPI95/networks/", "/Users/tho/Desktop/STRING/PrePPI95/");
		ppi = new PPIN("mixed_data/human_ppi_0.95.tsv");
		STRING_network = ppi.getAsSTRINGWeighted(false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/PrePPI95/unpruned.tsv");
		
		System.out.println("BioGRID");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/BioGRID/networks/", "/Users/tho/Desktop/STRING/BioGRID/");
		ppi = new PPIN("mixed_data/human_biogrid.tsv");
		STRING_network = ppi.getAsSTRINGWeighted(false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/BioGRID/unpruned.tsv");
		
		System.out.println("HIPPIE");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/HIPPIE/networks/", "/Users/tho/Desktop/STRING/HIPPIE/");
		ppi = new PPIN("mixed_data/human_hippie.tsv");
		STRING_network = ppi.getAsSTRINGWeighted(false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/HIPPIE/unpruned.tsv");
		
		System.out.println("HPRD");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/HPRD/networks/", "/Users/tho/Desktop/STRING/HPRD/");
		ppi = new PPIN("mixed_data/human_hprd.tsv");
		STRING_network = ppi.getAsSTRINGWeighted(false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/HPRD/unpruned.tsv");
	}
}
