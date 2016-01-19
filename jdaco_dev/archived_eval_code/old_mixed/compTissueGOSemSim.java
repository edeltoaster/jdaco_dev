package old_mixed;

import java.io.File;
import framework.PPIN;
import framework.Utilities;

public class compTissueGOSemSim {
	
	// STRING with weights WITHOUT coexpression incorporated
	public static PPIN STRING_ref = new PPIN("/Users/tho/Dropbox/Work/tissue specificity/STRING/GOSemSim/pGOSemSim_intact.tsv");
	
	public static void computeForNetwork(String input_folder, String results_out_folder) {
		
		// make output folder
		new File(results_out_folder).mkdir();

		System.out.println("Reading and processing samples ...");
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders(input_folder, "ppin.tsv")) {
			
			String sample_path = file.getAbsolutePath();
			String file_name = file.getName();
			
			System.out.println("Processing "+ file_name);
			PPIN ppi = new PPIN(sample_path);
			PPIN STRING_network = new PPIN(ppi, STRING_ref, false);
			STRING_network.writePPIN(results_out_folder + file_name);
		}
	}
	
	public static void main(String[] args) {
		new File("/Users/tho/Desktop/STRING/").mkdir();
		
		System.out.println("IntAct");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/IntAct/networks/", "/Users/tho/Desktop/STRING/IntAct/");
		System.out.println("PrePPI");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/PrePPI/networks/", "/Users/tho/Desktop/STRING/PrePPI/");
		
		System.out.println("PrePPI95");
		computeForNetwork("/Users/tho/Dropbox/Work/tissue specificity/PrePPI95/networks/", "/Users/tho/Desktop/STRING/PrePPI95/");
		
		// unpruned ones
		PPIN ppi = new PPIN("mixed_data/human_intact.tsv");
		PPIN STRING_network = new PPIN(ppi, STRING_ref, false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/intact_STRING.tsv");
		
		ppi = new PPIN("mixed_data/human_ppi.tsv");
		STRING_network = new PPIN(ppi, STRING_ref, false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/preppi_STRING.tsv");
		
		ppi = new PPIN("mixed_data/human_ppi_0.95.tsv");
		STRING_network = new PPIN(ppi, STRING_ref, false);
		STRING_network.writePPIN("/Users/tho/Desktop/STRING/preppi95_STRING.tsv");
		
	}
}
