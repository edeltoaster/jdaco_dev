package stem_cell_complexomes_older;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import framework.GOAnnotator;
import framework.Utilities;

public class definitions {
	static String daco_results_folder = "/Users/tho/GDrive/Work/projects/stem_cell_complexome/DACO_PrePPIhc_TPMgene/res5/";
	static String networks_folder = "/Users/tho/Desktop/PrePPIhc_TPMgene_networks/";
	
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	static Set<String> pluri_factors = new HashSet<>(Arrays.asList("Q01860", "P48431", "Q9H9S0"));
	
	static String binding_data = "/Users/tho/GDrive/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz";
	static GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags_retrieved.txt.gz");
	
	public static void main(String[] args) {
		// running updates the GO annotations
		GOAnnotator goa = new GOAnnotator("9606", true, "/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags.txt");
		goa.writeRetrievedData("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags_retrieved.txt.gz");
		goa.printTagInformation();
	}
}
