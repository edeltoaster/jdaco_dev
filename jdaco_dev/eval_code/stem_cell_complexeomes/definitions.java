package stem_cell_complexeomes;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import framework.GOAnnotator;
import framework.Utilities;

public class definitions {
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_PrePPIhc_TPMgene/res5/";
	static String networks_folder = "/Users/tho/Desktop/PrePPIhc_TPMgene_networks/";
	
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	static Set<String> pluri_factors = new HashSet<>(Arrays.asList("Q01860", "P48431", "Q9H9S0"));
	
	static String binding_data = "/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10/hocomoco_v10_EPD_2.5k.txt.gz";
	static GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/simple_tags_retrieved.txt.gz");
}