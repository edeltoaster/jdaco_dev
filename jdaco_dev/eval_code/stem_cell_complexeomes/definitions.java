package stem_cell_complexeomes;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import framework.GOAnnotator;
import framework.Utilities;

public class definitions {
	static String daco_results_folder = "res_95_5/";
	static double pvalue = 0.01;
	static String networks_folder = "ENCODE_networks/";
	static boolean check_supersets = false;
	
	static Set<String> seed = Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	static Set<String> pluri_factors = new HashSet<>(Arrays.asList("Q01860", "P48431", "Q9H9S0"));
	
	static String binding_data = "mixed_data/hocomoco_v10_EPD_v4_5k.txt.gz";
	static GOAnnotator goa = new GOAnnotator("mixed_data/stem_tags_retrieved.txt.gz");
	
	static int no_threads = 40;
	
	public static void main(String[] args) {
		// running updates the GO annotations
		GOAnnotator goa = new GOAnnotator("9606", true, "/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags.txt");
		goa.writeRetrievedData("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags_retrieved.txt.gz");
		goa.printTagInformation();
	}
}
