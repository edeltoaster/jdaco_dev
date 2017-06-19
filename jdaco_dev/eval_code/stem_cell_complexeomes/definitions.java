package stem_cell_complexeomes;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import framework.GOAnnotator;
import framework.Utilities;

public class definitions {
	static String daco_results_folder = "res_95_5/";
	static String networks_folder = "ENCODE_networks/";
	
	static double pvalue = 0.001;
	static boolean parametric = false;
	static boolean check_supersets = false;
	
	static int no_threads = 48;
	
	static String seed_file = "mixed_data/hocomoco_human_TFs_v10.txt.gz";
	static Set<String> seed = Utilities.readEntryFile(seed_file);
	
	static Set<String> pluri_factors = new HashSet<>(Arrays.asList("Q01860", "P48431", "Q9H9S0"));
	
	static String binding_data = "mixed_data/hocomoco_v10_EPD_v4_5k.txt.gz";
	static int d_min = -30;
	static int d_max = 30;
	
	static String GOA_def_file = "mixed_data/stem_tags_retrieved.txt.gz";
	static GOAnnotator goa = new GOAnnotator(GOA_def_file);
	
	public static void printParameters() {
		System.out.println("DACO results folder : " + daco_results_folder);
		System.out.println("networks folder : " + networks_folder);
		
		System.out.println("p-value : " + pvalue);
		System.out.println("parametric : " + parametric);
		System.out.println("checking supersets : " + check_supersets);
		
		System.out.println("no_threads : " + no_threads);
		
		System.out.println("seed file : " + seed_file);
		
		System.out.println("considered as pluri factors : " + pluri_factors);
		
		System.out.println("binding data file : " + binding_data);
		System.out.println("d_min : " + d_min);
		System.out.println("d_max : " + d_max);
		
		System.out.println("GOA definition file : " + GOA_def_file);
	}
	
	public static void main(String[] args) {
		// running updates the GO annotations
		GOAnnotator goa = new GOAnnotator("9606", true, "/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags.txt");
		goa.writeRetrievedData("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags_retrieved.txt.gz");
		goa.printTagInformation();
	}
}
