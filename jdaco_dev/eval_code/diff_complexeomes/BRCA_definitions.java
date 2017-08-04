package diff_complexeomes;

import java.util.Set;

import framework.GOAnnotator;
import framework.Utilities;

public class BRCA_definitions {
	static String daco_results_folder = "res_99_5/";
	static String networks_folder = "BRCA_networks_0/";
	static String diff_compl_output_folder = "diffnet_results_99_5_-25-25/"; // -15-10 should cover most known dimers, up to 25 even more; few may go higher; -25 as 24 is largest motif and overlap should be allowed
	
	static double qvalue = 0.05;
	static boolean parametric = false;
	static boolean paired = true;
	static boolean check_supersets = false;
	static double min_variant_fraction = 0.5;
	
	static int no_threads = 48;
	
	static String seed_file = "mixed_data/hocomoco_human_TFs_v10.txt.gz";
	static Set<String> seed = Utilities.readEntryFile(seed_file);
	
	static String binding_data = "mixed_data/hocomoco_v10_EPD_v4_5k.txt.gz";
	static int d_min = -25;
	static int d_max = 25;
	
	static String GOA_def_file = "mixed_data/stem_tags_retrieved.txt.gz";
	static GOAnnotator goa = new GOAnnotator(GOA_def_file);
	
	static int SPEnrich_iterations = 10000;
	static int SPEnrich_compl_part_threshold = 10;
	
	public static void printParameters() {
		System.out.println("DACO results folder : " + daco_results_folder);
		System.out.println("networks folder : " + networks_folder);
		System.out.println("results output folder : " + diff_compl_output_folder);
		
		System.out.println("q-value : " + qvalue);
		System.out.println("parametric : " + parametric);
		System.out.println("paired : " + paired);
		System.out.println("checking supersets : " + check_supersets);
		System.out.println("min variant fraction : " + min_variant_fraction);
		
		System.out.println("no_threads : " + no_threads);
		
		System.out.println("seed file : " + seed_file);
		
		System.out.println("binding data file : " + binding_data);
		System.out.println("d_min : " + d_min);
		System.out.println("d_max : " + d_max);
		
		//System.out.println("GOA definition file : " + GOA_def_file);
		
		System.out.println("SPE iterations : " + SPEnrich_iterations);
		System.out.println("SPE compl. part. threshold : " + SPEnrich_compl_part_threshold);
	}
	
	public static void main(String[] args) {
		// running updates the GO annotations
		GOAnnotator goa = new GOAnnotator("9606", false, "/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags.txt");
		goa.writeRetrievedData("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags_retrieved.txt.gz");
		goa.printTagInformation();
	}
}
