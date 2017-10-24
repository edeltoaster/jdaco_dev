package diff_compl_mono_geu;

import java.util.Set;

import framework.Utilities;

public class definitions {
	static String daco_results_folder = "res_95_5/";
	static String networks_folder = "networks/";
	static String diff_complex_output_folder = "diffcompl_results_95_5_-25-25/";
	static String diff_tfc_output_folder = "difftfc_results_95_5_-25-25/"; // -15-10 should cover most known dimers, up to 25 even more; few may go higher; -25 as 24 is largest motif and overlap should be allowed
	static String qr_output_folder = "q_results_95_5/";
	
	static double qvalue = 0.05;
	static boolean parametric = false;
	static boolean paired = false;
	static boolean check_supersets = false;
	static double min_variant_fraction = 0.75;
	
	static int no_threads = 48;
	
	static String seed_file = "mixed_data/hocomoco_human_TFs_v10.txt.gz";
	static Set<String> seed = Utilities.readEntryFile(seed_file);
	
	static int SPCEnrich_iterations = 10000;
	static int SPCEnrich_compl_part_threshold = 5;
	static int SPEnrich_iterations = 10000;
	static int SPEnrich_compl_part_threshold = 10;
	
	static String binding_data = "/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10_EPD_v5_2k.txt.gz";
	static int d_min = -25;
	static int d_max = 25;
	
	public static void printInitParameters() {
		System.out.println("DACO results folder : " + daco_results_folder);
		System.out.println("networks folder : " + networks_folder);
		System.out.println("diff compl results output folder : " + diff_complex_output_folder);
		System.out.println("diff tfc results output folder : " + diff_tfc_output_folder);
		System.out.println("quantified results output folder : " + qr_output_folder);
		
		System.out.println("q-value : " + qvalue);
		System.out.println("parametric : " + parametric);
		System.out.println("paired : " + paired);
		System.out.println("checking supersets : " + check_supersets);
		System.out.println("min variant fraction : " + min_variant_fraction);
		
		System.out.println("no_threads : " + no_threads);
		
		System.out.println("seed file : " + seed_file);
		
		System.out.println("SPC iterations : " + SPCEnrich_iterations);
		System.out.println("SPC compl. part. threshold : " + SPCEnrich_compl_part_threshold);
		System.out.println("SPE iterations : " + SPEnrich_iterations);
		System.out.println("SPE compl. part. threshold : " + SPEnrich_compl_part_threshold);
	}
	
	public static void printBindingDataParameters() {
		System.out.println("Binding data : " + binding_data);
		System.out.println("d_min/d_max : " + d_min + "-" + d_max);
	}
}
