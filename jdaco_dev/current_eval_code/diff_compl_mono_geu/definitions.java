package diff_compl_mono_geu;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.Utilities;

public class definitions {
	static String daco_results_folder = "res_95_5/";
	static String local_daco_results_folder = "/Users/tho/GDrive/Work/projects/CompleXChange/results/res_95_5/";
	static String networks_folder = "networks/";
	static String diff_out_folder = "diff_results_95_5/";
	static String qr_output_folder = "q_results_95_5/";
	
	static double qvalue = 0.05;
	static boolean parametric = false;
	static boolean paired = false;
	static boolean check_supersets = false;
	static double min_variant_fraction = 0.75;
	
	static int no_threads = 64;
	
	static String seed_file = "mixed_data/hocomoco_human_TFs_v10.txt.gz";
	static Set<String> seed = Utilities.readEntryFile(seed_file);
	
	static int SPCEnrich_iterations = 10000;
	static int SPCEnrich_compl_part_threshold = 5;
	static int SPEnrich_iterations = 10000;
	static int SPEnrich_compl_part_threshold = 10;
	
	static String binding_data = "/Users/tho/GDrive/Work/data_general/binding_sites/hocomoco_v10_EPD_v5_2k.txt.gz";
	static int d_min = -15; // -15-10 should cover most known dimers, up to 25 even more; few may go higher; -25 as 24 is largest motif and overlap should be allowed
	static int d_max = 10; // 10?
	
	static String compl_results = "/Users/tho/GDrive/Work/projects/CompleXChange/results/diff_results_95_5/unpaired_nonparametric/mono_dcd_compl.txt";
	static List<String> markers = Arrays.asList("P08637", "O75015", "P08571", "P41597", "P49238", "P14151", "P20701"); // CD16, CD16, CD14, CCR2, CX3CR1, SELL/CD62L, ITGAL
	static Map<String, String> marker_dir = new HashMap<String, String>() {
		private static final long serialVersionUID = 1L;
		{
			put("P08637", "+"); // CD16
			put("O75015", "+"); // CD16
			put("P08571", "-"); // CD14
			put("P41597", "-"); // CCR2 -> not sufficiently expressed in expr data
			put("P49238", "+"); // CX3CR1
			put("P14151", "-"); // SELL / CD62L
			put("P20701", "+"); // ITGAL
		}
	};
	
	public static void printInitParameters() {
		System.out.println("DACO results folder : " + daco_results_folder);
		System.out.println("networks folder : " + networks_folder);
		System.out.println("diff results output folder : " + diff_out_folder);
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

		System.out.println("Markers: " + markers);
	}
	
	public static void printBindingDataParameters() {
		System.out.println("Binding data : " + binding_data);
		System.out.println("d_min/d_max : " + d_min + "-" + d_max);
	}
}
