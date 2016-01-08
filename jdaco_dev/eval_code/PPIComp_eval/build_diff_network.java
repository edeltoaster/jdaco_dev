package PPIComp_eval;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.Utilities;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.StrPair;

public class build_diff_network {
	
	static int no_threads = 48;
	static double FDR = 0.05;
	
	// needs to be run on server
	static String network_folder = "BRCA_networks/";
	static String results_root = "BRCA_diffnet/";
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder);
		System.out.println();
		
		// define relations
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"normal", "tumor"});
		
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			System.out.println("read data");
			Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + state1 + "/");
			Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + state2 + "/");
			
			System.out.println("start processing");
			RewiringDetector rd = new RewiringDetector(g1, g2, FDR, no_threads, System.out);
			
			System.out.println("start output");
			rd.writeDiffnet(results_folder + state1 + "_" + state2 + ".txt");
			
			Map<String, List<StrPair>> major_alt_splice_switches = rd.determineAltSplicingSwitches(true, false);
			Utilities.writeEntries(major_alt_splice_switches.keySet(), results_folder + state1 + "_" + state2 + "_major_AS_proteins.txt");
			
			Map<String, List<StrPair>> all_alt_splice_switches = rd.determineAltSplicingSwitches(false, true);
			Utilities.writeEntries(all_alt_splice_switches.keySet(), results_folder + state1 + "_" + state2 + "_contributing_AS_proteins.txt");
			
			
			double P_rew_rounded = (double) Math.round( rd.getP_rew() * 1000d) / 1000d;
			System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getSignificantlyRewiredInteractions().size() + " dIAs" );
			
			System.out.println(major_alt_splice_switches.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches).size() + " diff. interactions.");
			System.out.println(all_alt_splice_switches.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches).size() + " diff. interactions that are mainly driven by AS events.");
			
			List<String> minReasons = rd.getMinMostLikelyReasons();
			System.out.println(minReasons.size() + " alterations can explain all significant changes.");
			Utilities.writeEntries(minReasons, results_folder + state1 + "_" + state2 + "_min_reasons.txt");
			rd.writeProteinAttributes(results_folder + state1 + "_" + state2 + "_protein_attributes.txt");
			System.out.println();

		}
		
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		new File(results_root).mkdir();
		
		for (File f:Utilities.listDirectoriesAndFilesWithinFolder(new File(network_folder))) {
			
			if (!f.isDirectory())
				continue;
			
			String threshold_results = f.getName();
			process(f.getAbsolutePath() + "/", results_root + threshold_results + "/");
			
		}
		
	}
}
