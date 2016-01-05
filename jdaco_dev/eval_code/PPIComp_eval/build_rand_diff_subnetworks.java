package PPIComp_eval;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.Utilities;
import framework.RewiringDetector;
import framework.StrPair;

public class build_rand_diff_subnetworks {
	
	static double FDR = 0.05;
	static double fraction = 0.5;
	static int min_size = 3;
	static int iterations = 5;
	
	// needs to be run on server
	static String network_folder = "BRCA_networks/";
	static String results_root = "BRCA_diffnets/";
	
	public static Map<String, ConstructedNetworks> getRandomSubset(Map<String, ConstructedNetworks> data, int min_size, double fraction) {
		Map<String, ConstructedNetworks> subset = new HashMap<String, ConstructedNetworks>();
		int sample_size = Math.max(min_size, (int) fraction * data.size());
		
		List<String> samples = new ArrayList<>(data.keySet());
		Collections.shuffle(samples);
		
		for (int i=0;i<sample_size;i++) {
			String sample = samples.get(i);
			subset.put(sample, data.get(sample));
		}
		
		return subset;
	}
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder + ", " + fraction + " fraction used.");
		System.out.println();
		
		// define relations
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"normal", "tumor"});
		
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			System.out.println("read data");
			Map<String, ConstructedNetworks> g1 = ConstructedNetworks.readNetworks(network_folder + state1 + "/");
			Map<String, ConstructedNetworks> g2 = ConstructedNetworks.readNetworks(network_folder + state2 + "/");
			
			for (int i=1;i<=iterations;i++) {
				String run_id = state1 + "_" + state2 + "_" + i;
				System.out.println("start iteration " + i);
				
				System.out.println("start sampling for " + i);
				Map<String, ConstructedNetworks> g1s = getRandomSubset(g1, min_size, fraction);
				Map<String, ConstructedNetworks> g2s = getRandomSubset(g2, min_size, fraction);
				
				System.out.println("start RD calculations for " + i);
				RewiringDetector rd = new RewiringDetector(g1s, g2s, FDR, 32, System.out);
				
				System.out.println("start output of " + i);
				rd.writeDiffnet(results_folder + run_id + ".txt");
				
				Map<String, List<StrPair>> major_alt_splice_switches = rd.determineAltSplicingSwitches(true, false);
				Utilities.writeEntries(major_alt_splice_switches.keySet(), results_folder + run_id + "_major_AS_proteins.txt");
				
				Map<String, List<StrPair>> all_alt_splice_switches = rd.determineAltSplicingSwitches(false, true);
				Utilities.writeEntries(all_alt_splice_switches.keySet(), results_folder + run_id + "_contributing_AS_proteins.txt");
				
				
				double P_rew_rounded = (double) Math.round( rd.getP_rew() * 1000d) / 1000d;
				System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getInteractionReasonsMap().size() + " dIAs" );
				
				System.out.println(major_alt_splice_switches.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches).size() + " diff. interactions.");
				System.out.println(all_alt_splice_switches.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches).size() + " diff. interactions that are mainly driven by AS events.");
				
				List<String> minReasons = rd.getMinMostLikelyReasons();
				System.out.println(minReasons.size() + " alterations can explain all significant changes.");
				Utilities.writeEntries(minReasons, results_folder + run_id + "_min_reasons.txt");
				System.out.println();
			}

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
			process(f.getAbsolutePath() + "/", results_root + threshold_results + "_rand/");
			
		}
		
	}
}
