package CD8_subtypes_public_and_SFB;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.Utilities;
import framework.DataQuery;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.StrPair;

public class build_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/CD8_subsets_public/CD8_networks_0.0/";
	static String results_root = "/Users/tho/Desktop/N_MNP_diffnets/";
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder);
		System.out.println();
		
		// define the developmental relations that are tested
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"N", "TMNP"});
		
		for (String[] s:relations) {
			
			String out_path = results_folder + String.join("_", s) + "/";
			new File(out_path).mkdir();
			
			String state1 = s[0];
			String state2 = s[1];
			Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + state1 + "/");
			Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + state2 + "/");
			
			RewiringDetector rd1 = new RewiringDetector(g1, g2, FDR);
			System.out.print("Processing " + state1 + " (" + g1.keySet().size() + ") vs " + state2 + " (" + g2.keySet().size() + ") : ");
			rd1.writeDiffnet(out_path + state1 + "_" + state2 + ".txt");
			
			Map<String, List<StrPair>> major_alt_splice_switches = rd1.determineAltSplicingSwitches(true, false);
			Utilities.writeEntries(major_alt_splice_switches.keySet(), out_path + state1 + "_" + state2 + "_major_AS_proteins.txt");
			
			Map<String, List<StrPair>> all_alt_splice_switches = rd1.determineAltSplicingSwitches(false, true);
			Utilities.writeEntries(all_alt_splice_switches.keySet(), out_path + state1 + "_" + state2 + "_contributing_AS_proteins.txt");
			
			
			double P_rew_rounded = (double) Math.round( rd1.getP_rew() * 1000d) / 1000d;
			System.out.println(rd1.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd1.getInteractionReasonsCountMap().size() + " dIAs" );
			
			System.out.println(major_alt_splice_switches.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches).size() + " diff. interactions.");
			System.out.println(all_alt_splice_switches.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches).size() + " diff. interactions that are mainly driven by AS events.");
			
			List<String> minReasons = rd1.getMinMostLikelyReasons();
			System.out.println(minReasons.size() + " alterations can explain all significant changes.");
			Utilities.writeEntries(minReasons, out_path + state1 + "_" + state2 + "_min_reasons.txt");
			rd1.writeProteinAttributes(out_path + state1 + "_" + state2 + "_protein_attributes.txt");
			
			System.out.println();
		}
		
		System.out.println();
	}
	
	public static void main(String[] args) {
		DataQuery.enforceSpecificEnsemblRelease("84");
		process(network_folder, results_root);
	}
}
