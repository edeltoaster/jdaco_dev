package PPIComp_hemato;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.Utilities;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.StrPair;

public class build_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/0.31/";
	static String results_root = "/Users/tho/Desktop/BLUEPRINT_diffnets/";
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder);
		System.out.println();
		
		// define the developmental relations that are tested
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"HSC", "MPP"});
		relations.add(new String[]{"MPP", "CMP", "CLP"});
		relations.add(new String[]{"CMP", "MEP", "GMP"});
		relations.add(new String[]{"MEP", "MK", "EB"});
		relations.add(new String[]{"GMP", "N", "M"});
		relations.add(new String[]{"CLP", "CD4"});
		
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
			
			// do more stuff if lineage branch
			if (s.length == 3) {
				System.out.println();
				String state3 = s[2];
				Map<String, RewiringDetectorSample> g3 = RewiringDetectorSample.readNetworks(network_folder + state3 + "/");
				
				// same as before
				RewiringDetector rd2 = new RewiringDetector(g1, g3, FDR);
				System.out.print("Processing " + state1 + " (" + g1.keySet().size() + ") vs " + state3 + " (" + g3.keySet().size() + ") : ");
				rd2.writeDiffnet(out_path + state1 + "_" + state3 + ".txt");
				
				Map<String, List<StrPair>> major_alt_splice_switches2 = rd2.determineAltSplicingSwitches(true, false);
				Utilities.writeEntries(major_alt_splice_switches2.keySet(), out_path + state1 + "_" + state3 + "_major_AS_proteins.txt");
				
				Map<String, List<StrPair>> all_alt_splice_switches2 = rd2.determineAltSplicingSwitches(false, true);
				Utilities.writeEntries(all_alt_splice_switches2.keySet(), out_path + state1 + "_" + state3 + "_contributing_AS_proteins.txt");
				
				
				double P_rew_rounded2 = (double) Math.round( rd2.getP_rew() * 1000d) / 1000d;
				System.out.println(rd2.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded2 + ", " + rd2.getInteractionReasonsCountMap().size() + " dIAs" );
				
				System.out.println(major_alt_splice_switches2.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches2).size() + " diff. interactions.");
				System.out.println(all_alt_splice_switches2.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches2).size() + " diff. interactions that are mainly driven by AS events.");
				
				List<String> minReasons2 = rd2.getMinMostLikelyReasons();
				System.out.println(minReasons2.size() + " alterations can explain all significant changes.");
				Utilities.writeEntries(minReasons2, out_path + state1 + "_" + state3 + "_min_reasons.txt");
				rd2.writeProteinAttributes(out_path + state1 + "_" + state3 + "_protein_attributes.txt");
				
				// look into shared/exclusive changes
				Set<StrPair> sign_IA_state2 = new HashSet<>(rd1.getSignificantlyRewiredInteractions());
				Set<StrPair> sign_IA_state3 = new HashSet<>(rd2.getSignificantlyRewiredInteractions());
				
				Set<StrPair> sign_IA_shared = new HashSet<>(sign_IA_state2);
				sign_IA_shared.retainAll(sign_IA_state3);
				
				Set<StrPair> sign_IA_excl_state2 = new HashSet<>(sign_IA_state2);
				sign_IA_excl_state2.removeAll(sign_IA_state3);
				Set<StrPair> sign_IA_excl_state3 = new HashSet<>(sign_IA_state3);
				sign_IA_excl_state3.removeAll(sign_IA_state2);
				
				// output to files
				Utilities.writeEntries(sign_IA_shared, out_path + state2 + "_" + state3 + "_shared.txt");
				Utilities.writeEntries(sign_IA_excl_state2, out_path + state1 + "_" + state2 + "_excl.txt");
				Utilities.writeEntries(sign_IA_excl_state3, out_path + state1 + "_" + state3 + "_excl.txt");
				
				// also compute and store minReasons of subsets
				Utilities.writeEntries(rd1.getMinMostLikelyReasons(sign_IA_shared), out_path + state1 + "_" + state2 + "_shared_min_reasons.txt");
				Utilities.writeEntries(rd2.getMinMostLikelyReasons(sign_IA_shared), out_path + state1 + "_" + state3 + "_shared_min_reasons.txt");
				Utilities.writeEntries(rd1.getMinMostLikelyReasons(sign_IA_excl_state2), out_path + state1 + "_" + state2 + "_excl_min_reasons.txt");
				Utilities.writeEntries(rd2.getMinMostLikelyReasons(sign_IA_excl_state3), out_path + state1 + "_" + state3 + "_excl_min_reasons.txt");
			}
			
			System.out.println();
		}
		
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		process(network_folder, results_root);
	}
}
