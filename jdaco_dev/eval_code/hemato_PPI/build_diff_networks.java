package hemato_PPI;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.Utilities;
import framework.RewiringDetector;
import framework.StrPair;

public class build_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/";
	//static String results_root = "/Users/tho/Desktop/BLUEPRINT_diffnets/";
	static String results_root = "/Users/tho/Desktop/BLUEPRINT_diffnets_filtered/";
	
	public static Map<String, ConstructedNetworks> filterVenous(Map<String, ConstructedNetworks> input) {
		Map<String, ConstructedNetworks> filtered = new HashMap<String, ConstructedNetworks>();
		for (String sample:input.keySet()) {
			if (sample.startsWith("Venous"))
				continue;
			filtered.put(sample, input.get(sample));
		}
		return filtered;
	}
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder);
		System.out.println();
		
		// define relations
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"HSC", "MPP"});
		
		relations.add(new String[]{"MPP", "CMP"});
		relations.add(new String[]{"MPP", "CLP"});
		
		relations.add(new String[]{"CMP", "MEP"});
		relations.add(new String[]{"CMP", "GMP"});
		
		relations.add(new String[]{"MEP", "MK"});
		relations.add(new String[]{"MEP", "EB"});
		
		relations.add(new String[]{"GMP", "N"});
		relations.add(new String[]{"GMP", "M"});
		
		relations.add(new String[]{"CLP", "CD4"});
		
		// additional relations
		relations.add(new String[]{"CMP", "CLP"});
		relations.add(new String[]{"MEP", "GMP"});
		relations.add(new String[]{"MK", "EB"});
		relations.add(new String[]{"N", "M"});
		
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			//Map<String, ConstructedNetworks> g1 = ConstructedNetworks.readNetworks(network_folder + state1 + "/");
			//Map<String, ConstructedNetworks> g2 = ConstructedNetworks.readNetworks(network_folder + state2 + "/");
			Map<String, ConstructedNetworks> g1 = filterVenous(ConstructedNetworks.readNetworks(network_folder + state1 + "/"));
			Map<String, ConstructedNetworks> g2 = filterVenous(ConstructedNetworks.readNetworks(network_folder + state2 + "/"));
			
			RewiringDetector rd = new RewiringDetector(g1, g2, FDR, 4);
			System.out.print("Processing " + state1 + " (" + g1.keySet().size() + ") vs " + state2 + " (" + g2.keySet().size() + ") : ");
			rd.writeDiffnet(results_folder + state1 + "_" + state2 + ".txt");
			
			Map<String, List<StrPair>> major_alt_splice_switches = rd.determineAltSplicingSwitches(true, false);
			Utilities.writeEntries(major_alt_splice_switches.keySet(), results_folder + state1 + "_" + state2 + "_major_AS_proteins.txt");
			
			Map<String, List<StrPair>> all_alt_splice_switches = rd.determineAltSplicingSwitches(false, true);
			Utilities.writeEntries(all_alt_splice_switches.keySet(), results_folder + state1 + "_" + state2 + "_contributing_AS_proteins.txt");
			
			
			double P_rew_rounded = (double) Math.round( Utilities.getMean(rd.getP_rews().values() ) * 1000d) / 1000d;
			System.out.println(rd.getP_rews().size()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getInteractionReasonsMap().size() + " dIAs" );
			
			System.out.println(major_alt_splice_switches.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches).size() + " diff. interactions.");
			System.out.println(all_alt_splice_switches.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches).size() + " diff. interactions that are mainly driven by AS events.");
			System.out.println();
			
//			// distribution of unfiltered changes
//			List<String> count_distr = new LinkedList<>();
//			for (Double d:rd.getDifferentialNetwork().values()) {
//				double ad = Math.abs(d);
//				count_distr.add( Double.toString(ad) );
//			}
//			
//			Utilities.writeEntries(count_distr, results_folder + state1 + "_" + state2 + "_raw_count_distr.txt");
		}
		
		System.out.println();
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		new File(results_root).mkdir();
		
		for (File f:new File(network_folder).listFiles()) {
			
			if (!f.isDirectory())
				continue;
			
			String threshold_results = f.getName();
			process(f.getAbsolutePath() + "/", results_root + threshold_results + "/");
			
		}
		
	}
}
