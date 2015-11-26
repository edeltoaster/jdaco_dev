package hemato_PPI;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.ConstructedNetworks;
import framework.Utilities;
import framework.RewiringDetector;
import framework.StrPair;

public class build_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/";
	static String results_root = "/Users/tho/Desktop/BLUEPRINT_diffnets/";
	
	public static void process(String network_folder, String results_folder) {
		
		new File(results_folder).mkdir();
		
		System.out.println("Analysis for " + network_folder + ", writing to " + results_folder);
		
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
		
		Set<String> AS_proteins = new HashSet<>();
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			
			Map<String, ConstructedNetworks> g1 = ConstructedNetworks.readNetworks(network_folder + state1 + "/");
			Map<String, ConstructedNetworks> g2 = ConstructedNetworks.readNetworks(network_folder + state2 + "/");
			
			RewiringDetector rd = new RewiringDetector(g1, g2, FDR);
			System.out.print("Processing " + state1 + " vs " + state2 + " : ");
			rd.writeDiffnet(results_folder + state1 + "_" + state2 + ".txt");
			
			Map<String, List<StrPair>> alt_splice_switches = rd.determineAltSplicingSwitches(true);
			Utilities.writeEntries(alt_splice_switches.keySet(), results_folder + state1 + "_" + state2 + "_AS_proteins.txt");
			AS_proteins.addAll(alt_splice_switches.keySet());
			double P_rew_rounded = (double) Math.round( Utilities.getMean(rd.getP_rews().values() ) * 1000d) / 1000d;
			System.out.println(rd.getP_rews().size()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getInteractionReasonsMap().size() + " dIAs" );
			
			System.out.println(alt_splice_switches.keySet().size() + " alt. spliced proteins that affect " + alt_splice_switches.values().stream().mapToInt(e->e.size()).sum() + " interactions.");
		}
		
		Utilities.writeEntries(AS_proteins, results_folder + "all_AS_proteins.txt");
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
