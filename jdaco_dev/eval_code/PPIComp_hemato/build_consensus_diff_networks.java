package PPIComp_hemato;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.RewiringDetectorSample;
import framework.StrPair;

public class build_consensus_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BLUEPRINT_networks/0.31/";
	
	public static void process(String network_folder) {

		System.out.println("Consensus-Analysis for " + network_folder);
		System.out.println();
		
		// define the developmental relations that are tested
		List<String[]> relations = new LinkedList<String[]>();
		relations.add(new String[]{"HSC", "MPP"});
		relations.add(new String[]{"MPP", "CMP"});
		relations.add(new String[]{"MPP", "CLP"});
		relations.add(new String[]{"CMP", "MEP"});
		relations.add(new String[]{"CMP", "GMP"});
		relations.add(new String[]{"MEP", "EB"});
		relations.add(new String[]{"MEP", "MK"});
		relations.add(new String[]{"GMP", "N"});
		relations.add(new String[]{"GMP", "M"});
		relations.add(new String[]{"CLP", "CD4"});
		
		Set<StrPair> all_interactions = null;
		for (String[] s:relations) {
			String state1 = s[0];
			String state2 = s[1];
			System.out.println(state1 + " -> " + state2);
			Map<String, RewiringDetectorSample> g1 = RewiringDetectorSample.readNetworks(network_folder + state1 + "/");
			Map<String, RewiringDetectorSample> g2 = RewiringDetectorSample.readNetworks(network_folder + state2 + "/");
			
			Set<StrPair> g1_interactions = null;
			for (RewiringDetectorSample rws:g1.values()) {
				if (g1_interactions == null)
					g1_interactions = new HashSet<>(rws.getInteractions());
				else
					g1_interactions.retainAll(rws.getInteractions());
				
				if (all_interactions == null)
					all_interactions = new HashSet<>(rws.getInteractions());
				else
					all_interactions.retainAll(rws.getInteractions());
			}
			
			Set<StrPair> g2_interactions = null;
			for (RewiringDetectorSample rws:g2.values()) {
				if (g2_interactions == null)
					g2_interactions = new HashSet<>(rws.getInteractions());
				else
					g2_interactions.retainAll(rws.getInteractions());
				
				if (all_interactions == null)
					all_interactions = new HashSet<>(rws.getInteractions());
				else
					all_interactions.retainAll(rws.getInteractions());
			}
			
			System.out.println(g1_interactions.size() + " & " + g2_interactions.size() + " -> " + (g2_interactions.size() - g1_interactions.size()) );
		}
		
		System.out.println();
		System.out.println("conserved in all: " + all_interactions.size());
	}
	
	public static void main(String[] args) {
		process(network_folder);
	}
}
