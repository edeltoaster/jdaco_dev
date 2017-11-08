package diff_compl_mono_geu;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DiffComplexDetector;
import framework.RegulatoryNetwork;

public class check_mono_complexes_GRN {
	static String compl_results = "/Users/tho/Dropbox/Work/projects/diff_complexome/results/diff_results_95_5_-25-25/mono_dcd_compl.txt";
	static int no_threads = 4;
	static List<String> markers = Arrays.asList("P08637", "O75015", "P08571");
	static String out_folder = "/Users/tho/Desktop/";
	
	public static void main(String[] args) {
		
		DiffComplexDetector.SignSortedComplexesResult res = DiffComplexDetector.readSignSortedComplexResult(compl_results);
		Set<String> relevant_proteins = new HashSet<>(markers);
		res.getMemberSeedComb().values().forEach(c -> relevant_proteins.addAll(c));
		
		System.out.println("Read binding data ...");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, relevant_proteins, relevant_proteins);
		
		System.out.println("Building RegNet ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(res.getSignificanceSortedComplexes(), bdh, definitions.d_min, definitions.d_max, no_threads, 1);
		
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(out_folder + "regnet.txt", 2);
		regnet.writeNodeTable(out_folder + "regnet_nodes.txt");
	}
}
