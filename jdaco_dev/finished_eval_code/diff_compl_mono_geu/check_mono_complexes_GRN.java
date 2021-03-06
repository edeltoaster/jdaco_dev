package diff_compl_mono_geu;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.GOAnnotator;
import framework.RegulatoryNetwork;

public class check_mono_complexes_GRN {
	static int no_threads = 4;
	static String out_folder = "/Users/tho/Desktop/";
	static Map<String, String> up_name = DataQuery.getUniprotToGeneNameMap("homo_sapiens_core_90_38");
	static GOAnnotator goa = new GOAnnotator("9606", "mixed_data/mono_tags.txt.gz");
	
	static double q_threshold = 0.001;
	
	public static void main(String[] args) {
		
		DiffComplexDetector.SignSortedComplexesResult res = DiffComplexDetector.readSignSortedComplexResult(definitions.compl_results, q_threshold);
		
		Set<String> relevant_TFs = new HashSet<>();
		Set<String> relevant_targets = new HashSet<>(definitions.markers);
		res.getSignificanceSortedComplexes().forEach(c -> relevant_TFs.addAll(res.getMemberSeedComb().get(c)));
		
		System.out.println("Read binding data ...");
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data, relevant_TFs, relevant_targets);
		
		System.out.println("Building RegNet ...");
		RegulatoryNetwork regnet = new RegulatoryNetwork(res.getSignificanceSortedComplexes(), relevant_TFs, bdh, definitions.d_min, definitions.d_max, no_threads, 1);
		
		Map<String, Map<String,String>> annotational_data = new HashMap<String, Map<String,String>>();
		
		// relate TFCs to sign. dereg. complexes
		Map<HashSet<String>, List<HashSet<String>>> tfc_complex_map = new HashMap<>();
		for (HashSet<String> complex : res.getMemberSeedComb().keySet()) {
			HashSet<String> tfc = res.getMemberSeedComb().get(complex);
			if (!tfc_complex_map.containsKey(tfc))
				tfc_complex_map.put(tfc, new LinkedList<HashSet<String>>());
			tfc_complex_map.get(tfc).add(complex);
		}
		
		// convert and store complex sets to strings of their gene names
		Map<HashSet<String>, String> complex_genes_map = new HashMap<>();
		res.getSignificanceSortedComplexes().forEach(c -> complex_genes_map.put(c, String.join("/", c.stream().map(p -> up_name.getOrDefault(p, p)).collect(Collectors.toList()))));
		
		// store directions
		Map<HashSet<String>, String> complex_dir_map = new HashMap<>();
		Map<String, String> dir_map = new HashMap<>();
		res.getSignificanceSortedComplexes().forEach(c -> complex_dir_map.put(c, complex_genes_map.get(c) + ":" + res.getDirections().get(c)));
		tfc_complex_map.keySet().forEach(tfc -> dir_map.put(tfc.toString(), String.join(",", tfc_complex_map.get(tfc).stream().map(c -> complex_dir_map.get(c)).collect(Collectors.toList()))));
		
		// store fold-change
		Map<HashSet<String>, String> complex_fc_map = new HashMap<>();
		Map<String, String> fc_map = new HashMap<>();
		res.getSignificanceSortedComplexes().forEach(c -> complex_fc_map.put(c, complex_genes_map.get(c) + ":" + String.format(Locale.US, "%.3g", res.getFoldChange().get(c))));
		tfc_complex_map.keySet().forEach(tfc -> fc_map.put(tfc.toString(), String.join(",", tfc_complex_map.get(tfc).stream().map(c -> complex_fc_map.get(c)).collect(Collectors.toList()))));
		
		// store median-change
		Map<HashSet<String>, String> complex_med_map = new HashMap<>();
		Map<String, String> med_map = new HashMap<>();
		res.getSignificanceSortedComplexes().forEach(c -> complex_med_map.put(c, complex_genes_map.get(c) + ":" + String.format(Locale.US, "%.3g", res.getMeanChange().get(c))));
		tfc_complex_map.keySet().forEach(tfc -> med_map.put(tfc.toString(), String.join(",", tfc_complex_map.get(tfc).stream().map(c -> complex_med_map.get(c)).collect(Collectors.toList()))));
				

		// store GOA tags
		Map<HashSet<String>, String> complex_GOA_map = new HashMap<>();
		Map<String, String> goa_map = new HashMap<>();
		res.getSignificanceSortedComplexes().forEach(c -> complex_GOA_map.put(c, complex_genes_map.get(c) + ":" + goa.rateProteins(c)));
		tfc_complex_map.keySet().forEach(tfc -> goa_map.put(tfc.toString(), String.join(",", tfc_complex_map.get(tfc).stream().map(c -> complex_GOA_map.get(c)).collect(Collectors.toList()))));
		
		annotational_data.put("Directions", dir_map);
		annotational_data.put("Fold-Changes", fc_map);
		annotational_data.put("Median-Changes", med_map);
		annotational_data.put("GOA", goa_map);
		annotational_data.put("Marker-Directions", definitions.marker_dir);
		
		
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(out_folder + "regnet.txt", 2);
		regnet.writeRegulatoryNetwork(out_folder + "regnet_all.txt", 1);
		regnet.writeNodeTable(out_folder + "regnet_nodes.txt", annotational_data);
		
		
		/*
		 * incorporating interactions among TFs
		 */
		
		// TFs incorporated as targets
		res.getSignificanceSortedComplexes().forEach(c -> relevant_targets.addAll(res.getMemberSeedComb().get(c)));
		
		System.out.println("Read binding data 2 ...");
		bdh = new BindingDataHandler(definitions.binding_data, relevant_TFs, relevant_targets);
		
		System.out.println("Building RegNet 2...");
		regnet = new RegulatoryNetwork(res.getSignificanceSortedComplexes(), relevant_TFs, bdh, definitions.d_min, definitions.d_max, no_threads, 1);

		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(out_folder + "regnet2.txt", 2);
		regnet.writeRegulatoryNetwork(out_folder + "regnet2_all.txt", 1);
		regnet.writeNodeTable(out_folder + "regnet_nodes2.txt");
		
		// prune network
		regnet.pruneToLargestSCCs();
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork(out_folder + "regnet2p.txt", 2);
		regnet.writeRegulatoryNetwork(out_folder + "regnet2p_all.txt", 1);
		regnet.writeNodeTable(out_folder + "regnet2p_nodes.txt");
		
	}
}
