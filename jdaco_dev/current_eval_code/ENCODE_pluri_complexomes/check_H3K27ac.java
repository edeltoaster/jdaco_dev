package ENCODE_pluri_complexomes;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.Utilities;

public class check_H3K27ac {
	static String diff_compl_file = "/Users/tho/GDrive/Work/projects/stem_cell_complexome/diff_compl_0.75/preppi_sign.txt";
	
	public static void main(String[] args) {
		
		// read diff complexes
		List<Set<String>> tfcs = new LinkedList<>();
		Map<Set<String>, String> dir_map = new HashMap<>();
		for (String s : Utilities.readFile(diff_compl_file)) {
			if (s.startsWith("(sub)"))
					continue;
			Set<String> tfc = new HashSet<String>(Arrays.asList(s.split(" ")[0].split("/")));
			dir_map.put(tfc, s.split(" ")[1]);
			tfcs.add(tfc);
		}
		
		Set<String> Hac_proteins = new HashSet<>();
		for (String s:Utilities.readFile("mixed_data/H_ac.tsv")) {
			if (s.startsWith("GENE"))
					continue;
			Hac_proteins.add(s.split("\t")[1]);
		}
		
		Map<String, Integer> n = new HashMap<>();
		n.put("+", 0);
		n.put("-", 0);
		Map<String, Integer> n_ep300 = new HashMap<>();
		n_ep300.put("+", 0);
		n_ep300.put("-", 0);
		for (Set<String> tfc : tfcs) {
			Set<String> ov = new HashSet<>(tfc);
			ov.retainAll(Hac_proteins);
			if (ov.size() > 0) {
				n.put(dir_map.get(tfc), n.get(dir_map.get(tfc)) + 1);
			}
			if (tfc.contains("Q09472"))
				n_ep300.put(dir_map.get(tfc), n_ep300.get(dir_map.get(tfc)) + 1);
		}
		System.out.println(n + " " + n_ep300);
	}

}
