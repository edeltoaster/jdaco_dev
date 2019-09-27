package KDEL;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class retrieve_annotation_data {

	static int min_annotated_genes = 4;
	static List<String> go_terms = new LinkedList<>();
	
	public static void writeOut(String file, Map<String, Set<String>> data) {
		List<String> lines = new LinkedList<String>();
		for (String term : go_terms) {
			if (data.get(term).size() < min_annotated_genes)
				continue;
			String temp = term + "\t"+ String.join(",", data.get(term));
			lines.add(temp);
		}
		Utilities.writeEntries(lines, file);
	}
	public static void main(String[] args) {
		// read Björn's GO terms of interest
		
		Map<String, String> go_tois = new HashMap<>();
		for (String s : Utilities.readFile("/Users/tho/GDrive/Work/projects/KDEL/new_stuff/enrichment/GO_TOI.csv")) {
			if (!s.contains(":"))
				continue;
			String[] spl = s.trim().split(";");
			go_tois.put(spl[0], spl[2]);
			go_terms.add(spl[2]);
		}
		
		Map<String, Set<String>> annotation_map = new HashMap<String, Set<String>>();
		String human = "9606";
		String out_file = "/Users/tho/Desktop/annotation_data.tsv";
		
		for (String go_term : go_tois.keySet() ) {
			String go_descr = go_tois.get(go_term);
			Set<String> temp = DataQuery.getGenesWithGO(go_term, human);
			System.out.println(go_descr + " : " + temp.size());
			annotation_map.put(go_descr, temp);
		}
		
		writeOut(out_file, annotation_map);
		
		
		System.out.println("Only experimental data:");
		annotation_map.clear();
		out_file = "/Users/tho/Desktop/annotation_data_exp.tsv";
		
		for (String go_term : go_tois.keySet() ) {
			String go_descr = go_tois.get(go_term);
			Set<String> temp = DataQuery.getGenesWithGO(go_term, human, false, true);
			System.out.println(go_descr + " : " + temp.size());
			annotation_map.put(go_descr, temp);
		}
		
		writeOut(out_file, annotation_map);
	}

}
