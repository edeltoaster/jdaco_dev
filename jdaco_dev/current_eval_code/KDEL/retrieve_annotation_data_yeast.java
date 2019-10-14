package KDEL;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class retrieve_annotation_data_yeast {

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
		
		
		// read Björn's yeast GO terms of interest
		Map<String, Set<String>> annotation_map = new HashMap<String, Set<String>>();
		String yeast = "4932";
		String out_file = "/Users/tho/Desktop/yeast_annotation_data_exp.tsv";

		Set<String> temp = DataQuery.getGenesWithGO("GO:0015031", yeast);
		System.out.println("protein transport : " + temp.size());
		annotation_map.put("protein transport", temp);
		temp = DataQuery.getGenesWithGO("GO:0006886", yeast);
		System.out.println("intracellular protein transport : " + temp.size());
		annotation_map.put("intracellular protein transport", temp);
		
		writeOut(out_file, annotation_map);
	}

}
