package KDEL;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class retrieve_annotation_data_yeast {
	
	public static void writeOut(String file, Map<String, Set<String>> data) {
		List<String> lines = new LinkedList<String>();
		for (String term : data.keySet()) {
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
		annotation_map.put("protein transport", temp);
		
		temp = DataQuery.getGenesWithGO("GO:0006886", yeast);
		annotation_map.put("intracellular protein transport", temp);
		
		writeOut(out_file, annotation_map);
	}

}
