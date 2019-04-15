package KDEL;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class retrieve_annotation_data {

	public static void writeOut(String file, Map<String, Set<String>> data) {
		List<String> lines = new LinkedList<String>();
		for (String term:data.keySet()) {
			String temp = term + "\t"+ String.join(",", data.get(term));
			lines.add(temp);
		}
		Utilities.writeEntries(lines, file);
	}
	
	public static void main(String[] args) {
		Map<String, Set<String>> annotation_map = new HashMap<String, Set<String>>();
		String human = "9606";
		Set<String> temp;
		String out_file = "/Users/tho/Desktop/annotation_data.tsv";
		
		temp = DataQuery.getProteinsWithGO("GO:0045095", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0031424", human));
		System.out.println("keratinizaion:" + temp.size());
		annotation_map.put("keratinization", temp);
		
		temp = DataQuery.getProteinsWithGO("GO:0005840", human);
		System.out.println("ribosome:" + temp.size());
		annotation_map.put("ribosome", temp);
		
		temp = DataQuery.getProteinsWithGO("GO:0031982", human);
		System.out.println("vesicle:" + temp.size());
		annotation_map.put("vesicle", temp);
		
		temp = DataQuery.getProteinsWithGO("GO:0005829", human);
		System.out.println("cytosol:" + temp.size());
		annotation_map.put("cytosol", temp);
		
		temp = DataQuery.getProteinsWithGO("GO:0005634", human);
		System.out.println("nucleus:" + temp.size());
		annotation_map.put("nucleus", temp);
		
		writeOut(out_file, annotation_map);
	}

}
