package PPIComp_hemato;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.PPIN;
import framework.Utilities;

public class hematoGO {

	public static void writeOut(String file, Map<String, Set<String>> data) {
		List<String> lines = new LinkedList<String>();
		for (String term:data.keySet()) {
			String temp = term + "\t"+ String.join(",", data.get(term));
			lines.add(temp);
		}
		Utilities.writeEntries(lines, file);
	}
	
	public static void main(String[] args) {
		Set<String> ref_proteins = new PPIN("mixed_data/human_mentha_17_jan.txt.gz").getProteins();
		Map<String, Set<String>> hallmark_prot_map = new HashMap<String, Set<String>>();
		String human = "9606";
		Set<String> temp;
		String out_file = "/Users/tho/Desktop/hemato_proteins_30_05_16.tsv";
		
		temp = DataQuery.getProteinsWithGO("GO:0030097", human);
		temp.retainAll(ref_proteins);
		System.out.println(temp.size());
		hallmark_prot_map.put("Hemopoiesis", temp);
		
		writeOut(out_file, hallmark_prot_map);
	}

}
