package eval_TCGA;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.PPIN;
import framework.Utilities;

public class locGO {

	public static void writeOut(String file, Map<String, Set<String>> data, Set<String> ref_proteins) {
		List<String> lines = new LinkedList<String>();
		for (String protein:ref_proteins) {
			String temp = protein + "\t";
			for (String loc:data.keySet()) {
				if (data.get(loc).contains(protein))
					temp += loc;
			}
			
			if (temp.equals(protein + "\t"))
				temp += "none";
			
			lines.add(temp);
		}
		Utilities.writeEntries(lines, file);
	}

	public static void main(String[] args) {
		Set<String> ref_proteins = new PPIN("mixed_data/human_merged.tsv.gz").getProteins();
		Map<String, Set<String>> prot_map = new HashMap<String, Set<String>>();
		String human = "9606";
		Set<String> temp;
		String out_file = "/Users/tho/Desktop/loc_GO_040815.tsv";
		
		temp = DataQuery.getProteinsWithGO("GO:0005886", human);// plasma membrane
		temp.retainAll(ref_proteins);
		prot_map.put("membrane", temp);
		
		temp = DataQuery.getProteinsWithGO("GO:GO:0005829", human);
		temp.retainAll(ref_proteins);
		prot_map.put("cytosol", temp);
		
		temp = DataQuery.getProteinsWithGO("GO:0005634", human);
		temp.retainAll(ref_proteins);
		prot_map.put("nucleus", temp);
		
		
		writeOut(out_file, prot_map, ref_proteins);
	}

}
