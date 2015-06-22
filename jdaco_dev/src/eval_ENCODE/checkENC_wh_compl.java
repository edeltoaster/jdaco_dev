package eval_ENCODE;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DACOResultSet;
import framework.DataQuery;
import framework.Utilities;

public class checkENC_wh_compl {
	
	public static void main(String[] args) {
		Set<String> seed = Utilities.readEntryFile("mixed_data/hocomoco_tfs.txt");
		
		//System.out.println("sample #complexes #variants");
		Map<String, DACOResultSet> data_map = new HashMap<String, DACOResultSet>();
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/ENCODE/DACO/ana_10/results10/", ".csv")) {
			String sample = file.getName().split("\\.")[0];
			DACOResultSet res = new DACOResultSet(file.getAbsolutePath(), seed);
			data_map.put(sample, res);
			//System.out.println(sample + " " + res.getResult().size() + " " + res.getSeedToComplexMap().size());
		}
		
		
		// get sets found in all H1ESC
		HashSet<HashSet<String>> ESC_tfc = new HashSet<HashSet<String>>(data_map.get("CT_H1_1").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_2").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_3").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_4").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CS_H1").getSeedToComplexMap().keySet());
		
		System.out.println("in all H1:" + ESC_tfc.size());
		
		for (String s:data_map.keySet()) {
			if (!s.startsWith("BM"))
				continue;
			ESC_tfc.removeAll(data_map.get(s).getSeedToComplexMap().keySet());
		}
		
		/**
		// filter single TFs
		List<HashSet<String>> to_del = new LinkedList<HashSet<String>>();
		for (HashSet<String> s:ESC_tfc) {
			if (s.size() == 1)
				to_del.add(s);
		}
		ESC_tfc.removeAll(to_del);
		**/
		System.out.println("all H1 only:" + ESC_tfc.size());
		
		Map<HashSet<String>, HashSet<String>> TF_to_all = new HashMap<HashSet<String>, HashSet<String>>();
		
		for (HashSet<String> tfv:ESC_tfc) {
			if (!TF_to_all.containsKey(tfv))
				TF_to_all.put(tfv, new HashSet<String>());
			for (String sample:new String[]{"CT_H1_1", "CT_H1_2", "CT_H1_3", "CT_H1_4", "CS_H1"}) {
				for (HashSet<String> complex:data_map.get(sample).getSeedToComplexMap().get(tfv)) {
					TF_to_all.get(tfv).addAll(complex);
				}
			}
		}
		
		// writeout
		
		List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
		Map<String, String> up_hgnc = new HashMap<String, String>();
		for (String[] s:HGNC_map) {
			up_hgnc.put(s[1], s[0]);
		}
		
		List<String> to_write = new LinkedList<String>();
		int i = 1;
		for (HashSet<String> tfc:ESC_tfc) {
			String temp = (i++) + ": ";
			boolean first = true;
			for (String s:tfc) {
				if (first) {
					temp += up_hgnc.get(s);
					first = false;
				}
				else 
					temp += "," + up_hgnc.get(s);
			}
			to_write.add(temp);
		}
		
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/complex_map.txt");
		
		new File("/Users/tho/Desktop/complexes/").mkdir();
		i = 1;
		for (HashSet<String> tfc:ESC_tfc) {
			Utilities.writeEntries(TF_to_all.get(tfc), "/Users/tho/Desktop/complexes/" + (i++) + ".txt");
		}
	}
}
