package sc_bonn_TCGA;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DACOResultSet;
import framework.DataQuery;
import framework.RegulatoryNetwork;
import framework.Utilities;

public class checkENC_among_H1_exkl_only {
	
	public static double jaccard(Set<String> a, Set<String> b) {
		Set<String> intersection = new HashSet<String>(a);
		intersection.retainAll(b);
		Set<String> union = new HashSet<String>(a);
		union.addAll(b);
		return ((double)intersection.size()) / union.size();
	}
	
	public static double set_jaccard(DACOResultSet a, DACOResultSet b) {
		Set<HashSet<String>> tfc_a = a.getSeedToComplexMap().keySet();
		Set<HashSet<String>> tfc_b = b.getSeedToComplexMap().keySet();
		//Set<HashSet<String>> tfc_a = a.getResult();
		//Set<HashSet<String>> tfc_b = b.getResult();

		Set<HashSet<String>> intersection = new HashSet<HashSet<String>>(tfc_a);
		intersection.retainAll(tfc_b);
		
		double sum = 0;
		for (HashSet<String> set_a:tfc_a) {
			double max = 0;
			for (HashSet<String> set_b:tfc_b) {
				double temp = jaccard(set_a, set_b);
				if (temp > max)
					max = temp;
			}
			sum += max;
		}
		
		Set<HashSet<String>> union = new HashSet<HashSet<String>>(tfc_a);
		union.addAll(tfc_b);
		
		return sum / union.size();
	}
	
	public static void main(String[] args) {
		Set<String> seed = Utilities.readEntryFile("mixed_data/hocomoco_tfs.txt.gz");
		
		//System.out.println("sample #complexes #variants");
		Map<String, DACOResultSet> data_map = new HashMap<String, DACOResultSet>();
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/ENCODE/DACO/ana_10/results10/", ".csv")) {
			String sample = file.getName().split("\\.")[0];
			DACOResultSet res = new DACOResultSet(file.getAbsolutePath(), seed);
			data_map.put(sample, res);
			//System.out.println(sample + " " + res.getResult().size() + " " + res.getSeedToComplexMap().size());
		}
		
		// check similarities
//		List<String> lines = new LinkedList<String>();
//		for (String a:data_map.keySet()) {
//			for (String b:data_map.keySet()) {
//				if (a.compareTo(b) <= 0)
//					continue;
//				String t = a+ "/" +b + " : " + set_jaccard(data_map.get(a), data_map.get(b));
//				lines.add(t);
//			}
//		}
//		Utilities.writeEntries(lines, "/Users/tho/Documents/workspace/ENCODE_check/compl_similarities.txt");
//		System.exit(0);
		
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
		
		List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
		Map<String, String> up_hgnc = new HashMap<String, String>();
		for (String[] s:HGNC_map) {
			up_hgnc.put(s[1], s[0]);
		}
		
		List<String> to_write = new LinkedList<String>();
		for (HashSet<String> tfc:ESC_tfc) {
			String temp = "";
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
		
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/stem_tfcs_only.txt");
		
		//System.exit(0);
		
		Set<String> involved_tfs = new HashSet<String>();
		for (HashSet<String> s:ESC_tfc)
			involved_tfs.addAll(s);
		
		Utilities.writeEntries(involved_tfs, "/Users/tho/Desktop/stem_tfs_only.txt");
		
		// more
		//System.exit(0);
		
		/////////////////////////////////
		// GRN STUFF
		/////////////////////////////////
		
		
		Set<String> relevant_targets = Utilities.readEntryFile("mixed_data/H1_acc_promotors_2k.txt.gz");
		// only data of relevant tfs and reachable of those
		
		
		
		// get sets found in all H1ESC
		ESC_tfc = new HashSet<HashSet<String>>(data_map.get("CT_H1_1").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_2").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_3").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_4").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CS_H1").getSeedToComplexMap().keySet());
		
		to_write = new LinkedList<String>();
		for (HashSet<String> tfc:ESC_tfc) {
			String temp = "";
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
		
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/stem_tfcs_all.txt");
		
		
		// restrict to the ones only found in H1 and write-out
		for (String s:data_map.keySet()) {
			if (!s.startsWith("BM"))
				continue;
			ESC_tfc.removeAll(data_map.get(s).getSeedToComplexMap().keySet());
		}
	
		// read binding data and restrict to reachable ones in DNAse data
		System.out.println("bdh reading");
		involved_tfs.retainAll(relevant_targets);
		//BindingDataHandler bdh = new BindingDataHandler("mixed_data/human_fimo_2k.txt.gz", seed, 00001, seed);
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/binding_sites/human_fimo_2k.txt.gz", involved_tfs, 0.0001, involved_tfs);
		
		RegulatoryNetwork regnet = new RegulatoryNetwork(ESC_tfc, bdh, -50, 50, 3, 1);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only.txt");
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only_min2.txt", 2);
		regnet.writeNodeTable("/Users/tho/Desktop/node_table_only.txt");
		
	}
}
