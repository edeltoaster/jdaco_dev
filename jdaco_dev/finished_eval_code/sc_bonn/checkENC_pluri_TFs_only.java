package sc_bonn;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DACOResultSet;
import framework.DataQuery;
import framework.Utilities;

public class checkENC_pluri_TFs_only {
	
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
		
		System.out.println("all H1 only:" + ESC_tfc.size());
		
		Set<String> all_needed = new HashSet<String>();
		HashSet<HashSet<String>> pluri_tfcs = new HashSet<HashSet<String>>();
		for (HashSet<String> c:ESC_tfc) {
			if (c.contains("Q01860") || c.contains("P48431") || c.contains("Q9H9S0")) {
				pluri_tfcs.add(c);
				all_needed.addAll(c);
			}
		}
		
		System.out.println("pluri: " + pluri_tfcs.size());
		List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
		Map<String, String> up_hgnc = new HashMap<String, String>();
		for (String[] s:HGNC_map) {
			up_hgnc.put(s[1], s[0]);
		}
		
		// add everything
		System.out.println("");
		Map<HashSet<String>, Set<HashSet<String>>>tfc_compl = new HashMap<HashSet<String>, Set<HashSet<String>>>();
		for (HashSet<String> tfc:ESC_tfc) {
			tfc_compl.put(tfc, new HashSet<HashSet<String>>());
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
			System.out.println("TFC: " + temp);
			for (String sample:new String[]{"CT_H1_1", "CT_H1_2", "CT_H1_3", "CT_H1_4", "CS_H1"}) {
				tfc_compl.get(tfc).addAll(data_map.get(sample).getSeedToComplexMap().get(tfc));
			}
			
			for (HashSet<String> c:tfc_compl.get(tfc)) {
				temp = "";
				first = true;
				for (String s:c) {
				if (first) {
					temp += up_hgnc.get(s);				
					first = false;
				}
				else 
					temp += "," + up_hgnc.get(s);
				}
				
				System.out.println("   " + temp);
			}
		}
		
		Set<String> relevant_targets = Utilities.readEntryFile("mixed_data/H1_acc_promotors_2k.txt.gz");
		// only data of relevant tfs and reachable of those
		System.out.println("DNAse allows for " + relevant_targets.size());
		
		// read binding data and restrict to reachable ones in DNAse data
		System.out.println("bdh reading");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/binding_sites/human_fimo_2k.txt.gz", all_needed, 0.0001, relevant_targets);//NANOG relatively "weak"
		System.out.println("reachable targets: " + bdh.getTargetsToTFsMap().keySet().size());
		System.out.println("Settings: 5k, DNAse targets, adj/med like before; rough_adj: -5-15");
		for (HashSet<String> tfc:pluri_tfcs) {
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
			
			System.out.println("Processing: " + temp);
			int common = bdh.getCommonTargets(tfc).size();
			int adj = bdh.getAdjacencyPossibilities(tfc, 0, 10, false).size();
			int rough_adj = bdh.getAdjacencyPossibilities(tfc, -5, 15, false).size();
			int mediated = bdh.getAdjacencyPossibilities(tfc, 10, 50, false).size();
			int coloc = bdh.getAdjacencyPossibilities(tfc, -50, 50, false).size();
			System.out.println("    Common: " + common);
			System.out.println("    Adj: " + adj);
			System.out.println("    R. adj: " + rough_adj);
			System.out.println("    Med: " + mediated);
			System.out.println("    Coloc: " + coloc);
		}
	}
}
