package sc_bonn;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;

public class checkENC_all {
	
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
		
		Map<String, DACOResultSet> data_map = new HashMap<String, DACOResultSet>();
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/ENCODE/DACO/ana_10/results10/", ".csv")) {
			String sample = file.getName().split("\\.")[0];
			DACOResultSet res = new DACOResultSet(file.getAbsolutePath(), seed);
			data_map.put(sample, res);
		}

		Set<String> relevant_targets = Utilities.readEntryFile("mixed_data/H1_acc_promotors_2k.txt.gz");
		// only data of relevant tfs and reachable of those
		
		// get sets found in all H1ESC
		HashSet<HashSet<String>> ESC_tfc = new HashSet<HashSet<String>>(data_map.get("CT_H1_1").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_2").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_3").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CT_H1_4").getSeedToComplexMap().keySet());
		ESC_tfc.retainAll(data_map.get("CS_H1").getSeedToComplexMap().keySet());
		
		// read binding data and restrict to reachable ones in DNAse data
		System.out.println("bdh reading");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/binding_sites/human_fimo_2k.txt.gz", seed, 0.00001, relevant_targets);
		System.out.println("TFs: " + bdh.getTFsWithBindingData().size());
		System.out.println("Possible targets: " + bdh.getTargetsToTFsMap().keySet().size());
		
		// first build for all
		System.out.println("regnet building (all hESC)");
		RegulatoryNetwork regnet = new RegulatoryNetwork(ESC_tfc, bdh, -50, 50, 3, 1);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_all.txt");
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_all_min2.txt", 2);
		regnet.writeNodeTable("/Users/tho/Desktop/node_table_all.txt");
		System.out.println("sink proteins: " + regnet.getSinkProteins().size());
		Utilities.writeEntries(regnet.getSinkProteins(),"/Users/tho/Desktop/sink_all.txt");
		System.out.println("sink proteins (2): " + regnet.getSinkProteins(2).size());
		Utilities.writeEntries(regnet.getSinkProteins(2),"/Users/tho/Desktop/sink_all_min2.txt");
		
		// restrict to the ones only found in H1 and write-out
		for (String s:data_map.keySet()) {
			if (!s.startsWith("BM"))
				continue;
			ESC_tfc.removeAll(data_map.get(s).getSeedToComplexMap().keySet());
		}
		
		System.out.println("regnet building (only hESC)");
		regnet = new RegulatoryNetwork(ESC_tfc, bdh, -50, 50, 3, 1);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only.txt");
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/regnet_only_min2.txt", 2);
		regnet.writeNodeTable("/Users/tho/Desktop/node_table_only.txt");
		System.out.println("sink proteins: " + regnet.getSinkProteins().size());
		Utilities.writeEntries(regnet.getSinkProteins(),"/Users/tho/Desktop/sink_only.txt");
		System.out.println("sink proteins (2): " + regnet.getSinkProteins(2).size());
		Utilities.writeEntries(regnet.getSinkProteins(2),"/Users/tho/Desktop/sink_only_min2.txt");
		
	}
}
