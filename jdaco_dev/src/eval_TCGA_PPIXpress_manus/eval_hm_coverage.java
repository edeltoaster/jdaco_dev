package eval_TCGA_PPIXpress_manus;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.Utilities;

public class eval_hm_coverage {
	
	public static Map<String,Set<String>> getHMSet(String file) {
		
		Map<String, Set<String>> map = new HashMap<String, Set<String>>();
		for (String line:Utilities.readEntryFile(file)) {
			String term = line.split("\\t")[0];
			Set<String> proteins = new HashSet<String>();
			
			for (String p:line.trim().split("\\t")[1].split(","))
				proteins.add(p);
			map.put(term, proteins);
		}
		
		return map;
	}
	
	public static void process(Set<String> proteins, String analysis, NetworkBuilder nb) {
		
		System.out.println(analysis);
		System.out.println("Set-size: " + proteins.size());
		System.out.println("IAs: " + nb.getMappingPercentage(proteins));
		System.out.println("P.: " + nb.getMappingDomainPercentage(proteins));
	}
	
	public static void main(String[] args) {
		DataQuery.enforceSpecificEnsemblRelease("79");
		
		PPIN intact = new PPIN("mixed_data/human_intact.txt.gz");
		PPIN biogrid = new PPIN("mixed_data/human_biogrid.txt.gz");
		
		NetworkBuilder nb_intact = new NetworkBuilder(intact);
		System.out.println(nb_intact.getDB());
		NetworkBuilder nb_bg = new NetworkBuilder(biogrid);
		
		
		process(biogrid.getProteins(), "BioGRID", nb_bg);
		Map<String,Set<String>> hm_map = getHMSet("/Users/tho/Dropbox/manuscripts/PPIXpress/eval/hallmarks/hallmarks_biogrid_05_May_15_no_IEA.tsv");
		Set<String> all_hm = new HashSet<String>();
		for (String term:hm_map.keySet()) {
			process(hm_map.get(term), "BioGRID ("+term+")", nb_bg);
			all_hm.addAll(hm_map.get(term));
		}
		process(all_hm, "BioGRID (all HM)", nb_bg);
		Set<String> non_hm = new HashSet<String>(biogrid.getProteins());
		non_hm.removeAll(all_hm);
		process(non_hm, "BioGRID (non HM)", nb_bg);
		
		System.out.println("");
		
		process(intact.getProteins(), "IntAct", nb_intact);
		hm_map = getHMSet("/Users/tho/Dropbox/manuscripts/PPIXpress/eval/hallmarks/hallmarks_intact_05_May_15_no_IEA.tsv");
		all_hm.clear();
		for (String term:hm_map.keySet()) {
			process(hm_map.get(term), "IntAct ("+term+")", nb_intact);
			all_hm.addAll(hm_map.get(term));
		}
		process(all_hm, "IntAct (all HM)", nb_intact);
		non_hm = new HashSet<String>(intact.getProteins());
		non_hm.removeAll(all_hm);
		process(non_hm, "IntAct (non HM)", nb_intact);
	}
}
