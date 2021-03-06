package PPIXpress_BRCA;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class long_transcr_per_prot {
	static Map<String, LinkedList<String>> transcript_to_proteins = new HashMap<String, LinkedList<String>>();
	static Map<String, List<String>> transcr_dom = new HashMap<String, List<String>>();
	static Map<String, String> prot_isoform_map;
			
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
	
	public static double[] test(Map<String, Float> tumor, Map<String, Float> normal, Set<String> proteins, int i) {
		
		Map<String, LinkedList<Integer>> isoforms_tumor = new HashMap<String, LinkedList<Integer>>();
		for (String transcript:tumor.keySet()) {
			if (!transcript_to_proteins.containsKey(transcript)) // if no protein coding transcript
				continue;
			for (String protein:transcript_to_proteins.get(transcript)) {
				if (proteins != null && !proteins.contains(protein))
					continue;
				if (!isoforms_tumor.containsKey(protein)) {
					isoforms_tumor.put(protein, new LinkedList<Integer>());
				}
				
				Set<String> here = new HashSet<String>(transcr_dom.getOrDefault(transcript, new LinkedList<String>()));
				Set<String> longest = new HashSet<String>(transcr_dom.getOrDefault(prot_isoform_map.get(protein), new LinkedList<String>()));
				
				if (here.equals(longest))
					isoforms_tumor.get(protein).add(1);
				else
					isoforms_tumor.get(protein).add(0);
			}
		}
		
		double t_p_t = 
			(isoforms_tumor.values().stream().mapToDouble(
					s -> s.stream().mapToInt(t->t).sum() / (double) s.size() ).sum()) 
				/ isoforms_tumor.keySet().size();
		
		
		Map<String, LinkedList<Integer>> isoforms_normal = new HashMap<String, LinkedList<Integer>>();
		for (String transcript:normal.keySet()) {
			if (!transcript_to_proteins.containsKey(transcript)) // if no protein coding transcript
				continue;
			for (String protein:transcript_to_proteins.get(transcript)) {
				if (proteins != null && !proteins.contains(protein))
					continue;
				if (!isoforms_normal.containsKey(protein)) {
					isoforms_normal.put(protein, new LinkedList<Integer>());
				}
				
				Set<String> here = new HashSet<String>(transcr_dom.getOrDefault(transcript, new LinkedList<String>()));
				Set<String> longest = new HashSet<String>(transcr_dom.getOrDefault(prot_isoform_map.get(protein), new LinkedList<String>()));
				
				if (here.equals(longest))
					isoforms_normal.get(protein).add(1);
				else
					isoforms_normal.get(protein).add(0);
			}
		}
		double t_p_n = 
				(isoforms_normal.values().stream().mapToDouble(
						s -> s.stream().mapToInt(t->t).sum() / (double) s.size() ).sum()) 
					/ isoforms_normal.keySet().size();
		
		return new double[]{t_p_t, t_p_n};
	}
	
	public static void process(Set<String> proteins, String analysis) {
		
		System.out.println(analysis);
		if (proteins != null)
			System.out.println("Set-size: " + proteins.size());
		
		// split data
		List<String> tumor_data = new ArrayList<String>();
		List<String> normal_data = new ArrayList<String>();
		
		for (File file:Utilities.getAllPrefixMatchingFilesInSubfolders("/Users/tho/GDrive/Work/tissue_spec/cancer_nets/TCGA/", "BRCA")) {
			String file_name = file.getAbsolutePath();
			if (file_name.endsWith("tumor.txt.gz"))
				tumor_data.add(file_name);
			else
				normal_data.add(file_name);
		}
		
		List<Double> val_t = new LinkedList<Double>();
		List<Double> val_n = new LinkedList<Double>();
		List<Double> val = new LinkedList<Double>();
		
		for (int i = 0; i < tumor_data.size(); i++) {
			
			String tumor_file = tumor_data.get(i);
			String normal_file = normal_data.get(i);
			
			double[] bla = test(TranscriptAbundanceReader.readTCGAIsoformRSEM(tumor_file, 1), TranscriptAbundanceReader.readTCGAIsoformRSEM(normal_file, 1), proteins, i);
			val_t.add(bla[0]);
			val_n.add(bla[1]);
			val.add(bla[0]);
			val.add(bla[1]);
		}
		
		System.out.println("T: " + Utilities.getMean(val_t) + " +- " + Utilities.getStd(val_t));
		System.out.println("N: " + Utilities.getMean(val_n) + " +- " + Utilities.getStd(val_n));
		System.out.println("M: " + Utilities.getMean(val) + " +- " + Utilities.getStd(val));
		System.out.println("d: " + (Utilities.getMean(val_t)-Utilities.getMean(val_n)));
	
		System.out.println("");
	}
	
	public static void main(String[] args) {
		System.out.println("Fraction of expr transcripts per protein with same domain comp as longest isoform:");
		
		transcr_dom = DataQuery.getTranscriptsDomains("homo_sapiens_core_79_38");
		prot_isoform_map = DataQuery.getIsoformTranscriptsOfProteins("homo_sapiens_core_79_38");
		
		for (String[] naming:DataQuery.getGenesTranscriptsProteins("homo_sapiens_core_79_38")) {
			String transcript = naming[1];
			String protein = naming[2];
			if (!transcript_to_proteins.containsKey(transcript))
				transcript_to_proteins.put(transcript, new LinkedList<String>());
			transcript_to_proteins.get(transcript).add(protein);
		}
		
		// get all prefixes
		Set<String> cancers = new HashSet<String>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/GDrive/Work/tissue_spec/cancer_nets/TCGA/", ".txt.gz")) {
			cancers.add(f.getName().split("_")[0]);
		}
		

		
		PPIN intact = new PPIN("mixed_data/human_intact.txt.gz");
		PPIN biogrid = new PPIN("mixed_data/human_biogrid.txt.gz");
		
		
		process(biogrid.getProteins(), "BioGRID");
		Map<String,Set<String>> hm_map = getHMSet("/Users/tho/GDrive/manuscripts/PPIXpress/eval/hallmarks/hallmarks_biogrid_05_May_15_no_IEA.tsv");
		Set<String> all_hm = new HashSet<String>();
		for (String term:hm_map.keySet()) {
			process(hm_map.get(term), "BioGRID ("+term+")");
			all_hm.addAll(hm_map.get(term));
		}
		process(all_hm, "BioGRID (all HM)");
		Set<String> non_hm = new HashSet<String>(biogrid.getProteins());
		non_hm.removeAll(all_hm);
		process(non_hm, "BioGRID (non HM)");
		
		System.out.println("");
		
		process(intact.getProteins(), "IntAct");
		hm_map = getHMSet("/Users/tho/GDrive/manuscripts/PPIXpress/eval/hallmarks/hallmarks_intact_05_May_15_no_IEA.tsv");
		all_hm.clear();
		for (String term:hm_map.keySet()) {
			process(hm_map.get(term), "IntAct ("+term+")");
			all_hm.addAll(hm_map.get(term));
		}
		process(all_hm, "IntAct (all HM)");
		non_hm = new HashSet<String>(intact.getProteins());
		non_hm.removeAll(all_hm);
		process(non_hm, "IntAct (non HM)");
	}

}
