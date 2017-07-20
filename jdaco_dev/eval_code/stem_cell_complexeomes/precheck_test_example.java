package stem_cell_complexeomes;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class precheck_test_example {

	public static String expr_data = "/Users/tho/Desktop/quantified_samples/";
	
	public static void main(String[] args) {
		Set<String> relevant_tfs = new HashSet<>();
		relevant_tfs.add("O43474"); // KLF4
		relevant_tfs.add("Q01860"); // POU5F1
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		Map<String, List<String>> transcr_prot_map = new HashMap<>();
		for (String[] genes_transcr_prot:DataQuery.getGenesTranscriptsProteins(db)) {
			String transcript = genes_transcr_prot[1];
			String protein = genes_transcr_prot[2];
			if (!transcr_prot_map.containsKey(transcript))
				transcr_prot_map.put(transcript, new LinkedList<String>());
			transcr_prot_map.get(transcript).add(protein);
		}
		
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz", 0.0001, relevant_tfs);
		Set<String> targets = bdh.getAdjacencyPossibilities(relevant_tfs, definitions.d_min, definitions.d_max, false);
		Map<String, Map<String, Float>> protein_abundance_per_sample = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_data, ".tsv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			if (!sample.contains("H1-hESC") && !sample.contains("H7-hESC") && !sample.contains("induced-pluripotent-stem-cell"))
				continue;
			
			Map<String, Float> expr_per_prot = TranscriptAbundanceReader.readSample(f.getAbsolutePath(), -1, true, "");
			for (String transcr:expr_per_prot.keySet()) {
				float transcr_expr = expr_per_prot.get(transcr);
				if (!transcr_prot_map.containsKey(transcr))
					continue;
				for (String prot:transcr_prot_map.get(transcr))
					expr_per_prot.put(prot, expr_per_prot.getOrDefault(prot, 0f) + transcr_expr);
			}
			protein_abundance_per_sample.put(sample, expr_per_prot);
		}
		
		Set<String> all_proteins = new HashSet<>();
		protein_abundance_per_sample.values().forEach(s->all_proteins.addAll(s.keySet()));
		List<String> to_write = new LinkedList<>();
		to_write.add("protein expression target");
		for (String protein:all_proteins) {
			List<Double> all_data = new LinkedList<>();
			for (String sample:protein_abundance_per_sample.keySet()) {
				all_data.add((double)protein_abundance_per_sample.get(sample).getOrDefault(protein, 0f));
			}
			
			double expr = Utilities.getMean(all_data);
			String in_target = "no";
			if (targets.contains(protein))
				in_target = "yes";
			
			to_write.add(protein + " " + expr + " " + in_target); 
		}
		
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/oct_klf_test.txt");
		
	}
}
