package breastcancer_proteome_complexome;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPEnrichment;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class diff_complexes_rel_mentha {
	static String net_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/mentha_networks/";
	static String compl_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/res/";
	static String seed_file = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/JDACO_run_stuff/hocomoco_human_core_TFs_v11.txt.gz";
	static String qr_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/relquant_compl/";
	static String out_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/rel_log2_diffout/";
	
	static Map<String, QuantDACOResultSet> all_data = new HashMap<>();
	static Map<String, QuantDACOResultSet> Basal_data = new HashMap<>();
	static Map<String, QuantDACOResultSet> HER2_data = new HashMap<>();
	static Map<String, QuantDACOResultSet> LumA_data = new HashMap<>();
	static Map<String, QuantDACOResultSet> LumB_data = new HashMap<>();
	static Map<String, QuantDACOResultSet> Normal_data = new HashMap<>();
	
	public static Map<String, QuantDACOResultSet> getDataset(String subtype) {
		switch (subtype) {
			case "Basal":
				return Basal_data;
			case "HER2":
				return HER2_data;
			case "LumA":
				return LumA_data;
			case "LumB":
				return LumB_data;
			case "Normal":
				return Normal_data;
		}
		return Normal_data;
	}
	
	public static void load_data(boolean write) {
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(compl_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			System.out.println("Processing " + sample);
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed_file, net_folder + sample + "_major-transcripts_rel.txt.gz");
			
			qdr.convertToLog2();
			
			if (write)
				qdr.writeQuantifiedResult(qr_folder + sample + ".txt.gz");
			
			all_data.put(sample, qdr);
			
			getDataset(sample.split("_")[1]).put(sample, qdr);
		}
	}
	
	public static void main(String[] args) {
		load_data(false);
		
		List<String> subtypes = new LinkedList<>();
		subtypes.add("Basal");
		subtypes.add("HER2");
		subtypes.add("LumA");
		subtypes.add("LumB");
		subtypes.add("Normal");
		
		System.out.println();
		for (String dataset:subtypes) {
			System.out.println("Checking " + dataset);
			Map<String, QuantDACOResultSet> target_data = getDataset(dataset); 
			Map<String, QuantDACOResultSet> all_other_data = new HashMap<>(all_data);
			target_data.keySet().stream().forEach(s -> all_other_data.remove(s));
			System.out.println(target_data.size() + " to " + all_other_data.size());
			
			DiffComplexDetector dcd = new DiffComplexDetector(all_other_data, target_data, 0.05, false, false, false, 0.75, Runtime.getRuntime().availableProcessors());
			System.out.println(dcd.getRawPValues().size() + " complexes tested, " + dcd.getSignificanceSortedComplexes().size() + " significant.");
			
			if (dcd.getSignificanceSortedComplexes().size() > 0)
				dcd.writeSignSortedComplexes(out_folder + dataset + "_hr.txt", true);
			
			SPEnrichment spe = dcd.calculateSPEnrichment(0.05, 10000, 10);
			System.out.println(spe.getSignificanceSortedSeedProteins().size() + " enriched TFs.");
			
			if (spe.getSignificanceSortedSeedProteins().size() > 0)
				spe.writeSignificantSeedProteins(out_folder + dataset + "_spe.txt");
			
		}
	}
}
