package CD8_subtypes_public;


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
import framework.GOAnnotator;
import framework.Utilities;


public class check_complexes {
	
	static String results_folder = "/Users/tho/Dropbox/Work/projects/CD8_subsets_public/CD8_DACO_0.0/res7/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	
	public static void all_TFcombinations() {
		System.out.println("all cell types, TF combinations only");
		Map<String, DACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			//String cell_type = sample.split("_")[1];
			
			//if (cell_type.equals("N") || cell_type.equals("TMNP"))
			results.put(sample, new DACOResultSet(f.getAbsolutePath(), seed));
		}
		
//		for (String sample1:results.keySet())
//			for (String sample2:results.keySet()) {
//				System.out.println(sample1 + " / " + sample2 + " : " + results.get(sample1).getComplexSimilarity(results.get(sample2)) + " " + results.get(sample1).getSeedVariantSimilarity(results.get(sample2)));
//			}
		
		// gather what is found in all samples of a cell type
		Map<String, Set<Set<String>>> exclusive_TFComb = new HashMap<>();
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			if (exclusive_TFComb.containsKey(cell_type))
				exclusive_TFComb.get(cell_type).retainAll(results.get(sample).getGeneralSeedToComplexMap().keySet());
			else
				exclusive_TFComb.put(cell_type, new HashSet<Set<String>>(results.get(sample).getGeneralSeedToComplexMap().keySet()));
		}
		
		// remove what is found in any other cell type
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			for (String cell_type2:exclusive_TFComb.keySet()) {
				if (cell_type.equals(cell_type2))
					continue;
				exclusive_TFComb.get(cell_type2).removeAll(results.get(sample).getGeneralSeedToComplexMap().keySet());
			}
		}
		
		// prune down to keep binding data smaller
		Set<String> TFC_left = new HashSet<>();
		for (Set<Set<String>> excl_TFCs:exclusive_TFComb.values()) {
			for (Set<String> TFC: excl_TFCs)
				TFC_left.addAll(TFC);
		}
		
		// check exclusive complexes
		GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/simple_tags_retrieved.txt.gz");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10/hocomoco_v10_EPD_2.5k.txt.gz", 1.0, TFC_left);
		for (String cell_type:exclusive_TFComb.keySet()) {
			System.out.println(cell_type + " " + exclusive_TFComb.get(cell_type).size());
			
			for (Set<String> excl_TFC:exclusive_TFComb.get(cell_type)) {
				
				// gather relevant TF combinations from all samples
				List<Set<String>> complexes = new LinkedList<>();
				for (String sample:results.keySet()) {
					String cell_type2 = sample.split("_")[1];
					if (cell_type2.equals(cell_type)) 
						complexes.addAll(results.get(sample).getSeedToComplexMap().get(excl_TFC));
				}
				
				// annotate with function
				List<List<String>> hgnc_complexes = new LinkedList<>();
				for (Set<String> compl:complexes)
					hgnc_complexes.add(DataQuery.batchHGNCProteinsGenes(compl));
				System.out.println(DataQuery.batchHGNCProteinsGenes(excl_TFC) + " -> " +  hgnc_complexes.size() + " : " + goa.rateListsOfProteins(complexes));
				
				// annotate targets
				Set<String> targets = bdh.getAdjacencyPossibilities(excl_TFC, -50, 50, false);
				System.out.println("shared targets: " + bdh.getCommonTargets(excl_TFC).size() + " restr. targets: " +targets.size());
				//Utilities.writeEntries(targets, "/Users/tho/Desktop/allcomb_" + cell_type + "_" + String.join("-", DataQuery.batchHGNCProteinsGenes(excl_TFC)) + "_targets.txt");
			}
		}
	}
	
	
	public static void N_MNP_TFcombinations() {
		System.out.println("N_MNP only, TF combinations only");
		Map<String, DACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			String cell_type = sample.split("_")[1];
			
			if (cell_type.equals("N") || cell_type.equals("TMNP"))
				results.put(sample, new DACOResultSet(f.getAbsolutePath(), seed));
		}
		
		// gather what is found in all samples of a cell type
		Map<String, Set<Set<String>>> exclusive_TFComb = new HashMap<>();
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			if (exclusive_TFComb.containsKey(cell_type))
				exclusive_TFComb.get(cell_type).retainAll(results.get(sample).getGeneralSeedToComplexMap().keySet());
			else
				exclusive_TFComb.put(cell_type, new HashSet<Set<String>>(results.get(sample).getGeneralSeedToComplexMap().keySet()));
		}
		
		// remove what is found in any other cell type
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			for (String cell_type2:exclusive_TFComb.keySet()) {
				if (cell_type.equals(cell_type2))
					continue;
				exclusive_TFComb.get(cell_type2).removeAll(results.get(sample).getGeneralSeedToComplexMap().keySet());
			}
		}
		
		// prune down to keep binding data smaller
		Set<String> TFC_left = new HashSet<>();
		for (Set<Set<String>> excl_TFCs:exclusive_TFComb.values()) {
			for (Set<String> TFC: excl_TFCs)
				TFC_left.addAll(TFC);
		}
		
		// check exclusive complexes
		GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/simple_tags_retrieved.txt.gz");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10/hocomoco_v10_EPD_2.5k.txt.gz", 1.0, TFC_left);
		for (String cell_type:exclusive_TFComb.keySet()) {
			System.out.println(cell_type + " " + exclusive_TFComb.get(cell_type).size());
			
			for (Set<String> excl_TFC:exclusive_TFComb.get(cell_type)) {
				
				// gather relevant TF combinations from all samples
				List<Set<String>> complexes = new LinkedList<>();
				for (String sample:results.keySet()) {
					String cell_type2 = sample.split("_")[1];
					if (cell_type2.equals(cell_type)) 
						complexes.addAll(results.get(sample).getSeedToComplexMap().get(excl_TFC));
				}
				
				// annotate with function
				List<List<String>> hgnc_complexes = new LinkedList<>();
				for (Set<String> compl:complexes)
					hgnc_complexes.add(DataQuery.batchHGNCProteinsGenes(compl));
				System.out.println(DataQuery.batchHGNCProteinsGenes(excl_TFC) + " -> " +  hgnc_complexes.size() + " : " + goa.rateListsOfProteins(complexes));
				
				// annotate targets
				Set<String> targets = bdh.getAdjacencyPossibilities(excl_TFC, -50, 50, false);
				System.out.println("shared targets: " + bdh.getCommonTargets(excl_TFC).size() + " restr. targets: " +targets.size());
				//Utilities.writeEntries(targets, "/Users/tho/Desktop/NMNP_" + cell_type + "_" + String.join("-", DataQuery.batchHGNCProteinsGenes(excl_TFC)) + "_targets.txt");
			}
		}
	}
	
	public static void main(String[] args) {
		N_MNP_TFcombinations();
		System.out.println();
		all_TFcombinations();
	}
}
