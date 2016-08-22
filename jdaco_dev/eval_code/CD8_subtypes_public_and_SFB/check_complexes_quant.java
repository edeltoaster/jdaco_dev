package CD8_subtypes_public_and_SFB;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import framework.DataQuery;
import framework.GOAnnotator;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/CD8_subsets_public/CD8_DACO_0.0/res7/";
	static String networks_folder_pre = "/Users/tho/Dropbox/Work/projects/CD8_subsets_public/CD8_networks_0.0/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	static GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/simple_tags_retrieved.txt.gz");
	
	public static double[] getDoubleArray(List<Double> list) {
		if (list.size() == 0)
			return new double[]{0.0};
		return list.stream().mapToDouble(d->d).toArray();
	}
	
	
	public static void N_MNP_TFcombinations() {
		System.out.println("N_MNP_TFcomb");
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			String cell_type = sample.split("_")[1];
			
			if (cell_type.equals("N") || cell_type.equals("TMNP"))
				results.put(sample, new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder_pre + cell_type + "/" + sample + "_major-transcripts.txt.gz"));
		}
		
		Set<HashSet<String>> TFvariants = new HashSet<>();
		for (QuantDACOResultSet qdr:results.values())
			TFvariants.addAll(qdr.getSeedToComplexMap().keySet());
		
		Map<HashSet<String>, LinkedList<Double>> N_TFV_abundance = new HashMap<>();
		Map<HashSet<String>, LinkedList<Double>> TMNP_TFV_abundance = new HashMap<>();
		
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			Map<HashSet<String>, Double> sample_abundances = results.get(sample).getAbundanceOfSeedVariantsComplexes();
			for (HashSet<String> TFvariant:TFvariants) {
				double abundance = sample_abundances.getOrDefault(TFvariant, 0.0);
				if (cell_type.equals("N")) {
					if (!N_TFV_abundance.containsKey(TFvariant))
						N_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					N_TFV_abundance.get(TFvariant).add(abundance);
				} 
				else { // TMNP
					if (!TMNP_TFV_abundance.containsKey(TFvariant))
						TMNP_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					TMNP_TFV_abundance.get(TFvariant).add(abundance);
				}
			}
		}
		
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(N_TFV_abundance.get(TFvariant)), getDoubleArray(TMNP_TFV_abundance.get(TFvariant)));
			
			if (pm < 0.05)
				System.out.println(DataQuery.batchHGNCProteinsGenes(TFvariant) + " : " + pm + "  -> " + N_TFV_abundance.get(TFvariant) + " vs " + TMNP_TFV_abundance.get(TFvariant));
		}
	}
	
	
	public static void MNP_ALL_TFcombinations() {
		System.out.println("MNP_vs_ALL_TFcomb");
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			String cell_type = sample.split("_")[1];
			
			results.put(sample, new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder_pre + cell_type + "/" + sample + "_major-transcripts.txt.gz"));
		}
		
		Set<HashSet<String>> TFvariants = new HashSet<>();
		for (QuantDACOResultSet qdr:results.values())
			TFvariants.addAll(qdr.getSeedToComplexMap().keySet());
		
		Map<HashSet<String>, LinkedList<Double>> other_TFV_abundance = new HashMap<>();
		Map<HashSet<String>, LinkedList<Double>> TMNP_TFV_abundance = new HashMap<>();
		
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			Map<HashSet<String>, Double> sample_abundances = results.get(sample).getAbundanceOfSeedVariantsComplexes();
			for (HashSet<String> TFvariant:TFvariants) {
				double abundance = sample_abundances.getOrDefault(TFvariant, 0.0);
				if (cell_type.equals("TMNP")) {
					if (!TMNP_TFV_abundance.containsKey(TFvariant))
						TMNP_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					TMNP_TFV_abundance.get(TFvariant).add(abundance);

				} 
				else { // all others
					if (!other_TFV_abundance.containsKey(TFvariant))
						other_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					other_TFV_abundance.get(TFvariant).add(abundance);
				}
			}
		}
		
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(TMNP_TFV_abundance.get(TFvariant)), getDoubleArray(other_TFV_abundance.get(TFvariant)));
			
			if (pm < 0.05) {
				List<Set<String>> complexes = new LinkedList<>();
				for (String sample:results.keySet()) {
					String cell_type2 = sample.split("_")[1];
					if (cell_type2.equals("TMNP")) 
						if (results.get(sample).getSeedToComplexMap().containsKey(TFvariant)) // could also be not found in MNP cells
							complexes.addAll(results.get(sample).getSeedToComplexMap().get(TFvariant));
				}
				System.out.println(DataQuery.batchHGNCProteinsGenes(TFvariant) + " : " + pm + "  -> " + goa.rateListsOfProteins(complexes) + " " + TMNP_TFV_abundance.get(TFvariant) + " vs " + other_TFV_abundance.get(TFvariant));
			}
		}
	}
	
	
	public static void MNP_ALL_TFcombinationsCleaned() {
		System.out.println("MNP_vs_ALL_TFcomb");
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			String cell_type = sample.split("_")[1];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder_pre + cell_type + "/" + sample + "_major-transcripts.txt.gz");
			System.out.println(sample + " " + qdr.getResult().size());
			qdr.removeOpposinglyAnnotatedComplexes(goa);
			System.out.println(sample + " " + qdr.getResult().size());
			results.put(sample, qdr);
		}
		
		Set<HashSet<String>> TFvariants = new HashSet<>();
		for (QuantDACOResultSet qdr:results.values())
			TFvariants.addAll(qdr.getSeedToComplexMap().keySet());
		
		Map<HashSet<String>, LinkedList<Double>> other_TFV_abundance = new HashMap<>();
		Map<HashSet<String>, LinkedList<Double>> TMNP_TFV_abundance = new HashMap<>();
		
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[1];
			Map<HashSet<String>, Double> sample_abundances = results.get(sample).getAbundanceOfSeedVariantsComplexes();
			for (HashSet<String> TFvariant:TFvariants) {
				double abundance = sample_abundances.getOrDefault(TFvariant, 0.0);
				if (cell_type.equals("TMNP")) {
					if (!TMNP_TFV_abundance.containsKey(TFvariant))
						TMNP_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					TMNP_TFV_abundance.get(TFvariant).add(abundance);

				} 
				else { // all others
					if (!other_TFV_abundance.containsKey(TFvariant))
						other_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					other_TFV_abundance.get(TFvariant).add(abundance);
				}
			}
		}
		
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(TMNP_TFV_abundance.get(TFvariant)), getDoubleArray(other_TFV_abundance.get(TFvariant)));
			
			if (pm < 0.05) {
				List<Set<String>> complexes = new LinkedList<>();
				for (String sample:results.keySet()) {
					String cell_type2 = sample.split("_")[1];
					if (cell_type2.equals("TMNP")) 
						if (results.get(sample).getSeedToComplexMap().containsKey(TFvariant)) // could also be not found in MNP cells
							complexes.addAll(results.get(sample).getSeedToComplexMap().get(TFvariant));
				}
				System.out.println(DataQuery.batchHGNCProteinsGenes(TFvariant) + " : " + pm + "  -> " + goa.rateListsOfProteins(complexes) + " " + TMNP_TFV_abundance.get(TFvariant) + " vs " + other_TFV_abundance.get(TFvariant));
			}
		}
	}
	
	
	public static void main(String[] args) {
		MNP_ALL_TFcombinationsCleaned();
		//MNP_ALL_TFcombinations();
		//System.out.println();
		//System.out.println();
		//N_MNP_TFcombinations();
	}
}
