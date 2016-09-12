package CD8_subtypes_public_and_SFB;


import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.GOAnnotator;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/CD8_DACO_0.0/res7/";
	static String networks_folder_pre = "/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/CD8_networks_0.0/";
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
		
		//increased: Plaur, Ilirn, Anpep, Ccl3, Mmp19, Serpb2, Cxcl3, Il1b and Sppi -> Q03405, P18510, P15144, P10147, Q99542, P05120, P19876, P01584, P10451
		Set<String> increased = new HashSet<>(Arrays.asList("Q03405","P18510", "P15144", "P10147", "Q99542", "P05120", "P19876", "P01584", "P10451"));
		// expressed effectors: Ifng, Gzmb, Il1b and Cxcl3 -> P01579, P10144, P01584, P19876
		Set<String> expressed_effectors = new HashSet<>(Arrays.asList("P01579", "P10144", "P01584", "P19876"));
		// not expressed: Nkg7, Fasl and Il22 -> Q16617, P48023, Q9GZX6
		Set<String> not_expressed = new HashSet<>(Arrays.asList("Q16617", "P48023", "Q9GZX6"));
		
		Map<String, Set<String>> annotated_targets = new HashMap<>();
		annotated_targets.put("increased", increased);
		annotated_targets.put("effectors", expressed_effectors);
		annotated_targets.put("not_expressed", not_expressed);
		Set<String> of_interest = new HashSet<>(increased);
		of_interest.addAll(expressed_effectors);
		of_interest.addAll(not_expressed);
		annotated_targets.put("of_interest", of_interest);
		
		System.out.println("MNP_vs_ALL_TFcomb");
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			String cell_type = sample.split("_")[1];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder_pre + cell_type + "/" + sample + "_major-transcripts.txt.gz");
			qdr.removeOpposinglyAnnotatedComplexes(goa);
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
		
		Set<String> allow_in_binding_data = new HashSet<>();
		allow_in_binding_data.addAll(increased);
		allow_in_binding_data.addAll(expressed_effectors);
		allow_in_binding_data.addAll(not_expressed);
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		Set<HashSet<String>> signTFvariants = new HashSet<>();
		Map<String, String> effect = new HashMap<>();
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(TMNP_TFV_abundance.get(TFvariant)), getDoubleArray(other_TFV_abundance.get(TFvariant)));
			
			if (pm < 0.05) {
				List<Set<String>> complexes = new LinkedList<>();
				signTFvariants.add(TFvariant);
				for (String sample:results.keySet()) {
					String cell_type2 = sample.split("_")[1];
					if (cell_type2.equals("TMNP")) 
						if (results.get(sample).getSeedToComplexMap().containsKey(TFvariant)) { // could also be not found in MNP cells
							complexes.addAll(results.get(sample).getSeedToComplexMap().get(TFvariant));
							allow_in_binding_data.addAll(TFvariant);
						}
				}
				
				String direction = "+";
				double test_median = Utilities.getMedian(TMNP_TFV_abundance.get(TFvariant));
				double other_median = Utilities.getMedian(other_TFV_abundance.get(TFvariant));
				if (test_median < other_median)
					direction = "-";
				
				System.out.println(direction + " " + DataQuery.batchHGNCProteinsGenes(TFvariant) + " : " + pm + "  -> " + goa.rateListsOfProteins(complexes) + 
						" -> " + test_median + " vs " + other_median + " , " + 
						complexes.stream().map(s->DataQuery.batchHGNCProteinsGenes(s)).collect(Collectors.toList()) + " , " + 
						TMNP_TFV_abundance.get(TFvariant) + " vs " + other_TFV_abundance.get(TFvariant));
				effect.put(TFvariant.toString(), goa.rateListsOfProteins(complexes));
			}
		}
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Regulatory_effect", effect);
		
		System.out.println("reading binding data ...");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10/hocomoco_v10_EPD_2.5k.txt.gz", allow_in_binding_data, allow_in_binding_data);
		System.out.println("done.");
		
		RegulatoryNetwork regnet = new RegulatoryNetwork(signTFvariants, bdh);
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/" + "TMNP_netout.txt", 2);
		regnet.writeNodeTable("/Users/tho/Desktop/" + "TMNP_newnodeout.txt", annotational_data);
		
	}
	
	
	public static void main(String[] args) {
		MNP_ALL_TFcombinationsCleaned();
		//MNP_ALL_TFcombinations();
		//System.out.println();
		//System.out.println();
		//N_MNP_TFcombinations();
	}
}
