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


public class check_hemato_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/hemato_DACO_0.0/res7/";
	static String networks_folder_pre = "/Users/tho/Dropbox/Work/projects/CD8_subtypes_public_and_SFB/quant_hemato_networks_0.0/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	static GOAnnotator goa = new GOAnnotator("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/simple_tags_retrieved.txt.gz");
	
	public static double[] getDoubleArray(List<Double> list) {
		if (list.size() == 0)
			return new double[]{0.0};
		return list.stream().mapToDouble(d->d).toArray();
	}
	
	
	public static void CLP_CD4_TFcombinations() {
		System.out.println("CLP_CD4_TFcomb");
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder_pre + sample + "_major-transcripts.txt.gz");
			//qdr.removeOpposinglyAnnotatedComplexes(goa);
			results.put(sample, qdr);
		}
		
		Set<HashSet<String>> TFvariants = new HashSet<>();
		for (QuantDACOResultSet qdr:results.values())
			TFvariants.addAll(qdr.getSeedToComplexMap().keySet());
		
		Map<HashSet<String>, LinkedList<Double>> CLP_TFV_abundance = new HashMap<>();
		Map<HashSet<String>, LinkedList<Double>> CD4_TFV_abundance = new HashMap<>();
		
		for (String sample:results.keySet()) {
			String cell_type = sample.split("_")[0];
			Map<HashSet<String>, Double> sample_abundances = results.get(sample).getAbundanceOfSeedVariantsComplexes();
			for (HashSet<String> TFvariant:TFvariants) {
				double abundance = sample_abundances.getOrDefault(TFvariant, 0.0);
				if (cell_type.equals("CLP")) {
					if (!CLP_TFV_abundance.containsKey(TFvariant))
						CLP_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					CLP_TFV_abundance.get(TFvariant).add(abundance);
				} 
				else {
					if (!CD4_TFV_abundance.containsKey(TFvariant))
						CD4_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					CD4_TFV_abundance.get(TFvariant).add(abundance);
				}
			}
		}
		
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		System.out.println(TFvariants.size() + " tests.");
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(CLP_TFV_abundance.get(TFvariant)), getDoubleArray(CD4_TFV_abundance.get(TFvariant)));
			
			if (pm < 0.05) {
				System.out.println(DataQuery.batchHGNCProteinsGenes(TFvariant) + " : " + pm + "  -> " + CLP_TFV_abundance.get(TFvariant) + " vs " + CD4_TFV_abundance.get(TFvariant));
			}
		}
	}
	
	public static void main(String[] args) {
		CLP_CD4_TFcombinations();
	}
}
