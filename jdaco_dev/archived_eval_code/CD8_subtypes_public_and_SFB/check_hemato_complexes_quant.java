package CD8_subtypes_public_and_SFB;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.GOAnnotator;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_hemato_complexes_quant {
	
	static String daco_results_folder = "/Users/tho/GDrive/Work/projects/CD8_subtypes_public_and_SFB/hemato_DACO_0.0/res7/";
	static String networks_folder_pre = "/Users/tho/GDrive/Work/projects/CD8_subtypes_public_and_SFB/quant_hemo_networks_0.0/";
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
				else if (cell_type.equals("CD4")){
					if (!CD4_TFV_abundance.containsKey(TFvariant))
						CD4_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					CD4_TFV_abundance.get(TFvariant).add(abundance);
				}
			}
		}
		
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		TTest tt = new TTest();
		List<String> mwu_out = new LinkedList<>();
		List<String> tt_out = new LinkedList<>();
		System.out.println(TFvariants.size() + " tests.");
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(CLP_TFV_abundance.get(TFvariant)), getDoubleArray(CD4_TFV_abundance.get(TFvariant)));
			double pt = tt.tTest(getDoubleArray(CLP_TFV_abundance.get(TFvariant)), getDoubleArray(CD4_TFV_abundance.get(TFvariant)));
			
			String sig_m = "-";
			if (pm <0.05)
				sig_m = "+";
			
			String sig_t = "-";
			if (pt <0.05)
				sig_t = "+";
			
			mwu_out.add(sig_m + " " + DataQuery.batchHGNCNamesFromProteins(TFvariant) + " : " + pm + "  -> " + CLP_TFV_abundance.get(TFvariant) + " vs " + CD4_TFV_abundance.get(TFvariant));
			tt_out.add(sig_t + " " + DataQuery.batchHGNCNamesFromProteins(TFvariant) + " : " + pt + "  -> " + CLP_TFV_abundance.get(TFvariant) + " vs " + CD4_TFV_abundance.get(TFvariant));
		}
		
		Utilities.writeEntries(mwu_out, "/Users/tho/Desktop/mwu_out.txt");
		Utilities.writeEntries(tt_out, "/Users/tho/Desktop/tt_out.txt");
	}
	
	public static void hemato_all_vs_all() {
		System.out.println("all_TFcomb");
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		Set<String> cell_types = new HashSet<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			String cell_type = sample.split("_")[0];
			cell_types.add(cell_type);
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder_pre + sample + "_major-transcripts.txt.gz");
			qdr.removeOpposinglyAnnotatedComplexes(goa);
			results.put(sample, qdr);
		}
		
		Set<HashSet<String>> TFvariants = new HashSet<>();
		for (QuantDACOResultSet qdr:results.values())
			TFvariants.addAll(qdr.getSeedToComplexMap().keySet());
		
		// prune down to keep binding data smaller
		Set<String> TFC_left = new HashSet<>();
		for (HashSet<String> TFC:TFvariants)
				TFC_left.addAll(TFC);
		
		// only among TFs
		System.out.println("reading binding data ...");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/GDrive/Work/data_general/binding_sites/hocomoco_v10/hocomoco_v10_EPD_2.5k.txt.gz", TFC_left, TFC_left);
		System.out.println("done.");
		
		for (String test_cell_type:cell_types) {
			System.out.println("Checking " + test_cell_type);
			
			Map<HashSet<String>, LinkedList<Double>> test_TFV_abundance = new HashMap<>();
			Map<HashSet<String>, LinkedList<Double>> other_TFV_abundance = new HashMap<>();
			
			for (String sample:results.keySet()) {
				String cell_type = sample.split("_")[0];
				Map<HashSet<String>, Double> sample_abundances = results.get(sample).getAbundanceOfSeedVariantsComplexes();
				for (HashSet<String> TFvariant:TFvariants) {
					double abundance = sample_abundances.getOrDefault(TFvariant, 0.0);
					if (cell_type.equals(test_cell_type)) {
						if (!test_TFV_abundance.containsKey(TFvariant))
							test_TFV_abundance.put(TFvariant, new LinkedList<Double>());
						test_TFV_abundance.get(TFvariant).add(abundance);
					} 
					else {
						if (!other_TFV_abundance.containsKey(TFvariant))
							other_TFV_abundance.put(TFvariant, new LinkedList<Double>());
						other_TFV_abundance.get(TFvariant).add(abundance);
					}
				}
			}
			
			// check significance
			MannWhitneyUTest mwu = new MannWhitneyUTest();
			List<String> out = new LinkedList<>();
			Map<HashSet<String>, Double> test_results = new HashMap<>();
			for (HashSet<String> TFvariant:TFvariants) {
				double pm = mwu.mannWhitneyUTest(getDoubleArray(test_TFV_abundance.get(TFvariant)), getDoubleArray(other_TFV_abundance.get(TFvariant)));
				test_results.put(TFvariant, pm);
			}
			
			Map<String, String> effect = new HashMap<>();
			Map<HashSet<String>, Double> adj_test_results = Utilities.convertRawPValuesToBHFDR(test_results, 0.05);
			for (HashSet<String> TFvariant:adj_test_results.keySet()) {
				String direction = "+";
				double test_median = Utilities.getMedian(test_TFV_abundance.get(TFvariant));
				double other_median = Utilities.getMedian(other_TFV_abundance.get(TFvariant));
				if (test_median < other_median)
					direction = "-";
				List<Set<String>> complexes = new LinkedList<>();
				for (String sample:results.keySet()) {
					String cell_type2 = sample.split("_")[0];
					if (cell_type2.equals(test_cell_type)) 
						if (results.get(sample).getSeedToComplexMap().containsKey(TFvariant))
							complexes.addAll(results.get(sample).getSeedToComplexMap().get(TFvariant));
				}
				
				out.add(direction + " " + DataQuery.batchHGNCNamesFromProteins(TFvariant) + " : " + goa.rateCollectionOfProteins(complexes) + ", "+ adj_test_results.get(TFvariant) + "  -> " + test_median + " vs " + other_median + " -> " + test_TFV_abundance.get(TFvariant) + " vs " + other_TFV_abundance.get(TFvariant));
				effect.put(TFvariant.toString(), goa.rateCollectionOfProteins(complexes));
			}		
			
			// adding annotational data
			Map<String, Map<String,String>> annotational_data = new HashMap<>();
			annotational_data.put("Regulatory_effect", effect);
			
			if (out.size() > 0) {
				Utilities.writeEntries(out, "/Users/tho/Desktop/" + test_cell_type + "_textout.txt");
				RegulatoryNetwork regnet = new RegulatoryNetwork(adj_test_results.keySet(), bdh);
				regnet.writeRegulatoryNetwork("/Users/tho/Desktop/" + test_cell_type + "_netout.txt", 2);
				regnet.writeNodeTable("/Users/tho/Desktop/" + test_cell_type + "_nodeout.txt", annotational_data);
			}
		}
		
	}
	
	public static void main(String[] args) {
		//CLP_CD4_TFcombinations();
		hemato_all_vs_all();
	}
}
