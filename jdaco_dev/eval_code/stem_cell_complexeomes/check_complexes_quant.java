package stem_cell_complexeomes;


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
	
	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_results/preppi95_TPM/res7/";
	static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/preppi_TPM_networks/";
//	static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/DACO_results/mentha_TPM/res5/";
//	static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/mentha_TPM_networks/";
	static Set<String> seed = Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/hocomoco_human_TFs_v10.txt.gz");
	static GOAnnotator goa = new GOAnnotator("9606", false, "/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/stem_tags.txt");
	
	public static double[] getDoubleArray(List<Double> list) {
		if (list.size() == 0)
			return new double[]{0.0};
		return list.stream().mapToDouble(d->d).toArray();
	}
	
	public static void main(String[] args) {
		Map<String, QuantDACOResultSet> results = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), seed, networks_folder + sample + "_major-transcripts.txt.gz");
			System.out.println(sample + " " + qdr.getResult().size());
			qdr.removeOpposinglyAnnotatedComplexes(goa);
			System.out.println(qdr.getResult().size());
			results.put(sample, qdr);
		}
		
		Set<HashSet<String>> TFvariants = new HashSet<>();
		for (QuantDACOResultSet qdr:results.values())
			TFvariants.addAll(qdr.getSeedToComplexMap().keySet());
		
		// prune down to keep binding data smaller
		Set<String> TFC_left = new HashSet<>();
		for (HashSet<String> TFC:TFvariants)
				TFC_left.addAll(TFC);
		
		// collect all abundance values
		Map<HashSet<String>, LinkedList<Double>> hESC_TFV_abundance = new HashMap<>();
		Map<HashSet<String>, LinkedList<Double>> BM_TFV_abundance = new HashMap<>();
		for (String sample:results.keySet()) {
			Map<HashSet<String>, Double> sample_abundances = results.get(sample).getAbundanceOfSeedVariantsComplexes();
			for (HashSet<String> TFvariant:TFvariants) {
				double abundance = sample_abundances.getOrDefault(TFvariant, 0.0);
				
				if (sample.startsWith("BM")) {
					if (!BM_TFV_abundance.containsKey(TFvariant))
						BM_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					BM_TFV_abundance.get(TFvariant).add(abundance);
				}
				else {
					if (!hESC_TFV_abundance.containsKey(TFvariant))
						hESC_TFV_abundance.put(TFvariant, new LinkedList<Double>());
					hESC_TFV_abundance.get(TFvariant).add(abundance);
				}
			}
		}
		
		// test
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		List<String> pos_out = new LinkedList<>();
		List<String> neg_out = new LinkedList<>();
		System.out.println(TFvariants.size() + " tests.");
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> TFvariant:TFvariants) {
			double pm = mwu.mannWhitneyUTest(getDoubleArray(hESC_TFV_abundance.get(TFvariant)), getDoubleArray(BM_TFV_abundance.get(TFvariant)));
			test_results.put(TFvariant, pm);
		}
		
		// mult. hypo. adjust and report
		Map<HashSet<String>, Double> adj_test_results = Utilities.convertRawPValuesToBHFDR(test_results, 0.05);
		
		for (HashSet<String> TFvariant:adj_test_results.keySet()) {
			double hESC_median = Utilities.getMedian(hESC_TFV_abundance.get(TFvariant));
			double BM_median = Utilities.getMedian(BM_TFV_abundance.get(TFvariant));
			
			String sign = "+";
			if (hESC_median < BM_median)
				sign = "-";
			
			Set<Set<String>> complexes = new HashSet<>();
			for (String sample:results.keySet()) {
				if (sample.contains("hESC")) 
					if (results.get(sample).getSeedToComplexMap().containsKey(TFvariant))
						complexes.addAll(results.get(sample).getSeedToComplexMap().get(TFvariant));
			}
			
			String out = sign + " " + DataQuery.batchHGNCProteinsGenes(TFvariant) + " : " + adj_test_results.get(TFvariant) + " (" + test_results.get(TFvariant) + ")  -> " + hESC_median + " vs " + BM_median + " , " + goa.rateCollectionOfProteins(complexes) + ", " + TFvariant + ", "+ hESC_TFV_abundance.get(TFvariant) + " vs " + BM_TFV_abundance.get(TFvariant);
			
			if (sign.equals("+"))
				pos_out.add(out);
			else
				neg_out.add(out);
		
		}
		
		Utilities.writeEntries(pos_out, "/Users/tho/Desktop/stem_excl.txt");
		Utilities.writeEntries(neg_out, "/Users/tho/Desktop/non_stem.txt");
	}
}
