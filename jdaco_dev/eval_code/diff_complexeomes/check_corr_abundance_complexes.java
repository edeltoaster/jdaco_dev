package diff_complexeomes;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import diff_complexeomes.definitions;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class check_corr_abundance_complexes {
	
	public static Map<String, Double> count_complex_participation(QuantDACOResultSet qdr) {
		Map<String, Double> protein_complex_association = new HashMap<>();
		for (HashSet<String> complex: qdr.getResult()) {
			complex.stream().forEach(p->protein_complex_association.put(p, 1 + protein_complex_association.getOrDefault(p, 0.0)));
		}
		return protein_complex_association;
	}
	
	public static void main(String[] args) {
		PearsonsCorrelation pcorr = new PearsonsCorrelation();
		
		List<Double> corrs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			Map<String, Double> protein_complex_association = count_complex_participation(qdr);
			
			List<Double> counts = new LinkedList<>();
			List<Double> abundances = new LinkedList<>();
			for (String protein:protein_complex_association.keySet()) {
				double count = protein_complex_association.get(protein);
				double abundance = qdr.getProteinAbundance(protein);
				counts.add(count);
				abundances.add(abundance);
			}

			double corr = pcorr.correlation(Utilities.getDoubleArray(counts), Utilities.getDoubleArray(abundances));
			System.out.println(sample + " " + corr);
			corrs.add(corr);
		}
		System.out.println();
		System.out.println("overall corr: " + Utilities.getMean(corrs) + "+-" + Utilities.getStd(corrs));
	}
}
