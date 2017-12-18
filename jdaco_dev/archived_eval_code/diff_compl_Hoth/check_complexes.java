package diff_compl_Hoth;


import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;


public class check_complexes {

	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> c1 = new HashMap<>();
		Map<String, QuantDACOResultSet> c2 = new HashMap<>();
		
		Map<String, String> sample_condition_map = new HashMap<>();
		Map<String, String> sample_cell_map = new HashMap<>();
		for (String line:Utilities.readFile("sample_association.tsv")) {
			String[] spl = line.trim().split("\\s+");
			sample_condition_map.put(spl[3], spl[1]);
			sample_cell_map.put(spl[3], spl[2]);
		}
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String sample_number = f.getName().split("_")[2];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			if (sample_condition_map.get(sample_number).equals("CTRL"))
				c1.put(sample_number, qdr);
			else if (sample_condition_map.get(sample_number).equals("CoCult"))
				c2.put(sample_number, qdr);
		}
		
		System.out.println("Determine diff. complexome");
		DiffComplexDetector dcd = new DiffComplexDetector(c1, c2, definitions.qvalue, definitions.parametric, definitions.paired, definitions.check_supersets, definitions.min_variant_fraction, definitions.no_threads);
		dcd.writeSignSortedComplexes("CTRL_CoCult_sign.txt", false);
		dcd.writeSignSortedComplexes("CTRL_CoCult_signh.txt", true);
		dcd.writeSignSortedVariants("CTRL_CoCult_tfcs_sign.txt", false);
		dcd.writeSignSortedVariants("CTRL_CoCult_tfcs_signh.txt", true);
		
		System.out.println("Determine seed protein comb. enrichment");
		dcd.calculateSPCEnrichment(definitions.qvalue, definitions.SPCEnrich_iterations, definitions.SPCEnrich_compl_part_threshold).writeSignificantSeedProteinCombinations("CTRL_CoCult_SPC_sign.txt");
		
		System.out.println("Determine seed protein enrichment");
		dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold).writeSignificantSeedProteins("CTRL_CoCult_SP_sign.txt");
	}
}
