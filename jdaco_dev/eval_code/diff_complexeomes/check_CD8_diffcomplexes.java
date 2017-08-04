package diff_complexeomes;


import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import diff_complexeomes.definitions;
import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffSeedVarDetector;
import framework.DiffSeedVarDetector.SPEnrichment;
import framework.GOAnnotator;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_CD8_diffcomplexes {
	
	public static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/CD8_subtypes/res_99_5/";
	public static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/CD8_subtypes/CD8_networks/";
	public static String diff_compl_output_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/CD8_subtypes/diffnet_results/";
	public static int no_threads = 4;
	
	//increased: Plaur, Ilirn, Anpep, Ccl3, Mmp19, Serpb2, Cxcl3, Il1b and Sppi -> Q03405, P18510, P15144, P10147, Q99542, P05120, P19876, P01584, P10451
	public static Set<String> increased = new HashSet<>(Arrays.asList("Q03405","P18510", "P15144", "P10147", "Q99542", "P05120", "P19876", "P01584", "P10451"));
	// expressed effectors: Ifng, Gzmb, Il1b and Cxcl3 -> P01579, P10144, P01584, P19876
	public static Set<String> expressed_effectors = new HashSet<>(Arrays.asList("P01579", "P10144", "P01584", "P19876"));
	// not expressed: Nkg7, Fasl and Il22 -> Q16617, P48023, Q9GZX6
	public static Set<String> not_expressed = new HashSet<>(Arrays.asList("Q16617", "P48023", "Q9GZX6"));
	
	public static GOAnnotator goa = new GOAnnotator("9606", false, "mixed_data/stem_tags.txt");
	
	
	public static void check_diff_compl(String result_string, Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, boolean parametric) {
		System.out.println("Running " + result_string);
		
		String output_folder = diff_compl_output_folder + result_string + "/";
		
		// for annotations
		Map<String, String> proteins_of_interest = new HashMap<>();
		increased.stream().forEach(p -> proteins_of_interest.put(p, "increased"));
		expressed_effectors.stream().forEach(p -> proteins_of_interest.put(p, "expr_effectors"));
		not_expressed.stream().forEach(p -> proteins_of_interest.put(p, "not_expr"));
		Set<String> interesting_targets = new HashSet<>(increased);
		interesting_targets.addAll(expressed_effectors);
		interesting_targets.addAll(not_expressed);
		
		System.out.println("Determine differential complexomes ...");
		Set<String> involved_tfs = new HashSet<>();
		DiffSeedVarDetector dcd = new DiffSeedVarDetector(group1, group2, definitions.qvalue, parametric, false, definitions.check_supersets, 0.0, no_threads);
		
		List<HashSet<String>> tf_variants = new LinkedList<>();
		Map<String, String> effect = new HashMap<>();
		Map<String, String> summarized_effect = new HashMap<>();
		Map<String, String> abundances = new HashMap<>();
		Map<String, String> actual_complexes = new HashMap<>();
		Map<String, String> directions = new HashMap<>();
		List<String> res_pos_all = new LinkedList<>();
		List<String> res_neg_all = new LinkedList<>();
		
		System.out.println(dcd.getNumberOfTests() + " TF variants tested.");
		System.out.println(dcd.getSignificanceSortedVariants().size() + " diff. TF variants.");
		
		if (dcd.getSignificanceSortedVariants().size() == 0) {
			System.out.println();
			return;
		}
		
		for (HashSet<String> variant:dcd.getSignificanceSortedVariants()) {
			
			String sign = dcd.getSignificantVariantsDirections().get(variant);
			
			String hgncs = DataQuery.batchHGNCNamesFromProteins(variant).toString();
			double pval = dcd.getSignificantVariantsQValues().get(variant);
			involved_tfs.addAll(variant);
			String out_string = sign + " " + hgncs + ", " + pval;
			
			// distinguish between increased/positive abundance and diminishing/negative abundance
			if (sign.equals("-")) {
				res_neg_all.add(out_string);
			} else {
				// everything in network
				res_pos_all.add(out_string);
			}
			
			tf_variants.add(variant);
			directions.put(variant.toString(), sign);
			
			String[] annotation_data = DiffSeedVarDetector.getSortedComplexesAnnotations(variant, sign, goa, group1, group2);
			effect.put(variant.toString(), annotation_data[1]);
			summarized_effect.put(variant.toString(), annotation_data[2]);
			actual_complexes.put(variant.toString(), annotation_data[0]);
			abundances.put(variant.toString(), annotation_data[3]);
		}
		
		System.out.println(res_pos_all.size() +"+, " + res_neg_all.size() + "-");
		
		new File(output_folder).mkdir();
		
		// write results
		Utilities.writeEntries(res_pos_all, output_folder + "res_pos_all.txt");
		Utilities.writeEntries(res_neg_all, output_folder + "res_neg_all.txt");
		
		interesting_targets.addAll(involved_tfs);
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz", involved_tfs, 0.0001, interesting_targets);
		
		/**
		 *  writing network data
		 */
		
		System.out.println("Building ...");
		RegulatoryNetwork pluri_regnet = new RegulatoryNetwork(tf_variants, bdh, definitions.d_min, definitions.d_max, no_threads, 1);
		System.out.println(pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(output_folder + "regnet.txt");
		Map<String, Map<String,String>> annotational_data = new HashMap<>();
		annotational_data.put("Epi_effect", summarized_effect);
		annotational_data.put("Epi_effect_details", effect);
		annotational_data.put("Actual_complexes", actual_complexes);
		annotational_data.put("Mean_abundances", abundances);
		annotational_data.put("T_context", proteins_of_interest);
		annotational_data.put("Direction", directions);
		pluri_regnet.writeNodeTable(output_folder + "nodetable.txt", annotational_data);
		// pruning
		pluri_regnet.pruneToLargestSCCs();
		System.out.println("SCC: " + pluri_regnet.getSizesStr());
		pluri_regnet.writeRegulatoryNetwork(output_folder + "regnet_pruned.txt");
		pluri_regnet.writeNodeTable(output_folder + "nodetable_pruned.txt", annotational_data);
		
		/**
		 * Playing around with seed protein enrichment
		 */
		
		System.out.println("Calculating TF enrichment ...");
		SPEnrichment tf_enrich = dcd.calculateTFEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
		List<String> pos_tf_enrich_out = new LinkedList<>();
		List<String> neg_tf_enrich_out = new LinkedList<>();
		for (String tf:tf_enrich.getSignificanceSortedSeedProteins()) {
			String dir = tf_enrich.getSignificantSeedProteinDirections().get(tf);
			String out = tf + " " + DataQuery.getHGNCNameFromProtein(tf) + " " + tf_enrich.getSignificantSeedProteinQvalues().get(tf);
			
			if (dir.equals("+"))
				pos_tf_enrich_out.add(out);
			else
				neg_tf_enrich_out.add(out);
		}
		Utilities.writeEntries(pos_tf_enrich_out, output_folder + "tf_enrich_pos.txt");
		Utilities.writeEntries(neg_tf_enrich_out, output_folder + "tf_enrich_neg.txt");
		System.out.println();
	}

	public static void main(String[] args) {
		definitions.printParameters();
		System.out.println("folders and more overwritten for CD8 stuff");
		System.out.println();
		
		new File(diff_compl_output_folder).mkdir();
		
		// TMNP_vs_all
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("TMNP"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other T samples : " + group1.size());
		System.out.println("T_MNP samples : " + group2.size());
		
		check_diff_compl("TMNP_vs_all_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// TMNP_vs_N
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("N_"))
				group1.put(sample, qdr);
			else if (sample.contains("TMNP")) {
				group2.put(sample, qdr);
			}
		}
		System.out.println("T_N samples : " + group1.size());
		System.out.println("T_MNP samples : " + group2.size());
		
		check_diff_compl("TMNP_vs_N_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// TMNP_vs_Eff
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("EM"))
				group1.put(sample, qdr);
			else if (sample.contains("TMNP")) {
				group2.put(sample, qdr);
			}
		}
		System.out.println("T_EM(RA) samples : " + group1.size());
		System.out.println("T_MNP samples : " + group2.size());
		
		check_diff_compl("TMNP_vs_Eff_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// TMNP_vs_Mem
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("CM"))
				group1.put(sample, qdr);
			else if (sample.contains("TMNP")) {
				group2.put(sample, qdr);
			}
		}
		System.out.println("T_CM samples : " + group1.size());
		System.out.println("T_MNP samples : " + group2.size());
		
		check_diff_compl("TMNP_vs_Mem_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// N_to_Eff
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("N_"))
				group1.put(sample, qdr);
			else if (sample.contains("EM")) {
				group2.put(sample, qdr);
			}
		}
		System.out.println("T_N samples : " + group1.size());
		System.out.println("T_EM(RA) samples : " + group2.size());
		
		check_diff_compl("N_to_Eff_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// N_to_Mem
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("N_"))
				group1.put(sample, qdr);
			else if (sample.contains("CM")) {
				group2.put(sample, qdr);
			}
		}
		System.out.println("T_N samples : " + group1.size());
		System.out.println("T_CM samples : " + group2.size());
		
		check_diff_compl("N_to_Mem_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// TMNP vs all non-para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("TMNP"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_TMNP samples : " + group2.size());
		
		check_diff_compl("TMNP_vs_all_non-parametric", group1, group2, false);
		group1.clear();
		group2.clear();
		
		// Eff vs all non-para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("EM"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_EM(RA) samples : " + group2.size());
		
		check_diff_compl("Eff_vs_all_non-parametric", group1, group2, false);
		group1.clear();
		group2.clear();
		
		// Mem vs all non-para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("CM"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_CM samples : " + group2.size());
		
		check_diff_compl("Mem_vs_all_non-parametric", group1, group2, false);
		group1.clear();
		group2.clear();
		
		// N vs all non-para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("N_"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_N samples : " + group2.size());
		
		check_diff_compl("N_vs_all_non-parametric", group1, group2, false);
		group1.clear();
		group2.clear();
		
		// Eff vs all para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("EM"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_EM(RA) samples : " + group2.size());
		
		check_diff_compl("Eff_vs_all_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// Mem vs all para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("CM"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_CM samples : " + group2.size());
		
		check_diff_compl("Mem_vs_all_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// N vs all para
		group1 = new HashMap<>();
		group2 = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (!sample.contains("N_"))
				group1.put(sample, qdr);
			else {
				group2.put(sample, qdr);
			}
		}
		System.out.println("all other samples : " + group1.size());
		System.out.println("T_N samples : " + group2.size());
		
		check_diff_compl("N_vs_all_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
	}
}
