package diff_compl_blood;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.DiffSeedVarDetector;
import framework.DiffSeedVarDetector.SPEnrichment;
import framework.GOAnnotator;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;
import framework.Utilities;


public class check_hemato_diffcomplexes {
	
	public static String daco_results_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/hemato_check/res_99_5/";
	public static String networks_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/hemato_check/BP_networks/";
	public static String diff_compl_output_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/hemato_check/diffnet_results/";
	public static int no_threads = 4;

	public static GOAnnotator goa = new GOAnnotator("9606", false, "mixed_data/stem_tags.txt");
		
	public static void check_diff_compl(String result_string, Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, boolean parametric) {
		System.out.println("Running " + result_string);
		
		String output_folder = diff_compl_output_folder + result_string + "/";
		
		// for annotations
		Map<String, String> proteins_of_interest = new HashMap<>(); // markers of interest?
//		increased.stream().forEach(p -> proteins_of_interest.put(p, "increased"));
//		expressed_effectors.stream().forEach(p -> proteins_of_interest.put(p, "expr_effectors"));
//		not_expressed.stream().forEach(p -> proteins_of_interest.put(p, "not_expr"));
		Set<String> interesting_targets = new HashSet<>();
//		interesting_targets.addAll(expressed_effectors);
//		interesting_targets.addAll(not_expressed);
		
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
		SPEnrichment tf_enrich = dcd.calculateSPEnrichment(definitions.qvalue, definitions.SPEnrich_iterations, definitions.SPEnrich_compl_part_threshold);
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
		
		// HSC -> MPP
		Map<String, QuantDACOResultSet> group1 = new HashMap<>();
		Map<String, QuantDACOResultSet> group2 = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("HSC"))
				group1.put(sample, qdr);
			else if (sample.contains("MPP")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("HSC samples : " + group1.size());
		System.out.println("MPP samples : " + group2.size());
		
		check_diff_compl("HSC_to_MPP_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// MPP -> CMP
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("MPP"))
				group1.put(sample, qdr);
			else if (sample.contains("CMP")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("MPP samples : " + group1.size());
		System.out.println("CMP samples : " + group2.size());
		
		check_diff_compl("MPP_to_CMP_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// MPP -> CLP
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("MPP"))
				group1.put(sample, qdr);
			else if (sample.contains("CLP")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("MPP samples : " + group1.size());
		System.out.println("CLP samples : " + group2.size());
		
		check_diff_compl("MPP_to_CLP_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// CMP -> MEP
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("CMP"))
				group1.put(sample, qdr);
			else if (sample.contains("MEP")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("CMP samples : " + group1.size());
		System.out.println("MEP samples : " + group2.size());
		
		check_diff_compl("CMP_to_MEP_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// CMP -> GMP
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("CMP"))
				group1.put(sample, qdr);
			else if (sample.contains("GMP")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("CMP samples : " + group1.size());
		System.out.println("GMP samples : " + group2.size());
		
		check_diff_compl("CMP_to_GMP_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// MEP -> EB
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("MEP"))
				group1.put(sample, qdr);
			else if (sample.contains("EB")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("MEP samples : " + group1.size());
		System.out.println("EB samples : " + group2.size());
		
		check_diff_compl("MEP_to_EB_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// MEP -> MK
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("MEP"))
				group1.put(sample, qdr);
			else if (sample.contains("MK")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("MEP samples : " + group1.size());
		System.out.println("MK samples : " + group2.size());
		
		check_diff_compl("MEP_to_MK_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// GMP -> N
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("GMP"))
				group1.put(sample, qdr);
			else if (sample.contains("N-")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("GMP samples : " + group1.size());
		System.out.println("N samples : " + group2.size());
		
		check_diff_compl("GMP_to_N_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// GMP -> M
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("GMP"))
				group1.put(sample, qdr);
			else if (sample.contains("M-v")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("GMP samples : " + group1.size());
		System.out.println("M samples : " + group2.size());
		
		check_diff_compl("GMP_to_M_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
		
		// CLP -> CD4
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.contains("CLP"))
				group1.put(sample, qdr);
			else if (sample.contains("CD4")){
				group2.put(sample, qdr);
			}
		}
		System.out.println("CLP samples : " + group1.size());
		System.out.println("CD4 samples : " + group2.size());
		
		check_diff_compl("CLP_to_CD4_parametric", group1, group2, true);
		group1.clear();
		group2.clear();
	}
}
