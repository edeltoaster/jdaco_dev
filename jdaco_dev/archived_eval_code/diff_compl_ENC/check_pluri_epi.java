package diff_compl_ENC;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.HypergeometricDistribution;

import framework.BindingDataHandler;
import framework.DataQuery;
import framework.Utilities;

public class check_pluri_epi {
	public static String pluri_epi_folder = "/Users/tho/GDrive/Work/projects/stem_cell_complexome/pluri_epi_check/";
	public static String binding_data = "/Users/tho/GDrive/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz";
	public static String nodetable_file = "/Users/tho/Desktop/pluri/diffcompl_results_99_5_-25-25/nodetable_pruned.txt";
	
	public static void main(String[] args) {
		
		// load ENCODE ChIP-seq data on histone marks in H1
		Map<String, Set<String>> covered_promoters = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(pluri_epi_folder + "covered_promoters/", ".bed.gz")) {
			String modification = f.getName().split("\\.")[0];
			Set<String> covered = new HashSet<>();
			for (String line:Utilities.readFile(f.getAbsolutePath()))
				covered.add(line.split("\\s+")[3]);
			covered_promoters.put(modification, covered);
		}
		
		// load and parse data to check
		Map<Set<String>, List<String>> to_check = new HashMap<>();
		Map<Set<String>, String> actualcompl_epi_map = new HashMap<>();
		Map<Set<String>, String> directions_map = new HashMap<>();
		
		Set<String> TFs_seen = new HashSet<>();
		for (String line:Utilities.readFile(nodetable_file)) {
			if (line.startsWith("Node"))
				continue;
			String[] spl = line.trim().split("\\s+");
			// Node Nodetype Gene Direction GO_details GO_overall Med_abundance_G1 Med_abundance_G2
			String overall_epi_effect = spl[5];
			
			if (overall_epi_effect.equals("/"))
				continue;

			Set<String> tf_combinations = new HashSet<>(Arrays.asList(spl[0].split("/")));
			TFs_seen.addAll(tf_combinations);
			
			List<String> modifications = Arrays.asList(overall_epi_effect.split(",")).stream().filter(m->m.length()>4 && !m.contains("!")).collect(Collectors.toList());
			to_check.put(tf_combinations, modifications);
			actualcompl_epi_map.put(tf_combinations, spl[4]);
			directions_map.put(tf_combinations, spl[3]);
		}
		
		// check
		System.out.println("Read binding data for " + TFs_seen.size() + " TFs.");
		BindingDataHandler bnd = new BindingDataHandler(binding_data, TFs_seen);
		List<String> all_targets = new ArrayList<>(bnd.getTargetsToTFsMap().keySet());
		for (Set<String> tf_combinations:to_check.keySet()) {
			Set<String> targets = bnd.getAdjacencyPossibilities(tf_combinations, pluri_definitions.d_min, pluri_definitions.d_max, false);
			
			// pre-outputs
			System.out.println(tf_combinations + " <-> " + DataQuery.batchHGNCNamesFromProteins(tf_combinations));
			String epi_string = actualcompl_epi_map.get(tf_combinations);
			System.out.println("epi: " + epi_string);
			Map<String, String> compl_epi = Arrays.asList(actualcompl_epi_map.get(tf_combinations).split(",")).stream().filter(s -> s.split(":").length == 2).map(s->s.split(":")).collect(Collectors.toMap(s -> s[0], s -> s[1]));
			Map<String, String> compl_dir = Arrays.asList(directions_map.get(tf_combinations).split(",")).stream().map(s->s.split(":")).collect(Collectors.toMap(s -> s[0], s -> s[1]));
			for (String modification_string:to_check.get(tf_combinations)) {
				System.out.println(modification_string);
				if (!modification_string.contains("!")) {
					String directions = String.join(",", compl_epi.keySet().stream().filter(c -> compl_epi.get(c).equals(modification_string)).map(c -> c + ":" + compl_dir.get(c)).collect(Collectors.toList()));
					System.out.println("Directions: " + directions);
				}
				Set<String> overlap_set = new HashSet<>(targets);
				Set<String> complete_set = new HashSet<>(all_targets);
				for (String modification:modification_string.split("\\|")) {
					String sign = modification.substring(modification.length()-1);
					String mark = modification.substring(0, modification.length()-1);
					Set<String> mark_covered = covered_promoters.get(mark);
					
					if (sign.equals("+")) {
						overlap_set.retainAll(mark_covered);
						complete_set.retainAll(mark_covered);
					}
					else {
						overlap_set.removeAll(mark_covered);
						complete_set.removeAll(mark_covered);
					}
				}
				// hypergeometric test
				HypergeometricDistribution hyper = new HypergeometricDistribution(all_targets.size(), complete_set.size(), targets.size());
				// calculates P(X >= overlap)
				double hyper_p = hyper.upperCumulativeProbability(overlap_set.size());
				double ref_overlap_fraction = overlap_set.size() / (double) targets.size();

				System.out.println(modification_string + " " + ref_overlap_fraction + "% -> " + hyper_p);
				
				if (modification_string.contains("|")) {
					for (String modification:modification_string.split("\\|")) {
						overlap_set = new HashSet<>(targets);
						complete_set = new HashSet<>(all_targets);
						String sign = modification.substring(modification.length()-1);
						String mark = modification.substring(0, modification.length()-1);
						Set<String> mark_covered = covered_promoters.get(mark);
						
						if (sign.equals("+")) {
							overlap_set.retainAll(mark_covered);
							complete_set.retainAll(mark_covered);
						}
						else {
							overlap_set.removeAll(mark_covered);
							complete_set.removeAll(mark_covered);
						}
						
						// hypergeometric test
						hyper = new HypergeometricDistribution(all_targets.size(), complete_set.size(), targets.size());
						// calculates P(X >= overlap)
						hyper_p = hyper.upperCumulativeProbability(overlap_set.size());
						ref_overlap_fraction = overlap_set.size() / (double) targets.size();

						System.out.println("  " + modification + " " + ref_overlap_fraction + "% -> " + hyper_p);
					}

				}		
			}
			System.out.println();
		}
	}
}
