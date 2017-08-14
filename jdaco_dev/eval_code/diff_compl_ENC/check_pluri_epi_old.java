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

public class check_pluri_epi_old {
	public static String pluri_epi_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/pluri_epi_check/";
	public static String binding_data = "/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz";
	public static String nodetable_file = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/diff_results_99_5_-25-25/pluri_nodetable_pruned.txt";
	
	public static void main(String[] args) {
		
		// TODO: determine abundances, quantiles ...
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
		Map<Set<String>, String> actualcompl_meanabun_map = new HashMap<>();
		Map<String, String> meanabun_map = new HashMap<>();
		Map<Set<String>, String> actualcompl_sampleabun_map = new HashMap<>();
		Set<String> TFs_seen = new HashSet<>();
		for (String line:Utilities.readFile(nodetable_file)) {
			if (line.startsWith("Node"))
				continue;
			String[] spl = line.trim().split("\\s+");
			// Node Nodetype Gene Mean_abundances Epi_effect_details Epi_effect Actual_complexes
			String epi_effect = spl[5];
			
			if (epi_effect.equals("/"))
				continue;

			Set<String> tf_combinations = new HashSet<>(Arrays.asList(spl[0].split("/")));
			TFs_seen.addAll(tf_combinations);
			
			List<String> modifications = Arrays.asList(epi_effect.split(","));
			to_check.put(tf_combinations, modifications);
			actualcompl_epi_map.put(tf_combinations, spl[4]);
			String mean_abundances = spl[3];
			actualcompl_meanabun_map.put(tf_combinations, mean_abundances);
			actualcompl_sampleabun_map.put(tf_combinations, spl[6]);
			
			// parse mean abundances to gain a map in the form of complex->complex:abundance
			Arrays.asList(mean_abundances.split(",")).stream().forEach(s -> meanabun_map.put(s.split(":")[0], s));
		}
		
		// check
		System.out.println("Read binding data for " + TFs_seen.size() + " TFs.");
		BindingDataHandler bnd = new BindingDataHandler(binding_data, TFs_seen);
		List<String> all_targets = new ArrayList<>(bnd.getTargetsToTFsMap().keySet());
		for (Set<String> tf_combinations:to_check.keySet()) {
			Set<String> targets = bnd.getAdjacencyPossibilities(tf_combinations, pluri_definitions.d_min, pluri_definitions.d_max, false);
			
			// pre-outputs
			System.out.println(tf_combinations + " <-> " + DataQuery.batchHGNCNamesFromProteins(tf_combinations));
			System.out.println("socc: " + actualcompl_sampleabun_map.get(tf_combinations));
			System.out.println("mabun: " + actualcompl_meanabun_map.get(tf_combinations));
			String epi_string = actualcompl_epi_map.get(tf_combinations);
			System.out.println("epi: " + epi_string);
			String epi_abun_string = String.join(",", Arrays.asList(epi_string.split(",")).stream().map(s -> meanabun_map.get(s.split(":")[0])).collect(Collectors.toList()));
			System.out.println("epi_ab: " + epi_abun_string);
			
			for (String modification_string:to_check.get(tf_combinations)) {
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
