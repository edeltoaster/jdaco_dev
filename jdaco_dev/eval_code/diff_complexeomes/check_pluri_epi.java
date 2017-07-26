package diff_complexeomes;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.Utilities;

public class check_pluri_epi {
	public static String pluri_epi_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/pluri_epi_check/";
	public static String binding_data = "/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz";
	static int iterations = 10000;
	
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
		Set<String> TFs_seen = new HashSet<>();
		for (String line:Utilities.readFile(pluri_epi_folder + "scripts/to_check.txt")) {
			String[] spl = line.trim().split("\\s+");
			Set<String> tf_combinations = new HashSet<>(Arrays.asList(spl[0].split("/")));
			TFs_seen.addAll(tf_combinations);
			List<String> modifications = Arrays.asList(spl[1].split(","));
			to_check.put(tf_combinations, modifications);
		}
		
		// check
		System.out.println("Read binding data for " + TFs_seen.size() + " TFs.");
		BindingDataHandler bnd = new BindingDataHandler(binding_data, TFs_seen);
		List<String> all_targets = new ArrayList<>(bnd.getTargetsToTFsMap().keySet());
		for (Set<String> tf_combinations:to_check.keySet()) {
			Set<String> targets = bnd.getAdjacencyPossibilities(tf_combinations, definitions.d_min, definitions.d_max, false);
			System.out.println(tf_combinations);
			for (String modification_string:to_check.get(tf_combinations)) {
				Set<String> overlap_set = new HashSet<>(targets);
				for (String modification:modification_string.split("\\|")) {
					String sign = modification.substring(modification.length()-1);
					String mark = modification.substring(0, modification.length()-1);
					Set<String> mark_covered = covered_promoters.get(mark);
					
					if (sign.equals("+"))
						overlap_set.retainAll(mark_covered);
					else
						overlap_set.removeAll(mark_covered);
				}
				
				double ref_overlap_fraction = overlap_set.size() / (double) targets.size();
				
				// see if by chance
				double count = 0.0;
				for (int i = 0; i < iterations; i++) {
					for (String modification:modification_string.split("\\|")) {
						String sign = modification.substring(modification.length()-1);
						String mark = modification.substring(0, modification.length()-1);
						Set<String> mark_covered = covered_promoters.get(mark);
						Collections.shuffle(all_targets);
						overlap_set = new HashSet<>(all_targets.subList(0, targets.size()));
						if (sign.equals("+"))
							overlap_set.retainAll(mark_covered);
						else
							overlap_set.removeAll(mark_covered);
					}
					double overlap_fraction = overlap_set.size() / (double) targets.size();
					if (overlap_fraction > ref_overlap_fraction)
						count += 1;
				}
				
				double p = Math.max(1.0 / iterations, count / iterations);
				System.out.println(modification_string + " " + ref_overlap_fraction + " " + p);
			}
			
		}
	}
}
