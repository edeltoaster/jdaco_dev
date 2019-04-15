package diff_compl_mono_geu;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.BindingDataHandler;
import framework.DACOResultSet;
import framework.DiffComplexDetector;
import framework.RegulatoryNetwork;
import framework.Utilities;

public class check_mono_complexes_marker_target_enr {
	static int no_threads = 4;
	static int no_iterations = 1000;
	static BindingDataHandler bdh;
	
	public static int getNumberOfComplexesTargetingMarkers(List<HashSet<String>> complexes) {
		RegulatoryNetwork regnet = new RegulatoryNetwork(complexes, definitions.seed, bdh, definitions.d_min, definitions.d_max, no_threads, 2);
		return regnet.getComplexToTargets().size();
	}
	
	public static void main(String[] args) {
		
		DiffComplexDetector.SignSortedComplexesResult sign_complexes = DiffComplexDetector.readSignSortedComplexResult(definitions.compl_results);
		int sign_complexes_size = sign_complexes.getSignificanceSortedComplexes().size();
		
		Set<String> relevant_targets = new HashSet<>(definitions.markers);
		
		System.out.println("Read binding data with all TFs but only markers as the targets...");
		bdh = new BindingDataHandler(definitions.binding_data, definitions.seed, relevant_targets);
		
		System.out.println("from " + sign_complexes_size + " sign.: " + getNumberOfComplexesTargetingMarkers(sign_complexes.getSignificanceSortedComplexes()));
		
		System.out.println("Reading all complexes");
		int ncms = 0;
		int cms = 0;
		Map<HashSet<String>, Integer> count_map_cm = new HashMap<>();
		Map<HashSet<String>, Integer> count_map_ncm = new HashMap<>();
		Set<HashSet<String>> unfiltered_complexes = new HashSet<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.local_daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			if (sample.startsWith("HG"))
				continue;
			
			DACOResultSet qdr = new DACOResultSet(f.getAbsolutePath(), definitions.seed);
			unfiltered_complexes.addAll(qdr.getResult());
			
			if (sample.startsWith("non")) {
				qdr.getResult().stream().forEach(c -> count_map_ncm.put(c, count_map_ncm.getOrDefault(c, 0) + 1));
				ncms++;
			}
			else {
				qdr.getResult().stream().forEach(c -> count_map_cm.put(c, count_map_cm.getOrDefault(c, 0) + 1));
				cms++;
			}
			
		}
		
		List<HashSet<String>> all_complexes = new ArrayList<HashSet<String>>(unfiltered_complexes);
		System.out.println("pre-pruning: " + all_complexes.size());
		double min_count_cms = definitions.min_variant_fraction * cms;
		double min_count_ncms = definitions.min_variant_fraction * ncms;
		for (HashSet<String> complex:unfiltered_complexes) {
			// filter
			if (count_map_ncm.getOrDefault(complex, 0) < min_count_ncms && count_map_cm.getOrDefault(complex, 0) < min_count_cms)
				all_complexes.remove(complex);
		}
		System.out.println("post-pruning: " + all_complexes.size());
		
		for (int i = 0; i < no_iterations; i++) {
			Collections.shuffle(all_complexes);
			List<HashSet<String>> test_complexes = new ArrayList<>(all_complexes.subList(0, sign_complexes_size));
			int hits = getNumberOfComplexesTargetingMarkers(test_complexes);
			System.out.println("from " + test_complexes.size() + " test.: " + hits);
		}
		
	}
}
