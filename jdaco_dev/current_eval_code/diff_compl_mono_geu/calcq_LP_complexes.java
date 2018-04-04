package diff_compl_mono_geu;


import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.TimeUnit;

import framework.QuantDACOResultSet;
import framework.Utilities;


public class calcq_LP_complexes {

	public static void main(String[] args) {
		definitions.printInitParameters();
	
		System.out.println();

		Map<String, QuantDACOResultSet> data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			
			QuantDACOResultSet qdr = new QuantDACOResultSet(f.getAbsolutePath(), definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");
			
			if (sample.startsWith("HG") || sample.startsWith("non") || sample.startsWith("class"))
				data.put(sample, qdr);
		}
		
		for (Entry<String, QuantDACOResultSet> e:data.entrySet()) {
			System.out.println("calculating " + e.getKey());
			long start = System.currentTimeMillis();
			LP_algo lp = new LP_algo(e.getValue());
			long end = System.currentTimeMillis();
			
			long duration_LP = TimeUnit.MILLISECONDS.toSeconds(end - start);
			
			start = System.currentTimeMillis();
			e.getValue().getAbundanceOfComplexes();
			end = System.currentTimeMillis();
			
			long duration_approx = TimeUnit.MILLISECONDS.toSeconds(end - start);
			System.out.println(duration_LP + " " + duration_approx);
			
			Map<HashSet<String>, Double> abundances = lp.getAbundanceOfComplexes();
			List<String> to_write = new LinkedList<>();
			
			for (HashSet<String> cluster : abundances.keySet())
				to_write.add( String.join(",", cluster) + " " + abundances.get(cluster));
			
			Utilities.writeEntries(to_write, "lp_q_results_95_5/" + e.getValue() + ".txt.gz");
		}
		
	}
}
