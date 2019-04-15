package diff_compl_mono_geu;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import framework.QuantDACOResultSet;
import framework.Utilities;

public class calcq_LP_complexes_single {

	public static void main(String[] args) {

		String sample = args[0];
		QuantDACOResultSet qdr = new QuantDACOResultSet(definitions.daco_results_folder + sample + ".csv.gz", definitions.seed, definitions.networks_folder + sample + "_major-transcripts.txt.gz");

		
		System.out.println("calculating " + sample);
		long start = System.currentTimeMillis();
		LP_algo lp = new LP_algo(qdr);
		long end = System.currentTimeMillis();
		
		long duration_LP = TimeUnit.MILLISECONDS.toSeconds(end - start);
		
		start = System.currentTimeMillis();
		qdr.getAbundanceOfComplexes();
		end = System.currentTimeMillis();
		
		long duration_approx = TimeUnit.MILLISECONDS.toSeconds(end - start);
		System.out.println(duration_LP + " " + duration_approx);
		
		Map<HashSet<String>, Double> abundances = lp.getAbundanceOfComplexes();
		List<String> to_write = new LinkedList<>();
		
		for (HashSet<String> cluster : abundances.keySet())
			to_write.add( String.join(",", cluster) + " " + abundances.get(cluster));
		
		Utilities.writeEntries(to_write, "lp_q_results_95_5/" + sample + ".txt.gz");
		
	}
}
