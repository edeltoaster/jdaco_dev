package diff_compl_mono_geu;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import framework.QuantDACOResultSet;
import framework.Utilities;

public class collect_abundance_estimation_LP {
	static String tmp_folder = "tmp/";
	
	public static double[] calcResults(String art_compl, String art_prot_rem, String LP_compl, String mt) {
		// read artificial complexes
		Map<HashSet<String>, Double> complab_art = QuantDACOResultSet.readQuantifiedResult(art_compl);
		
		// read artificial remaining abundances
		Map<String, Double> remaining_art = new HashMap<>();
		for (String s:Utilities.readEntryFile(art_prot_rem)) {
			String[] spl = s.trim().split("\\s+");
			remaining_art.put(spl[0], Double.parseDouble(spl[1]));
		}
		
		// read LP predicted complex abundances
		Map<HashSet<String>, Double> complab_LP = QuantDACOResultSet.readQuantifiedResult(LP_compl);
		
		// read major-transcripts file
		Map<String, Double> prot_abundance_LP = new HashMap<>();
		for (String s:Utilities.readEntryFile(mt)) {
			String[] spl = s.trim().split("\\s+");
			prot_abundance_LP.put(spl[0], Double.parseDouble(spl[2]));
		}
		
		// set remaining abundance
		Map<String, Double> remaining_LP = new HashMap<>();
		Map<String, Double> used_abundance = new HashMap<>();
		for (Entry<HashSet<String>, Double> e:complab_LP.entrySet()) {
			for (String p:e.getKey())
				used_abundance.put(p, used_abundance.getOrDefault(p, 0.0) + e.getValue());
		}
		for (String p:prot_abundance_LP.keySet()) {
			remaining_LP.put(p, prot_abundance_LP.get(p) - used_abundance.getOrDefault(p, 0.0));
		}
		remaining_LP.keySet().removeIf(k -> remaining_LP.get(k) == 0.0);
		
		
		// build evaluation data structures
		List<Double> artificial = new LinkedList<>();
		List<Double> evaluated = new LinkedList<>();
		for (HashSet<String> complex: complab_art.keySet()) {
			artificial.add(complab_art.get(complex));
			evaluated.add(complab_LP.get(complex));
		}
		List<Double> rem_artificial = new LinkedList<>();
		List<Double> rem_evaluated = new LinkedList<>();
		for (String protein:remaining_art.keySet()) { // artificial data structure involves all proteins that are in complexes
			rem_artificial.add(remaining_art.get(protein));
			rem_evaluated.add(remaining_LP.getOrDefault(protein, 0.0)); // here, limiting proteins are NOT listed
		}
		
		// evaluate
		PearsonsCorrelation pcorr = new PearsonsCorrelation();
		double[] results = new double[4];
		results[0] = pcorr.correlation(Utilities.getDoubleArray(artificial), Utilities.getDoubleArray(evaluated));
		results[1] = Utilities.getRMSD(artificial, evaluated);
		results[2] = pcorr.correlation(Utilities.getDoubleArray(rem_artificial), Utilities.getDoubleArray(rem_evaluated));
		results[3] = Utilities.getRMSD(rem_artificial, rem_evaluated);
		
		return results;
	}
	
	public static void main(String[] args) {
		
		List<String> all_iterations = new LinkedList<>();
		all_iterations.add("sample std prefactor iter corr_compl rmsd_compl corr_rem rmsd_rem");
		
		int i = 0;
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(tmp_folder, "_LP.txt.gz")) {
			
			// load data
			// sample + "-" + std_factor + "-" + remaining_prefactor + "-" + iteration;
			String run_id = f.getName().split("_LP")[0];
			String[] run_id_spl = run_id.split("-");
			String sample = run_id_spl[0];
			String std = run_id_spl[1];
			String prefactor = run_id_spl[2];
			String n = run_id_spl[3];
			
			String LP_results_file = f.getAbsolutePath();
			String mt = tmp_folder + run_id + "_mt.txt.gz";
			String art_compl_ab = tmp_folder + run_id + "_artcomplab.txt.gz";
			String art_prot_ab = tmp_folder + run_id + "_artpab.txt.gz";
			
			// add to output
			double result[] = calcResults(art_compl_ab, art_prot_ab, LP_results_file, mt);
			all_iterations.add(sample + " " + std + " " + prefactor + " " + n + " " + result[0] + " " + result[1] + " " + result[2] + " " + result[3]);
			i++;
			
			if (i % 500 == 0) {
				System.out.println(i + "th file.");
				Utilities.writeEntries(all_iterations, "perf_all_iterations_LP.txt.gz");
			}
		}
		
		Utilities.writeEntries(all_iterations, "perf_all_iterations_LP.txt.gz");
	}
}
