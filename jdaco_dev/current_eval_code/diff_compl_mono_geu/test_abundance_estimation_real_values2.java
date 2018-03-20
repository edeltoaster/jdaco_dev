package diff_compl_mono_geu;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import framework.Utilities;

public class test_abundance_estimation_real_values2 {

	static List<Double> rem_data_total = new LinkedList<>();
	static List<Double> rel_distr_data_total = new LinkedList<>();
	static List<String> distr_data_total = new LinkedList<>();
	
	public static Map<Set<String>, Double> readQDR(String path) {
		Map<Set<String>, Double> qdr = new HashMap<>();
		
		for (String s:Utilities.readFile(path)) {
			String[] line = s.trim().split("\\s+");
			qdr.put(new HashSet<String>(Arrays.asList(line[0].split(","))), Double.parseDouble(line[1]));
		}
		
		return qdr;
	}
	
	public static Map<String, Double> readProtAbun(String path) {
		Map<String, Double> res = new HashMap<>();
		
		for (String s:Utilities.readFile(path)) {
			String[] line = s.trim().split("\\s+");
			res.put(line[0], Double.parseDouble(line[2]));
		}
		
		return res;
	}
	
	public static void checkSample(String sample, String qdr_file, String major_tr_file) {
		
		// get data
		Map<Set<String>, Double> qdr = readQDR(qdr_file);
		Map<String, Double> prot_abun = readProtAbun(major_tr_file);

		
		// determine limiting proteins -> lowest abundance
		Set<String> limiting_proteins = new HashSet<>();
		Map<String, Double> sum_in_complexes = new HashMap<>();
		Map<String, List<Set<String>>> prot_compl_map = new HashMap<>();
		for (Set<String> compl:qdr.keySet()) {
			double compl_abun = qdr.get(compl);
			for (String p:compl) {
				sum_in_complexes.put(p, sum_in_complexes.getOrDefault(p, 0.0) + compl_abun);
				if (!prot_compl_map.containsKey(p))
					prot_compl_map.put(p, new LinkedList<Set<String>>());
				prot_compl_map.get(p).add(compl);
			}
		}
		
		Set<String> prots_in_compl = new HashSet<>(sum_in_complexes.keySet());
		Map<String, Double> remaining_abun = new HashMap<>();
		for (String p:prots_in_compl)
			remaining_abun.put(p, prot_abun.get(p) - sum_in_complexes.get(p));
		
		for (Set<String> compl:qdr.keySet()) {
			double lowest = 999999;
			String lowest_prot = null;
			for (String p:compl) {
				double p_rem = remaining_abun.get(p);
				if (p_rem < lowest) {
					lowest = p_rem;
					lowest_prot = p;
				}
			}
			limiting_proteins.add(lowest_prot); // no change if already in set
		}
		
		List<Double> rem_data = new ArrayList<>(prots_in_compl.size());
		List<Double> distr_data = new ArrayList<>(prots_in_compl.size());
		for (String prot:prots_in_compl) {
			
			if (limiting_proteins.contains(prot)) {
				
				List<Double> rel_compl_distr = new LinkedList<>();
				List<String> compl_distr = new LinkedList<>();
				for (Set<String> compl:prot_compl_map.get(prot)) {
						rel_compl_distr.add(qdr.get(compl));
						compl_distr.add(Double.toString(qdr.get(compl) / sum_in_complexes.get(prot)));
				}
				
				// determine avg distr of limiting proteins
				final double equal = sum_in_complexes.get(prot) / compl_distr.size();
				distr_data.add(Utilities.getMean(rel_compl_distr.stream().map(c -> (Math.abs(c-equal) / equal) ).collect(Collectors.toList())));
				distr_data_total.add(String.join(",", compl_distr));
				
			} else { // if not limiting
				// determine remaining prefactor
				double factor = remaining_abun.get(prot) / sum_in_complexes.get(prot);
				rem_data.add(factor);
			}

		}
		
		System.out.println(sample + " " + Utilities.getMedian(rem_data) + " " + Utilities.getMedian(distr_data));
		
		rem_data_total.add(Utilities.getMedian(rem_data));
		rel_distr_data_total.add(Utilities.getMedian(distr_data));
		
	}
	
	public static void main(String[] args) {
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.qr_output_folder, ".txt.gz")) {
			String sample = f.getName().split("\\.")[0];
			String qdr = f.getAbsolutePath();
			String major_transcript_file = definitions.networks_folder + sample + "_major-transcripts.txt.gz";
			checkSample(sample, qdr, major_transcript_file);
		}
		
		System.out.println();
		System.out.println("mean remaining of medians: " + Utilities.getMean(rem_data_total));
		System.out.println("mean distr of medians (of means): " + Utilities.getMean(rel_distr_data_total));
		System.out.println("limiting proteins determined as lowest remaining abundance in complex");
		Utilities.writeEntries(distr_data_total, "real_distr_data.csv.gz");
	}
}
