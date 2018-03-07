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

public class test_abundance_estimation_real_values {

	static double limiting_prot_threshold = 1.0E-7;
	static List<Double> rem_data_total = new LinkedList<>();
	static List<Double> distr_data_total = new LinkedList<>();
	
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
		Map<Set<String>, Double> qdr = readQDR(qdr_file);
		Map<String, Double> prot_abun = readProtAbun(major_tr_file);
		
		Set<String> prots_in_compl = new HashSet<>();
		qdr.keySet().stream().forEach(compl -> prots_in_compl.addAll(compl));
		
		List<Double> rem_data = new ArrayList<>(prots_in_compl.size());
		List<Double> distr_data = new ArrayList<>(prots_in_compl.size());
		
		for (String prot:prots_in_compl) {
			double sum_in_compl = 0.0;
			List<Double> compl_distr = new ArrayList<>(prots_in_compl.size());
			for (Set<String> compl:qdr.keySet()) {
				if (compl.contains(prot)) {
					sum_in_compl += qdr.get(compl);
					compl_distr.add(qdr.get(compl));
				}
			}
			
			// determine remaining prefactor
			double remaining = prot_abun.get(prot) - sum_in_compl;
			double factor = remaining / sum_in_compl;
			
			if (remaining >= limiting_prot_threshold)
				rem_data.add(factor);
			
			// determine avg distr
			final double equal = sum_in_compl / compl_distr.size();
			
			if (remaining < limiting_prot_threshold)
				distr_data.add(Utilities.getMean(compl_distr.stream().map(c -> (Math.abs(c-equal) / equal) ).collect(Collectors.toList())));
		}
		
		System.out.println(sample + " " + Utilities.getMedian(rem_data) + " " + Utilities.getMedian(distr_data));
		
		rem_data_total.add(Utilities.getMedian(rem_data));
		distr_data_total.add(Utilities.getMedian(distr_data));
		
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
		System.out.println("mean distr of medians (of means): " + Utilities.getMean(distr_data_total));
		System.out.println("limiting threshold: " + limiting_prot_threshold);
	}
}
