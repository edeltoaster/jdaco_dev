package diff_compl_mono_geu;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import framework.QuantDACOResultSet;
import framework.Utilities;

public class prep_abundance_estimation_LP {
	static int no_iterations = 20;
	static double[] stds = new double[]{0.2, 1.0};
	static double[] prefactors = new double[]{5, 30, 60};
	static String tmp_folder = "tmp/";
	
	static Random rnd = new Random(System.currentTimeMillis());
	static Map<Integer, List<ArrayList<Double>>> real_distr = new HashMap<>();

	/**
	 * Realistic simulation of the "equal distribution"-model on the basis of real data and noise regarding the equality of the distribution of abundance values and 
	 * the portion of remaining protein abundances.
	 * @param daco_outfile
	 * @param major_transcript_file
	 * @param std_factor
	 * @param remaining_prefactor
	 * @return
	 */
	public static Object[] simulate_sample_model(String daco_outfile, String major_transcript_file, double std_factor, double remaining_prefactor) {
		// load real complex result
		QuantDACOResultSet dr = new QuantDACOResultSet(daco_outfile, "mixed_data/hocomoco_human_TFs_v10.txt.gz", major_transcript_file);
		List<HashSet<String>> complexes = new ArrayList<>(dr.getResult());
		List<Float> tr_abundance_values = new ArrayList<>(dr.getTranscriptAbundance().values());
		
		// construct appropriate protein->transcript data structure
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		for (HashSet<String> complex:complexes)
			for (String protein:complex)
				if (!protein_to_assumed_transcript.containsKey(protein))
					protein_to_assumed_transcript.put(protein, "T" + protein);
		
		// select limiting proteins
		HashMap<String, LinkedList<HashSet<String>>> limiting_protein_to_complex = new HashMap<>();
		Collections.shuffle(complexes, rnd); // order is shuffled to ensure very different choices of limiting proteins
		for (HashSet<String> complex:complexes) {
			Set<String> temp_set = new HashSet<>(limiting_protein_to_complex.keySet());
			temp_set.retainAll(complex);
			if (temp_set.size() == 0) { // if there is no limiting protein yet, choose one
				int chosen_index = rnd.nextInt(complex.size());
				String chosen = (String) complex.toArray()[chosen_index];
				limiting_protein_to_complex.put(chosen, new LinkedList<HashSet<String>>());
			}
		}
		
		// note in which complexes limiting proteins are found
		for (HashSet<String> complex:complexes) {
			List<String> temp_list = new ArrayList<>(complex);
			temp_list.retainAll(limiting_protein_to_complex.keySet());
			
			if (temp_list.size() > 1)// there may be one or more proteins limiting, chose one randomly and link that
				Collections.shuffle(temp_list, rnd);
			for (String limiting_protein:temp_list) { 
				limiting_protein_to_complex.get(limiting_protein).add(complex);
				break;
			}
		}
		limiting_protein_to_complex.keySet().removeIf(e->limiting_protein_to_complex.get(e).size() == 0);
		
		// distribute limiting proteins roughly equally using uniform distribution
		Map<HashSet<String>, Double> artificial_complex_abundance = new HashMap<>();
		Map<String, List<Double>> limiting_protein_distribution = new HashMap<>();
		for (String limiting_protein:limiting_protein_to_complex.keySet()) {
			int chosen_index = rnd.nextInt(tr_abundance_values.size());
			double prot_abundance = tr_abundance_values.get(chosen_index); // will be fast as it is an array list
			final int size = limiting_protein_to_complex.get(limiting_protein).size();
			double complex_abundance_mean = prot_abundance / size;
			List<Double> abundances = new ArrayList<>(size);
			
			// introduce the noise as a factor of the mean, correct the sum afterwards and ensure that there are only abundances > 0
			do {
				abundances.clear();
				double sum = 0;
				for (int i = 0; i < size; i++) {
					double alteration_amount = rnd.nextDouble() * std_factor * complex_abundance_mean;
					if (rnd.nextBoolean())
						alteration_amount *= -1;
					abundances.add(complex_abundance_mean + alteration_amount);
					sum += alteration_amount;
				}
				sum /= size;
				sum *= -1;
				final double fsum = sum;
				abundances = abundances.stream().map(d->d + fsum).collect(Collectors.toList());
			} while (abundances.stream().anyMatch( d -> d <= 0));
			
			// collect for statistics
			List<Double> relative_deviation_distribution = abundances.stream().map(a -> Math.abs(a-complex_abundance_mean)/complex_abundance_mean).collect(Collectors.toList());
			limiting_protein_distribution.put(limiting_protein, relative_deviation_distribution);
			
			// shuffle and associate as many as needed values with complexes that are limited by that protein
			Collections.shuffle(abundances, rnd);
			int i = 0;
			for (HashSet<String> complex:limiting_protein_to_complex.get(limiting_protein)) {
				artificial_complex_abundance.put(complex, abundances.get(i));
				++i;
			}
		}
		
		// sum up transcript abundances
		Map<String, Float> transcript_abundance = new HashMap<>();
		for (HashSet<String> complex:artificial_complex_abundance.keySet()) {
			double abundance = artificial_complex_abundance.get(complex);
			for (String protein:complex) {
				String transcript_string = "T" + protein;
				transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + (float) abundance);
			}
		}
		
		// add to non-limiting proteins
		Map<String, Double> remaining_protein_abundance = new HashMap<>();
		Map<String, Double> rel_remaining_protein_abundance = new HashMap<>();
		for (String transcript:transcript_abundance.keySet()) {
			String prot = transcript.substring(1);
			if (limiting_protein_to_complex.containsKey(prot)) {
				remaining_protein_abundance.put(prot, 0.0);
				continue;
			}
			float compl_abundance = rnd.nextFloat() * ((float) remaining_prefactor) * transcript_abundance.get(transcript);
			remaining_protein_abundance.put(prot, (double) compl_abundance);
			rel_remaining_protein_abundance.put(prot, (double) compl_abundance / transcript_abundance.get(transcript));
			transcript_abundance.put(transcript, transcript_abundance.getOrDefault(transcript, 0f) + compl_abundance);
		}
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		
		Object[] output = new Object[5];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		output[2] = remaining_protein_abundance; // only remaining abundance of all proteins that are in complexes, limiting proteins have rem_abundance 0.0
		output[3] = limiting_protein_distribution; // for statistics on the model
		output[4] = rel_remaining_protein_abundance; // for statistics on the model
		
		return output;
	}
	
	/**
	 * Helper function to run analyses on on a given sample in a batch fashion.
	 * @param std_factor
	 * @param remaining_prefactor
	 * @return
	 */
	public static void prep_sample_model_run(String daco_result_file, String major_transcripts_file, double std_factor, double remaining_prefactor, int iteration) {
		
		// get results of simulation
		Object[] simulation = simulate_sample_model(daco_result_file, major_transcripts_file, std_factor, remaining_prefactor);
		QuantDACOResultSet qdr = (QuantDACOResultSet) simulation[0];
		@SuppressWarnings("unchecked")
		Map<HashSet<String>, Double> artificial_complex_abundance = (Map<HashSet<String>, Double>) simulation[1];
		@SuppressWarnings("unchecked")
		Map<String, Double> art_remaining_protein_abundance = (Map<String, Double>) simulation[2];
		
		// write files
		String sample = daco_result_file.split("/")[daco_result_file.split("/").length-1].split("\\.")[0];
		String run_id = sample + "-" + std_factor + "-" + remaining_prefactor + "-" + iteration;
		String compl_res = tmp_folder + run_id + "_daco.txt.gz";
		String mt = tmp_folder + run_id + "_mt.txt.gz";
		String art_compl_ab = tmp_folder + run_id + "_artcomplab.txt.gz";
		String art_prot_ab = tmp_folder + run_id + "_artpab.txt.gz";
		qdr.writeResult(compl_res);
		qdr.writeMajorTranscriptFile(mt);
		
		// write artificial data
		List<String> to_write = new LinkedList<>();
		for (HashSet<String> cluster : artificial_complex_abundance.keySet())
			to_write.add( String.join(",", cluster) + " " + artificial_complex_abundance.get(cluster));
		Utilities.writeEntries(to_write, art_compl_ab);
		to_write.clear();
		for (String p : art_remaining_protein_abundance.keySet())
			to_write.add( p + " " + art_remaining_protein_abundance.get(p));
		Utilities.writeEntries(to_write, art_prot_ab);
	}
	
	/**
	 * Tests prediction on real data
	 */
	public static void benchmark() {
		
		// read all data
		Map<String, String> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String daco_result = f.getAbsolutePath();
			String major_transcript_file = definitions.networks_folder + sample + "_major-transcripts.txt.gz";
			data.put(daco_result, major_transcript_file);
		}
		
		List<String> all_iterations = new LinkedList<>();
		all_iterations.add("sample std prefactor iter corr_compl rmsd_compl corr_rem rmsd_rem");
		
		System.out.println("Running LP benchmark on M/GEU");
		System.out.println(no_iterations + " iterations for each sample (" + data.size() + " samples)");
		System.out.println("STDevs: " + Arrays.toString(stds));
		System.out.println("Prefactors: " + Arrays.toString(prefactors));

		for (double std:stds)
			for (double prefactor:prefactors) {
				System.out.println("Running " + std + "/" + prefactor);
				System.out.flush();
				
				for (Entry<String, String> sample : data.entrySet()) {
					IntStream.range(0, no_iterations).parallel().boxed().forEach(d -> prep_sample_model_run(sample.getKey(), sample.getValue(), std, prefactor, d));
				}
			}
		
		System.out.println("Finished!");
	}
	
	
	/*
	 * Analyses using real distribution values and variable remaining values.
	 */
	
	/**
	 * Samples from a real distribution of distributions
	 * @param size
	 * @return
	 */
	public static ArrayList<Double> getDistr(int size) {
		ArrayList<Double> distr = null;
		
		if (real_distr.containsKey(size))
			distr = real_distr.get(size).get(rnd.nextInt(real_distr.get(size).size()));
		else {
			int size1 = size / 2;
			int size2 = size - size1;
			distr = new ArrayList<Double>(getDistr(size1));
			distr.addAll(getDistr(size2));
			distr = new ArrayList<Double>(distr.stream().map(f -> (f/2.0)).collect(Collectors.toList()));
		}
		
		return distr;
	}
	
	/**
	 * Realistic simulation of the "equal distribution"-model on the basis of real data and noise regarding the portion of remaining protein abundances.
	 * @param daco_outfile
	 * @param major_transcript_file
	 * @param remaining_prefactor
	 * @return
	 */
	public static Object[] simulate_sample_model2(String daco_outfile, String major_transcript_file, double remaining_prefactor) {
		// load real complex result
		QuantDACOResultSet dr = new QuantDACOResultSet(daco_outfile, "mixed_data/hocomoco_human_TFs_v10.txt.gz", major_transcript_file);
		List<HashSet<String>> complexes = new ArrayList<>(dr.getResult());
		List<Float> tr_abundance_values = new ArrayList<>(dr.getTranscriptAbundance().values());
		
		// construct appropriate protein->transcript data structure
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		for (HashSet<String> complex:complexes)
			for (String protein:complex)
				if (!protein_to_assumed_transcript.containsKey(protein))
					protein_to_assumed_transcript.put(protein, "T" + protein);
		
		// select limiting proteins
		HashMap<String, LinkedList<HashSet<String>>> limiting_protein_to_complex = new HashMap<>();
		Collections.shuffle(complexes, rnd); // order is shuffled to ensure very different choices of limiting proteins
		for (HashSet<String> complex:complexes) {
			Set<String> temp_set = new HashSet<>(limiting_protein_to_complex.keySet());
			temp_set.retainAll(complex);
			if (temp_set.size() == 0) { // if there is no limiting protein yet, choose one
				int chosen_index = rnd.nextInt(complex.size());
				String chosen = (String) complex.toArray()[chosen_index];
				limiting_protein_to_complex.put(chosen, new LinkedList<HashSet<String>>());
			}
		}
		
		// note in which complexes limiting proteins are found
		for (HashSet<String> complex:complexes) {
			List<String> temp_list = new ArrayList<>(complex);
			temp_list.retainAll(limiting_protein_to_complex.keySet());
			
			if (temp_list.size() > 1)// there may be one or more proteins limiting, chose one randomly and link that
				Collections.shuffle(temp_list, rnd);
			for (String limiting_protein:temp_list) { 
				limiting_protein_to_complex.get(limiting_protein).add(complex);
				break;
			}
		}
		limiting_protein_to_complex.keySet().removeIf(e->limiting_protein_to_complex.get(e).size() == 0);
		
		
		Map<HashSet<String>, Double> artificial_complex_abundance = new HashMap<>();
		Map<String, List<Double>> limiting_protein_distribution = new HashMap<>();
		for (String limiting_protein:limiting_protein_to_complex.keySet()) {
			int chosen_index = rnd.nextInt(tr_abundance_values.size());
			int size = limiting_protein_to_complex.get(limiting_protein).size();
			double prot_abundance = tr_abundance_values.get(chosen_index); // will be fast as it is an array list
			double complex_abundance_mean = prot_abundance / size;
			
			// determine distribution from real values
			List<Double> abundances = new ArrayList<>(size);
			ArrayList<Double> distr = getDistr(size);
			for (int i = 0; i < size; i++) {
				abundances.add(prot_abundance * distr.get(i));
			}
			
			// collect for statistics
			List<Double> relative_deviation_distribution = abundances.stream().map(a -> Math.abs(a-complex_abundance_mean)/complex_abundance_mean).collect(Collectors.toList());
			limiting_protein_distribution.put(limiting_protein, relative_deviation_distribution);
			
			// shuffle and associate as many as needed values with complexes that are limited by that protein
			Collections.shuffle(abundances, rnd);
			int i = 0;
			for (HashSet<String> complex:limiting_protein_to_complex.get(limiting_protein)) {
				artificial_complex_abundance.put(complex, abundances.get(i));
				++i;
			}
		}
		
		// sum up transcript abundances
		Map<String, Float> transcript_abundance = new HashMap<>();
		for (HashSet<String> complex:artificial_complex_abundance.keySet()) {
			double abundance = artificial_complex_abundance.get(complex);
			for (String protein:complex) {
				String transcript_string = "T" + protein;
				transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + (float) abundance);
			}
		}
		
		// add to non-limiting proteins
		Map<String, Double> remaining_protein_abundance = new HashMap<>();
		Map<String, Double> rel_remaining_protein_abundance = new HashMap<>();
		for (String transcript:transcript_abundance.keySet()) {
			String prot = transcript.substring(1);
			if (limiting_protein_to_complex.containsKey(prot)) {
				remaining_protein_abundance.put(prot, 0.0);
				continue;
			}
			float compl_abundance = rnd.nextFloat() * ((float) remaining_prefactor) * transcript_abundance.get(transcript);
			remaining_protein_abundance.put(prot, (double) compl_abundance);
			rel_remaining_protein_abundance.put(prot, (double) compl_abundance / transcript_abundance.get(transcript));
			transcript_abundance.put(transcript, transcript_abundance.getOrDefault(transcript, 0f) + compl_abundance);
		}
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		
		Object[] output = new Object[5];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		output[2] = remaining_protein_abundance; // only remaining abundance of all proteins that are in complexes, limiting proteins have rem_abundance 0.0
		output[3] = limiting_protein_distribution; // for statistics on the model
		output[4] = rel_remaining_protein_abundance; // for statistics on the model
		
		return output;
	}
	
	/**
	 * Helper function to run analyses on on a given sample in a batch fashion. Version for real distr values.
	 * @param std_factor
	 * @param remaining_prefactor
	 * @return
	 */
	public static void prep_sample_model_run2(String daco_result_file, String major_transcripts_file, double remaining_prefactor, int iteration) {
		
		// get results of simulation
		Object[] simulation = simulate_sample_model2(daco_result_file, major_transcripts_file, remaining_prefactor);
		QuantDACOResultSet qdr = (QuantDACOResultSet) simulation[0];
		@SuppressWarnings("unchecked")
		Map<HashSet<String>, Double> artificial_complex_abundance = (Map<HashSet<String>, Double>) simulation[1];
		@SuppressWarnings("unchecked")
		Map<String, Double> art_remaining_protein_abundance = (Map<String, Double>) simulation[2];
		
		// prepare evaluation
		// write files
		String sample = daco_result_file.split("/")[daco_result_file.split("/").length-1].split("\\.")[0];
		String run_id = sample + "-" + "rd" + "-" + remaining_prefactor + "-" + iteration;
		String compl_res = tmp_folder + run_id + "_daco.txt.gz";
		String mt = tmp_folder + run_id + "_mt.txt.gz";
		String art_compl_ab = tmp_folder + run_id + "_artcomplab.txt.gz";
		String art_prot_ab = tmp_folder + run_id + "_artpab.txt.gz";
		qdr.writeResult(compl_res);
		qdr.writeMajorTranscriptFile(mt);
		
		// write artificial data
		List<String> to_write = new LinkedList<>();
		for (HashSet<String> cluster : artificial_complex_abundance.keySet())
			to_write.add( String.join(",", cluster) + " " + artificial_complex_abundance.get(cluster));
		Utilities.writeEntries(to_write, art_compl_ab);
		to_write.clear();
		for (String p : art_remaining_protein_abundance.keySet())
			to_write.add( p + " " + art_remaining_protein_abundance.get(p));
		Utilities.writeEntries(to_write, art_prot_ab);
	}
	
	/**
	 * Tests prediction on real data, uses real distributions
	 */
	public static void benchmark2() {
		
		// read all data
		Map<String, String> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String daco_result = f.getAbsolutePath();
			String major_transcript_file = definitions.networks_folder + sample + "_major-transcripts.txt.gz";
			data.put(daco_result, major_transcript_file);
		}
		
		List<String> all_iterations = new LinkedList<>();
		all_iterations.add("sample std prefactor iter corr_compl rmsd_compl corr_rem rmsd_rem");
		
		System.out.println("Running LP benchmark on M/GEU using real distributions");
		System.out.println(no_iterations + " iterations for each sample (" + data.size() + " samples)");
		System.out.println("real distributions: rd");
		System.out.println("Prefactors: " + Arrays.toString(prefactors));
		
		String std = "rd";
		for (double prefactor:prefactors) {
			System.out.println("Running " + std + "/" + prefactor);
			System.out.flush();
			
			for (Entry<String, String> sample : data.entrySet()) {
				IntStream.range(0, no_iterations).parallel().boxed().forEach(d -> prep_sample_model_run2(sample.getKey(), sample.getValue(), prefactor, d));
			}
		}
		
		System.out.println("Finished!");
	}
	
	
	/**
	 * Realistic simulation of the "rnd distribution"-model on the basis of rnd data and noise regarding the portion of remaining protein abundances.
	 * @param daco_outfile
	 * @param major_transcript_file
	 * @param remaining_prefactor
	 * @return
	 */
	public static Object[] simulate_sample_model3(String daco_outfile, String major_transcript_file, double remaining_prefactor) {
		// load real complex result
		QuantDACOResultSet dr = new QuantDACOResultSet(daco_outfile, "mixed_data/hocomoco_human_TFs_v10.txt.gz", major_transcript_file);
		List<HashSet<String>> complexes = new ArrayList<>(dr.getResult());
		List<Float> tr_abundance_values = new ArrayList<>(dr.getTranscriptAbundance().values());
		
		// construct appropriate protein->transcript data structure
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		for (HashSet<String> complex:complexes)
			for (String protein:complex)
				if (!protein_to_assumed_transcript.containsKey(protein))
					protein_to_assumed_transcript.put(protein, "T" + protein);
		
		// select limiting proteins
		HashMap<String, LinkedList<HashSet<String>>> limiting_protein_to_complex = new HashMap<>();
		Collections.shuffle(complexes, rnd); // order is shuffled to ensure very different choices of limiting proteins
		for (HashSet<String> complex:complexes) {
			Set<String> temp_set = new HashSet<>(limiting_protein_to_complex.keySet());
			temp_set.retainAll(complex);
			if (temp_set.size() == 0) { // if there is no limiting protein yet, choose one
				int chosen_index = rnd.nextInt(complex.size());
				String chosen = (String) complex.toArray()[chosen_index];
				limiting_protein_to_complex.put(chosen, new LinkedList<HashSet<String>>());
			}
		}
		
		// note in which complexes limiting proteins are found
		for (HashSet<String> complex:complexes) {
			List<String> temp_list = new ArrayList<>(complex);
			temp_list.retainAll(limiting_protein_to_complex.keySet());
			
			if (temp_list.size() > 1)// there may be one or more proteins limiting, chose one randomly and link that
				Collections.shuffle(temp_list, rnd);
			for (String limiting_protein:temp_list) { 
				limiting_protein_to_complex.get(limiting_protein).add(complex);
				break;
			}
		}
		limiting_protein_to_complex.keySet().removeIf(e->limiting_protein_to_complex.get(e).size() == 0);
		
		
		Map<HashSet<String>, Double> artificial_complex_abundance = new HashMap<>();
		Map<String, List<Double>> limiting_protein_distribution = new HashMap<>();
		for (String limiting_protein:limiting_protein_to_complex.keySet()) {
			int chosen_index = rnd.nextInt(tr_abundance_values.size());
			int size = limiting_protein_to_complex.get(limiting_protein).size();
			double prot_abundance = tr_abundance_values.get(chosen_index); // will be fast as it is an array list
			double complex_abundance_mean = prot_abundance / size;
			
			// determine distribution from real values
			List<Double> abundances = new ArrayList<>(size);
			ArrayList<Double> distr = test_abundance_estimation.getRndDistr(size);
			for (int i = 0; i < size; i++) {
				abundances.add(prot_abundance * distr.get(i));
			}
			
			// collect for statistics
			List<Double> relative_deviation_distribution = abundances.stream().map(a -> Math.abs(a-complex_abundance_mean)/complex_abundance_mean).collect(Collectors.toList());
			limiting_protein_distribution.put(limiting_protein, relative_deviation_distribution);
			
			// shuffle and associate as many as needed values with complexes that are limited by that protein
			Collections.shuffle(abundances, rnd);
			int i = 0;
			for (HashSet<String> complex:limiting_protein_to_complex.get(limiting_protein)) {
				artificial_complex_abundance.put(complex, abundances.get(i));
				++i;
			}
		}
		
		// sum up transcript abundances
		Map<String, Float> transcript_abundance = new HashMap<>();
		for (HashSet<String> complex:artificial_complex_abundance.keySet()) {
			double abundance = artificial_complex_abundance.get(complex);
			for (String protein:complex) {
				String transcript_string = "T" + protein;
				transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + (float) abundance);
			}
		}
		
		// add to non-limiting proteins
		Map<String, Double> remaining_protein_abundance = new HashMap<>();
		Map<String, Double> rel_remaining_protein_abundance = new HashMap<>();
		for (String transcript:transcript_abundance.keySet()) {
			String prot = transcript.substring(1);
			if (limiting_protein_to_complex.containsKey(prot)) {
				remaining_protein_abundance.put(prot, 0.0);
				continue;
			}
			float compl_abundance = rnd.nextFloat() * ((float) remaining_prefactor) * transcript_abundance.get(transcript);
			remaining_protein_abundance.put(prot, (double) compl_abundance);
			rel_remaining_protein_abundance.put(prot, (double) compl_abundance / transcript_abundance.get(transcript));
			transcript_abundance.put(transcript, transcript_abundance.getOrDefault(transcript, 0f) + compl_abundance);
		}
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		
		Object[] output = new Object[5];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		output[2] = remaining_protein_abundance; // only remaining abundance of all proteins that are in complexes, limiting proteins have rem_abundance 0.0
		output[3] = limiting_protein_distribution; // for statistics on the model
		output[4] = rel_remaining_protein_abundance; // for statistics on the model
		
		return output;
	}
	
	/**
	 * Helper function to run analyses on on a given sample in a batch fashion. Version for rnd distr values.
	 * @param std_factor
	 * @param remaining_prefactor
	 * @return
	 */
	public static void prep_sample_model_run3(String daco_result_file, String major_transcripts_file, double remaining_prefactor, int iteration) {
		
		// get results of simulation
		Object[] simulation = simulate_sample_model3(daco_result_file, major_transcripts_file, remaining_prefactor);
		QuantDACOResultSet qdr = (QuantDACOResultSet) simulation[0];
		@SuppressWarnings("unchecked")
		Map<HashSet<String>, Double> artificial_complex_abundance = (Map<HashSet<String>, Double>) simulation[1];
		@SuppressWarnings("unchecked")
		Map<String, Double> art_remaining_protein_abundance = (Map<String, Double>) simulation[2];
		
		// prepare evaluation
		// write files
		String sample = daco_result_file.split("/")[daco_result_file.split("/").length-1].split("\\.")[0];
		String run_id = sample + "-" + "rnd" + "-" + remaining_prefactor + "-" + iteration;
		String compl_res = tmp_folder + run_id + "_daco.txt.gz";
		String mt = tmp_folder + run_id + "_mt.txt.gz";
		String art_compl_ab = tmp_folder + run_id + "_artcomplab.txt.gz";
		String art_prot_ab = tmp_folder + run_id + "_artpab.txt.gz";
		qdr.writeResult(compl_res);
		qdr.writeMajorTranscriptFile(mt);
		
		// write artificial data
		List<String> to_write = new LinkedList<>();
		for (HashSet<String> cluster : artificial_complex_abundance.keySet())
			to_write.add( String.join(",", cluster) + " " + artificial_complex_abundance.get(cluster));
		Utilities.writeEntries(to_write, art_compl_ab);
		to_write.clear();
		for (String p : art_remaining_protein_abundance.keySet())
			to_write.add( p + " " + art_remaining_protein_abundance.get(p));
		Utilities.writeEntries(to_write, art_prot_ab);
	}
	
	/**
	 * Tests prediction on real data, uses real distributions
	 */
	public static void benchmark3() {
		
		// read all data
		Map<String, String> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(definitions.daco_results_folder, ".csv.gz")) {
			String sample = f.getName().split("\\.")[0];
			String daco_result = f.getAbsolutePath();
			String major_transcript_file = definitions.networks_folder + sample + "_major-transcripts.txt.gz";
			data.put(daco_result, major_transcript_file);
		}
		
		List<String> all_iterations = new LinkedList<>();
		all_iterations.add("sample std prefactor iter corr_compl rmsd_compl corr_rem rmsd_rem");
		
		System.out.println("Running LP benchmark on M/GEU using real distributions");
		System.out.println(no_iterations + " iterations for each sample (" + data.size() + " samples)");
		System.out.println("rnd distributions");
		System.out.println("Prefactors: " + Arrays.toString(prefactors));
		
		String std = "rnd";
		for (double prefactor:prefactors) {
			System.out.println("Running " + std + "/" + prefactor);
			System.out.flush();
			
			for (Entry<String, String> sample : data.entrySet()) {
				IntStream.range(0, no_iterations).parallel().boxed().forEach(d -> prep_sample_model_run3(sample.getKey(), sample.getValue(), prefactor, d));
			}
		}
		
		System.out.println("Finished!");
	}
	
	
	public static void main(String[] args) {
		// read realistic data
		for (String s:Utilities.readFile("real_distr_data.csv.gz")) {
			String[] spl = s.trim().split(",");
			int len = spl.length;
			ArrayList<Double> dspl = new ArrayList<>(len);
			for (int i = 0; i < len; i++)
				dspl.add(Double.parseDouble(spl[i]));
			if (!real_distr.containsKey(len))
				real_distr.put(len, new ArrayList<ArrayList<Double>>());
			real_distr.get(len).add(dspl);
		}
		
		try {
			
			benchmark2();
			
			System.out.println();
			
			benchmark();
			
			System.out.println();
			
			benchmark3();
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
}
