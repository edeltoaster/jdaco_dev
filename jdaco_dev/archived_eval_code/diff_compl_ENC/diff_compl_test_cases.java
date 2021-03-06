package diff_compl_ENC;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import framework.DACOResultSet;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class diff_compl_test_cases {
	public static Random rnd = new Random(System.currentTimeMillis());
	
	public static Map<HashSet<String>, Double> bio_example() {
		QuantDACOResultSet qdr = new QuantDACOResultSet("mixed_data/A172_1_1_ENCSR580GSX.csv.gz", "mixed_data/hocomoco_human_TFs_v10.txt.gz", "mixed_data/A172_1_1_ENCSR580GSX_major-transcripts.txt.gz");
		return qdr.getAbundanceOfComplexes();
	}
	
	/**
	 * Given a DACO result and some parameters, returns [QuantDACOResultSet qdr, Map<HashSet<String>, Double> artificial_complex_abundance] as an Object[]
	 * @param daco_outfile
	 * @param min_abundance
	 * @param max_abundance
	 * @return
	 */
	public static Object[] simulate_sample(String daco_outfile, double min_abundance, double max_abundance, double remaining_prefactor) {
		// load real complex result
		DACOResultSet dr = new DACOResultSet(daco_outfile, "mixed_data/hocomoco_human_TFs_v10.txt.gz");
		
		// construct appropriate protein->transcript data structure
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		for (HashSet<String> complex:dr.getResult())
			for (String protein:complex)
				if (!protein_to_assumed_transcript.containsKey(protein))
					protein_to_assumed_transcript.put(protein, "T" + protein);
		
		// assign random abundance values to complexes
		final double scaled_max_abundance = max_abundance - min_abundance;
		Map<HashSet<String>, Double> artificial_complex_abundance = dr.getResult().stream().collect(Collectors.toMap(e -> e, e -> (rnd.nextDouble() * scaled_max_abundance) + min_abundance));
		
		// summarize transcript abundances
		Map<String, Float> transcript_abundance = new HashMap<>();
		for (HashSet<String> complex:artificial_complex_abundance.keySet()) {
			float compl_abundance = (float) ((double) artificial_complex_abundance.get(complex));
			for (String protein:complex) {
				String transcript_string = "T" + protein;
				transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + compl_abundance);
			}
		}
		
		// select limiting proteins
		Set<String> limiting_proteins = new HashSet<>();
		for (HashSet<String> complex:artificial_complex_abundance.keySet()) {
			Set<String> temp_set = new HashSet<>(limiting_proteins);
			temp_set.retainAll(complex);
			if (temp_set.size() == 0) { // if there is no limiting protein yet, choose one
				int chosen_index = rnd.nextInt(complex.size());
				String chosen = (String) complex.toArray()[chosen_index];
				limiting_proteins.add(chosen);
			}
		}
		
		// factor in remaining abundances
		for (String protein:protein_to_assumed_transcript.keySet()) {
			// limiting proteins do not have remaining abundance
			if (limiting_proteins.contains(protein))
				continue;
			String transcript_string = "T" + protein;
			float compl_abundance = rnd.nextFloat() * ((float) remaining_prefactor) * transcript_abundance.getOrDefault(transcript_string, 0f);
			transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + compl_abundance);
		}
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		Object[] output = new Object[2];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		return output;
	}
	
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
		HashMap<String, LinkedList<HashSet<String>>> limiting_protein_to_all_complexes = new HashMap<>();
		Collections.shuffle(complexes, rnd); // order is shuffled to ensure very different choices of limiting proteins
		for (HashSet<String> complex:complexes) {
			Set<String> temp_set = new HashSet<>(limiting_protein_to_complex.keySet());
			temp_set.retainAll(complex);
			if (temp_set.size() == 0) { // if there is no limiting protein yet, choose one
				int chosen_index = rnd.nextInt(complex.size());
				String chosen = (String) complex.toArray()[chosen_index];
				limiting_protein_to_complex.put(chosen, new LinkedList<HashSet<String>>());
				limiting_protein_to_all_complexes.put(chosen, new LinkedList<HashSet<String>>());
			}
		}
		
		// note in which complexes limiting proteins are found
		for (HashSet<String> complex:complexes) {
			List<String> temp_list = new ArrayList<>(complex);
			temp_list.retainAll(limiting_protein_to_complex.keySet());
			for (String limiting_protein:temp_list)
				limiting_protein_to_all_complexes.get(limiting_protein).add(complex);
			
			if (temp_list.size() > 0)// there may be one or more proteins limiting, chose one randomly and link that
				Collections.shuffle(temp_list, rnd);
			for (String limiting_protein:temp_list) { 
				limiting_protein_to_complex.get(limiting_protein).add(complex);
				break;
			}
		}
		limiting_protein_to_complex.keySet().removeIf(e->limiting_protein_to_complex.get(e).size() == 0);
		
		// distribute limiting proteins roughly equally using normal distribution
		Map<HashSet<String>, Double> artificial_complex_abundance = new HashMap<>();
		Map<String, List<Double>> limiting_protein_distribution = new HashMap<>();
		for (String limiting_protein:limiting_protein_to_complex.keySet()) {
			int chosen_index = rnd.nextInt(tr_abundance_values.size());
			double prot_abundance = tr_abundance_values.get(chosen_index); // will be fast as it is an array list
			double complex_abundance_mean = prot_abundance / limiting_protein_to_all_complexes.get(limiting_protein).size();
			List<Double> abundances = new ArrayList<>(limiting_protein_to_all_complexes.get(limiting_protein).size());
			
			// introduce the noise as a factor of the mean, correct the sum afterwards and ensure that there are only abundances > 0
			do {
				abundances.clear();
				double sum = 0;
				for (int i = 0; i < limiting_protein_to_all_complexes.get(limiting_protein).size(); i++) {
					double alteration_amount = rnd.nextDouble() * std_factor * complex_abundance_mean;
					if (rnd.nextBoolean())
						alteration_amount *= -1;
					abundances.add(complex_abundance_mean + alteration_amount);
					sum += alteration_amount;
				}
				sum /= limiting_protein_to_all_complexes.get(limiting_protein).size();
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
		for (String transcript:transcript_abundance.keySet()) {
			String prot = transcript.substring(1);
			if (limiting_protein_to_complex.containsKey(prot)) {
				remaining_protein_abundance.put(prot, 0.0);
				continue;
			}
			float compl_abundance = rnd.nextFloat() * ((float) remaining_prefactor) * transcript_abundance.get(transcript);
			remaining_protein_abundance.put(prot, (double) compl_abundance);
			transcript_abundance.put(transcript, transcript_abundance.getOrDefault(transcript, 0f) + compl_abundance);
		}
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		
		Object[] output = new Object[4];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		output[2] = remaining_protein_abundance; // only remaining abundance of all proteins that are in complexes, limiting proteins have rem_abundance 0.0
		output[3] = limiting_protein_distribution; // for statistics on the model
		
		return output;
	}
	
	/**
	 * Helper function to run analyses on sample A172_1_1_ENCSR580GSX (simply the very first one) in a batch fashion.
	 * @param std_factor
	 * @param remaining_prefactor
	 * @return
	 */
	public static double[] simulate_sample_model_run(double std_factor, double remaining_prefactor, int iteration, String[] sample_construction_outputs) {
		
		// get results of simulation
		Object[] simulation = simulate_sample_model("mixed_data/A172_1_1_ENCSR580GSX.csv.gz", "mixed_data/A172_1_1_ENCSR580GSX_major-transcripts.txt.gz", std_factor, remaining_prefactor);
		QuantDACOResultSet qdr = (QuantDACOResultSet) simulation[0];
		@SuppressWarnings("unchecked")
		Map<HashSet<String>, Double> artificial_complex_abundance = (Map<HashSet<String>, Double>) simulation[1];
		@SuppressWarnings("unchecked")
		Map<String, Double> art_remaining_protein_abundance = (Map<String, Double>) simulation[2];
		@SuppressWarnings("unchecked")
		Map<String, List<Double>> limiting_protein_distribution = (Map<String, List<Double>>) simulation[3];
		
		// get some sample construction values
		double complex_abundance_min = artificial_complex_abundance.values().stream().min(Double::compareTo).get();
		double complex_abundance_mean = Utilities.getMean(artificial_complex_abundance.values());
		double complex_abundance_std = Utilities.getStd(artificial_complex_abundance.values());
		double complex_abundance_max = artificial_complex_abundance.values().stream().max(Double::compareTo).get();
		
		List<Double> means = limiting_protein_distribution.values().stream().map(l->Utilities.getMean(l)).collect(Collectors.toList());
		List<Double> mins = limiting_protein_distribution.values().stream().map(l->l.stream().min(Double::compareTo).get()).collect(Collectors.toList());
		List<Double> maxs = limiting_protein_distribution.values().stream().map(l->l.stream().max(Double::compareTo).get()).collect(Collectors.toList());
		double lim_distr_min = Utilities.getMean(mins);
		double lim_distr_mean = Utilities.getMean(means);
		double lim_distr_std = Utilities.getStd(means);
		double lim_distr_max = Utilities.getMean(maxs);
		
		double rem_abundance_min = art_remaining_protein_abundance.values().stream().min(Double::compareTo).get();
		double rem_abundance_mean = Utilities.getMean(art_remaining_protein_abundance.values());
		double rem_abundance_std = Utilities.getStd(art_remaining_protein_abundance.values());
		double rem_abundance_max = art_remaining_protein_abundance.values().stream().max(Double::compareTo).get();
		
		// prepare evaluation
		Map<HashSet<String>, Double> eval_complex_abundance = qdr.getAbundanceOfComplexes();
		Map<String, Double> eval_remaining_protein_abundance = qdr.getRemainingAbundanceOfProteins();
		
		List<Double> artificial = new LinkedList<>();
		List<Double> evaluated = new LinkedList<>();
		for (HashSet<String> complex: artificial_complex_abundance.keySet()) {
			artificial.add(artificial_complex_abundance.get(complex));
			evaluated.add(eval_complex_abundance.get(complex));
		}
		List<Double> rem_artificial = new LinkedList<>();
		List<Double> rem_evaluated = new LinkedList<>();
		for (String protein:art_remaining_protein_abundance.keySet()) { // artificial datastructure involves all proteins that are in complexes
			rem_artificial.add(art_remaining_protein_abundance.get(protein));
			rem_evaluated.add(eval_remaining_protein_abundance.getOrDefault(protein, 0.0)); // here, limiting proteins are NOT listed
		}
		
		// evaluate
		PearsonsCorrelation pcorr = new PearsonsCorrelation();
		double[] results = new double[4];
		results[0] = pcorr.correlation(Utilities.getDoubleArray(artificial), Utilities.getDoubleArray(evaluated));
		results[1] = Utilities.getRMSD(artificial, evaluated);
		results[2] = pcorr.correlation(Utilities.getDoubleArray(rem_artificial), Utilities.getDoubleArray(rem_evaluated));
		results[3] = Utilities.getRMSD(rem_artificial, rem_evaluated);
		
		if (sample_construction_outputs != null) {
			String out = std_factor + " " + remaining_prefactor + " " + (iteration+1) + " " + complex_abundance_min + " " + complex_abundance_mean + " " + complex_abundance_std + " " + complex_abundance_max + " " + lim_distr_min + " " + lim_distr_mean + " " + lim_distr_std + " " + lim_distr_max + " " + rem_abundance_min + " " + rem_abundance_mean + " " +rem_abundance_std + " " + rem_abundance_max;
			sample_construction_outputs[iteration] = out;
		}
		
		return results;
	}
	
	/**
	 * Tests prediction to simulated data
	 */
	public static void benchmark() {
		
		List<String> all_iterations = new LinkedList<>();
		List<String> averaged = new LinkedList<>();
		List<String> sample_construction = new LinkedList<>();
		averaged.add("std prefactor corr_compl rmsd_compl corr_rem rmsd_rem");
		all_iterations.add("std prefactor iter corr_compl rmsd_compl corr_rem rmsd_rem");
		sample_construction.add("std prefactor iter complex_abundance_min complex_abundance_mean complex_abundance_std complex_abundance_max lim_distr_min lim_distr_mean lim_distr_std lim_distr_max rem_abundance_min rem_abundance_mean rem_abundance_std rem_abundance_max");
		
		int no_iterations = 100;
		double[] stds = new double[]{0.1, 0.25, 0.5, 0.75, 1.0};
		double[] prefactors = new double[]{0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5};
		
		for (double std:stds) 
			for (double prefactor:prefactors) {
				String[] sample_construction_outputs = new String[no_iterations];
				List<double[]> results = IntStream.range(0, no_iterations).boxed().parallel().map(d->simulate_sample_model_run(std, prefactor, d, sample_construction_outputs)).collect(Collectors.toList());
				List<Double> corrs = new LinkedList<>();
				List<Double> rmsds = new LinkedList<>();
				List<Double> rem_corrs = new LinkedList<>();
				List<Double> rem_rmsds = new LinkedList<>();
				int n = 1;
				for (double[] result:results) {
					all_iterations.add(std + " " + prefactor + " " + n + " " + result[0] + " " + result[1] + " " + result[2] + " " + result[3]);
					sample_construction.add(sample_construction_outputs[n-1]);
					corrs.add(result[0]);
					rmsds.add(result[1]);
					rem_corrs.add(result[2]);
					rem_rmsds.add(result[3]);
					++n;
				}
				averaged.add(std + " " + prefactor + " " + Utilities.getMean(corrs) + "+-" + Utilities.getStd(corrs) + " " + Utilities.getMean(rmsds) + "+-" + Utilities.getStd(rmsds)+ " " + Utilities.getMedian(rem_corrs) + "+-" + Utilities.getStd(rem_corrs) + " " + Utilities.getMedian(rem_rmsds) + "+-" + Utilities.getStd(rem_rmsds));
			}
		
		Utilities.writeEntries(all_iterations, "perf_all_iterations.txt");
		Utilities.writeEntries(averaged, "perf_averaged.txt");
		Utilities.writeEntries(sample_construction, "sample_construction.txt");
	}
	
	public static void main(String[] args) {
		benchmark();
	}
}
