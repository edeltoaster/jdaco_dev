package mixed;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import framework.DACOResultSet;
import framework.QuantDACOResultSet;
import framework.Utilities;

public class diff_compl_test_cases {
	
	// example why complexes that consist only of one protein do not work well here; gladly this cannot happen with DACO results!
	public static void single_prot_complex_test() {
		HashSet<HashSet<String>> results = new HashSet<>();
		results.add(new HashSet<>(Arrays.asList("PA", "PB")));
		results.add(new HashSet<>(Arrays.asList("PC", "PD")));
		results.add(new HashSet<>(Arrays.asList("PD", "PE")));
		results.add(new HashSet<>(Arrays.asList("PD", "PE", "PF")));
		results.add(new HashSet<>(Arrays.asList("PF")));
		
		Set<String> seed = new HashSet<>();
		seed.add("PB");
		
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		protein_to_assumed_transcript.put("PA", "TA");
		protein_to_assumed_transcript.put("PB", "TB");
		protein_to_assumed_transcript.put("PC", "TC");
		protein_to_assumed_transcript.put("PD", "TD");
		protein_to_assumed_transcript.put("PE", "TE");
		protein_to_assumed_transcript.put("PF", "TF");
		
		Map<String, Float> transcript_abundance = new HashMap<>();
		transcript_abundance.put("TA", 0.2f);
		transcript_abundance.put("TB", 1f);
		transcript_abundance.put("TC", 0.2f);
		transcript_abundance.put("TD", 0.5f);
		transcript_abundance.put("TE", 1f);
		transcript_abundance.put("TF", 1f);
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(results, seed, protein_to_assumed_transcript, transcript_abundance);
		System.out.println("out: " + qdr.getAbundanceOfComplexes());
	}
	
	// example where only a part of the remaining_abundance can be directly distributed to complexes that are constrained by the particular abundance
	public static void alterating_example() {
		HashSet<HashSet<String>> results = new HashSet<>();
		results.add(new HashSet<>(Arrays.asList("PA", "PB")));
		results.add(new HashSet<>(Arrays.asList("PC", "PD")));
		results.add(new HashSet<>(Arrays.asList("PD", "PE")));
		results.add(new HashSet<>(Arrays.asList("PE", "PF", "PG")));
		
		Set<String> seed = new HashSet<>();
		seed.add("PB");
		
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		protein_to_assumed_transcript.put("PA", "TA");
		protein_to_assumed_transcript.put("PB", "TB");
		protein_to_assumed_transcript.put("PC", "TC");
		protein_to_assumed_transcript.put("PD", "TD");
		protein_to_assumed_transcript.put("PE", "TE");
		protein_to_assumed_transcript.put("PF", "TF");
		protein_to_assumed_transcript.put("PG", "TG");
		
		Map<String, Float> transcript_abundance = new HashMap<>();
		transcript_abundance.put("TA", 0.2f);
		transcript_abundance.put("TB", 1.0f);
		transcript_abundance.put("TC", 0.2f);
		transcript_abundance.put("TD", 0.5f);
		transcript_abundance.put("TE", 1.5f);
		transcript_abundance.put("TF", 1.0f);
		transcript_abundance.put("TG", 1.0f);
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(results, seed, protein_to_assumed_transcript, transcript_abundance);
		System.out.println("out: " + qdr.getAbundanceOfComplexes());
	}
	
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
		Random rnd = new Random(System.currentTimeMillis());
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
	
	public static Object[] simulate_sample_model(String daco_outfile, double min_abundance, double max_abundance, double std_factor, double remaining_prefactor) {
		double scaled_max_abundance = max_abundance - min_abundance;
		
		// load real complex result
		DACOResultSet dr = new DACOResultSet(daco_outfile, "mixed_data/hocomoco_human_TFs_v10.txt.gz");
		List<HashSet<String>> complexes = new ArrayList<>(dr.getResult());
		
		// construct appropriate protein->transcript data structure
		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
		for (HashSet<String> complex:complexes)
			for (String protein:complex)
				if (!protein_to_assumed_transcript.containsKey(protein))
					protein_to_assumed_transcript.put(protein, "T" + protein);
		
		// select limiting proteins
		Random rnd = new Random(System.currentTimeMillis());
		HashMap<String, LinkedList<HashSet<String>>> limiting_protein_to_complex = new HashMap<>();
		Collections.shuffle(complexes, rnd); // order is shuffled to ensure very different choices of limiting proteins
		for (HashSet<String> complex:complexes) {
			Set<String> temp_set = new HashSet<>(limiting_protein_to_complex.keySet());
			temp_set.retainAll(complex);
			if (temp_set.size() == 0) { // if there is no limiting protein yet, choose one
				int chosen_index = rnd.nextInt(complex.size());
				String chosen = (String) complex.toArray()[chosen_index];
				limiting_protein_to_complex.put(chosen, new LinkedList<HashSet<String>>());;
			}
		}
		
		// note in which complexes limiting proteins are found
		for (HashSet<String> complex:complexes) {
			List<String> temp_list = new ArrayList<>(complex);
			temp_list.retainAll(limiting_protein_to_complex.keySet());
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
		for (String limiting_protein:limiting_protein_to_complex.keySet()) {
			double prot_abundance = (rnd.nextDouble() * scaled_max_abundance) + min_abundance;
			double complex_abundance_mean = prot_abundance / limiting_protein_to_complex.get(limiting_protein).size();
			double complex_abundance_std = complex_abundance_mean * std_factor;
			NormalDistribution distr = new NormalDistribution(complex_abundance_mean, complex_abundance_std);
			for (HashSet<String> complex:limiting_protein_to_complex.get(limiting_protein)) {
				double abundance = -1;
				while (abundance <= 0.0) // TODO: other distribution?
					abundance = distr.sample();
				artificial_complex_abundance.put(complex, abundance);
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
		for (String transcript:transcript_abundance.keySet()) {
			String prot = transcript.substring(1);
			if (limiting_protein_to_complex.containsKey(prot))
				continue;
			float compl_abundance = rnd.nextFloat() * ((float) remaining_prefactor) * transcript_abundance.get(transcript);
			transcript_abundance.put(transcript, transcript_abundance.getOrDefault(transcript, 0f) + compl_abundance);
		}
		
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		
		Object[] output = new Object[3];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		
		return output;
	}
	
	// for testing purposes
	public static void main(String[] args) {
		double[] stds = new double[]{0.01, 0.05, 0.1, 0.2};
		for (double std:stds) {
			for (int i = 0; i < 3; i++) {
				// get reference data
				System.out.println("std: " + std + ", construct ...");
				Object[] simulation = simulate_sample_model("mixed_data/A172_1_1_ENCSR580GSX.csv.gz", 0.0001, 10, std, 2);
				QuantDACOResultSet qdr = (QuantDACOResultSet) simulation[0];
				@SuppressWarnings("unchecked")
				Map<HashSet<String>, Double> artificial_complex_abundance = (Map<HashSet<String>, Double>) simulation[1];
				
				System.out.println("evaluate ...");
				Map<HashSet<String>, Double> eval_complex_abundance = qdr.getAbundanceOfComplexes();
				
				List<Double> artificial = new LinkedList<>();
				List<Double> evaluated = new LinkedList<>();
				for (HashSet<String> complex: artificial_complex_abundance.keySet()) {
					artificial.add(artificial_complex_abundance.get(complex));
					evaluated.add(eval_complex_abundance.get(complex));
				}
				
				PearsonsCorrelation pcorr = new PearsonsCorrelation();
				double corr = pcorr.correlation(Utilities.getDoubleArray(artificial), Utilities.getDoubleArray(evaluated)); // TODO: correlation a good measure here?
				System.out.println("std: " + std + ", corr: " + corr);
			}
		}
	}
}
