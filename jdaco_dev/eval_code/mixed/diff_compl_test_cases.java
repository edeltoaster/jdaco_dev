package mixed;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

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
	
	/**
	 * Given a DACO result and some parameters, returns [QuantDACOResultSet qdr, Map<HashSet<String>, Double> artificial_complex_abundance, Map<String, Double> remaining_abundances] as an Object[]. 
	 * Here, the approximation algorithm is run in between to adjust for the model.
	 * @param daco_outfile
	 * @param min_abundance
	 * @param max_abundance
	 * @return
	 */
	public static Object[] simulate_sample_model(String daco_outfile, double min_abundance, double max_abundance) {// TODO: think about model, likely too optimistic
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
		
		// build qdr object from the artificial data and compute abundances using model
		QuantDACOResultSet qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		artificial_complex_abundance = qdr.getAbundanceOfComplexes();
		Map<String, Double> remaining_abundances = qdr.getRemainingAbundanceOfProteins();
		
		// integrate into new artificial data that better matches the model
		transcript_abundance.clear();
		for (HashSet<String> complex:artificial_complex_abundance.keySet()) {
			float compl_abundance = (float) ((double) artificial_complex_abundance.get(complex));
			for (String protein:complex) {
				String transcript_string = "T" + protein;
				transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + compl_abundance);
			}
		}
		
		for (String protein:remaining_abundances.keySet()) {
			String transcript_string = "T" + protein;
			transcript_abundance.put(transcript_string, transcript_abundance.getOrDefault(transcript_string, 0f) + (float) (double)remaining_abundances.get(protein));
		}
			
		qdr = new QuantDACOResultSet(new HashSet<HashSet<String>>(dr.getResult()), Utilities.readEntryFile("mixed_data/hocomoco_human_TFs_v10.txt.gz"), protein_to_assumed_transcript, transcript_abundance);
		
		Object[] output = new Object[3];
		output[0] = qdr;
		output[1] = artificial_complex_abundance;
		output[2] = remaining_abundances;
		return output;
	}
	
	// for testing purposes
	public static void main(String[] args) {
		
		for (int i=0; i< 5; i++) {
			// get reference data
			System.out.println("construct ...");
			Object[] simulation = simulate_sample_model("mixed_data/A172_1_1_ENCSR580GSX.csv.gz", 0.00001, 10);
			QuantDACOResultSet qdr = (QuantDACOResultSet) simulation[0];
			@SuppressWarnings("unchecked")
			Map<HashSet<String>, Double> artificial_complex_abundance = (Map<HashSet<String>, Double>) simulation[1];
			@SuppressWarnings("unchecked")
			Map<String, Double> remaining_abundances = (Map<String, Double>) simulation[2];
			
			System.out.println("evaluate ...");
			Map<HashSet<String>, Double> eval_complex_abundance = qdr.getAbundanceOfComplexes();
			

			System.out.println(qdr.getRemainingAbundanceOfProteins().values().stream().reduce(0.0, Double::sum));
			System.out.println(remaining_abundances.values().stream().reduce(0.0, Double::sum));
			
			List<Double> artificial = new LinkedList<>();
			List<Double> evaluated = new LinkedList<>();
			for (HashSet<String> complex: artificial_complex_abundance.keySet()) {
				artificial.add(artificial_complex_abundance.get(complex));
				evaluated.add(eval_complex_abundance.get(complex));
			}
			
			PearsonsCorrelation pcorr = new PearsonsCorrelation();
			double corr = pcorr.correlation(Utilities.getDoubleArray(artificial), Utilities.getDoubleArray(evaluated));
			System.out.println("corr: " + corr);
		}
	}
}
