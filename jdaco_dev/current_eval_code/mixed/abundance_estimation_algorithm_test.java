package mixed;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import framework.QuantDACOResultSet;

public class abundance_estimation_algorithm_test {
	public static Random rnd = new Random(System.currentTimeMillis());
	
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
	
	public static void main(String[] args) {
		// to play around and test
	}
}
