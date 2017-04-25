package framework;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class QuantDACOResultSet extends DACOResultSet {
	private  Map<String, String> protein_to_assumed_transcript;
	private  Map<String, Double> transcript_abundance; // transcript data actually only given in float precision, but double simplifies the usage in statistical functions
	private Map<HashSet<String>, Double> cached_abundance_of_complexes; // convenient storage for complex quantification results of the non-simple method
	
	/*
	 * diverse constructorsand necessities
	 */
	
	public QuantDACOResultSet(String daco_out_file, String seed_file, String protein_to_assumed_transcript_file) {
		super(daco_out_file, seed_file);
		this.readProteinTranscriptFile(protein_to_assumed_transcript_file);
	}

	public QuantDACOResultSet(String daco_out_file, Set<String> seed, String protein_to_assumed_transcript_file) {
		super(daco_out_file, seed);
		this.readProteinTranscriptFile(protein_to_assumed_transcript_file);
	}
	
	public QuantDACOResultSet(DACOResultSet daco_result, String protein_to_assumed_transcript_file) {
		super(daco_result.getResult(), daco_result.getAbundantSeedProteins());
		this.readProteinTranscriptFile(protein_to_assumed_transcript_file);
	}
	
	public QuantDACOResultSet(DACOResultSet daco_result, ConstructedNetworks constructed_network) {
		super(daco_result.getResult(), daco_result.getAbundantSeedProteins());
		
		this.protein_to_assumed_transcript = constructed_network.getProteinToAssumedTranscriptMap();
		this.transcript_abundance = new HashMap<String, Double>(1024);
		for (String transcript:constructed_network.getTranscriptAbundanceMap().keySet())
			this.transcript_abundance.put(transcript, (double) constructed_network.getTranscriptAbundanceMap().get(transcript));
	}
	
	public QuantDACOResultSet(HashSet<HashSet<String>> results, Set<String> seed, Map<String, String> protein_to_assumed_transcript, Map<String, Float> transcript_abundance) {
		super(results, seed);
		this.protein_to_assumed_transcript = protein_to_assumed_transcript;
		this.transcript_abundance = new HashMap<String, Double>(1024);
		for (String transcript:transcript_abundance.keySet())
			this.transcript_abundance.put(transcript, (double) transcript_abundance.get(transcript));
	}
	
	private void readProteinTranscriptFile(String protein_to_assumed_transcript_file) {
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		this.transcript_abundance = new HashMap<String, Double>(1024);
		
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			protein_to_assumed_transcript.put(spl[0], spl[1]);
			transcript_abundance.put(spl[1], Double.parseDouble(spl[2])); // assumes file that includes abundance values, error otherwise
		}
	}
	
	/**
	 * Rebuilds useful data-structures on the basis of relevant seed, removes other results and resets cached results
	 * @param refined_seed
	 */
	public void rebuildData(Set<String> refined_seed) {
		super.rebuildData(refined_seed);
		this.cached_abundance_of_complexes = null; // additionally reset the cache as results will differ if changed occurred
	}
	/*
	 * quantification functions
	 */
	
	/**
	 * Quantify each complex with the abundance of its least abundant member
	 * @return
	 */
	public Map<HashSet<String>, Double> getSimpleAbundanceOfComplexes() {
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		
		for (HashSet<String> complex:this.getResult()) {
			double min_abundance = Double.MAX_VALUE;
			for (String protein:complex) {
				double abundance = this.getProteinAbundance(protein);
				if (abundance < min_abundance)
					min_abundance = abundance;
			}
			quantification_result.put(complex, min_abundance);
		}
		
		return quantification_result;
	}
	
	/**
	 * Quantify each complex with the abundance of its least abundant member, but take overall amount of into account
	 * @return
	 */
	public Map<HashSet<String>, Double> getAbundanceOfComplexes() {
		
		// get cached results if already computed
		if (this.cached_abundance_of_complexes != null)
			return this.cached_abundance_of_complexes;
		
		// count occurrence of each protein in complexes
		Map<String, Integer> protein_in_complexes_count = new HashMap<>();
		Map<String, List<HashSet<String>>> protein_in_complexes = new HashMap<>();
		
		for (HashSet<String> complex:this.getResult())
			for (String protein:complex) {
				if (!protein_in_complexes.containsKey(protein))
					protein_in_complexes.put(protein, new LinkedList<HashSet<String>>());
				protein_in_complexes.get(protein).add(complex);
				protein_in_complexes_count.put(protein, protein_in_complexes_count.getOrDefault(protein, 0) + 1);
			}
		
		// filter proteins that are only in one complex
		protein_in_complexes_count.keySet().removeIf(d -> protein_in_complexes_count.get(d) == 1);
		
		// set initial "effective" abundance per protein in each complex by equal distribution
		Map<HashSet<String>, Map<String, Double>> protein_abundance_per_complex = new HashMap<>();
		for (HashSet<String> complex:this.getResult()) {
			Map<String, Double> effective_protein_in_complex = new HashMap<>();
			for (String protein:complex)
				effective_protein_in_complex.put(protein, this.getProteinAbundance(protein) / protein_in_complexes_count.getOrDefault(protein, 1));
			protein_abundance_per_complex.put(complex, effective_protein_in_complex);
		}
		
		// iterative
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		boolean iterate = false;
		double last_to_distr = Double.MAX_VALUE;
		
		do {
			iterate = false;
			Map<String, Double> remaining_amount = new HashMap<>();
			Map<String, List<HashSet<String>>> distribute_to = new HashMap<>();
			for (HashSet<String> complex:this.getResult()) {
				double min_abundance = protein_abundance_per_complex.get(complex).values().stream().min(Double::compare).get();
				quantification_result.put(complex, min_abundance);
				
				// calculate remaining amount of proteins and add it up to overall remaining
				for (String protein:complex) {
					// only relevant for those in several complexes -> do not touch proteins that are only in one complex
					if (!protein_in_complexes_count.containsKey(protein))
						continue;
					double distance_to_min = protein_abundance_per_complex.get(complex).get(protein) - min_abundance;
					remaining_amount.put(protein, remaining_amount.getOrDefault(protein, 0.0) + distance_to_min);
					
					// if protein is the limiting protein, note that it should be increased in the complex here
					if (distance_to_min == 0.0) {
						if (!distribute_to.containsKey(protein))
							distribute_to.put(protein, new LinkedList<HashSet<String>>());
						distribute_to.get(protein).add(complex);
					}
					
					// limit to maximal value, remaining amount of protein (actually transcript) as was added to remaining_amount
					protein_abundance_per_complex.get(complex).put(protein, min_abundance);
				}
			}
			remaining_amount.keySet().removeIf(d -> remaining_amount.get(d) == 0.0); // can happen if protein is always the limiting member
			distribute_to.keySet().removeIf(d -> !remaining_amount.containsKey(d)); // don't distribute when there's nothing to distribute
			
			// do another iteration if there is still a decrease in the overall amount of re-distributable abundance values
			double current_to_distr = remaining_amount.values().stream().reduce(0.0, Double::sum);
			if ( (last_to_distr - current_to_distr) > 0.0001)
				iterate = true;
			last_to_distr = current_to_distr;
			
			// distribute remaining account equally on the limiting proteins
			for (String protein:distribute_to.keySet()) {
				double distr_amount = remaining_amount.get(protein) / distribute_to.get(protein).size();
				for (HashSet<String> complex:distribute_to.get(protein))
					protein_abundance_per_complex.get(complex).put(protein, protein_abundance_per_complex.get(complex).get(protein) + distr_amount);
				
				remaining_amount.remove(protein);
			}
			distribute_to.clear();
			
			// distribute for all other proteins
			for (String protein:remaining_amount.keySet()) {
				double distr_amount = remaining_amount.get(protein) / protein_in_complexes_count.get(protein);
				for (HashSet<String> complex:protein_in_complexes.get(protein))
					protein_abundance_per_complex.get(complex).put(protein, protein_abundance_per_complex.get(complex).get(protein) + distr_amount);
			}
			remaining_amount.clear();
			
		} while (iterate);
		
		// cache results
		this.cached_abundance_of_complexes = quantification_result;
		
		return quantification_result;
	}
	
	/**
	 * Quantify each seed variant with the overall (simple) abundance of complexes containing it
	 * @return
	 */
	public Map<HashSet<String>, Double> getSimpleAbundanceOfSeedVariantsComplexes() {
		Map<HashSet<String>, Double> individual_quantification_result = this.getSimpleAbundanceOfComplexes();
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		
		for (HashSet<String> seed_variant:this.getSeedToComplexMap().keySet()) {
			List<Double> abundance_values = new LinkedList<>();
			for (HashSet<String> complex:this.getSeedToComplexMap().get(seed_variant))
				abundance_values.add(individual_quantification_result.get(complex));
			quantification_result.put(seed_variant, abundance_values.stream().reduce(0.0, Double::sum));
		}
		
		return quantification_result;
	}
	
	/**
	 * Quantify each seed variant with the overall abundance of complexes containing it
	 * @return
	 */
	public Map<HashSet<String>, Double> getAbundanceOfSeedVariantsComplexes() {
		Map<HashSet<String>, Double> individual_quantification_result = this.getAbundanceOfComplexes();
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		
		for (HashSet<String> seed_variant:this.getSeedToComplexMap().keySet()) {
			double abundance_values = 0.0;
			for (HashSet<String> complex:this.getSeedToComplexMap().get(seed_variant))
				abundance_values += individual_quantification_result.get(complex);
			quantification_result.put(seed_variant, abundance_values);
		}
		
		return quantification_result;
	}
	
	
	/*
	 * distance function
	 */
	
	/**
	 * Euclidean distance on abundance values for all reference seed variants
	 */
	public double getSeedVariantSetsAbundanceDistance(Set<HashSet<String>> reference_universe, QuantDACOResultSet result_set2) {
		Map<HashSet<String>, Double> result_abundances1 = this.getAbundanceOfSeedVariantsComplexes();
		Map<HashSet<String>, Double> result_abundances2 = result_set2.getAbundanceOfSeedVariantsComplexes();
		
		double sum = 0.0;
		for (HashSet<String> seed_variant:reference_universe) {
			sum += Math.pow(result_abundances1.getOrDefault(seed_variant, 0.0) - result_abundances2.getOrDefault(seed_variant, 0.0), 2);
		}
		
		return Math.sqrt(sum);
	}
	
	
	/*
	 * getters
	 */
	
	/**
	 * Returns a map of proteins in the input network and their most abundant transcript
	 * @return
	 */
	public Map<String, String> getProteinToAssumedTranscript() {
		return protein_to_assumed_transcript;
	}

	/**
	 * Returns a map of most abundant transcripts per protein and their abundance
	 * @return
	 */
	public Map<String, Double> getTranscriptAbundance() {
		return transcript_abundance;
	}
	
	/**
	 * Returns the expression value of the most abundant transcript associated with this protein and 0.0 if it cannot be found
	 * @param protein
	 * @return
	 */
	public double getProteinAbundance(String protein) {
		return this.transcript_abundance.getOrDefault(this.protein_to_assumed_transcript.get(protein), 0.0);
	}
	
	// for testing purposes
//	public static void main(String[] args) {
//		HashSet<HashSet<String>> results = new HashSet<>();
//		results.add(new HashSet<>(Arrays.asList("PA", "PB")));
//		results.add(new HashSet<>(Arrays.asList("PB", "PC")));
//		results.add(new HashSet<>(Arrays.asList("PC", "PD")));
//		results.add(new HashSet<>(Arrays.asList("PC", "PE")));
//		
//		Set<String> seed = new HashSet<>();
//		seed.add("PB");
//		
//		Map<String, String> protein_to_assumed_transcript = new HashMap<>();
//		protein_to_assumed_transcript.put("PA", "TA");
//		protein_to_assumed_transcript.put("PB", "TB");
//		protein_to_assumed_transcript.put("PC", "TC");
//		protein_to_assumed_transcript.put("PD", "TD");
//		protein_to_assumed_transcript.put("PE", "TE");
//		
//		Map<String, Float> transcript_abundance = new HashMap<>();
//		transcript_abundance.put("TA", 0.2f);
//		transcript_abundance.put("TB", 1f);
//		transcript_abundance.put("TC", 1f);
//		transcript_abundance.put("TD", 0.5f);
//		transcript_abundance.put("TE", 1f);
//		
//		QuantDACOResultSet qdr = new QuantDACOResultSet(results, seed, protein_to_assumed_transcript, transcript_abundance);
//		System.out.println("out: " + qdr.getAbundanceOfComplexes());
//	}
}
