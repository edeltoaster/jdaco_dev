package framework;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

public class QuantDACOResultSet extends DACOResultSet {
	private  Map<String, String> protein_to_assumed_transcript;
	private  Map<String, Float> transcript_abundance; // transcript data only given in float precision, but double later simplifies the usage in statistical functions in calculations
	private Map<HashSet<String>, Double> cached_abundance_of_complexes; // convenient storage for complex quantification results of the non-simple method
	private Map<String, Double> cached_remaining_abundance_of_proteins; // convenient storage for advanced complex quantification results of the non-simple method
	
	/*
	 * diverse constructors and necessities
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
		this.transcript_abundance = new HashMap<String, Float>(1024);
		for (String transcript:constructed_network.getTranscriptAbundanceMap().keySet())
			this.transcript_abundance.put(transcript, constructed_network.getTranscriptAbundanceMap().get(transcript));
	}
	
	public QuantDACOResultSet(HashSet<HashSet<String>> results, Set<String> seed, Map<String, String> protein_to_assumed_transcript, Map<String, Float> transcript_abundance) {
		super(results, seed);
		this.protein_to_assumed_transcript = protein_to_assumed_transcript;
		this.transcript_abundance = new HashMap<String, Float>(1024);
		for (String transcript:transcript_abundance.keySet())
			this.transcript_abundance.put(transcript, transcript_abundance.get(transcript));
	}
	
	private void readProteinTranscriptFile(String protein_to_assumed_transcript_file) {
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		this.transcript_abundance = new HashMap<String, Float>(1024);
		
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			protein_to_assumed_transcript.put(spl[0], spl[1]);
			transcript_abundance.put(spl[1], Float.parseFloat(spl[2])); // assumes file that includes abundance values, error otherwise
		}
	}
	
	/**
	 * (Re-)Builds useful data-structures on the basis of relevant seed, removes other results and resets cached results
	 * @param seed
	 */
	@Override
	public void buildData(Set<String> seed) {
		super.buildData(seed);
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
	 * Quantify each complex with the abundance of its least abundant member, but take overall amount of into account.
	 * Returns cached results if default parameters have already been used before.
	 * @return
	 */
	public Map<HashSet<String>, Double> getAbundanceOfComplexes() {
		// get cached results if already computed
		if (this.cached_abundance_of_complexes != null)
			return this.cached_abundance_of_complexes;
		
		double convergence_limit = 1.0E-7 * this.protein_to_assumed_transcript.size();
		
		// cache results
		this.cached_abundance_of_complexes = getAbundanceOfComplexes(convergence_limit, 10000, false);
		
		return this.cached_abundance_of_complexes;
	}
	
	/**
	 * Quantify each complex with the abundance of its least abundant member, but take overall amount of into account (detailed)
	 * @return
	 */
	public Map<HashSet<String>, Double> getAbundanceOfComplexes(final double convergence_limit, final int max_iterations, boolean verbose) {

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
		
		
		// some pre-allocated objects and variables, to limit the amount of object generation
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		Map<String, Double> remaining_amount = new HashMap<>();
		Map<String, List<HashSet<String>>> distribute_to = new HashMap<>();
		Set<HashSet<String>> unsaturated_complexes = new HashSet<>(this.getResult());
		Set<HashSet<String>> saturated_complexes = new HashSet<>();
		boolean limiting_protein_changeable = false;
		double last_to_distr = Double.MAX_VALUE;
		double min_abundance;
		double distance_to_min;
		double current_to_distr;
		double distr_factor;
		double distr_amount;
		int iteration_no = 1;
		
		// iterative refinement
		do {
			// determine limiting proteins and remaining abundances
			for (HashSet<String> complex:unsaturated_complexes) {
				
				// complex abundance is limited by least abundant protein
				min_abundance = protein_abundance_per_complex.get(complex).values().stream().min(Double::compare).get();
				quantification_result.put(complex, min_abundance);
				
				// calculate remaining amount of proteins and add it up to overall remaining
				limiting_protein_changeable = false;
				for (String protein:complex) {
					// only relevant for those in several complexes -> do not touch proteins that are only in one complex
					if (!protein_in_complexes_count.containsKey(protein))
						continue;
					// update used amount of protein and limit to maximal value, remaining amount of protein (actually transcript) as was added to remaining_amount
					distance_to_min = protein_abundance_per_complex.get(complex).get(protein) - min_abundance;
					remaining_amount.put(protein, remaining_amount.getOrDefault(protein, 0.0) + distance_to_min);
					protein_abundance_per_complex.get(complex).put(protein, min_abundance);
					
					// if protein is the limiting protein, note that it should be increased in the complex here
					if (distance_to_min == 0.0) {
						limiting_protein_changeable = true;
						if (!distribute_to.containsKey(protein))
							distribute_to.put(protein, new LinkedList<HashSet<String>>());
						distribute_to.get(protein).add(complex);
					}
				}
				
				// if the limiting protein cannot be changed anyhow, take this complex as saturated and don't touch it anymore, adjust proteins
				if (!limiting_protein_changeable) {
					saturated_complexes.add(complex);
					for (String protein:complex) {
						if (protein_in_complexes_count.containsKey(protein))
							protein_in_complexes_count.put(protein, protein_in_complexes_count.get(protein) - 1);
						protein_in_complexes.get(protein).remove(complex);
					}
				}
			}
			
			// some output for testing
			if (verbose) {
				System.out.println("#iteration: " + iteration_no);
				System.out.println("    Current remaining amount: " + remaining_amount);
				System.out.println("    Current complex: " + quantification_result);
			}
			
			// remove saturated complexes that cannot be changed anyhow
			unsaturated_complexes.removeAll(saturated_complexes);
			saturated_complexes.clear();
			
			remaining_amount.keySet().removeIf(d -> remaining_amount.get(d) == 0.0); // limiting proteins
			distribute_to.keySet().removeIf(d -> !remaining_amount.containsKey(d)); // don't distribute when there's nothing to distribute
			
			// do another iteration if there is still a decrease in the overall amount of re-distributable abundance values
			current_to_distr = remaining_amount.values().stream().reduce(0.0, Double::sum);
			if ((last_to_distr - current_to_distr) < convergence_limit) //approx match to float accuracy
				break;
			last_to_distr = current_to_distr;
			
			// stop algorithm if max_iterations were done
			if (iteration_no > max_iterations) {
				System.err.println("max_iterations reached, system did not converge.");
				break;
			}
			
			// saturation function that goes from 0.99 to 0.09
			distr_factor = 1.89 - (1.8 / (1+Math.exp(-0.05* (iteration_no-1))));
			
			// distribute remaining account equally on the limiting proteins if worthwhile, algorithm stopped otherwise
			for (String protein:distribute_to.keySet()) {
				distr_amount = (distr_factor * remaining_amount.get(protein)) / distribute_to.get(protein).size();
				for (HashSet<String> complex:distribute_to.get(protein))
					protein_abundance_per_complex.get(complex).put(protein, protein_abundance_per_complex.get(complex).get(protein) + distr_amount);
				remaining_amount.put(protein, (1.0-distr_factor) * remaining_amount.get(protein));
			}
			distribute_to.clear();
			
			// distribute for all other proteins
			for (String protein:remaining_amount.keySet()) {
				distr_amount = remaining_amount.get(protein) / protein_in_complexes_count.get(protein);
				for (HashSet<String> complex:protein_in_complexes.get(protein))
					protein_abundance_per_complex.get(complex).put(protein, protein_abundance_per_complex.get(complex).get(protein) + distr_amount);
			}
			remaining_amount.clear();
			
			++iteration_no;
		} while (true); // while(true) actually nicest form to implement that! :-P
		
		// set remaining abundance
		Map<String, Double> used_abundance = new HashMap<>();
		for (Entry<HashSet<String>, Double> e:quantification_result.entrySet()) {
			for (String p:e.getKey())
				used_abundance.put(p, used_abundance.getOrDefault(p, 0.0) + e.getValue());
		}
		this.cached_remaining_abundance_of_proteins = new HashMap<>();
		for (String p:this.protein_to_assumed_transcript.keySet()) {
			this.cached_remaining_abundance_of_proteins.put(p, this.getProteinAbundance(p) - used_abundance.getOrDefault(p, 0.0));
		}
		this.cached_remaining_abundance_of_proteins.keySet().removeIf(k -> this.cached_remaining_abundance_of_proteins.get(k) == 0.0);
		
		return quantification_result;
	}
	
	/**
	 * Returns the remaining abundance value of the last abundance computation for each protein with > 0 remaining abundance value.
	 * @return
	 */
	public Map<String, Double> getRemainingAbundanceOfProteins() {
		
		if (this.cached_remaining_abundance_of_proteins == null)
			this.getAbundanceOfComplexes();
		return this.cached_remaining_abundance_of_proteins;
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
	 * Utilities
	 */
	
	/**
	 * Writes output to a file in the format:
	 * prot_A,prot_B abundance.
	 * @param out_file
	 */
	public void writeQuantifiedResult(String out_file) {
		Map<HashSet<String>, Double> abundances = getAbundanceOfComplexes();
		List<String> to_write = new LinkedList<>();
		
		for (HashSet<String> cluster : abundances.keySet())
			to_write.add( String.join(",", cluster) + " " + abundances.get(cluster));
		
		Utilities.writeEntries(to_write, out_file);
	}
	
	/**
	 * Reads written quantified output from file
	 * @param path
	 * @return
	 */
	public static Map<HashSet<String>, Double> readQuantifiedResult(String path) {
		Map<HashSet<String>, Double> abundances = new HashMap<>();
		
		for (String line:Utilities.readFile(path)) {
			String[] line_spl = line.trim().split("\\s+");
			abundances.put(new HashSet<String>(Arrays.asList(line_spl[0].split(","))), Double.parseDouble(line_spl[1]));
		}
		
		return abundances;
	}
	
	/**
	 * Reads all files from a certain folder as quantified results: filename-prefix -> quantified results map.
	 * Assumes the standard file ending is "_qr.txt(.gz).
	 * @param folder
	 * @return
	 */
	public static Map<String, Map<HashSet<String>, Double>> readQuantifiedResults(String folder) {
		Map<String, Map<HashSet<String>, Double>> data = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(folder, "_qr.txt")) {
			String sample = f.getName().split("_ppin")[0];
			data.put(sample, readQuantifiedResult(f.getAbsolutePath()));
		}

		return data;
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
	public Map<String, Float> getTranscriptAbundance() {
		return transcript_abundance;
	}
	
	/**
	 * Returns the expression value associated with this protein (most abundant transcript or sum of all transcripts) and 0.0 if it cannot be found
	 * @param protein
	 * @return
	 */
	public double getProteinAbundance(String protein) {
		return this.transcript_abundance.getOrDefault(this.protein_to_assumed_transcript.get(protein), 0f);
	}
}