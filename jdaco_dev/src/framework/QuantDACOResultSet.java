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
	
	/*
	 * diverse constructors
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
	
	private void readProteinTranscriptFile(String protein_to_assumed_transcript_file) {
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		this.transcript_abundance = new HashMap<String, Double>(1024);
		
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			protein_to_assumed_transcript.put(spl[0], spl[1]);
			transcript_abundance.put(spl[1], Double.parseDouble(spl[2])); // assumes file that includes abundance values, error otherwise
		}
	}

	
	/*
	 * functions
	 */
	
	/**
	 * Quantify each complex with the abundance of its least abundant member
	 * @return
	 */
	public Map<HashSet<String>, Double> getSimpleAbundanceOfComplexes() {
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		
		for (LinkedList<HashSet<String>> complexes:this.getSeedToComplexMap().values())
			for (HashSet<String> complex:complexes) {
				double min_abundance = Double.MAX_VALUE;
				for (String protein:complex) {
					double abundance = this.transcript_abundance.get(this.protein_to_assumed_transcript.get(protein));
					if (abundance < min_abundance)
						min_abundance = abundance;
				}
				quantification_result.put(complex, min_abundance);
			}
		
		return quantification_result;
	}
	
	/**
	 * Quantify each seed variant with the median (simple) abundance of complexes containing it
	 * @return
	 */
	public Map<HashSet<String>, Double> getSimpleAbundanceOfSeedVariantsComplexes() {
		Map<HashSet<String>, Double> individual_quantification_result = this.getSimpleAbundanceOfComplexes();
		Map<HashSet<String>, Double> quantification_result = new HashMap<>();
		
		for (HashSet<String> seed_variant:this.getSeedToComplexMap().keySet()) {
			List<Double> abundance_values = new LinkedList<>();
			for (HashSet<String> complex:this.getSeedToComplexMap().get(seed_variant))
				abundance_values.add(individual_quantification_result.get(complex));
			quantification_result.put(seed_variant, Utilities.getMedian(abundance_values));
		}
		
		return quantification_result;
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
	
}
