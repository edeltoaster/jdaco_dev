package framework;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;

/**
 * Utility class to estimate abundances for DACO results
 * @author Thorsten Will
 */
public class DACOResultQuant {
	
	private final DACOResultSet daco_result;
	private final Map<String, String> protein_to_assumed_transcript;
	private final Map<String, Float> transcript_abundance;
	
	/**
	 * Constructor for object input
	 * @param daco_result
	 * @param constructed_networks (needs to have non-null transcript_abundance data)
	 */
	public DACOResultQuant(DACOResultSet daco_result, ConstructedNetworks constructed_networks) {
		this.daco_result = daco_result;
		this.protein_to_assumed_transcript = constructed_networks.getProteinToAssumedTranscriptMap();
		this.transcript_abundance = constructed_networks.getTranscriptAbundanceMap();
	}
	
	/**
	 * Constructor for file input
	 * @param daco_result_file
	 * @param seed_file
	 * @param protein_to_assumed_transcript_file (file needs to include abundance values)
	 */
	public DACOResultQuant(String daco_result_file, String seed_file, String protein_to_assumed_transcript_file) {
		this.daco_result = new DACOResultSet(daco_result_file, seed_file);
		
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		this.transcript_abundance = new HashMap<String, Float>(1024);
		
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			this.protein_to_assumed_transcript.put(spl[0], spl[1]);
			this.transcript_abundance.put(spl[1], Float.parseFloat(spl[2])); // assumes file that includes abundance values, error otherwise
		}
	}
	
	/**
	 * Quantify each complex with the abundance of its least abundant member
	 * @return
	 */
	public Map<HashSet<String>, Float> getSimpleAbundanceOfComplexes() {
		Map<HashSet<String>, Float> quantification_result = new HashMap<>();
		
		for (LinkedList<HashSet<String>> complexes:this.daco_result.getSeedToComplexMap().values())
			for (HashSet<String> complex:complexes) {
				float min_abundance = Float.MAX_VALUE;
				for (String protein:complex) {
					float abundance = this.transcript_abundance.get(this.protein_to_assumed_transcript.get(protein));
					if (abundance < min_abundance)
						min_abundance = abundance;
				}
				quantification_result.put(complex, min_abundance);
			}
		
		return quantification_result;
	}

	
	/*
	 * getters
	 */
	
	public DACOResultSet getDACOResultSet() {
		return daco_result;
	}

	public Map<String, String> getProteinToAssumedTranscriptMap() {
		return protein_to_assumed_transcript;
	}

	public Map<String, Float> getTranscriptAbundanceMap() {
		return transcript_abundance;
	}
	
	
}
