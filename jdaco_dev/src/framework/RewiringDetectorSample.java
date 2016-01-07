package framework;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Stores interactions and protein_to_assumed_isoform map for a RewiringDetector sample,
 * similar to - but much slicker than - a ConstructedNetworks object
 * @author Thorsten Will
 */
public class RewiringDetectorSample {
	private final Set<StrPair> interactions;
	private final Map<String, String> protein_to_assumed_transcript;
	
	/**
	 * Constructor, builds objects from files
	 * @param ppin_file
	 * @param protein_to_assumed_transcript_file
	 */
	public RewiringDetectorSample(String ppin_file, String protein_to_assumed_transcript_file) {
		this.interactions = PPIN.readInteractionsFromPPINFile(ppin_file);
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			this.protein_to_assumed_transcript.put(spl[0], spl[1]);
		}
	}
	
	/**
	 * Returns the set of interactions in the network
	 * @return
	 */
	public Set<StrPair> getInteractions() {
		return interactions;
	}

	/**
	 * Returns a map of proteins in the network and their most abundant transcript
	 * @return
	 */
	public Map<String, String> getProteinToAssumedTranscriptMap() {
		return protein_to_assumed_transcript;
	}

	
	
	/*
	 * General purpose helpers
	 */
	
	/**
	 * Reads all matching PPIXpress output pairs of PPINs/major transcripts from a certain folder: 
	 * [folder]/[sample-string]_ppin.txt(.gz) and [folder]/[sample-string]_major-transcripts.txt(.gz)
	 * Assumes the standard file endings "_ppin.txt(.gz) and _major-transcripts.txt(.gz)".
	 * @param folder
	 * @return
	 */
	public static Map<String, RewiringDetectorSample> readNetworks(String folder) {
		Map<String, RewiringDetectorSample> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(folder, "_ppin.txt")) {
			String gz = "";
			if (f.getName().endsWith(".gz"))
				gz = ".gz";
			String pre = f.getAbsolutePath().split("_ppin")[0];
			String sample = f.getName().split("_ppin")[0];
			RewiringDetectorSample rds = new RewiringDetectorSample(pre + "_ppin.txt" + gz, pre + "_major-transcripts.txt" + gz);
			
			data.put(sample, rds);
		}

		return data;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((interactions == null) ? 0 : interactions.hashCode());
		result = prime * result
				+ ((protein_to_assumed_transcript == null) ? 0 : protein_to_assumed_transcript.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		RewiringDetectorSample other = (RewiringDetectorSample) obj;
		if (interactions == null) {
			if (other.interactions != null)
				return false;
		} else if (!interactions.equals(other.interactions))
			return false;
		if (protein_to_assumed_transcript == null) {
			if (other.protein_to_assumed_transcript != null)
				return false;
		} else if (!protein_to_assumed_transcript.equals(other.protein_to_assumed_transcript))
			return false;
		return true;
	}
}
