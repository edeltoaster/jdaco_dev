package framework;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Stores PPI and DDI interactions and protein_to_assumed_isoform map for a RewiringDetector sample,
 * similar to a ConstructedNetworks object but much slicker and tuned towards PPICompare
 * @author Thorsten Will
 */
public class RewiringDetectorSample {
	private final Set<StrPair> interactions;
	private final Map<String, String> protein_to_assumed_transcript;
	
	// sophisticated part of domain interaction data
	private final Map<String, Set<String>> protein_to_domains; // set of non-FB domains per protein in this network, not in there if no non-FB domain
	private final Map<String, Map<String, Set<String>>> protein_to_protein_by_domains; // one-sided reachable protein2s over given domains in protein1, not in there if no non-FB domain
	private final Set<String> decay_proteins;
	
	/**
	 * Constructor, builds objects from files
	 * @param ppin_file
	 * @param ddin_file
	 * @param protein_to_assumed_transcript_file
	 */
	public RewiringDetectorSample(String ppin_file, String ddin_file, String protein_to_assumed_transcript_file) {
		// set interactome outcome
		this.interactions = PPIN.readInteractionsFromPPINFile(ppin_file);
		
		// build AS-relevant domain info
		this.protein_to_domains = new HashMap<String, Set<String>>(1024);
		this.protein_to_protein_by_domains = new HashMap<String,  Map<String, Set<String>>>(1024);
		DDIN ddin = new DDIN(ddin_file);
		
		for (String protein1:ddin.getProtein_to_domains().keySet()) {
			for (String domain1:ddin.getProtein_to_domains().get(protein1)) {
				
				// domains have format 0|FB|Q9UKT9
				String domain_family = domain1.split("\\|")[1];
				
				// only non-FB domains
				if (domain_family.equals("FB"))
					continue;
				
				// ensure set is initialized and domain_family not yet processed for this protein
				if (!this.protein_to_domains.containsKey(protein1))
					this.protein_to_domains.put(protein1, new HashSet<String>(2));
				else if (this.protein_to_domains.get(protein1).contains(domain_family))
					continue;
				
				this.protein_to_domains.get(protein1).add(domain_family);
				
				for (String domain2:ddin.getDDIs().get(domain1)) {
					String protein2 = ddin.getDomain_to_protein().get(domain2);
					
					// ensure underlying datastructure is present, build if not
					if (!this.protein_to_protein_by_domains.containsKey(protein1))
						this.protein_to_protein_by_domains.put(protein1, new HashMap<String, Set<String>>(64));
					if (!this.protein_to_protein_by_domains.get(protein1).containsKey(protein2))
						this.protein_to_protein_by_domains.get(protein1).put(protein2, new HashSet<>(2));
					
					// note that protein2 is reachable from protein1 given that some domains of 1 are present 
					this.protein_to_protein_by_domains.get(protein1).get(protein2).add(domain_family);
				}
			}
		}
		
		// fill protein_to_assumed_transcript
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			
			if (spl.length >= 2)// to account for newer protein->assumed transcript files including expression data
				this.protein_to_assumed_transcript.put(spl[0], spl[1]); 
			else
				this.protein_to_assumed_transcript.put(spl[0], "unknown_transcript"); // NMD transcripts could bring up such artefacts with older versions of PPIXpress
		}
		
		// store decay transcripts: those that are in protein->transcript but have not even annotated FB domains
		this.decay_proteins = new HashSet<String>(this.protein_to_assumed_transcript.keySet());
		this.decay_proteins.removeAll(ddin.getProtein_to_domains().keySet());
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

	/**
	 * Returns a map of proteins in the network and the set of abundant non-FB domain types
	 * @return
	 */
	public Map<String, Set<String>> getProteinToUsedDomains() {
		return protein_to_domains;
	}

	/**
	 * Returns a map of proteins in the network and the proteins that can be connected on the DDIN using which non-F domain types
	 * @return
	 */
	public Map<String, Map<String, Set<String>>> getProteinToProteinByDomains() {
		return protein_to_protein_by_domains;
	}
	
	/**
	 * Returns all proteins that are not abundant due to known post-translational regulation on their most abundant transcript
	 * @return
	 */
	public Set<String> getDecayProteins() {
		return decay_proteins;
	}
	
	/*
	 * General purpose helpers
	 */

	/**
	 * Reads all matching PPIXpress output pairs of PPINs/major transcripts from a certain folder: 
	 * [folder]/[sample-string]_ppin.txt(.gz), [folder]/[sample-string]_ddin.txt(.gz) and [folder]/[sample-string]_major-transcripts.txt(.gz)
	 * Assumes the standard file endings "_ppin.txt(.gz), "_ddin.txt(.gz) and _major-transcripts.txt(.gz)".
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
			RewiringDetectorSample rds = new RewiringDetectorSample(pre + "_ppin.txt" + gz, pre + "_ddin.txt" + gz, pre + "_major-transcripts.txt" + gz);
			
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
