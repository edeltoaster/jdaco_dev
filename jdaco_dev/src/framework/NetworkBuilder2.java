package framework;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Implementation of the network construction method PPIXpress2
 * @author Thorsten Will
 */
public class NetworkBuilder2 {
	
	private PPIN original_ppi;
	private String organism_database;
	private boolean isoform_based;
	private Map<String, Set<String>> holistic_protein_domain_composition_map;
	private Map<String, List<String>> holistic_ddis = new HashMap<>();
	private HashMap<String, String> holistic_domain_to_protein = new HashMap<>();
	private Map<String, List<String>> holistic_protein_to_domains = new HashMap<>();
	private Map<String, LinkedList<String>> transcript_to_proteins = new HashMap<>();
	private Map<String, List<String>> transcript_to_domains;
	
	/**
	 * Constructs holistic network that can be used for many samples. Construction based on major transcripts.
	 * @param ppi
	 * @param isoform_based
	 */
	public NetworkBuilder2(PPIN ppi) {
		buildHolisticNetwork(ppi, true);
	}
	
	/**
	 * Constructs holistic network that can be used for many samples. Construction either based on major isoform or all transcripts.
	 * @param ppi
	 * @param isoform_based
	 */
	public NetworkBuilder2(PPIN ppi, boolean isoform_based) {
		buildHolisticNetwork(ppi, isoform_based);
	}
	
	/**
	 * Actually builds holistic network either based on major isoform or all transcripts
	 * @param ppi
	 * @param isoform_based
	 */
	private void buildHolisticNetwork(PPIN ppi, boolean isoform_based) {
		Set<String> proteins = ppi.getProteins();
		this.original_ppi = ppi;
		
		// retrieve annotational data, also for later usage
		this.organism_database = DataQuery.getEnsemblOrganismDatabaseFromProteins(proteins);
		this.isoform_based = isoform_based;
		if (this.isoform_based)
			this.holistic_protein_domain_composition_map = DataQuery.getIsoformProteinDomainMap(this.organism_database);
		else
			this.holistic_protein_domain_composition_map = DataQuery.getHolisticProteinDomainMap(this.organism_database);
		this.transcript_to_domains = DataQuery.getTranscriptsDomains(this.organism_database);
		for (String[] naming:DataQuery.getGenesTranscriptsProteins(this.organism_database)) {
			String transcript = naming[1];
			String protein = naming[2];
			if (!this.transcript_to_proteins.containsKey(transcript))
				this.transcript_to_proteins.put(transcript, new LinkedList<String>());
			this.transcript_to_proteins.get(transcript).add(protein);
		}
		
		// shrink to the ones of interest and adjust missing annotation
		this.holistic_protein_domain_composition_map.keySet().retainAll(proteins);
		for (String protein:proteins) {
			if (!this.holistic_protein_domain_composition_map.containsKey(protein))
				this.holistic_protein_domain_composition_map.put(protein, new HashSet<String>());
		}
		
		// temporary stored mapping
		Map<String, List<String>> domain_map = new HashMap<>();
		for (String protein:proteins) {
			holistic_protein_to_domains.put(protein, new LinkedList<String>());
			for (String domain_type:this.holistic_protein_domain_composition_map.get(protein)) {
				String domain_id = domain_type+"|"+protein;
				holistic_protein_to_domains.get(protein).add(domain_id);
				// add relations
				holistic_domain_to_protein.put(domain_id, protein);
				if (!domain_map.containsKey(domain_type))
					domain_map.put(domain_type, new LinkedList<String>());
				domain_map.get(domain_type).add(domain_id);
			}
		}
		
		// determine interactions found in the holistic network
		Map<String, List<String>> knownDDIs = DataQuery.getKnownDDIs();
		for (String domain_name:domain_map.keySet()) {
			if (!knownDDIs.containsKey(domain_name))
				continue;
			for (String target_domain_name:knownDDIs.get(domain_name)) {
				
				// no double occurance since undirected, don't think about domain types that are not there
				if (target_domain_name.compareTo(domain_name) < 0 || !domain_map.containsKey(target_domain_name))
					continue;
				
				for (String domain1:domain_map.get(domain_name))
					for (String domain2:domain_map.get(target_domain_name)) {
						String protein1 = holistic_domain_to_protein.get(domain1);
						String protein2 = holistic_domain_to_protein.get(domain2);
						
						// no self-interactions
						if (protein1.equals(protein2))
							continue;
						// no edge twice
						if (target_domain_name.equals(domain_name) && protein1.compareTo(protein2) < 0)
							continue;
						// not in PPI
						if (!ppi.getPartners().get(protein1).contains(protein2))
							continue;
						
						// if all checks passed: add DDI
						if (!holistic_ddis.containsKey(domain1))
							holistic_ddis.put(domain1, new LinkedList<String>());
						holistic_ddis.get(domain1).add(domain2);
						if (!holistic_ddis.containsKey(domain2))
							holistic_ddis.put(domain2, new LinkedList<String>());
						holistic_ddis.get(domain2).add(domain1);
					}
			}
		}
		
		// add fallback-interactions where it cannot be explained by any DDI
		for (String protein1:proteins)
			for (String protein2:ppi.getPartners().get(protein1)) {
				// do everything only once
				if (protein1.compareTo(protein2) < 0)
					continue;
				
				// is there a domain interaction supporting the interaction between the proteins?
				boolean support = false;
				for (String domain1:holistic_protein_to_domains.get(protein1)) {
					if (support) // break will only break the inner loop, thus support is caught in the outer one
						break;
					if (!holistic_ddis.containsKey(domain1))
						continue;
					for (String domain2:holistic_ddis.get(domain1))
						if (holistic_domain_to_protein.get(domain2).equals(protein2)) {
							support = true;
							break;
						}
				}
				
				// add fallback/artificial domains and interactions
				if (!support) {
					String fb1 = "FB|"+protein1;
					String fb2 = "FB|"+protein2;
					
					// add specific FB-domain to both, if needed
					if (!holistic_protein_to_domains.get(protein1).contains(fb1)) {
						holistic_protein_to_domains.get(protein1).add(fb1);
						holistic_domain_to_protein.put(fb1, protein1);
					}
					if (!holistic_protein_to_domains.get(protein2).contains(fb2)) {
						holistic_protein_to_domains.get(protein2).add(fb2);
						holistic_domain_to_protein.put(fb2, protein2);
					}
					
					// add FB<->FB interactions
					if (!holistic_ddis.containsKey(fb1))
						holistic_ddis.put(fb1, new LinkedList<String>());
					holistic_ddis.get(fb1).add(fb2);
					if (!holistic_ddis.containsKey(fb2))
						holistic_ddis.put(fb2, new LinkedList<String>());
					holistic_ddis.get(fb2).add(fb1);
				}
			}
		
		
		// remove unnecessary domains from data structures
		
		// determine domains without interactions
		List<String> removable_domains = new LinkedList<>();
		for (String domain:holistic_domain_to_protein.keySet())
			if (!holistic_ddis.containsKey(domain))
				removable_domains.add(domain);
		
		// remove them
		for (String domain:removable_domains) {
			// remove from protein -> domains
			holistic_protein_to_domains.get(holistic_domain_to_protein.get(domain)).remove(domain);
			// remove from domain -> protein
			holistic_domain_to_protein.remove(domain);
		}
	}
	
	/**
	 * Returns the underlying holistic DDIN
	 * @return
	 */
	public DDIN getHolisticDDIN() {
		return new DDIN(this.holistic_ddis, this.holistic_protein_to_domains, this.holistic_domain_to_protein);
	}
	
	/**
	 * Construct specific PPIN and DDIN from a map of abundant proteins and their expressed transcripts
	 * @param transcript_abundance
	 * @param remove_decayed
	 * @return
	 */
	public ConstructedNetworks constructAssociatedNetworksFromTranscriptAbundance(Map<String, Float> transcript_abundance, boolean remove_decayed) {
		// map abundant transcripts to proteins and those again to highest expressed transcript
		Map<String, String> major_transcript = new HashMap<>();
		for (String transcript:transcript_abundance.keySet()) {
			if (!this.transcript_to_proteins.containsKey(transcript)) // if no protein coding transcript
				continue;
			for (String protein:this.transcript_to_proteins.get(transcript)) {
				if (!major_transcript.containsKey(protein))
					major_transcript.put(protein, transcript);
				else { // compare case
					if (transcript_abundance.get(transcript) > transcript_abundance.get(major_transcript.get(protein)))
						major_transcript.put(protein, transcript);
				}
			}
		}
		
		// map that stores p1<->p2
		HashMap<String, Set<String>> ppi_partners = new HashMap<>(1024);
		HashMap<StrPair, Double> weights = new HashMap<>(8192);
		// maps that store DDI-stuff
		Map<String, List<String>> ddis = new HashMap<>(8192);// will later be converted to fixed data structure
		Map<String, List<String>> protein_to_domains = new HashMap<>(1024);// will later be converted to fixed data structure
		HashMap<String, String> domain_to_protein = new HashMap<>(8192);// HashMap since already final
		
		return null;
	}
	
	/**
	 * Returns the Ensembl DB that was used to gather the data
	 * @return
	 */
	public String getDB() {
		return this.organism_database;
	}
	
	/**
	 * Return fraction of PPIs associated to non-artificial domain interactions
	 * @return
	 */
	public double getMappingPercentage() {
		Set<StrPair> DDI_backed = new HashSet<>();
		
		for (String d1:this.holistic_ddis.keySet())
			for (String d2:this.holistic_ddis.get(d1)) {
				if (d1.compareTo(d2) <= 0)
					continue;
				
				String[] temp1 = d1.split("\\|");
				String[] temp2 = d2.split("\\|");
				
				if (temp1[0].equals("FB"))
					continue;
				
				DDI_backed.add(new StrPair(temp1[1], temp2[1]));
			}
		
		return DDI_backed.size() / (double) this.original_ppi.getInteractionsFast().size();
	}
	
	/**
	 * Return fraction of PPIs associated to non-artificial domain interactions (subset)
	 * @return
	 */
	public double getMappingPercentage(Set<String> proteins) {
		Set<StrPair> DDI_backed = new HashSet<>();
		Set<StrPair> all = new HashSet<>();
		
		// build reference
		for (String protein:proteins) {
			if (!this.original_ppi.contains(protein))
				continue;
			for (String protein2:this.original_ppi.getPartners().get(protein))
				all.add(new StrPair(protein, protein2));
		}
		
		// get backed portion
		for (String d1:this.holistic_ddis.keySet())
			for (String d2:this.holistic_ddis.get(d1)) {
				if (d1.compareTo(d2) <= 0)
					continue;
				
				String[] temp1 = d1.split("\\|");
				String[] temp2 = d2.split("\\|");
				
				if (temp1[0].equals("FB"))
					continue;
				
				if (proteins.contains(temp1[1]) || proteins.contains(temp2[1]))
					DDI_backed.add(new StrPair(temp1[1], temp2[1]));
			}
		
		return DDI_backed.size() / (double) all.size();
	}
	
	/**
	 * Return fraction of proteins with non-artificial domains that contribute to the mapping
	 * @return
	 */
	public double getMappingDomainPercentage() {
		int proteins_with_domain_annotation = 0;
		for (String protein:this.holistic_protein_to_domains.keySet())
			for (String domain:this.holistic_protein_to_domains.get(protein)) {
				String[] temp = domain.split("\\|");
				if (!temp[0].equals("FB")) {
					proteins_with_domain_annotation++;
					break;
				}
			}
		return proteins_with_domain_annotation / (double) this.original_ppi.getProteins().size();
	}
	
	/**
	 * Return fraction of proteins with non-artificial domains that contribute to the mapping (subset)
	 * @return
	 */
	public double getMappingDomainPercentage(Set<String> proteins) {
		int proteins_with_domain_annotation = 0;
	
		// define reference set
		Set<String> included = new HashSet<String>(this.original_ppi.getProteins());
		included.retainAll(proteins);
		
		for (String protein:this.holistic_protein_to_domains.keySet()) {
			if (!included.contains(protein))
				continue;
			for (String domain:this.holistic_protein_to_domains.get(protein)) {
				String[] temp = domain.split("\\|");
				if (!temp[0].equals("FB")) {
					proteins_with_domain_annotation++;
					break;
				}
			}
		}
		return proteins_with_domain_annotation / (double ) included.size();
	}
}
