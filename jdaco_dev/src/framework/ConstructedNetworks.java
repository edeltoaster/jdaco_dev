package framework;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Stores matching PPIN/DDIN and the exact Ensembl database used to retrieve the data
 * @author Thorsten Will
 */
public class ConstructedNetworks {
	
	private final PPIN ppi;
	private final DDIN ddi;
	private final Map<String, String> protein_to_assumed_transcript;
	private Map<String, Float> transcript_abundance = null;
	private final String db;
	private final boolean isoform_based;
	
	/**
	 * Constructor used by NetworkBuilder (no abundance information)
	 * @param ppi
	 * @param ddi
	 * @param protein_to_assumed_isoform
	 * @param db
	 * @param isoform_based
	 */
	public ConstructedNetworks(PPIN ppi, DDIN ddi, Map<String, String> protein_to_assumed_transcript, String db, boolean isoform_based) {
		this.ppi = ppi;
		this.ddi = ddi;
		this.protein_to_assumed_transcript = protein_to_assumed_transcript;
		this.db = db;
		this.isoform_based = isoform_based;
	}
	
	/**
	 * Constructor used by NetworkBuilder (with abundance information)
	 * @param ppi
	 * @param ddi
	 * @param protein_to_assumed_isoform
	 * @param transcript_abundance
	 * @param db
	 * @param isoform_based
	 */
	public ConstructedNetworks(PPIN ppi, DDIN ddi, Map<String, String> protein_to_assumed_transcript, Map<String, Float> transcript_abundance, String db, boolean isoform_based) {
		this.ppi = ppi;
		this.ddi = ddi;
		this.protein_to_assumed_transcript = protein_to_assumed_transcript;
		this.transcript_abundance = transcript_abundance;
		this.db = db;
		this.isoform_based = isoform_based;
	}
	
	/**
	 * Constructor to built object from files
	 * @param ppi_file
	 * @param ddi_file
	 * @param protein_to_assumed_isoform
	 * @param db
	 * @param isoform_based
	 */
	public ConstructedNetworks(String ppi_file, String ddi_file, String protein_to_assumed_transcript_file, String db, boolean isoform_based) {
		this.ppi = new PPIN(ppi_file);
		this.ddi = new DDIN(ddi_file);
		
		this.protein_to_assumed_transcript = new HashMap<String, String>(1024);
		boolean check = true;
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			this.protein_to_assumed_transcript.put(spl[0], spl[1]);
			
			// depending on information content of file
			if (check) {
				check = false;
				if (spl.length == 3)
					this.transcript_abundance = new HashMap<String, Float>(1024);
			}
			
			if (spl.length == 3)
				this.transcript_abundance.put(spl[1], Float.parseFloat(spl[2]));
		}
		
		this.db = db;
		this.isoform_based = isoform_based;
	}
	
	/**
	 * Returns constructed PPIN
	 * @return
	 */
	public PPIN getPPIN() {
		return this.ppi;
	}

	/**
	 * Returns the DDIN of the constructed PPIN
	 * @return
	 */
	public DDIN getDDIN() {
		return this.ddi;
	}
	
	/**
	 * Returns a map of proteins in the network and their most abundant transcript
	 * @return
	 */
	public Map<String, String> getProteinToAssumedTranscriptMap() {
		return this.protein_to_assumed_transcript;
	}
	
	/**
	 * Returns a map of most abundant transcripts per protein and their abundance/summarized abundance of all transcripts associated with the protein (depending on building process), null if not stored
	 * @return
	 */
	public Map<String, Float> getTranscriptAbundanceMap() {
		return this.transcript_abundance;
	}
	
	/**
	 * Writes the protein->transcript mapping to a file for later usage
	 * @param out_file
	 */
	public void writeProteinToAssumedTranscriptMap(String out_file) {
		List<String> to_write = new LinkedList<>();
		
		for (String protein:this.protein_to_assumed_transcript.keySet()) {
			String transcript = this.protein_to_assumed_transcript.get(protein);
			
			if (this.transcript_abundance != null)
				to_write.add( protein + " " + transcript + " " + this.transcript_abundance.get(transcript) );
			else
				to_write.add( protein + " " + transcript );
		}
		
		Utilities.writeEntries(to_write, out_file);
	}
	
	/**
	 * Returns Ensembl database that was used
	 * @return
	 */
	public String getDB() {
		return this.db;
	}
	
	/**
	 * Returns if the initial/holistic network was based on the domain-annotation of the major isoform only
	 * @return
	 */
	public boolean getIsoformBased() {
		return this.isoform_based;
	}
	
	/**
	 * Returns Ensembl release of database used
	 * @return
	 */
	public int getEnsemblRelease() {
		String[] split = this.db.split("_");
		return Integer.parseInt(split[split.length-2]);
	}
	
	/**
	 * Returns organism's assembly version of the Ensembl release
	 * @return
	 */
	public String getAssembly() {
		String[] split = this.db.split("_");
		return split[split.length-1];
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((db == null) ? 0 : db.hashCode());
		result = prime * result + ((ddi == null) ? 0 : ddi.hashCode());
		result = prime * result + (isoform_based ? 1231 : 1237);
		result = prime * result + ((ppi == null) ? 0 : ppi.hashCode());
		result = prime * result
				+ ((protein_to_assumed_transcript == null) ? 0 : protein_to_assumed_transcript.hashCode());
		result = prime * result + ((transcript_abundance == null) ? 0 : transcript_abundance.hashCode());
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
		ConstructedNetworks other = (ConstructedNetworks) obj;
		if (db == null) {
			if (other.db != null)
				return false;
		} else if (!db.equals(other.db))
			return false;
		if (ddi == null) {
			if (other.ddi != null)
				return false;
		} else if (!ddi.equals(other.ddi))
			return false;
		if (isoform_based != other.isoform_based)
			return false;
		if (ppi == null) {
			if (other.ppi != null)
				return false;
		} else if (!ppi.equals(other.ppi))
			return false;
		if (protein_to_assumed_transcript == null) {
			if (other.protein_to_assumed_transcript != null)
				return false;
		} else if (!protein_to_assumed_transcript.equals(other.protein_to_assumed_transcript))
			return false;
		if (transcript_abundance == null) {
			if (other.transcript_abundance != null)
				return false;
		} else if (!transcript_abundance.equals(other.transcript_abundance))
			return false;
		return true;
	}
	
}
