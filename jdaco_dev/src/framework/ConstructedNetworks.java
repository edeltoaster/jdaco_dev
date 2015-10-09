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
	private final String db;
	private final boolean isoform_based;
	
	/**
	 * Constructor to built object from files, protein_to_assumed_isoform will be empty
	 * @param ppi_file
	 * @param ddi_file
	 * @param db
	 * @param isoform_based
	 */
	public ConstructedNetworks(String ppi_file, String ddi_file, String db, boolean isoform_based) {
		this.ppi = new PPIN(ppi_file);
		this.ddi = new DDIN(ddi_file);
		this.protein_to_assumed_transcript = new HashMap<>();
		this.db = db;
		this.isoform_based = isoform_based;
	}
	
	/**
	 * Constructor used by NetworkBuilder
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
		
		this.protein_to_assumed_transcript = new HashMap<String, String>();
		for (String s:Utilities.readEntryFile(protein_to_assumed_transcript_file)) {
			String[] spl = s.trim().split("\\s+");
			this.protein_to_assumed_transcript.put(spl[0], spl[0]);
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
	 * Writes the protein->transcript mapping to a file for later usage
	 * @param out_file
	 */
	public void writeProteinToAssumedTranscriptMap(String out_file) {
		List<String> to_write = new LinkedList<>();
		
		for (String protein:this.protein_to_assumed_transcript.keySet())
			to_write.add( protein + " " + this.protein_to_assumed_transcript.get(protein));
		
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

}
