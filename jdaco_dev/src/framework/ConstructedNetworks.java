package framework;
/**
 * Stores matching PPIN/DDIN and the exact Ensembl database used to retrieve the data
 * @author Thorsten Will
 */
public class ConstructedNetworks {
	
	private PPIN ppi;
	private DDIN ddi;
	private String db;
	private boolean isoform_based;
	
	public ConstructedNetworks(String ppi_file, String ddi_file, String db, boolean isoform_based) {
		this.ppi = new PPIN(ppi_file);
		this.ddi = new DDIN(ddi_file);
		this.db = db;
		this.isoform_based = isoform_based;
	}
	
	public ConstructedNetworks(PPIN ppi, DDIN ddi, String db, boolean isoform_based) {
		this.ppi = ppi;
		this.ddi = ddi;
		this.db = db;
		this.isoform_based = isoform_based;
	}

	public PPIN getPPIN() {
		return this.ppi;
	}

	public DDIN getDDIN() {
		return this.ddi;
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
