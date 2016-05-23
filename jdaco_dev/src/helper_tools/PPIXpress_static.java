package helper_tools;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;

/**
 * PPIXpress helper tool for Cecilia Hernandez Rivas :-)
 * @author Thorsten Will
 */
public class PPIXpress_static {

	
	public static void printHelp() {
		System.out.println("Static PPIXpress PPIN/DDIN construction");
		System.out.println("usage: java -jar PPIXpress_static.jar [INPUT-NETWORK]");
		System.out.println();
		System.out.println("[INPUT-NETWORK]:");
		System.out.println("    A suitable PPIN file using UniProt Accessions or taxon:[NCBI-taxon] (mentha retrieval),");
		System.out.println("    e.g. taxon:559292 for the yeast S.cerevisiae S288C");
		System.exit(0);
	}
	
	public static void main(String[] args) {
		
		if (args.length < 1)
			printHelp();
		
		String original_network_path = args[0];
		
		// load PPIN
		PPIN original_network = null;
		if (original_network_path.startsWith("taxon:")) {
			String taxon_id = original_network_path.split(":")[1];
			System.out.print("Retrieving mentha interaction data for taxon " + taxon_id + " ... ");
			original_network = DataQuery.getMenthaNetwork(taxon_id, System.out);
			
			if (original_network.getSizes()[0] == 0) {
				System.out.println("no interaction data for taxon " + taxon_id + " available in mentha.");
				System.out.print("Retrieving IntAct interaction data instead ... ");
				original_network = DataQuery.getIntActNetwork(taxon_id, System.out);
			}
			System.out.println("done.");
			
		} else {
			System.out.println("Reading " + original_network_path + " (may take some time if ID conversion is necessary) ... ");
			original_network = new PPIN(original_network_path);
		}
		
		
		// retrieve annotation data and output versions
		String organism_database = DataQuery.getEnsemblOrganismDatabaseFromProteins(original_network.getProteins());
		String ensembl_version = organism_database.split("_")[organism_database.split("_").length-2];
		System.out.flush();
		
		System.out.print("Retrieving ENSEMBL " + ensembl_version + " data from database " + organism_database + " (may take some minutes) ... ");
		DataQuery.getGenesTranscriptsProteins(organism_database);
		System.out.print("50% ... ");
		System.out.flush();
		DataQuery.getIsoformProteinDomainMap(organism_database);
		System.out.println("100%");
		System.out.flush();
		
		System.out.print("Retrieving current interaction data fom 3did (" + DataQuery.get3didVersion() + ") ...");
		DataQuery.getKnownDDIs();
		System.out.println("100%");
		System.out.flush();
		
		
		// build DDIN
		ConstructedNetworks cn = NetworkBuilder.constructAssociatedIsoformNetworks(original_network);
		
		// write PPIN and DDIN
		cn.getPPIN().writePPIN("ppin.tsv");
		cn.getDDIN().writeDDIN("ddin.tsv");
		
		System.out.println("Files written to ppin.tsv and ddin.tsv");
	}
}
