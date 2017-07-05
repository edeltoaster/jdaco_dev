package stem_cell_complexeomes;

import java.util.Set;

import framework.DataQuery;
import framework.RegulatoryNetwork;

public class check_pluri_diff_regnets {
	public static String path = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/";
	
	// testing and checking some stuff first
	public static void main(String[] args) {
		Set<String> allosome_proteins = DataQuery.getAllosomeProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println(allosome_proteins.size());
		
		System.out.println("pluri only");
		RegulatoryNetwork regnet = new RegulatoryNetwork(path + definitions.diff_compl_output_folder + "pluri_regnet_only.txt");
		System.out.println(regnet.getSizesStr());
		regnet.pruneToLargestSCCs();
		System.out.println(regnet.getSizesStr());
		regnet.removeProteinSet(allosome_proteins);
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/pluri_SCC.txt");
		
		System.out.println("plurisub only");
		RegulatoryNetwork subregnet = new RegulatoryNetwork(path + definitions.diff_compl_output_folder + "plurisub_regnet_only.txt");
		System.out.println(subregnet.getSizesStr());
		subregnet.pruneToLargestSCCs();
		System.out.println(subregnet.getSizesStr());
		subregnet.removeProteinSet(allosome_proteins);
		System.out.println(subregnet.getSizesStr());
		subregnet.writeRegulatoryNetwork("/Users/tho/Desktop/plurisub_SCC.txt");
		subregnet.writeNodeTable("/Users/tho/Desktop/plurisub_SCC_node.txt");
	}
}
