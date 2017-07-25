package diff_complexeomes;

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
		RegulatoryNetwork regnet = new RegulatoryNetwork(path + definitions.diff_compl_output_folder + "pluri_regnet.txt");
		System.out.println(regnet.getSizesStr());
		regnet.pruneToLargestSCCs();
		System.out.println(regnet.getSizesStr());
		regnet.removeProteinSet(allosome_proteins);
		System.out.println(regnet.getSizesStr());
//		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/pluri_SCC.txt");
//		regnet.writeNodeTable("/Users/tho/Desktop/pluri_SCC_node.txt");
//		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/pluri_SCC2.txt", 2);
		
		System.out.println("plurisub only");
		RegulatoryNetwork subregnet = new RegulatoryNetwork(path + definitions.diff_compl_output_folder + "plurisub_regnet.txt");
		System.out.println(subregnet.getSizesStr());
		subregnet.pruneToLargestSCCs();
		System.out.println(subregnet.getSizesStr());
		subregnet.removeProteinSet(allosome_proteins);
		System.out.println(subregnet.getSizesStr());
//		subregnet.writeRegulatoryNetwork("/Users/tho/Desktop/plurisub_SCC.txt");
//		subregnet.writeNodeTable("/Users/tho/Desktop/plurisub_SCC_node.txt");
		
		System.out.println("nonpluri only");
		RegulatoryNetwork nonregnet = new RegulatoryNetwork(path + definitions.diff_compl_output_folder + "nonpluri_regnet.txt");
		System.out.println(nonregnet.getSizesStr());
		nonregnet.pruneToLargestSCCs();
		System.out.println(nonregnet.getSizesStr());
		nonregnet.removeProteinSet(allosome_proteins);
		System.out.println(nonregnet.getSizesStr());
//		nonregnet.writeRegulatoryNetwork("/Users/tho/Desktop/nonpluri.txt");
//		nonregnet.writeNodeTable("/Users/tho/Desktop/nonpluri_node.txt");
	}
}
