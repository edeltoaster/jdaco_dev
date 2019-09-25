package KDEL;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class retrieve_annotation_data2 {

	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		Map<String, String> uniprot_genename = DataQuery.getUniprotToGeneNameMap(db);
		System.out.println(uniprot_genename.size());
		Set<String> genes = new HashSet<>(uniprot_genename.values());
		System.out.println(genes.size());
		Utilities.writeEntries(genes, "/Users/tho/Desktop/all_ensembl_uniprot_genes.txt");
	}

}
