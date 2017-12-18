package diff_compl_Hoth;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class get_surface_proteins_genes_transcripts {
	
	
	public static void main(String[] args) {
		Set<String> surface_proteins = DataQuery.getProteinsWithGO("GO:0009986", "9606"); // GO for location: cell surface
		
		String db = DataQuery.getEnsemblOrganismDatabaseFromProteins(surface_proteins);
		System.out.println(db);
		List<String[]> conv_data = DataQuery.getGenesTranscriptsProteins(db);
		
		Set<String> surface_genes_transcripts = new HashSet<>();
		for (String[] tripel : conv_data) {
			String g = tripel[0];
			String tr = tripel[1];
			String p = tripel[2];
			if (surface_proteins.contains(p)) {
				surface_genes_transcripts.add(tr);
				surface_genes_transcripts.add(g);
			}
		}
		
		Utilities.writeEntries(surface_genes_transcripts, "/Users/tho/GDrive/Work/projects/capello/diff-expr/surface_23_11_17.txt.gz");
	}
}
