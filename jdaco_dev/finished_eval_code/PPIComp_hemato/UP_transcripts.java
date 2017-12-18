package PPIComp_hemato;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import framework.DataQuery;
import framework.Utilities;

public class UP_transcripts {
	
	public static void main(String[] args) {
		// study data was build using v83
		DataQuery.enforceSpecificEnsemblRelease("83");
		
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		
		Set<String> protein_transcripts = new HashSet<>(DataQuery.getGenesTranscriptsProteins(db).stream().map(s->s[1]).collect(Collectors.toList()));
		System.out.println(protein_transcripts.size());
		
		Utilities.writeEntries(protein_transcripts, "/Users/tho/GDrive/Work/projects/hemato_rewiring/diffnet_analysis/scripts/quantitative/protein_transcripts.txt.gz");
	}
}
