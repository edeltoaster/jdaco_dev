package mixed;

import java.util.LinkedList;
import java.util.List;

import framework.DataQuery;
import framework.Utilities;

public class generateTranscriptGeneProteinMap {
	
	public static void main(String[] args) {
		DataQuery.enforceSpecificEnsemblRelease("84");
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		
		long start = System.currentTimeMillis();
		List<String[]> res = DataQuery.getGenesTranscriptsProteins(db);
		long end = System.currentTimeMillis();
		
		System.out.println((end-start)/1000 + "s");
		
		System.out.println(res.size());
		
		List<String> to_write = new LinkedList<>();
		for (String[] s:res) {
			to_write.add(s[0] + " " + s[1] + " " + s[2]);
		}
		
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/ensembl84_mapping.txt");
	}
}
