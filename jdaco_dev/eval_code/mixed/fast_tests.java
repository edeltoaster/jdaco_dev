package mixed;

import java.util.List;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		//String db = DataQuery.getEnsemblOrganismDatabaseFromName("mus musculus");
		//String db = DataQuery.getEnsemblOrganismDatabaseFromName("saccharomyces cerevisiae");
		
		System.out.println(db);
		
		long start = System.currentTimeMillis();
		List<String[]> res = DataQuery.getGenesTranscriptsProteins(db);
		long end = System.currentTimeMillis();
		
		System.out.println((end-start)/1000 + "s");
		
		System.out.println(res.size());
	}
}
