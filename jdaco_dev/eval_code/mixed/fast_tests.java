package mixed;

import java.util.List;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("mus musculus");
		System.out.println(db);
		
		long start = System.currentTimeMillis();
		List<String[]> res = DataQuery.getGenesTranscriptsProteins(db);
		long end = System.currentTimeMillis();
		
		System.out.println((end-start)/1000 + "s");
		
		System.out.println(res.size());
	}
}
