package mixed;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		System.out.println(DataQuery.getTranscriptsDomains(db).size());
		System.out.println(DataQuery.getDecayTranscripts(db).size());
	}
}
