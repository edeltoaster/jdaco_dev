package mixed;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		DataQuery.get(db);
		// TODO: testing ground for protein sequence retrieval
	}
}
