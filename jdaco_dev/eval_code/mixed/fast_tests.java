package mixed;

import java.util.Map;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println(db);
		Map<String, Integer> bla = DataQuery.getTranscriptsCDNALength(db);
		System.out.println(bla.get("ENST00000380152")); // ref: 11986bp
		System.out.println(bla.get("ENST00000544455")); // ref: 10984bp
		System.out.println(bla.get("ENST00000325404")); // ref: 2513bp
	}
}
