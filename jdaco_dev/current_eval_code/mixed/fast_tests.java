package mixed;

import java.util.List;
import java.util.Map;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		Map<String, List<String>> bla = DataQuery.getTranscriptsELMMotifs(db);
		
		System.out.println(bla.size());
		System.out.println(bla);
	}
}
