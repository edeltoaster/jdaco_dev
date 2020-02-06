package mixed;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		Map<String, String> ELM_motifs = DataQuery.getELMMotifs();
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		Map<String, String> trans_protseq = DataQuery.getTranscriptsProteinSeq(db);
		
		for (String regexp : ELM_motifs.values()) {
			System.out.println("checking " + regexp);
			
			Pattern p = Pattern.compile(regexp);
			for (String seq : trans_protseq.values()) {
				Matcher m = p.matcher(seq);
				List<Integer> hits = new LinkedList<>();
				while (m.find()) {
					hits.add(m.start());
				}
				
				if (hits.size() > 0) {
					System.out.println(regexp);
					System.out.println(seq);
					System.out.println(hits);
					System.exit(0);
				}

			}
			
		}
	}
}
