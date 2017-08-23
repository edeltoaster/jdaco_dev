package sc_bonn;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.PPIN;
import framework.Utilities;

public class scGO {

	public static void writeOut(String file, Map<String, Set<String>> data) {
		List<String> lines = new LinkedList<String>();
		for (String term:data.keySet()) {
			String temp = term + "\t"+ String.join(",", data.get(term));
			lines.add(temp);
		}
		Utilities.writeEntries(lines, file);
	}

	public static void main(String[] args) {
		Set<String> ref_proteins = new PPIN("mixed_data/human_merged.tsv.gz").getProteins();
		Map<String, Set<String>> prot_map = new HashMap<String, Set<String>>();
		String human = "9606";
		Set<String> temp;
		String out_file = "/Users/tho/Desktop/mixed_GO.tsv";
		
		// http://www.ncbi.nlm.nih.gov/pubmed/23715547
		
		// glycolytic process GO:0006096
		System.out.println("glycolytic process");
		temp = DataQuery.getProteinsWithGO("GO:0006096", human);
		temp.retainAll(ref_proteins);
		prot_map.put("glycolytic_process", temp);
		
		// oxidative phosphorylation GO:0006119
		System.out.println("oxidative phosphorylation");
		temp = DataQuery.getProteinsWithGO("GO:0006119", human);
		temp.retainAll(ref_proteins);
		prot_map.put("oxidative_phosphorylation", temp);
		
		// response to oxidative stress GO:0006979
		System.out.println("response to oxidative stress");
		temp = DataQuery.getProteinsWithGO("GO:0006979", human);
		temp.retainAll(ref_proteins);
		prot_map.put("response_to_oxidative_stress", temp);
		
		// others
		
		// cell-cycle GO:0007049
		System.out.println("cell-cycle");
		temp = DataQuery.getProteinsWithGO("GO:0007049", human);
		temp.retainAll(ref_proteins);
		prot_map.put("cell_cycle", temp);
		
		// sequence-specific DNA binding transcription factor activity GO:0003700
		System.out.println("TFs");
		temp = DataQuery.getProteinsWithGO("GO:0003700", human);
		temp.retainAll(ref_proteins);
		prot_map.put("TFs", temp);
		
		// GO:0007165 signal transduction
		System.out.println("signaling");
		temp = DataQuery.getProteinsWithGO("GO:0007165", human);
		temp.retainAll(ref_proteins);
		prot_map.put("signal_transduction", temp);
		
		// GO:0045254 pyruvate dehydrogenase complex -> connects glycolysis and TCA
		System.out.println("pyruvate dehydrogenase omplex");
		temp = DataQuery.getProteinsWithGO("GO:0045254", human);
		temp.retainAll(ref_proteins);
		prot_map.put("pyruvate_dehydrogenase_complex", temp);
		
		writeOut(out_file, prot_map);
	}

}
