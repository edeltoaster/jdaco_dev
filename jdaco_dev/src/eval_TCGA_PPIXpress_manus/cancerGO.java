package eval_TCGA_PPIXpress_manus;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.PPIN;
import framework.Utilities;

public class cancerGO {

	public static void writeOut(String file, Map<String, Set<String>> data) {
		List<String> lines = new LinkedList<String>();
		for (String term:data.keySet()) {
			String temp = term + "\t"+ String.join(",", data.get(term));
			lines.add(temp);
		}
		Utilities.writeEntries(lines, file);
	}
	
	// from http://nar.oxfordjournals.org/content/42/22/13557/suppl/DC1
	public static void main(String[] args) {
		Set<String> ref_proteins = new PPIN("mixed_data/human_biogrid.txt.gz").getProteins();
		Map<String, Set<String>> hallmark_prot_map = new HashMap<String, Set<String>>();
		String human = "9606";
		Set<String> temp;
		String out_file = "/Users/tho/Desktop/hallmarks_biogrid_05_May_15.tsv";
		
		// Inducing Angiogenesis	 GO:0001525	 angiogenesis
		System.out.println("Inducing Angiogenesis");
		temp = DataQuery.getProteinsWithGO("GO:0001525", human);
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Inducing Angiogenesis", temp);
		
		// Enabling Replicative Immortality	 GO:0032202	 telomere assembly
		// Enabling Replicative Immortality	 GO:0000723	 telomere maintenance
		// Enabling Replicative Immortality	 GO:0090398	 cellular senescence
		// Enabling Replicative Immortality	 GO:0090399	 replicative senescence
		System.out.println("Enabling Replicative Immortality");
		temp = DataQuery.getProteinsWithGO("GO:0032202", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0000723", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0090398", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0090399", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Enabling Replicative Immortality", temp);
		
//		 Activating Invasion and Metastasis	 GO:0045216	 cell-cell junction organization
//		 Activating Invasion and Metastasis	 GO:0034329	 cell junction assembly
//		 Activating Invasion and Metastasis	 GO:0045217	 cell-cell junction maintenance
//		 Activating Invasion and Metastasis	 GO:0034334	 adherens junction maintenance
//		 Activating Invasion and Metastasis	 GO:0016477	 cell migration
//		 Activating Invasion and Metastasis	 GO:0010718	 positive regulation of epithelial to mesenchymal transition
//		 Activating Invasion and Metastasis	 GO:0007155	 cell adhesion
		System.out.println("Activating Invasion and Metastasis");
		temp = DataQuery.getProteinsWithGO("GO:0045216", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0034329", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0045217", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0034334", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0016477", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0010718", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0007155", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Activating Invasion and Metastasis", temp);
		
//		 Genome Instability and Mutation	 GO:0006281	 DNA repair
//		 Genome Instability and Mutation	 GO:0051383	 kinetochore organization
//		 Genome Instability and Mutation	 GO:0007062	 sister chromatid cohesion
//		 Genome Instability and Mutation	 GO:0000819	 sister chromatid segregation
//		 Genome Instability and Mutation	 GO:0051988	 regulation of attachment of spindle microtubules to kinetochore
//		 Genome Instability and Mutation	 GO:0030997	 regulation of centriole-centriole cohesion
//		 Genome Instability and Mutation	 GO:0046605	 regulation of centrosome cycle
//		 Genome Instability and Mutation	 GO:0060236	 regulation of mitotic spindle organization
//		 Genome Instability and Mutation	 GO:0090169	 regulation of spindle assembly
//		 Genome Instability and Mutation	 GO:0043146	 spindle stabilization
//		 Genome Instability and Mutation	 GO:0031577	 spindle checkpoint
		System.out.println("Genome Instability and Mutation");
		temp = DataQuery.getProteinsWithGO("GO:0006281", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0051383", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0007062", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0000819", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0051988", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0030997", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0046605", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0060236", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0090169", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0043146", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0031577", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Genome Instability and Mutation", temp);
		
//		 Resisting Cell Death	 GO:0060548	 negative regulation of cell death
//		 Resisting Cell Death	 GO:0012501	 programmed cell death
//		 Resisting Cell Death	 GO:0010941	 regulation of cell death
		System.out.println("Resisting Cell Death");
		temp = DataQuery.getProteinsWithGO("GO:0060548", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0012501", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0010941", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Resisting Cell Death", temp);
		
//		 Deregulating Cellular Energetics	 GO:0006091	 generation of precursor metabolites and energy
		System.out.println("Deregulating Cellular Energetics");
		temp = DataQuery.getProteinsWithGO("GO:0006091", human);
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Deregulating Cellular Energetics", temp);
		
//		 Sustaining Proliferative Signaling	 GO:0007166	 cell surface receptor signaling pathway
//		 Sustaining Proliferative Signaling	 GO:0070848	 response to growth factor
		System.out.println("Sustaining Proliferative Signaling");
		temp = DataQuery.getProteinsWithGO("GO:0007166", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0070848", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Sustaining Proliferative Signaling", temp);
		
//		 Tumor-Promoting Inflammation	 GO:0006954	 inflammatory response
//		 Tumor-Promoting Inflammation	 GO:0045321	 leukocyte activation
		System.out.println("Tumor-Promoting Inflammation");
		temp = DataQuery.getProteinsWithGO("GO:0006954", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0045321", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Tumor-Promoting Inflammation", temp);
		
//		 Avoiding Immune Destruction	 GO:0002507	 tolerance induction
//		 Avoiding Immune Destruction	 GO:0001910	 regulation of leukocyte mediated cytotoxicity
//		 Avoiding Immune Destruction	 GO:0019882	 antigen processing and presentation
//		 Avoiding Immune Destruction	 GO:0002767	 immune response-inhibiting cell surface receptor signaling pathway
		System.out.println("Avoiding Immune Destruction");
		temp = DataQuery.getProteinsWithGO("GO:0002507", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0001910", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0019882", human));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0002767", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Avoiding Immune Destruction", temp);
		
//		 Evading Growth Suppressors	 GO:0007049	 cell cycle
//		 Evading Growth Suppressors	 GO:0008283	 cell proliferation
		System.out.println("Evading Growth Suppressors");
		temp = DataQuery.getProteinsWithGO("GO:0007049", human);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0008283", human));
		temp.retainAll(ref_proteins);
		hallmark_prot_map.put("Evading Growth Suppressors", temp);
		
		writeOut(out_file, hallmark_prot_map);
	}

}
