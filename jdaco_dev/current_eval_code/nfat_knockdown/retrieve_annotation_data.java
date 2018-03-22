package nfat_knockdown;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.Utilities;

public class retrieve_annotation_data {

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
		Map<String, Set<String>> annotation_map = new HashMap<String, Set<String>>();
		String human = "9606";
		Set<String> temp;
		String out_file = "/Users/tho/Desktop/annotation_data_iea.tsv";
		boolean include_IEA = true;
		
		// Inducing Angiogenesis	 GO:0001525	 angiogenesis
		System.out.println("Inducing Angiogenesis");
		temp = DataQuery.getProteinsWithGO("GO:0001525", human, include_IEA, false, true);
		annotation_map.put("Inducing Angiogenesis", temp);
		
		// Enabling Replicative Immortality	 GO:0032202	 telomere assembly
		// Enabling Replicative Immortality	 GO:0000723	 telomere maintenance
		// Enabling Replicative Immortality	 GO:0090398	 cellular senescence
		// Enabling Replicative Immortality	 GO:0090399	 replicative senescence
		System.out.println("Enabling Replicative Immortality");
		temp = DataQuery.getProteinsWithGO("GO:0032202", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0000723", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0090398", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0090399", human, include_IEA, false, true));
		annotation_map.put("Enabling Replicative Immortality", temp);
		
//		 Activating Invasion and Metastasis	 GO:0045216	 cell-cell junction organization
//		 Activating Invasion and Metastasis	 GO:0034329	 cell junction assembly
//		 Activating Invasion and Metastasis	 GO:0045217	 cell-cell junction maintenance
//		 Activating Invasion and Metastasis	 GO:0034334	 adherens junction maintenance
//		 Activating Invasion and Metastasis	 GO:0016477	 cell migration
//		 Activating Invasion and Metastasis	 GO:0010718	 positive regulation of epithelial to mesenchymal transition
//		 Activating Invasion and Metastasis	 GO:0007155	 cell adhesion
		System.out.println("Activating Invasion and Metastasis");
		temp = DataQuery.getProteinsWithGO("GO:0045216", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0034329", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0045217", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0034334", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0016477", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0010718", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0007155", human, include_IEA, false, true));
		annotation_map.put("Activating Invasion and Metastasis", temp);
		
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
		temp = DataQuery.getProteinsWithGO("GO:0006281", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0051383", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0007062", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0000819", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0051988", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0030997", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0046605", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0060236", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0090169", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0043146", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0031577", human, include_IEA, false, true));
		annotation_map.put("Genome Instability and Mutation", temp);
		
//		 Resisting Cell Death	 GO:0060548	 negative regulation of cell death
//		 Resisting Cell Death	 GO:0012501	 programmed cell death
//		 Resisting Cell Death	 GO:0010941	 regulation of cell death
		System.out.println("Resisting Cell Death");
		temp = DataQuery.getProteinsWithGO("GO:0060548", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0012501", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0010941", human, include_IEA, false, true));
		annotation_map.put("Resisting Cell Death", temp);
		
//		 Deregulating Cellular Energetics	 GO:0006091	 generation of precursor metabolites and energy
		System.out.println("Deregulating Cellular Energetics");
		temp = DataQuery.getProteinsWithGO("GO:0006091", human, include_IEA, false, true);
		annotation_map.put("Deregulating Cellular Energetics", temp);
		
//		 Sustaining Proliferative Signaling	 GO:0007166	 cell surface receptor signaling pathway
//		 Sustaining Proliferative Signaling	 GO:0070848	 response to growth factor
		System.out.println("Sustaining Proliferative Signaling");
		temp = DataQuery.getProteinsWithGO("GO:0007166", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0070848", human, include_IEA, false, true));
		annotation_map.put("Sustaining Proliferative Signaling", temp);
		
//		 Tumor-Promoting Inflammation	 GO:0006954	 inflammatory response
//		 Tumor-Promoting Inflammation	 GO:0045321	 leukocyte activation
		System.out.println("Tumor-Promoting Inflammation");
		temp = DataQuery.getProteinsWithGO("GO:0006954", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0045321", human, include_IEA, false, true));
		annotation_map.put("Tumor-Promoting Inflammation", temp);
		
//		 Avoiding Immune Destruction	 GO:0002507	 tolerance induction
//		 Avoiding Immune Destruction	 GO:0001910	 regulation of leukocyte mediated cytotoxicity
//		 Avoiding Immune Destruction	 GO:0019882	 antigen processing and presentation
//		 Avoiding Immune Destruction	 GO:0002767	 immune response-inhibiting cell surface receptor signaling pathway
		System.out.println("Avoiding Immune Destruction");
		temp = DataQuery.getProteinsWithGO("GO:0002507", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0001910", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0019882", human, include_IEA, false, true));
		temp.addAll(DataQuery.getProteinsWithGO("GO:0002767", human, include_IEA, false, true));
		annotation_map.put("Avoiding Immune Destruction", temp);
		
//		 Evading Growth Suppressors	 GO:0007049	 cell cycle
//		 Evading Growth Suppressors	 GO:0008283	 cell proliferation
		System.out.println("Evading Growth Suppressors");
		temp = DataQuery.getProteinsWithGO("GO:0007049", human, include_IEA, false, true);
		temp.addAll(DataQuery.getProteinsWithGO("GO:0008283", human, include_IEA, false, true));
		annotation_map.put("Evading Growth Suppressors", temp);
		
		
		// other terms
		System.out.println("other terms");
		temp = DataQuery.getProteinsWithGO("GO:0005739", human, include_IEA, false, true);
		annotation_map.put("Mitochondrion", temp);
		temp = DataQuery.getProteinsWithGO("GO:0016491", human, include_IEA, false, true);
		annotation_map.put("Oxidoreductase activity", temp);
		temp = DataQuery.getProteinsWithGO("GO:0016209", human, include_IEA, false, true);
		annotation_map.put("Antioxidant activity", temp);
		temp = DataQuery.getProteinsWithGO("GO:0015085", human, include_IEA, false, true);
		annotation_map.put("Calcium ion transmembrane transporter activity", temp);
		
		writeOut(out_file, annotation_map);
	}

}
