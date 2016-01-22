package PPIComp_hemato;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks {
	static String BLUEPRINT_expr_folder = "/Users/tho/Desktop/BLUEPRINT_expr/";
	static String network_folder_pre = "/Users/tho/Desktop/BLUEPRINT_networks/";
	static Map<String, String> folder_type_map = new HashMap<>();
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	static List<String> parameters = new LinkedList<>();
	
	public static void loadAndStoreReferenceNetwork(String network_out) {
		PPIN ppin = DataQuery.getMenthaNetwork("9606");
		System.out.println(ppin.getSizesStr());
		ppin = ppin.updateUniprotAccessions();
		ppin.writePPIN(network_out);
		System.out.println(ppin.getSizesStr());
	}
	
	public static Set<String> checkExpressedGenes(List<Map<String, Float>> gene_data) {
		
		// count across all samples
		Map<String, Integer> count_map = new HashMap<>();
		for (Map<String, Float> expr:gene_data) {
			for (String s:expr.keySet()) {
				if (!count_map.containsKey(s))
					count_map.put(s, 1);
				else
					count_map.put(s, count_map.get(s)+1);
			}
		}
		
		Set<String> result = new HashSet<>();
		
		// org. publication counts transcripts only if at least in 2 samples
		for (String s:count_map.keySet()) {
			int count = count_map.get(s);
			if (count > 1)
				result.add(s);
		}
		
		return result;
	}
	
	public static Set<String> checkExpressedProtCodingTranscripts(List<Map<String, Float>> transcr_data) {
		
		Set<String> pc_trans = new HashSet<>();
		for (String[] d:DataQuery.getGenesTranscriptsProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"))) {
			pc_trans.add(d[1]);
		}
		
		// count across all samples
		Map<String, Integer> count_map = new HashMap<>();
		for (Map<String, Float> expr:transcr_data) {
			for (String s:expr.keySet()) {
				if (!pc_trans.contains(s))
					continue;
				if (!count_map.containsKey(s))
					count_map.put(s, 1);
				else
					count_map.put(s, count_map.get(s)+1);
			}
		}
		
		Set<String> result = new HashSet<>();
		
		// org. publication counts transcripts only if at least in 2 samples
		for (String s:count_map.keySet()) {
			int count = count_map.get(s);
			if (count > 1)
				result.add(s);
		}
		
		return result;
	}
	
	public static Set<String> checkExpressedProtCodingGenes(List<Map<String, Float>> gene_data) {
		
		Set<String> pc_genes = new HashSet<>();
		for (String[] d:DataQuery.getGenesTranscriptsProteins(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"))) {
			pc_genes.add(d[0]);
		}
		
		// count across all samples
		Map<String, Integer> count_map = new HashMap<>();
		for (Map<String, Float> expr:gene_data) {
			for (String s:expr.keySet()) {
				if (!pc_genes.contains(s))
					continue;
				if (!count_map.containsKey(s))
					count_map.put(s, 1);
				else
					count_map.put(s, count_map.get(s)+1);
			}
		}
		
		Set<String> result = new HashSet<>();
		
		// org. publication counts transcripts only if at least in 2 samples
		for (String s:count_map.keySet()) {
			int count = count_map.get(s);
			if (count > 1)
				result.add(s);
		}
		
		return result;
	}
	
	public static void preprocess() {
		for (String s:Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/eval_code/PPIComp_hemato/cell_types.txt")) {
			if (s.startsWith("#"))
				continue;
			String[] spl = s.trim().split(" ");
			folder_type_map.put(spl[0], spl[1]);
		}
		
		System.out.println("Original PPIN: " + "mixed_data/human_mentha_17_jan.txt.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did:" + DataQuery.get3didVersion());
		System.out.println("iPfam:" + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/human_mentha_17_jan.txt.gz");
		builder = new NetworkBuilder(original_ppin);
	}
	
	public static void process(double TPM_threshold) {
		
		parameters.add(Double.toString(TPM_threshold));
		
		System.out.println("BLUEPRINT net-builder, TPM threshold: " + TPM_threshold);
		String network_folder = network_folder_pre + TPM_threshold + "/";
		
		new File(network_folder).mkdir();
		
		Map<String, List<String>> data_map = new HashMap<>();
		List<Map<String, Float>> gene_abundance_data = new LinkedList<>();
		List<Map<String, Float>> pc_transcr_abundance_data = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(BLUEPRINT_expr_folder, ".rsem.tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			
			String cell_type = path_split[path_split.length-2];
			
			if (!folder_type_map.containsKey(cell_type))
				continue;
			
			cell_type = folder_type_map.get(cell_type);
			
			// filter N, M, NK, CD4 to venous blood samples
			if (cell_type.equals("M") || cell_type.equals("N") || cell_type.equals("NK") || cell_type.equals("CD4"))
				if (file_name.startsWith("Cord"))
					continue;
			
			String out_path = network_folder + cell_type + "/";
			
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, TPM_threshold);
			Map<String, Float> gene_expr = TranscriptAbundanceReader.readRSEMGenesTPM(path, TPM_threshold);
			
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(transcr_expr);
			cn.getPPIN().writePPIN(out_path + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(out_path + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(out_path + file_name + "_major-transcripts.txt.gz");
			
			// write path of network to data_map
			if (!data_map.containsKey(cell_type)) {
				data_map.put(cell_type, new LinkedList<>());
			}
			
			data_map.get(cell_type).add(out_path + file_name + "_ppin.txt.gz");
			
			// store expression data for all in publication
			if (cell_type.equals("CD4") || cell_type.equals("N") || cell_type.equals("M"))
				continue;
			gene_abundance_data.add(gene_expr);
			pc_transcr_abundance_data.add(transcr_expr);
		}
		
		System.out.println();
		
		for (String cell_type:data_map.keySet()) {

			List<Double> no_proteins = new LinkedList<>();
			List<Double> no_interactions = new LinkedList<>();
			
			for (String path:data_map.get(cell_type)) {
				PPIN ppi = new PPIN(path);
				int[] sizes = ppi.getSizes();
				no_proteins.add( (double) sizes[0] );
				no_interactions.add( (double) sizes[1]);
			}
			
			int prots = (int) Utilities.getMean(no_proteins);
			System.out.println();
			System.out.println(cell_type + ": " + data_map.get(cell_type).size() + " samples");
			System.out.println("Size: " + prots + "+-" + (int) Utilities.getStd(no_proteins) + " / " + (int) Utilities.getMean(no_interactions) + "+-" + (int) Utilities.getStd(no_interactions));
		}
		
//		System.out.println();
//		System.out.println(checkExpressedGenes(gene_abundance_data).size() + " expressed genes according to counting as in BLUEPRINT paper.");
//		System.out.println(checkExpressedProtCodingGenes(gene_abundance_data).size() + " expressed prot-c. genes according to counting as in BLUEPRINT paper.");
//		System.out.println(checkExpressedGenes(pc_transcr_abundance_data).size() + " expressed prot-coding transcr. according to counting as in BLUEPRINT paper.");
//		System.out.println();
	}
	
	public static void main(String[] args) {
		//loadAndStoreReferenceNetwork("mixed_data/human_mentha_19_jan.txt.gz");
		//System.exit(0);
		
		preprocess();
		
		new File(network_folder_pre).mkdir();
		
		for (double thr = 0.0; thr <= 1.0; thr += 0.01) {
			process(thr);
		}
		
	}
}
