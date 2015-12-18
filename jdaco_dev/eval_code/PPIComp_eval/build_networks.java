package PPIComp_eval;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks {
	static String expr_folder = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BRCA_expr/";
	static String network_folder_pre = "/Users/tho/Dropbox/Work/projects/hemato_rewiring/BRCA_networks/";
	static Map<String, String> folder_type_map = new HashMap<>();
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void loadAndStoreReferenceNetwork(String network_out) {
		PPIN ppin = DataQuery.getIntActNetwork("9606");
		ppin = ppin.mergeAll(DataQuery.getIRefIndexNetwork("9606"));
		ppin = ppin.updateUniprotAccessions();
		ppin.writePPIN(network_out);
		System.out.println(ppin.getSizesStr());
	}
	
	public static void preprocess() {
		System.out.println("Original PPIN: " + "mixed_data/human_merged_dec_4.txt.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did:" + DataQuery.get3didVersion());
		System.out.println("iPfam:" + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("mixed_data/human_merged_dec_4.txt.gz");
		builder = new NetworkBuilder(original_ppin);
	}
	
	public static void process(double threshold) {
		
		System.out.println("BRCA net-builder, threshold: " + threshold);
		String network_folder = network_folder_pre + threshold + "/";
		
		new File(network_folder).mkdir();
		
		Map<String, List<String>> data_map = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".txt.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String cell_type = path_split[path_split.length-2];
			
			String out_path = network_folder + cell_type + "/";
			
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readTCGAIsoformRSEM(path, threshold);
			
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(transcr_expr);
			cn.getPPIN().writePPIN(out_path + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(out_path + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(out_path + file_name + "_major-transcripts.txt.gz");
			
			// write path of network to data_map
			if (!data_map.containsKey(cell_type)) {
				data_map.put(cell_type, new LinkedList<>());
			}
			
			data_map.get(cell_type).add(out_path + file_name + "_ppin.txt.gz");
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
			
			System.out.println();
			System.out.println(cell_type + ": " + data_map.get(cell_type).size() + " samples");
			System.out.println("Size: " + (int) Utilities.getMean(no_proteins) + "+-" + (int) Utilities.getStd(no_proteins) + " / " + (int) Utilities.getMean(no_interactions) + "+-" + (int) Utilities.getStd(no_interactions));
		}
		
		System.out.println();
	}
	
	public static void main(String[] args) {
		//loadAndStoreReferenceNetwork("mixed_data/human_merged_dec_4.txt.gz");
		//System.exit(0);
		
		preprocess();
		
		new File(network_folder_pre).mkdir();
		
		process(1.0);
	}
}
