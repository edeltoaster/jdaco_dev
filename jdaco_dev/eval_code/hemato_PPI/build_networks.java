package hemato_PPI;

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
	static String BLUEPRINT_expr_folder = "/Users/tho/Desktop/BLUEPRINT_expr/";
	static String network_folder_pre = "/Users/tho/Desktop/BLUEPRINT_networks/";
	static Map<String, String> folder_type_map = new HashMap<>();
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void preprocess() {
		for (String s:Utilities.readEntryFile("eval_code/hemato_PPI/cell_types.txt")) {
			if (s.startsWith("#"))
				continue;
			String[] spl = s.trim().split(" ");
			folder_type_map.put(spl[0], spl[1]);
		}
		
		System.out.println("Original PPIN: " + "mixed_data/human_merged_6_Oct_15.tsv.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did:" + DataQuery.get3didVersion());
		System.out.println("iPfam:" + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("mixed_data/human_merged_6_Oct_15.tsv.gz");
		builder = new NetworkBuilder(original_ppin);
	}
	
	public static void process(double TPM_threshold) {
		
		System.out.println("BLUEPRINT net-builder, TPM threshold: " + TPM_threshold);
		String network_folder = network_folder_pre + TPM_threshold + "/";
		
		new File(network_folder).mkdir();
		
		Map<String, List<String>> data_map = new HashMap<>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(BLUEPRINT_expr_folder, ".rsem.tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String cell_type = path_split[path_split.length-2];
			
			if (!folder_type_map.containsKey(cell_type))
				continue;
			
			cell_type = folder_type_map.get(cell_type);
			
			String out_path = network_folder + cell_type + "/";
			
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, TPM_threshold));
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
		preprocess();
		
		process(0.0);
		process(0.03125);
		process(0.1);
		process(0.3);
		process(0.5);
		process(1.0);
	}
}
