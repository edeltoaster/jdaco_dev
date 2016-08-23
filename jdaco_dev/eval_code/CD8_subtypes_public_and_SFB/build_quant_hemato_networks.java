package CD8_subtypes_public_and_SFB;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_quant_hemato_networks {
	static String expr_folder = "/Users/tho/Desktop/BLUEPRINT_expr/";
	static String network_folder_pre = "/Users/tho/Desktop/quant_hemo_networks_0.0/";
	static Map<String, String> folder_type_map = new HashMap<>();
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void loadAndStoreReferenceNetwork(String network_out) {
		PPIN ppin = DataQuery.getMenthaNetwork("9606");
		System.out.println(ppin.getSizesStr());
		ppin = ppin.updateUniprotAccessions();
		ppin.writePPIN(network_out);
		System.out.println(ppin.getSizesStr());
	}
	
	
	public static void preprocess() {
		for (String s:Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/eval_code/PPIComp_hemato/cell_types.txt")) {
			if (s.startsWith("#"))
				continue;
			String[] spl = s.trim().split(" ");
			folder_type_map.put(spl[0], spl[1]);
		}
		
		System.out.println("Original PPIN: " + "mixed_data/human_mentha_8_jul.txt.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/human_mentha_8_jul.txt.gz");
		builder = new NetworkBuilder(original_ppin);
		
		System.out.println("Proteins mapped: " + builder.getMappingDomainPercentage());
		System.out.println("PPIs mapped: " + builder.getMappingPercentage());
		
		System.out.println("");
		
	}
	
	public static void process(double TPM_threshold) {
		
		System.out.println("hemato net-builder, TPM threshold: " + TPM_threshold);
		String network_folder = network_folder_pre;
		
		new File(network_folder).mkdir();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String sample = path_split[path_split.length-1].split("\\.")[0].split("-")[2];
			String cell_type = path_split[path_split.length-2];
			
			if (!folder_type_map.containsKey(cell_type))
				continue;
			
			cell_type = folder_type_map.get(cell_type);
			sample = cell_type + "_" + sample;
			// filter N, M, NK, CD4 to venous blood samples
			if (cell_type.equals("M") || cell_type.equals("N") || cell_type.equals("NK") || cell_type.equals("CD4"))
				if (sample.startsWith("Cord"))
					continue;
			
			String out_path = network_folder;
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, TPM_threshold);
			
			System.out.println("Processing" + path);
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(transcr_expr, true);
			cn.getPPIN().writePPIN(out_path + sample + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(out_path + sample + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(out_path + sample + "_major-transcripts.txt.gz");
		}
		
	}
	
	public static void main(String[] args) {
		//loadAndStoreReferenceNetwork("mixed_data/human_mentha_8_jul.txt.gz");
		//System.exit(0);
		preprocess();
		
		new File(network_folder_pre).mkdir();
		
		process(0.0);
		
	}
}
