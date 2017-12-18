package PPIComp_hemato;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks2 {
	static String BLUEPRINT_expr_folder = "/Users/tho/GDrive/Work/projects/hemato_rewiring/BLUEPRINT_expr/";
	static String network_folder_pre = "/Users/tho/Desktop/BLUEPRINT_networks/";
	static Map<String, String> folder_type_map = new HashMap<>();
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void preprocess() {
		System.out.println("Original PPIN: " + "mixed_data/human_mentha_17_jan.txt.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/human_mentha_17_jan.txt.gz");
		builder = new NetworkBuilder(original_ppin);
		
		System.out.println("Proteins mapped: " + builder.getMappingDomainPercentage());
		System.out.println("PPIs mapped: " + builder.getMappingPercentage());
	}
	
	public static void process(double TPM_threshold) {
		
		System.out.println("BLUEPRINT net-builder, TPM threshold: " + TPM_threshold);
		String network_folder = network_folder_pre + TPM_threshold + "/";
		
		new File(network_folder).mkdir();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(BLUEPRINT_expr_folder, ".rsem.tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			
			String cell_type = path_split[path_split.length-2];
			
			String out_path = network_folder + cell_type + "/";
			
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, TPM_threshold);
			
			ConstructedNetworks cn = builder.constructAssociatedWeightedNetworksFromTranscriptAbundance(transcr_expr, true);
			cn.getPPIN().writePPIN(out_path + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(out_path + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(out_path + file_name + "_major-transcripts.txt.gz");
		}
	}
	
	public static void main(String[] args) {
		
		// study data was build using v83
		DataQuery.enforceSpecificEnsemblRelease("83");
		
		preprocess();
		
		new File(network_folder_pre).mkdir();
		process(0.0);
		
	}
}
