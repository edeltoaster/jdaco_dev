package CD8_subtypes_public_and_SFB;

import java.io.File;
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
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String cell_type = path_split[5];
			String sample = file_name.split("-")[2];
			String out_path = network_folder;
			
			if (cell_type.contains("CD4"))
				cell_type = "CD4";
			else
				cell_type = "CLP";
			
			file_name = cell_type + "_" + sample;
			
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, TPM_threshold);
			
			System.out.println("Processing" + path);
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(transcr_expr, true);
			cn.getPPIN().writePPIN(out_path + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(out_path + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(out_path + file_name + "_major-transcripts.txt.gz");
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
