package CD8_subtypes_public;

import java.io.File;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks {
	static String expr_folder = "/Users/tho/Dropbox/Work/projects/CD8_subsets_public/expr_data/";
	static String network_folder_pre = "/Users/tho/Desktop/CD8_networks/";
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
		
		System.out.println("CD8 net-builder, TPM threshold: " + TPM_threshold);
		String network_folder = network_folder_pre + TPM_threshold + "/";
		
		new File(network_folder).mkdir();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String cell_type = file_name.split("_")[1];
			
			String out_path = network_folder + cell_type + "/";
			
			if (!new File(out_path).exists())
				new File(out_path).mkdir();
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readKallistoFile(path, TPM_threshold);
			
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
		
		for (double thr = 0.0; thr <= 1.0; thr += 0.01) {
			thr = Math.round(thr * 100.0) / 100.0;
			process(thr);
		}
		
	}
}
