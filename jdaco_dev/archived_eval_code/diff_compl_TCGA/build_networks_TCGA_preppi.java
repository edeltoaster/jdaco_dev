package diff_compl_TCGA;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks_TCGA_preppi {
	static String expr_folder = "transcr_expr_data/"; // intended to be run on server
	static String network_folder = "TCGA_networks/";
	static PPIN original_ppin;
	static NetworkBuilder builder;
	static double transcr_threshold = 0.0;
	static String organism_database = "";
	
	public static void preprocess() {
		DataQuery.switchServerGRCh37();
		System.out.println("Original PPIN: " + "mixed_data/human_PrePPI_17_01_17_hc.txt.gz"); // PrePPI, jan 2017
		organism_database = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println("Ensembl version: " + organism_database);
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("mixed_data/human_PrePPI_17_01_17_hc.txt.gz");
		System.out.println(original_ppin.getSizesStr());
		original_ppin = original_ppin.updateUniprotAccessions();
		System.out.println("Updating Uniprot Accs with " + DataQuery.getUniprotRelease());
		System.out.println(original_ppin.getSizesStr());
		System.out.println("upper 5% cutoff: " + original_ppin.getPercentile(5));
		builder = new NetworkBuilder(original_ppin);
		
		System.out.println("Proteins mapped: " + builder.getMappingDomainPercentage());
		System.out.println("PPIs mapped: " + builder.getMappingPercentage());
		System.out.println();
	}
	
	public static void process() {
		
		Map<String, Integer> matched_sample_count = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".txt.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String[] sample_annotations = file_name.split("_");
			String cancer_type = sample_annotations[0];
			String patient = sample_annotations[1];
			//String quant_granularity = sample_annotations[2];
			String condition = sample_annotations[3];
			
			file_name = cancer_type + "_" + patient + "_" + condition;
			System.out.println("Processing " + file_name);
			
			// normalized counts, sum of all genes etc
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.normalizeByTranscriptLength(TranscriptAbundanceReader.readTCGAIsoformRSEMFile(path, transcr_threshold, false), organism_database), true, true); //returns abundance as gene abundance (sum of expressed transcripts of gene)
			cn.getPPIN().writePPIN(network_folder + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + file_name + "_major-transcripts.txt.gz");
			System.out.println(cn.getPPIN().getSizesStr());
			
			if (condition.equals("normal"))
				matched_sample_count.put(cancer_type, matched_sample_count.getOrDefault(cancer_type, 0) + 1);
		}
		
		System.out.println();
		
		int n_total = 0;
		for (String cancer_type:matched_sample_count.keySet()) {
			System.out.println(cancer_type + " : " + matched_sample_count.get(cancer_type));
			n_total += matched_sample_count.get(cancer_type);
		}
		
		System.out.println("total : " + n_total + " matching samples, " + (2*n_total) + " total samples.");
	}
	
	public static void main(String[] args) {

		preprocess();
		
		new File(network_folder).mkdir();

		process();

	}
}
