package diff_compl_network_construction;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks_BRCA_preppi {
	static String expr_folder = "expr_data/"; // intended to be run on server
	static String network_folder = "BRCA_networks_0/";
	static PPIN original_ppin;
	static NetworkBuilder builder;
	static double transcr_threshold = 0.0;
	static String organism_database;
	
	public static void preprocess() {
		DataQuery.switchServerGRCh37();
		organism_database = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		System.out.println("Original PPIN: " + "mixed_data/human_PrePPI_17_01_17_hc.txt.gz"); // PrePPI, jan 2017
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
	}
	
	public static void process() {
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".txt.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String[] sample_annotations = file_name.split("_");
			String cancer_type = sample_annotations[0];
			String patient = sample_annotations[1];
			String quant_granularity = sample_annotations[2];
			String condition = sample_annotations[3];
			
			if (!cancer_type.equals("BRCA"))
				continue;
			
			if (!quant_granularity.equals("transcripts"))
				continue;
			
			file_name = patient + "_" + condition;
			System.out.println("Processing " + file_name);
			
			// normalized counts, sum of all transcripts etc
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.normalizeByTranscriptLength(TranscriptAbundanceReader.readTCGAIsoformRSEMFile(path, transcr_threshold, false), organism_database), true, true); //returns abundance as gene abundance (sum of expressed transcripts of gene)
			cn.getPPIN().writePPIN(network_folder + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + file_name + "_major-transcripts.txt.gz");
			System.out.println(cn.getPPIN().getSizesStr());
		}
		
		System.out.println();
	}
	
	public static void main(String[] args) {

		DataQuery.enforceSpecificEnsemblRelease("89");
		preprocess();
		
		new File(network_folder).mkdir();

		process();

	}
}
