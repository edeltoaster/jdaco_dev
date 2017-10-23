package diff_compl_mono_geu;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks_preppi {
	static String expr_folder = "expr_data/"; // intended to be run on server
	static String network_folder = "networks/";
	static PPIN original_ppin;
	static NetworkBuilder builder;
	static double transcr_threshold = 0.0;
	static String organism_database = "";
	
	public static void preprocess() {
		DataQuery.enforceSpecificEnsemblRelease("90");
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
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			System.out.println("Processing " + file_name);
			
			// normalized counts, sum of all genes etc
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readKallistoFile(path, transcr_threshold), true, true);
			cn.getPPIN().writePPIN(network_folder + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + file_name + "_major-transcripts.txt.gz");
			System.out.println(cn.getPPIN().getSizesStr());
		}
	}
	
	public static void main(String[] args) {

		preprocess();
		
		new File(network_folder).mkdir();

		process();

	}
}
