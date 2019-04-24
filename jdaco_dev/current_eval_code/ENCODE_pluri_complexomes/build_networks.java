package ENCODE_pluri_complexomes;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks { // intended to be run on a server
	static String expr_folder = "quantified_samples/";
	
	static String preppi_in = "mixed_data/human_PrePPI_17_01_17_hc.txt.gz";
	static String preppi_network_folder = "preppi_networks/";
	
	static String mentha_in = "mixed_data/human_mentha_14_04_19.txt.gz";
	static String mentha_network_folder = "mentha_networks/";
	
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void preprocess(String in_file) {
		System.out.println("Original PPIN: " + in_file);
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN(in_file);
		System.out.println(original_ppin.getSizesStr());
		original_ppin = original_ppin.updateUniprotAccessions();
		System.out.println("Updating Uniprot Accs with " + DataQuery.getUniprotRelease());
		System.out.println(original_ppin.getSizesStr());
		System.out.println("upper 5% cutoff: " + original_ppin.getPercentile(5));
		builder = new NetworkBuilder(original_ppin);
		
		System.out.println("Proteins mapped: " + builder.getMappingDomainPercentage());
		System.out.println("PPIs mapped: " + builder.getMappingPercentage());
	}
	
	public static void process(String network_folder) {
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			System.out.println("Processing " + file_name);
			
			//returns abundance as gene abundance (sum of expressed transcripts of gene)
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readKallistoFile(path, 0.0), true, true); 
			
			cn.getPPIN().writePPIN(network_folder + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + file_name + "_major-transcripts.txt.gz");
			System.out.println(cn.getPPIN().getSizesStr());
		}
		
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		preprocess(mentha_in);
		new File(mentha_network_folder).mkdir();
		process(mentha_network_folder);
		
		System.out.println();
		
		preprocess(preppi_in);
		new File(preppi_network_folder).mkdir();
		process(preppi_network_folder);
		
	}
}
