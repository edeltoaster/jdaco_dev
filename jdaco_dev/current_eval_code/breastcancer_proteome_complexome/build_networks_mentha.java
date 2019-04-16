package breastcancer_proteome_complexome;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks_mentha {
	static String expr_folder = "/Users/tho/GDrive/Work/projects/breastcancer_proteome_complexomes/abundances/";
	static String network_folder = "/Users/tho/Desktop/mentha_networks/";
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void preprocess() {
		System.out.println("Original PPIN: " + "mixed_data/human_mentha_14_04_19.txt.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("mixed_data/human_mentha_14_04_19.txt.gz");
		System.out.println(original_ppin.getSizesStr());
		original_ppin = original_ppin.updateUniprotAccessions();
		System.out.println("Updating Uniprot Accs with " + DataQuery.getUniprotRelease());
		System.out.println(original_ppin.getSizesStr());
		System.out.println("upper 5% cutoff: " + original_ppin.getPercentile(5));
		System.out.println("upper 10% cutoff: " + original_ppin.getPercentile(10));
		builder = new NetworkBuilder(original_ppin);
		
		System.out.println("Proteins mapped: " + builder.getMappingDomainPercentage());
		System.out.println("PPIs mapped: " + builder.getMappingPercentage());
	}
	
	public static void process() {
		String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".txt.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0] + "." + path_split[path_split.length-1].split("\\.")[1];
			System.out.println("Processing " + file_name);
			
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromGeneAbundance(TranscriptAbundanceReader.readSample(path, 0.0, true, db), true);
			cn.getPPIN().writePPIN(network_folder + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + file_name + "_major-transcripts.txt.gz");
			System.out.println(cn.getPPIN().getSizesStr());
		}
		
		System.out.println();
	}
	
	public static void main(String[] args) {

		preprocess();
		
		new File(network_folder).mkdir();

		process();

	}
}
