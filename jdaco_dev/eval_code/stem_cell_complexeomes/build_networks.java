package stem_cell_complexeomes;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks {
	static String expr_folder = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/expr_data/";
	static String network_folder = "/Users/tho/Desktop/networks/";
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
		System.out.println("Original PPIN: " + "mixed_data/human_mentha_25_oct.txt.gz");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		
		original_ppin = new PPIN("/Users/tho/git/jdaco_dev/jdaco_dev/mixed_data/human_mentha_25_oct.txt.gz");
		builder = new NetworkBuilder(original_ppin);
		
		System.out.println("Proteins mapped: " + builder.getMappingDomainPercentage());
		System.out.println("PPIs mapped: " + builder.getMappingPercentage());
	}
	
	public static void process() {
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			System.out.println("Processing " + file_name);
			
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readKallistoFile(path, 0.0), true);
			cn.getPPIN().writePPIN(network_folder + file_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + file_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + file_name + "_major-transcripts.txt.gz");
			System.out.println(cn.getPPIN().getSizesStr());
		}
		
		System.out.println();
	}
	
	public static void main(String[] args) {
		//loadAndStoreReferenceNetwork("mixed_data/human_mentha_25_oct.txt.gz");
		//System.exit(0);
		
		preprocess();
		
		new File(network_folder).mkdir();

		process();

	}
}
