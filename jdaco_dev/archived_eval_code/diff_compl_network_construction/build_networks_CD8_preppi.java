package diff_compl_network_construction;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class build_networks_CD8_preppi {
	static String expr_folder = "/Users/tho/GDrive/Work/projects/stem_cell_complexome/CD8_subtypes/expr_data/";
	static String network_folder = "/Users/tho/Desktop/CD8_networks/";
	static PPIN original_ppin;
	static NetworkBuilder builder;
	
	public static void preprocess() {
		System.out.println("Original PPIN: " + "mixed_data/human_PrePPI_17_01_17_hc.txt.gz"); // PrePPI, jan 2017
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
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
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(expr_folder, ".tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			String sample_name = file_name;
			System.out.println("Processing " + sample_name);
			
			ConstructedNetworks cn = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readKallistoFile(path, 0.0), true, true); //returns abundance as gene abundance (sum of expressed transcripts of gene)
			cn.getPPIN().writePPIN(network_folder + sample_name + "_ppin.txt.gz");
			cn.getDDIN().writeDDIN(network_folder + sample_name + "_ddin.txt.gz");
			cn.writeProteinToAssumedTranscriptMap(network_folder + sample_name + "_major-transcripts.txt.gz");
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
