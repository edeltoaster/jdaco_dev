package sc_bonn_TCGA;

import java.io.File;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class builtIndCalTechNetworks {
	public static double expr_threshold = 0.3;
	public static String out_path = "/Users/tho/Desktop/CalTech_nets_ind/";
	
	public static void main(String[] args) {
		System.out.println("Build CalTech data networks");
		PPIN ppi = new PPIN("mixed_data/human_merged.tsv.gz");
		System.out.println("human preppi");
		System.out.println(ppi.getSizesStr());
		
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("expr-threshold: " + expr_threshold);
		new File(out_path).mkdir();
		
		NetworkBuilder builder = new NetworkBuilder(ppi);
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/ENCODE/CalTech/", "gtf.gz")) {
			String sample = f.getName().split("RnaSeq")[1].split("R2x75")[0].toUpperCase();
			int replicate = Integer.parseInt(f.getName().split("Rep")[1].split("V3")[0]);
			System.out.println("Processing " + sample + " (" + replicate +")");
			ConstructedNetworks net = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readGencodeGTFTranscripts(f.getAbsolutePath(), expr_threshold));
			
			net.getPPIN().writePPIN(out_path + sample + "_" + replicate +"_ppin.tsv");
			System.out.println(net.getPPIN().getSizesStr());
			
			net.getDDIN().writeDDIN(out_path + sample + "_" + replicate + "_ddin.tsv");
			System.out.println(net.getDDIN().getSizesStr());
		}
	}

}
