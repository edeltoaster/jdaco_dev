package eval_ENCODE;

import java.io.File;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class builtCSHLNetworks {
	//public static double expr_threshold = 0.01;
	public static String out_path = "/Users/tho/Desktop/CSHL_nets/";
	
	public static void main(String[] args) {
		System.out.println("Build CSHL data networks");
		PPIN ppi = new PPIN("mixed_data/human_merged.tsv.gz");
		System.out.println("human preppi");
		System.out.println(ppi.getSizesStr());
		
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		//System.out.println("expr-threshold: " + expr_threshold);
		System.out.println("expr-threshold: as in publication");
		new File(out_path).mkdir();
		
		NetworkBuilder builder = new NetworkBuilder(ppi);
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/ENCODE/CSHL/", "gtf.gz")) {
			String sample = f.getName().split("RnaSeq")[1].split("Cell")[0].toUpperCase();
			System.out.println("Processing " + sample);
			
			//Map<String, Float> abundances = TranscriptAbundanceReader.readGencodeGTFTranscripts(f.getAbsolutePath(), expr_threshold);
			Map<String, Float> abundances = TranscriptAbundanceReader.readCSHLData(f.getAbsolutePath());
			
			if (abundances.get("no_samples") < 2) {
				System.out.println("only 1 sample -> withdraw");
				continue;
			}
			abundances.remove("no_samples");
			
			ConstructedNetworks net = builder.constructAssociatedNetworksFromTranscriptAbundance(abundances);
			
			net.getPPIN().writePPIN(out_path + sample + "_ppin.tsv");
			System.out.println(net.getPPIN().getSizesStr());
			
			net.getDDIN().writeDDIN(out_path + sample + "_ddin.tsv");
			System.out.println(net.getDDIN().getSizesStr());
			
		}
	}

}
