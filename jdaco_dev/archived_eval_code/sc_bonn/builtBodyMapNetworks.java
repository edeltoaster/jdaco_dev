package sc_bonn;

import java.io.File;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class builtBodyMapNetworks {
	public static double expr_threshold = 0.3;
	public static String out_path = "/Users/tho/Desktop/bodymap_nets/";
	
	public static void main(String[] args) {
		System.out.println("Build bodymap data networks");
		PPIN ppi = new PPIN("mixed_data/human_merged.tsv.gz");
		System.out.println("human preppi");
		System.out.println(ppi.getSizesStr());
		
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("expr-threshold: " + expr_threshold);
		new File(out_path).mkdir();
		
		NetworkBuilder builder = new NetworkBuilder(ppi);
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/tissue_spec/bodymap/", ".gz")) {
			String sample = f.getName().split(".fpkm")[0];
			System.out.println("Processing " + sample);
			
			Map<String, Float> abundances = TranscriptAbundanceReader.readCufflinksTranscriptsFPKM(f.getAbsolutePath(), expr_threshold);
			ConstructedNetworks net = builder.constructAssociatedNetworksFromTranscriptAbundance(abundances, true);
			
			net.getPPIN().writePPIN(out_path + sample + "_ppin.tsv");
			System.out.println(net.getPPIN().getSizesStr());
			
			net.getDDIN().writeDDIN(out_path + sample + "_ddin.tsv");
			System.out.println(net.getDDIN().getSizesStr());
			
		}
	}

}
