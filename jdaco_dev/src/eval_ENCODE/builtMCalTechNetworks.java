package eval_ENCODE;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class builtMCalTechNetworks {
	public static double expr_threshold = 1;
	public static String out_path = "/Users/tho/Desktop/CalTech_nets_merged/";
	
	public static void main(String[] args) {
		System.out.println("Build CalTech data networks");
		PPIN ppi = new PPIN("mixed_data/human_ppi.tsv.gz");
		System.out.println("human preppi");
		System.out.println(ppi.getSizesStr());
		
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("expr-threshold: " + expr_threshold);
		new File(out_path).mkdir();
		
		NetworkBuilder builder = new NetworkBuilder(ppi);
		Map<String, LinkedList<Map<String, Float>>> replications = new HashMap<String, LinkedList<Map<String,Float>>>();
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/ENCODE/CalTech/", "gtf.gz")) {
			String sample = f.getName().split("RnaSeq")[1].split("R2x75")[0].toUpperCase();
			if (!replications.containsKey(sample))
				replications.put(sample, new LinkedList<Map<String,Float>>());
			replications.get(sample).add(TranscriptAbundanceReader.readGencodeGTFTranscripts(f.getAbsolutePath(), 0.0));
		}
		
		for (String sample:replications.keySet()) {
			System.out.println("Processing: " + sample);
			Map<String, Float> abundances = TranscriptAbundanceReader.averageAbundances(replications.get(sample), expr_threshold);
			ConstructedNetworks net = builder.constructAssociatedNetworksFromTranscriptAbundance(abundances);
			
			net.getPPIN().writePPIN(out_path + sample + "_ppin.tsv");
			System.out.println(net.getPPIN().getSizesStr());
			
			net.getDDIN().writeDDIN(out_path + sample + "_ddin.tsv");
			System.out.println(net.getDDIN().getSizesStr());
		}
	}

}
