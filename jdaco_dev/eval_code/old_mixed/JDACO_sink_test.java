package old_mixed;

import java.io.File;
import java.util.Set;

import framework.BindingDataHandler;
import framework.ConstructedNetworks;
import framework.DACO;
import framework.DACOResultSet;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.RegulatoryNetwork;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class JDACO_sink_test {
	
	static int no_threads = 40;
	static int depth = 6;
	static int percentiles = 5;
	static int d_min = -50;
	static int d_max = 50;
	static int min_TFs = 2;
	static NetworkBuilder network_builder = new NetworkBuilder(new PPIN("mixed_data/human_ppi.tsv.gz"));
	static BindingDataHandler bdh = new BindingDataHandler("/Users/tho/Dropbox/Work/binding_sites/human_fimo_2k.txt.gz");
	static String out_folder = "PrePPI_spl_d6_med_min2/";
	
	public static void process(String sample, String path) {
		System.gc();
		System.out.println("Processing " + sample);
		
		// build networks
		ConstructedNetworks networks = network_builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readCufflinksTranscriptsFPKM(path, 1.0));
				
		// run DACO
		double pair_thr = networks.getPPIN().getPercentile(percentiles);
		double prob_thr = networks.getPPIN().getSampledComplexCutoff(percentiles, depth);
				
		System.out.println("starting DACO");
		DACO daco_exec = new DACO(networks, no_threads, depth, pair_thr, prob_thr, null);
		System.gc();
		DACOResultSet res = daco_exec.batchGrowPairs(bdh.getTFsWithBindingData());
		System.out.println(res.getResult().size() + " complexes found.");
		res.writeCSV(out_folder + "DACO_res/" + sample +".csv");
		
		// get outcome in GRN
		System.out.println("building GRN");
		RegulatoryNetwork grn = new RegulatoryNetwork(res, bdh, d_min, d_max, no_threads);
		
		System.out.println("writing out");
		for (int i = 2; i <= grn.getBiggestTFComplexWithTargets(); i++) {
			Set<String> sink_prot = grn.getSinkProteins(i);
			Utilities.writeEntries(sink_prot, out_folder + "sink_proteins/" + sample + "_" + i +".txt");
			sink_prot.retainAll(networks.getPPIN().getProteins());
			Utilities.writeEntries(sink_prot, out_folder + "fsink_proteins/" + sample + "_" + i +".txt");
			grn.writeRegulatoryNetwork(out_folder + "GRNs/" + sample + "_" + i +".txt", i);
		}
		
	}
	
	public static void main(String[] args) {
		new File(out_folder).mkdir();
		new File(out_folder + "sink_proteins").mkdir();
		new File(out_folder + "fsink_proteins").mkdir();
		new File(out_folder + "GRNs").mkdir();
		new File(out_folder + "DACO_res").mkdir();
		for (File file:Utilities.getAllMatchingFilesInSubfolders("bodymap", "isoforms.fpkm_tracking.gz")) {
			String sample = file.getParentFile().getName();
			String sample_path = file.getParentFile().getAbsolutePath();
			
			//if (!sample.equals("testes") && !sample.equals("thyroid"))
				//continue;
			
			process(sample, sample_path+"/isoforms.fpkm_tracking.gz");
			
		}

	}
}
