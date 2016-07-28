package eval_TCGA;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.StrPair;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

// partially changed blablubb

public class TCGAdiff {
	static PPIN original_ppi;
	static double rsem_cutoff = 1.0;
	static BinomialTest test = new BinomialTest();
	static int min_samples = 15;
	
	public static HashSet<String> FDRwriteOutNet(Map<StrPair, Double> diff, int samples, String out_path, double P_rew, double FDR) {
		List<String> helper = new LinkedList<String>();
		HashSet<String> sign_net = new HashSet<String>();
		Map<StrPair, Double> test_map = new HashMap<StrPair, Double>();
		Map<Double, LinkedList<StrPair>> p2pair = new HashMap<Double, LinkedList<StrPair>>();
		for (StrPair pair:diff.keySet()) {
			int v = (int) Math.abs(diff.get(pair));
			
			if (v == 0)
				continue;
			double raw_p = test.binomialTest(samples, (int) Math.abs(v), P_rew, AlternativeHypothesis.GREATER_THAN);
			test_map.put(pair, raw_p);
			
			if (!p2pair.containsKey(raw_p))
				p2pair.put(raw_p, new LinkedList<StrPair>());
			p2pair.get(raw_p).add(pair);
		}
		
		List<Double> p_values = new ArrayList<Double>(test_map.values());
		int m = p_values.size();
		Collections.sort(p_values);
		int k = 1;
		int largest_k = -1;
		Map<Double, Double> rawp2adjp = new HashMap<Double, Double>();
		
		// find largest k
		for (double p:p_values) {
			if (p <= k * FDR / m)
				largest_k = k;
			rawp2adjp.put(p, (p* m) / k); // if multiple have the same rank, take the biggest
			k++;
		}
		
		helper.add("Protein1 Protein2 Type Count Probability p-val p-val_adj");
		k = 1;
		p_values = new LinkedList<Double>(new HashSet<Double>(p_values));
		Collections.sort(p_values);
		for (double p:p_values) {
			if (k > largest_k) { // remaining ones not deemed significant
				// reset counter for output
				k--;
				break;
			}
			for (StrPair pair:p2pair.get(p)) {
				double v = diff.get(pair);
				if (v == 0)
					continue;
				String sign = "-";
				if (Math.signum(v) == +1) 
					sign = "+";
				helper.add(pair.getL() + " " + pair.getR() + " " + sign + " " + (int) Math.abs(v) + " " + Math.abs(v / samples) + " " + p + " " + rawp2adjp.get(p));
				sign_net.add(pair.getL() + " " + pair.getR() + " " + sign);
				k++;
			}
		}
		System.out.println(k + " diff-IAs with FDR below " + FDR);
		Utilities.writeEntries(helper, out_path);
		
		return sign_net;
	}
	
	public static void process(NetworkBuilder builder, String out_folder, String TCGA_prefix, double FDR) {
		
		// split data
		List<String> tumor_data = new ArrayList<String>();
		List<String> normal_data = new ArrayList<String>();
		Map<StrPair, Integer> overall_added = new HashMap<StrPair, Integer>();
		Map<StrPair, Integer> overall_lost = new HashMap<StrPair, Integer>();
		
		for (File file:Utilities.getAllPrefixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/cancer/TCGA/", TCGA_prefix)) {
			String file_name = file.getAbsolutePath();
			if (file_name.endsWith("_transcripts_tumor.txt.gz"))
				tumor_data.add(file_name);
			else if (file_name.endsWith("_transcripts_normal.txt.gz")) {
				normal_data.add(file_name);
			}
		}
		
		System.out.println(TCGA_prefix.split("_")[0]+": " + tumor_data.size() + " matched samples.");
		
		if (tumor_data.size() < min_samples) {
			System.out.println("Too few samples.");
			return;
		}
		
		List<Double> normal_nodes = new LinkedList<Double>();
		List<Double> tumor_nodes = new LinkedList<Double>();
		List<Double> normal_edges = new LinkedList<Double>();
		List<Double> tumor_edges = new LinkedList<Double>();
		Map<String, Double> P_rew = new HashMap<String, Double>();
		
		Map<String, HashSet<String>> diffnets = new HashMap<String, HashSet<String>>();
		
		new File(out_folder).mkdir();
		
		// over all patients
		for (int i = 0; i < tumor_data.size(); i++) {
			
			String tumor_file = tumor_data.get(i);
			String normal_file = normal_data.get(i);
			String[] path = tumor_file.split("/");
			String patient_id = path[path.length-1].split("_")[1];
			
			// transcript-based
			PPIN tumor_net = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readTCGAIsoformRSEM(tumor_file, rsem_cutoff), true).getPPIN();
			tumor_nodes.add((double) tumor_net.getSizes()[0]);
			tumor_edges.add((double) tumor_net.getSizes()[1]);
			
			PPIN normal_net = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readTCGAIsoformRSEM(normal_file, rsem_cutoff), true).getPPIN();
			normal_nodes.add((double) normal_net.getSizes()[0]);
			normal_edges.add((double) normal_net.getSizes()[1]);
			
			Set<StrPair> added_interactions = tumor_net.removeAll(normal_net).getInteractions();
			Set<StrPair> lost_interactions = normal_net.removeAll(tumor_net).getInteractions();
			
			P_rew.put(patient_id, ((double)added_interactions.size()+ lost_interactions.size())/Math.min(tumor_net.getSizes()[1], normal_net.getSizes()[1]));
			
			// count
			HashSet<String> temp_net = new HashSet<String>();
			
			for (StrPair pair:added_interactions) {
				temp_net.add(pair.getL() + " " + pair.getR() + " +");
				if(!overall_added.containsKey(pair))
					overall_added.put(pair, 0);
				overall_added.put(pair, overall_added.get(pair)+1 );
			}
			for (StrPair pair:lost_interactions) {
				temp_net.add(pair.getL() + " " + pair.getR() + " -");
				if(!overall_lost.containsKey(pair))
					overall_lost.put(pair, 0);
				overall_lost.put(pair, overall_lost.get(pair)+1 );
			}
			
			// output diff-network
			Utilities.writeEntries(temp_net, out_folder + patient_id+"_diffnet.txt");
			normal_net.writePPIN(out_folder + patient_id+"_normalppi.txt");
			tumor_net.writePPIN(out_folder + patient_id+"_tumorppi.txt");
			diffnets.put(patient_id, temp_net);
		}
		
		// do diff.
		Map<StrPair, Double> diff = new HashMap<StrPair, Double>();
		for (StrPair pair:overall_added.keySet()) {
			if (!diff.containsKey(pair))
				diff.put(pair, 0.0);
			diff.put(pair, diff.get(pair) + overall_added.get(pair));
		}
		for (StrPair pair:overall_lost.keySet()) {
			if (!diff.containsKey(pair))
				diff.put(pair, 0.0);
			diff.put(pair, diff.get(pair) - overall_lost.get(pair));
		}
		
		// sizes
		System.out.println("sizes:");
		System.out.println("n: " + Utilities.getMean(normal_nodes) + " / " + Utilities.getMean(normal_edges) + " +- " + Utilities.getStd(normal_nodes) + " / " + Utilities.getStd(normal_edges));
		System.out.println("t: " + Utilities.getMean(tumor_nodes) + " / " + Utilities.getMean(tumor_edges) + " +- " + Utilities.getStd(tumor_nodes) + " / " + Utilities.getStd(tumor_edges));
		System.out.println("P_rew: " + Utilities.getMean(P_rew.values()) + " +- " + Utilities.getStd(P_rew.values()));
		
		HashSet<String> sign_net = FDRwriteOutNet(diff, tumor_data.size(), out_folder + "sign_diffnet.txt", Utilities.getMean(P_rew.values()), FDR);
		
		// write P_rews
		List<String> temp_list = new LinkedList<String>();
		for (String patient:P_rew.keySet())
			temp_list.add(patient + " " + Double.toString(P_rew.get(patient)));
		Utilities.writeEntries(temp_list, out_folder + "P_rew.txt");
		System.out.println();
		
		// write sign. patient subnetworks
		for (String patient_id:diffnets.keySet()) {
			// prune nets to sign. parts
			diffnets.get(patient_id).retainAll(sign_net);
			
			// output
			Utilities.writeEntries(diffnets.get(patient_id), out_folder + patient_id+"_signdiffnet.txt");
		}
	}
	
	public static void main(String[] args) {
		System.out.println("TCGA diff-builder, human merged, FDR: 0.05");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("iPfam version: " + DataQuery.getIPfamVersion());
		System.out.println("3did version: " + DataQuery.get3didVersion());
		System.out.println();
		
		original_ppi = new PPIN("mixed_data/human_merged.tsv.gz");
		NetworkBuilder builder = new NetworkBuilder(original_ppi);
		
		// get all prefixes
		Set<String> cancers = new HashSet<String>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/cancer/TCGA/", ".txt.gz")) {
			cancers.add(f.getName().split("_")[0]);
		}
		
		new File("/Users/tho/Desktop/merged_fdr5/").mkdir();
		for (String cancer:cancers) {
			process(builder, "/Users/tho/Desktop/merged_fdr5/" + cancer + "/", cancer, 0.05);
		}

	}

}
