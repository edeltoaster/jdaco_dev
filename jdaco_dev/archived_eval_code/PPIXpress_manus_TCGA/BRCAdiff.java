package PPIXpress_manus_TCGA;

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

public class BRCAdiff {
	static PPIN original_ppi;
	static double rsem_cutoff = 1.0;
	static BinomialTest test = new BinomialTest();

	public static void writeOutNet(Map<StrPair, Double> diff, int samples, String out_path, double P_rew, double P_cutoff) {
		List<String> helper = new LinkedList<String>();
		
		Map<Integer, Double> test_map = new HashMap<Integer, Double>();
		for (StrPair pair:diff.keySet()) {
			int v = (int) Math.abs(diff.get(pair));
			
			if (v == 0)
				continue;
			
			if (!test_map.containsKey(v))
				test_map.put(v, test.binomialTest(samples, (int) Math.abs(v), P_rew, AlternativeHypothesis.GREATER_THAN));
		}
		
		helper.add("Protein1 Protein2 Type Count Probability p-val");
		int n = 0;
		for (StrPair pair:diff.keySet()) {
			double v = diff.get(pair);
			
			if (v == 0)
				continue;
			
			String sign = "-";
			if (Math.signum(v) == +1) 
				sign = "+";
			double p = test_map.get((int) Math.abs(v));
			if (p >= P_cutoff)
				continue;
			helper.add(pair.getL() + " " + pair.getR() + " " + sign + " " + (int) Math.abs(v) + " " + Math.abs(v / samples) + " " + p);
			n++;
		}
		System.out.println(n + " diff-IAs with p below " + P_cutoff);
		Utilities.writeEntries(helper, out_path);
	}
	
	public static void FDRwriteOutNet(Map<StrPair, Double> diff, int samples, String out_path, double P_rew, double FDR) {
		List<String> helper = new LinkedList<String>();
		
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
				k++;
			}
		}
		System.out.println(k + " diff-IAs with FDR below " + FDR);
		Utilities.writeEntries(helper, out_path);
	}
	
	public static void process(NetworkBuilder builder, String out_folder, String TCGA_prefix, double FDR) {
		
		new File(out_folder).mkdir();
		
		// split data
		List<String> tumor_data = new ArrayList<String>();
		List<String> normal_data = new ArrayList<String>();
		Map<StrPair, Integer> overall_added = new HashMap<StrPair, Integer>();
		Map<StrPair, Integer> overall_lost = new HashMap<StrPair, Integer>();
		
		for (File file:Utilities.getAllPrefixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/tissue_spec/cancer_nets/TCGA/", TCGA_prefix)) {
			String file_name = file.getAbsolutePath();
			if (file_name.endsWith("tumor.txt.gz"))
				tumor_data.add(file_name);
			else
				normal_data.add(file_name);
		}
		
		System.out.println(TCGA_prefix.split("_")[0]+": " + tumor_data.size() + " matched samples.");
		List<Double> normal_nodes = new LinkedList<Double>();
		List<Double> tumor_nodes = new LinkedList<Double>();
		List<Double> normal_edges = new LinkedList<Double>();
		List<Double> tumor_edges = new LinkedList<Double>();
		List<Double> P_rew = new LinkedList<Double>();
		for (int i = 0; i < tumor_data.size(); i++) {
			
			String tumor_file = tumor_data.get(i);
			String normal_file = normal_data.get(i);
			
			
			// transcript-based
			PPIN tumor_net = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readTCGAIsoformRSEM(tumor_file, rsem_cutoff), true).getPPIN();
			tumor_nodes.add((double) tumor_net.getSizes()[0]);
			tumor_edges.add((double) tumor_net.getSizes()[1]);
			
			PPIN normal_net = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readTCGAIsoformRSEM(normal_file, rsem_cutoff), true).getPPIN();
			normal_nodes.add((double) normal_net.getSizes()[0]);
			normal_edges.add((double) normal_net.getSizes()[1]);
			
			Set<StrPair> added_interactions = tumor_net.removeAll(normal_net).getInteractions();
			Set<StrPair> lost_interactions = normal_net.removeAll(tumor_net).getInteractions();
			
			P_rew.add( ((double)added_interactions.size()+ lost_interactions.size())/Math.min(tumor_net.getSizes()[1], normal_net.getSizes()[1]));
			
			// count
			for (StrPair pair:added_interactions) {
				if(!overall_added.containsKey(pair))
					overall_added.put(pair, 0);
				overall_added.put(pair, overall_added.get(pair)+1 );
			}
			for (StrPair pair:lost_interactions) {
				if(!overall_lost.containsKey(pair))
					overall_lost.put(pair, 0);
				overall_lost.put(pair, overall_lost.get(pair)+1 );
			}
			
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
		System.out.println("t-based sizes:");
		System.out.println("n: " + Utilities.getMean(normal_nodes) + " / " + Utilities.getMean(normal_edges) + " +- " + Utilities.getStd(normal_nodes) + " / " + Utilities.getStd(normal_edges));
		System.out.println("t: " + Utilities.getMean(tumor_nodes) + " / " + Utilities.getMean(tumor_edges) + " +- " + Utilities.getStd(tumor_nodes) + " / " + Utilities.getStd(tumor_edges));
		System.out.println("P_rew: " + Utilities.getMean(P_rew) + " +- " + Utilities.getStd(P_rew));
		new File(out_folder).mkdir();
		FDRwriteOutNet(diff, tumor_data.size(), out_folder + "t_net.txt", Utilities.getMean(P_rew), FDR);
		
		// write stuff for statistics
		List<String> temp_list = new LinkedList<String>();
		for (double d:P_rew)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "t_P_rew.txt");
		
		temp_list.clear();
		for (double d:normal_nodes)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "t_normal_nodes.txt");
		
		temp_list.clear();
		for (double d:tumor_nodes)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "t_tumor_nodes.txt");
		
		temp_list.clear();
		for (double d:normal_edges)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "t_normal_edges.txt");
		
		temp_list.clear();
		for (double d:tumor_edges)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "t_tumor_edges.txt");
		
		/**
		 * Now the same based on genes
		 */
		
		if (!out_folder.contains("all_")) // only do that once
			return;
		
		// reset data, split data
		overall_added = new HashMap<StrPair, Integer>();
		overall_lost = new HashMap<StrPair, Integer>();
		
		normal_nodes.clear();
		tumor_nodes.clear();
		normal_edges.clear();
		tumor_edges.clear();
		P_rew.clear();
		temp_list.clear();
		for (int i = 0; i < tumor_data.size(); i++) {
			String tumor_file = tumor_data.get(i);
			String normal_file = normal_data.get(i);
			
			
			// gene-based
			PPIN tumor_net = builder.constructAssociatedNetworksFromGeneAbundance(TranscriptAbundanceReader.getGeneAbundanceFromTCGAIsoformRSEM(tumor_file, rsem_cutoff), false).getPPIN();
			tumor_nodes.add((double) tumor_net.getSizes()[0]);
			tumor_edges.add((double) tumor_net.getSizes()[1]);
			
			PPIN normal_net = builder.constructAssociatedNetworksFromGeneAbundance(TranscriptAbundanceReader.getGeneAbundanceFromTCGAIsoformRSEM(normal_file, rsem_cutoff), false).getPPIN();
			normal_nodes.add((double) normal_net.getSizes()[0]);
			normal_edges.add((double) normal_net.getSizes()[1]);
			
			Set<StrPair> added_interactions = tumor_net.removeAll(normal_net).getInteractions();
			Set<StrPair> lost_interactions = normal_net.removeAll(tumor_net).getInteractions();
			
			P_rew.add( ((double)added_interactions.size()+ lost_interactions.size())/Math.min(tumor_net.getSizes()[1], normal_net.getSizes()[1]));
			
			// count
			for (StrPair pair:added_interactions) {
				if(!overall_added.containsKey(pair))
					overall_added.put(pair, 0);
				overall_added.put(pair, overall_added.get(pair)+1 );
			}
			for (StrPair pair:lost_interactions) {
				if(!overall_lost.containsKey(pair))
					overall_lost.put(pair, 0);
				overall_lost.put(pair, overall_lost.get(pair)+1 );
			}
			
		}
		
		// do diff.
		diff = new HashMap<StrPair, Double>();
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
		System.out.println("g-based sizes:");
		System.out.println("n: " + Utilities.getMean(normal_nodes) + " / " + Utilities.getMean(normal_edges) + " +- " + Utilities.getStd(normal_nodes) + " / " + Utilities.getStd(normal_edges));
		System.out.println("t: " + Utilities.getMean(tumor_nodes) + " / " + Utilities.getMean(tumor_edges) + " +- " + Utilities.getStd(tumor_nodes) + " / " + Utilities.getStd(tumor_edges));
		System.out.println("P_rew: " + Utilities.getMean(P_rew) + " +- " + Utilities.getStd(P_rew));
		
		FDRwriteOutNet(diff, normal_data.size(), out_folder + "g_net.txt", Utilities.getMean(P_rew), FDR);
		
		// write stuff for statistics
		for (double d:P_rew)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "g_P_rew.txt");
			
		temp_list.clear();
		for (double d:normal_nodes)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "g_normal_nodes.txt");
		
		temp_list.clear();
		for (double d:tumor_nodes)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "g_tumor_nodes.txt");
		
		temp_list.clear();
		for (double d:normal_edges)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "g_normal_edges.txt");
		
		temp_list.clear();
		for (double d:tumor_edges)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "g_tumor_edges.txt");
		
		System.out.println("");
	}
	
	public static void main(String[] args) {
		DataQuery.switchServer("ensembldb.ensembl.org:3306");
		DataQuery.enforceSpecificEnsemblRelease("79");
		DataQuery.enforceSpecific3didRelease("2015_02");
		System.out.println("TCGA diff-builder, FDR: 0.05");
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did:" + DataQuery.get3didVersion());
		System.out.println("iPfam:" + DataQuery.getIPfamVersion());
		
		
		PPIN intact = new PPIN("mixed_data/human_intact.txt.gz");
		PPIN biogrid = new PPIN("mixed_data/human_biogrid.txt.gz");
		
		System.out.println("");
		
		System.out.println("ALL-DDI");
		
		original_ppi = intact;
		NetworkBuilder builder = new NetworkBuilder(original_ppi);
		
		System.out.println("IntAct");
		process(builder, "/Users/tho/Desktop/all_intact/", "BRCA", 0.05);
		
		original_ppi = biogrid;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("BioGRID");
		process(builder, "/Users/tho/Desktop/all_biogrid/", "BRCA", 0.05);
		
		DataQuery.localDDIsOnly();
		//DataQuery.stricterLocalDDIs();
		
		System.out.println("PRE_HC");
		
		original_ppi = intact;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("IntAct");
		process(builder, "/Users/tho/Desktop/local_intact/", "BRCA", 0.05);
		
		original_ppi = biogrid;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("BioGRID");
		process(builder, "/Users/tho/Desktop/local_biogrid/", "BRCA", 0.05);
		
		
		DataQuery.onlyRetrievedDDIs();
		
		System.out.println("RETRIEVED");
		
		original_ppi = intact;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("IntAct");
		process(builder, "/Users/tho/Desktop/retrieved_intact/", "BRCA", 0.05);
		
		original_ppi = biogrid;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("BioGRID");
		process(builder, "/Users/tho/Desktop/retrieved_biogrid/", "BRCA", 0.05);
		
		// not working like this, execute extra
		DataQuery.localDDIsOnly();
		DataQuery.stricterLocalDDIs();
		
		System.out.println("PRE_VHC");
		
		original_ppi = intact;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("IntAct");
		process(builder, "/Users/tho/Desktop/less_intact/", "BRCA", 0.05);
		
		original_ppi = biogrid;
		builder = new NetworkBuilder(original_ppi);
		
		System.out.println("BioGRID");
		process(builder, "/Users/tho/Desktop/less_biogrid/", "BRCA", 0.05);
		
	}

}
