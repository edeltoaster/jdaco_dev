package hemato_PPI;

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

import framework.PPIN;
import framework.StrPair;
import framework.Utilities;

public class build_diff_networks {
	
	static BinomialTest test = new BinomialTest();
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Desktop/BLUEPRINT_networks_0.03125/";
	static String results_root = "/Users/tho/Desktop/test/";
	
	
	public static HashSet<String> FDRwriteOutNet(Map<StrPair, Double> diff, int samples, String out_file, double P_rew, double FDR) {
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
		Utilities.writeEntries(helper, out_file);
		
		return sign_net;
	}
	
	public static void processPPINs(List<PPIN> group1, List<PPIN> group2, double FDR, String out_folder) {
		Map<StrPair, Integer> overall_added = new HashMap<StrPair, Integer>();
		Map<StrPair, Integer> overall_lost = new HashMap<StrPair, Integer>();
		
		System.out.println("Processing " + out_folder + ":");
		System.out.println("Comparing " + group1.size() + " with " + group2.size() + " samples.");

		List<Double> g1_nodes = new LinkedList<Double>();
		List<Double> g2_nodes = new LinkedList<Double>();
		List<Double> g1_edges = new LinkedList<Double>();
		List<Double> g2_edges = new LinkedList<Double>();
		List<Double> P_rew = new LinkedList<Double>();
		
		// sizes / edges
		for (PPIN ppin1:group1) {
			g1_nodes.add( (double) ppin1.getSizes()[0]);
			g1_edges.add((double) ppin1.getSizes()[1]);
		}
		
		for (PPIN ppin1:group2) {
			g2_nodes.add( (double) ppin1.getSizes()[0]);
			g2_edges.add((double) ppin1.getSizes()[1]);
		}
		
		// comparison
		for (PPIN ppin1:group1) {
			for (PPIN ppin2:group2) {
				
				Set<StrPair> added_interactions = ppin2.removeAllIAs(ppin1).getInteractions();
				Set<StrPair> lost_interactions = ppin1.removeAllIAs(ppin2).getInteractions();
				P_rew.add( ((double)added_interactions.size()+ lost_interactions.size())/Math.min(ppin1.getSizes()[1], ppin2.getSizes()[1]));
				
				
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
		
		System.out.println("g1: " + (int) Utilities.getMean(g1_nodes) + "+-" + (int) Utilities.getStd(g1_nodes) + " / " + (int) Utilities.getMean(g1_edges) + "+-" + (int) Utilities.getStd(g1_edges));
		System.out.println("g2: " + (int) Utilities.getMean(g2_nodes) + "+-" + (int) Utilities.getStd(g2_nodes) + " / " + (int) Utilities.getMean(g2_edges) + "+-" + (int) Utilities.getStd(g2_edges));
		System.out.println("P_rew: " + Utilities.getMean(P_rew) + " +- " + Utilities.getStd(P_rew));
		
		new File(out_folder).mkdir();
		FDRwriteOutNet(diff, P_rew.size(), out_folder + "diffnet.tsv", Utilities.getMean(P_rew), FDR);
		
		// write P_rews for statistics
		List<String> temp_list = new LinkedList<String>();
		for (double d:P_rew)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "P_rew_distr.txt");
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		// read all data
		
		List<PPIN> HSCs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "HSC/", ".tsv.gz"))
			HSCs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> MPPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "MPP/", ".tsv.gz"))
			MPPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> CLPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "CLP/", ".tsv.gz"))
			CLPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> CMPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "CMP/", ".tsv.gz"))
			CMPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> GMPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "GMP/", ".tsv.gz"))
			GMPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> MEPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "MEP/", ".tsv.gz"))
			MEPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> EBs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "EB/", ".tsv.gz"))
			EBs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> MKs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "MK/", ".tsv.gz"))
			MKs.add(new PPIN(f.getAbsolutePath()));
		
		// process data
		
		processPPINs(HSCs, MPPs, FDR, results_root + "HSC_MPP/");
		processPPINs(MPPs, CMPs, FDR, results_root + "MPP_CMP/");
		processPPINs(MPPs, CLPs, FDR, results_root + "MPP_CLP/");
		processPPINs(CMPs, MEPs, FDR, results_root + "CMP_MEP/");
		processPPINs(CMPs, GMPs, FDR, results_root + "CMP_GMP/");
		processPPINs(MEPs, MKs, FDR, results_root + "MEP_MK/");
		processPPINs(MEPs, EBs, FDR, results_root + "MEP_EB/");
	}
}
