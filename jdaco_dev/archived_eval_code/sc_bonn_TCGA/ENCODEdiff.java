package sc_bonn_TCGA;

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
import framework.PPIN;
import framework.StrPair;
import framework.Utilities;

public class ENCODEdiff {

	static BinomialTest test = new BinomialTest();
	static double FDR = 0.05;
	
	public static void FDRwriteOutNet(Map<StrPair, Double> diff, int samples, String out_folder, double P_rew, double FDR) {
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
		Set<String> proteins_add = new HashSet<String>();
		Set<String> proteins_del = new HashSet<String>();
		int max = -1;
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
				
				if ((int) Math.abs(v) > max)
					max = (int) Math.abs(v);
				
				helper.add(pair.getL() + " " + pair.getR() + " " + sign + " " + (int) Math.abs(v) + " " + Math.abs(v / samples) + " " + p + " " + rawp2adjp.get(p));
				
				if (sign.equals("-")) {
					proteins_del.add(pair.getL());
					proteins_del.add(pair.getR());
				} else {
					proteins_add.add(pair.getL());
					proteins_add.add(pair.getR());
				}
					
				k++;
			}
		}
		System.out.println(k + " diff-IAs with FDR below " + FDR);
		
		// all proteins
		Set<String> all_proteins = new HashSet<String>(proteins_add);
		all_proteins.addAll(proteins_del);
		System.out.println("#all proteins: " + all_proteins.size());
		
		// all "mixed" proteins
		Set<String> mixed_proteins = new HashSet<String>(proteins_add);
		mixed_proteins.retainAll(proteins_del);
		System.out.println("#mixed proteins: " + mixed_proteins.size());
		
		// write node attr. table
		List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
		Map<String, String> up_hgnc = new HashMap<String, String>();
		for (String[] s:HGNC_map) {
			up_hgnc.put(s[1], s[0]);
		}
		
		List<String> to_write = new LinkedList<String>();
		to_write.add("Protein HGNC_name mixed");
		for (String protein:all_proteins) {
			to_write.add(protein + " " + up_hgnc.get(protein) + " " + mixed_proteins.contains(protein));
		}
		Utilities.writeEntries(to_write, out_folder + "diff_net_node_attr.txt");
		Utilities.writeEntries(all_proteins, out_folder + "diff_net_all_proteins.txt");
		Utilities.writeEntries(mixed_proteins, out_folder + "diff_net_mixed_proteins.txt");
		Utilities.writeEntries(helper, out_folder + "diff_net.txt");
		
		System.out.println("only highest");
		helper.clear();
		helper.add("Protein1 Protein2 Type Count Probability p-val p-val_adj");
		k = 1;
		p_values = new LinkedList<Double>(new HashSet<Double>(p_values));
		Collections.sort(p_values);
		proteins_add = new HashSet<String>();
		proteins_del = new HashSet<String>();
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
				
				if ((int) Math.abs(v) != max)
					continue;
				
				helper.add(pair.getL() + " " + pair.getR() + " " + sign + " " + (int) Math.abs(v) + " " + Math.abs(v / samples) + " " + p + " " + rawp2adjp.get(p));
				
				if (sign.equals("-")) {
					proteins_del.add(pair.getL());
					proteins_del.add(pair.getR());
				} else {
					proteins_add.add(pair.getL());
					proteins_add.add(pair.getR());
				}
					
				k++;
			}
		}
		System.out.println(k + " diff-IAs with highest count");
		
		// all proteins
		all_proteins = new HashSet<String>(proteins_add);
		all_proteins.addAll(proteins_del);
		System.out.println("#all proteins: " + all_proteins.size());
		
		// all "mixed" proteins
		mixed_proteins = new HashSet<String>(proteins_add);
		mixed_proteins.retainAll(proteins_del);
		System.out.println("#mixed proteins: " + mixed_proteins.size());
		
		// write node attr. table
		to_write.clear();
		to_write.add("Protein HGNC_name mixed");
		for (String protein:all_proteins) {
			to_write.add(protein + " " + up_hgnc.get(protein) + " " + mixed_proteins.contains(protein));
		}
		Utilities.writeEntries(to_write, out_folder + "h_diff_net_node_attr.txt");
		Utilities.writeEntries(all_proteins, out_folder + "h_diff_net_all_proteins.txt");
		Utilities.writeEntries(mixed_proteins, out_folder + "h_diff_net_mixed_proteins.txt");
		Utilities.writeEntries(helper, out_folder + "h_diff_net.txt");
	}
	
	public static void processPPINs(List<PPIN> group1, List<PPIN> group2, String out_folder, double FDR) {
		Map<StrPair, Integer> overall_added = new HashMap<StrPair, Integer>();
		Map<StrPair, Integer> overall_lost = new HashMap<StrPair, Integer>();
		
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
				Set<StrPair> added_interactions = ppin2.removeAll(ppin1).getInteractions();
				Set<StrPair> lost_interactions = ppin1.removeAll(ppin2).getInteractions();
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
		
		System.out.println("g1: " + Utilities.getMean(g1_nodes) + " / " + Utilities.getMean(g1_edges) + " +- " + Utilities.getStd(g1_nodes) + " / " + Utilities.getStd(g1_edges));
		System.out.println("g2: " + Utilities.getMean(g2_nodes) + " / " + Utilities.getMean(g2_edges) + " +- " + Utilities.getStd(g2_nodes) + " / " + Utilities.getStd(g2_edges));
		System.out.println("P_rew: " + Utilities.getMean(P_rew) + " +- " + Utilities.getStd(P_rew));
		new File(out_folder).mkdir();
		FDRwriteOutNet(diff, P_rew.size(), out_folder, Utilities.getMean(P_rew), FDR);
		
		// write P_rews for statistics
		List<String> temp_list = new LinkedList<String>();
		for (double d:P_rew)
			temp_list.add(Double.toString(d));
		Utilities.writeEntries(temp_list, out_folder + "P_rew_distr.txt");
	}
	
	public static void main(String[] args) {
		System.out.println("CT ESCs vs BM tissues");
		List<PPIN> ESCs = new LinkedList<PPIN>();
		List<PPIN> tissues = new LinkedList<PPIN>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Documents/workspace/ENCODE_check/data/merged_networks/CalTech_nets_ind/", "_ppin.tsv")) {
			if (f.getName().startsWith("H1HESC"))
				ESCs.add(new PPIN(f.getAbsolutePath()));
		}
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Documents/workspace/ENCODE_check/data/merged_networks/CSHL_nets/", "_ppin.tsv")) {
			if (f.getName().startsWith("H1HESC"))
				ESCs.add(new PPIN(f.getAbsolutePath()));
		}
		
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Documents/workspace/ENCODE_check/data/merged_networks/bodymap_nets/", "_ppin.tsv")) {
			tissues.add(new PPIN(f.getAbsolutePath()));
		}
		
		processPPINs(ESCs, tissues, "/Users/tho/Desktop/ESCs_tissues/", FDR);
	}

}
