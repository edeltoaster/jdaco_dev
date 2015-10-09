package framework;

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

public class RewiringDetector {
	
	private static BinomialTest binom_test = new BinomialTest();
	private Map<String, ConstructedNetworks> group1;
	private Map<String, ConstructedNetworks> group2;
	private double FDR;
	private String out_folder;
	
	// TODO: do right
	
	private String checkReason(boolean addition, StrPair interaction, String sample1, String sample2) {
		String p1 = interaction.getL();
		String p2 = interaction.getR();
		Map<String, String> m1 = group1.get(sample1).getProteinToAssumedTranscriptMap();
		Map<String, String> m2 = group2.get(sample2).getProteinToAssumedTranscriptMap();
		List<String> reasons = new LinkedList<>();
		
		// check protein-level
		if (addition) { // interaction found in sample2 but not sample1 -> maybe one or both proteins not expressed in sample1?
			if (!m1.containsKey(p1))
				reasons.add(p1);
			if (!m1.containsKey(p2))
				reasons.add(p2);
		} else { // vice versa
			if (!m2.containsKey(p1))
				reasons.add(p1);
			if (!m2.containsKey(p2))
				reasons.add(p2);
		}
		
		if (reasons.size() != 0)
			return String.join("/", reasons);
		
		// else: check on transcript-level
		if (!m1.get(p1).equals(m2.get(p1))) {
			reasons.add( p1 + ":" + m1.get(p1) + "->" + m2.get(p1) );
		}
		
		if (!m1.get(p2).equals(m2.get(p2))) {
			reasons.add( p2 + ":" + m1.get(p2) + "->" + m2.get(p2) );
		}
		
		return String.join("/", reasons);
	}
	
	// stump
	public RewiringDetector(Map<String, ConstructedNetworks> group1, Map<String, ConstructedNetworks> group2, double FDR, String out_folder) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.out_folder = out_folder;
		
		process();
	}
	
	private void process() {
		Map<StrPair, Integer> overall_added = new HashMap<>();
		Map<StrPair, Integer> overall_lost = new HashMap<>();
		
		System.out.println("Processing " + this.out_folder + ":");
		System.out.println("Comparing " + this.group1.size() + " with " + this.group2.size() + " samples.");

		List<Double> g1_nodes = new LinkedList<>();
		List<Double> g2_nodes = new LinkedList<>();
		List<Double> g1_edges = new LinkedList<>();
		List<Double> g2_edges = new LinkedList<>();
		List<Double> P_rew = new LinkedList<>();
		
		// sizes / edges
		for (ConstructedNetworks cn:this.group1.values()) {
			PPIN ppin = cn.getPPIN();
			g1_nodes.add( (double) ppin.getSizes()[0]);
			g1_edges.add((double) ppin.getSizes()[1]);
		}
		
		for (ConstructedNetworks cn:this.group2.values()) {
			PPIN ppin = cn.getPPIN();
			g2_nodes.add( (double) ppin.getSizes()[0]);
			g2_edges.add((double) ppin.getSizes()[1]);
		}
		
		// comparison
		for (ConstructedNetworks cn1:group1.values()) {
			PPIN ppin1 = cn1.getPPIN();
			for (ConstructedNetworks cn2:group2.values()) {
				PPIN ppin2 = cn2.getPPIN();
				Set<StrPair> added_interactions = ppin2.removeAllIAs(ppin1).getInteractions();
				Set<StrPair> lost_interactions = ppin1.removeAllIAs(ppin2).getInteractions();
				P_rew.add( ((double)added_interactions.size()+ lost_interactions.size())/Math.min(ppin1.getSizes()[1], ppin2.getSizes()[1]));
				
				
				// count
				
				for (StrPair pair:added_interactions) {
					if(!overall_added.containsKey(pair))
						overall_added.put(pair, 0);
					overall_added.put(pair, overall_added.get(pair) + 1 );
				}
				
				for (StrPair pair:lost_interactions) {
					if(!overall_lost.containsKey(pair))
						overall_lost.put(pair, 0);
					overall_lost.put(pair, overall_lost.get(pair) + 1 );
				}
				
			}
		}
		
		// do diff.
		Map<StrPair, Double> diff = new HashMap<>();
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
		
		writeOut(diff, P_rew.size(), out_folder + "diffnet.tsv", Utilities.getMean(P_rew), FDR);
		
		// write P_rews for statistics
		List<String> temp_list = new LinkedList<>();
		for (double d:P_rew)
			temp_list.add(Double.toString(d));
		
		Utilities.writeEntries(temp_list, out_folder + "P_rew_distr.txt");
		
		System.out.println();
	}
	
	// stump
	private void writeOut(Map<StrPair, Double> diff, int samples, String out_file, double P_rew, double FDR) {
		List<String> helper = new LinkedList<>();
		Map<StrPair, Double> test_map = new HashMap<>();
		Map<Double, LinkedList<StrPair>> p2pair = new HashMap<>();
		for (StrPair pair:diff.keySet()) {
			int v = (int) Math.abs(diff.get(pair));
			
			if (v == 0)
				continue;
			
			double raw_p = binom_test.binomialTest(samples, Math.abs(v), P_rew, AlternativeHypothesis.GREATER_THAN);
			test_map.put(pair, raw_p);
			
			if (!p2pair.containsKey(raw_p))
				p2pair.put(raw_p, new LinkedList<>());
			p2pair.get(raw_p).add(pair);
		}
		
		List<Double> p_values = new ArrayList<>(test_map.values());
		int m = p_values.size();
		Collections.sort(p_values);
		int k = 1;
		int largest_k = -1;
		Map<Double, Double> rawp2adjp = new HashMap<>();
		
		// find largest k
		for (double p:p_values) {
			if (p <= k * FDR / m)
				largest_k = k;
			rawp2adjp.put(p, (p* m) / k); // if multiple have the same rank, take the biggest
			k++;
		}
		
		helper.add("Protein1 Protein2 Type Count Probability p-val p-val_adj reason");
		k = 1;
		p_values = new LinkedList<Double>(new HashSet<>(p_values));
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
				boolean addition = false;
				if (Math.signum(v) == +1) {
					sign = "+";
					addition = true;
				}
				
				// check for all samples
				List<String> reasons = new LinkedList<>();
				for (String sample1:this.group1.keySet())
					for (String sample2:this.group2.keySet()) {
						reasons.add(checkReason(addition, pair, sample1, sample2));
					}
				
				helper.add(pair.getL() + " " + pair.getR() + " " + sign + " " + (int) Math.abs(v) + " " + Math.abs(v / samples) + " " + p + " " + rawp2adjp.get(p) + " " + String.join(",", reasons));
				k++;
			}
		}
		
		System.out.println(k + " diff-IAs with FDR below " + FDR);
		Utilities.writeEntries(helper, out_file);
	}
}
