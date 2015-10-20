package framework;

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
	private boolean strict_denominator = false;
	
	private Map<String, Double> g1_nodes = new HashMap<>();
	private Map<String, Double> g2_nodes = new HashMap<>();
	private Map<String, Double> g1_edges = new HashMap<>();
	private Map<String, Double> g2_edges = new HashMap<>();
	private Map<String, Double> P_rews = new HashMap<>();
	
	// all detected changes
	private Map<StrPair, Double> differential_network = new HashMap<>();
	private Map<Double, Double> rawp2adjp = new HashMap<>();
	
	// only stored for significant detected changes
	private List<StrPair> significantly_rewired_interactions = new LinkedList<>();
	private Map<StrPair, Double> interaction_p_map = new HashMap<>();
	private Map<StrPair, Boolean> interaction_direction_map = new HashMap<>();
	private Map<StrPair, Map<String, List<String>>> interaction_reasons_map = new HashMap<>();
	private Map<StrPair, Map<String, Integer>> interaction_reasons_count_map = new HashMap<>();
	private Map<StrPair, List<String>> interaction_sorted_reasons_map = new HashMap<>();
	
	/**
	 * For a given change in an interaction, check for the reason in both samples and return a report
	 * @param addition
	 * @param interaction
	 * @param sample1
	 * @param sample2
	 * @return
	 */
	private String checkReason(boolean addition, StrPair interaction, String sample1, String sample2) {
		String p1 = interaction.getL();
		String p2 = interaction.getR();
		Map<String, String> m1 = group1.get(sample1).getProteinToAssumedTranscriptMap();
		Map<String, String> m2 = group2.get(sample2).getProteinToAssumedTranscriptMap();
		List<String> reasons = new LinkedList<>();
		
		// check protein-level
		if (addition) { // interaction found in sample2 but not sample1 -> maybe one or both proteins not expressed in sample1?
			if (!m1.containsKey(p1) && m2.containsKey(p1))
				reasons.add(p1 + "(gain)");
			if (!m1.containsKey(p2) && m2.containsKey(p2))
				reasons.add(p2 + "(gain)");
		} else { // vice versa
			if (!m2.containsKey(p1) && m1.containsKey(p1))
				reasons.add(p1 + "(loss)");
			if (!m2.containsKey(p2) && m1.containsKey(p2))
				reasons.add(p2 + "(loss)");
		}
		
		if (reasons.size() != 0)
			return String.join("/", reasons);
		
		// else: check on transcript-level
		if (m1.containsKey(p1) && m2.containsKey(p1) && !m1.get(p1).equals(m2.get(p1))) {
			reasons.add( p1 + "(" + m1.get(p1) + "->" + m2.get(p1) + ")");
		}
		
		if (m1.containsKey(p2) && m2.containsKey(p2) && !m1.get(p2).equals(m2.get(p2))) {
			reasons.add( p2 + "(" + m1.get(p2) + "->" + m2.get(p2) + ")");
		}
		
		// note that there was no change observed
		if (reasons.size() == 0)
			return "no_change";
		
		return String.join("/", reasons);
	}
	
	public RewiringDetector(Map<String, ConstructedNetworks> group1, Map<String, ConstructedNetworks> group2, double FDR) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		
		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	public RewiringDetector(Map<String, ConstructedNetworks> group1, Map<String, ConstructedNetworks> group2, double FDR, boolean strict_denominator) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.strict_denominator = strict_denominator;
		
		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	private void determineGroupwiseDifferences() {
		Map<StrPair, Integer> overall_added = new HashMap<>();
		Map<StrPair, Integer> overall_lost = new HashMap<>();
		
		// sizes / edges
		for (String sample:this.group1.keySet()) {
			PPIN ppin = this.group1.get(sample).getPPIN();
			this.g1_nodes.put(sample, (double) ppin.getSizes()[0]);
			this.g1_edges.put(sample, (double) ppin.getSizes()[1]);
		}
		
		for (String sample:this.group2.keySet()) {
			PPIN ppin = this.group2.get(sample).getPPIN();
			this.g2_nodes.put(sample, (double) ppin.getSizes()[0]);
			this.g2_edges.put(sample, (double) ppin.getSizes()[1]);
		}
		
		// comparison
		for (String sample1:this.group1.keySet()) {
			PPIN ppin1 = this.group1.get(sample1).getPPIN();
			for (String sample2:group2.keySet()) {
				PPIN ppin2 = this.group2.get(sample2).getPPIN();
				Set<StrPair> added_interactions = ppin2.removeAllIAs(ppin1).getInteractions();
				Set<StrPair> lost_interactions = ppin1.removeAllIAs(ppin2).getInteractions();
				double denominator = (double) ppin1.mergeAllIAs(ppin2).getSizes()[1];
				
				if (this.strict_denominator)
					denominator = Math.min(ppin1.getSizes()[1], ppin2.getSizes()[1]);
				
				this.P_rews.put(sample1 + "-" + sample2, ( (double) added_interactions.size() + lost_interactions.size() ) / denominator);
				
				
				// actual count
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
		
		// fill differential network
		for (StrPair pair:overall_added.keySet())
			this.differential_network.put(pair, this.differential_network.getOrDefault(pair, 0.0) + overall_added.get(pair));
		
		for (StrPair pair:overall_lost.keySet())
			this.differential_network.put(pair, this.differential_network.getOrDefault(pair, 0.0) - overall_lost.get(pair));
		
		// clean from zeros (results of +1 and -1 ...)
		this.differential_network.entrySet().removeIf( e -> e.getValue().equals(0.0));
	}
	
	private void assessRewiring() {
		Map<StrPair, Double> test_map = new HashMap<>();
		Map<Double, LinkedList<StrPair>> p2pair = new HashMap<>();
		double P_rew = Utilities.getMean(this.P_rews.values());
		
		// workaround: can happen in VERY FEW cases
		if (P_rew > 1.0) {
			System.err.println("P_rew too high, less strict denominator is advised.");
			return;
		}
		
		int groupwise_comparisons = this.group1.size() * this.group2.size();
		
		for (StrPair pair:this.differential_network.keySet()) {
			int v = (int) Math.abs(this.differential_network.get(pair));
			
			double raw_p = binom_test.binomialTest(groupwise_comparisons, v, P_rew, AlternativeHypothesis.GREATER_THAN);
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
		
		// find largest k
		for (double p:p_values) {
			if (p <= k * FDR / m)
				largest_k = k;
			this.rawp2adjp.put(p, (p* m) / k); // if multiple have the same rank, take the biggest
			k++;
		}
		
		k = 1;
		p_values = new LinkedList<>(new HashSet<>(p_values));
		Collections.sort(p_values);
		for (double p:p_values) {
			if (k > largest_k) { // remaining ones not deemed significant
				// reset counter for output
				k--;
				break;
			}

			for (StrPair pair:p2pair.get(p)) {
				
				this.interaction_p_map.put(pair, p);
				this.significantly_rewired_interactions.add(pair);
				
				double v = this.differential_network.get(pair);
				boolean addition = false;
				if (Math.signum(v) == +1) {
					addition = true;
				}
				this.interaction_direction_map.put(pair, addition);
				
				// check the exact reason for the difference between all samples
				Map<String, List<String>> reasons = new HashMap<>();
				for (String sample1:this.group1.keySet())
					for (String sample2:this.group2.keySet()) {
						String reason = checkReason(addition, pair, sample1, sample2);
						if (!reasons.containsKey(reason))
							reasons.put(reason, new LinkedList<>());
						reasons.get(reason).add(sample1 + "-" + sample2);
					}
				
				this.interaction_reasons_map.put(pair, reasons);
				
				// building sorted (descending) reasons-string for output
				Map<String, Integer> reasons_count = new HashMap<>();
				for (String reason: reasons.keySet())
					reasons_count.put(reason, reasons.get(reason).size());
				
				this.interaction_reasons_count_map.put(pair, reasons_count);
				
				List<String> sorted_reasons = new LinkedList<>(reasons.keySet());
				sorted_reasons.sort((e2,e1) -> reasons_count.get(e1).compareTo(reasons_count.get(e2)));
				sorted_reasons.replaceAll(e->e + ":" + reasons_count.get(e));
				
				this.interaction_sorted_reasons_map.put(pair, sorted_reasons);
				
				k++;
			}
		}
	}
	// initial preprocessing up to here

	public Map<String, ConstructedNetworks> getGroup1() {
		return group1;
	}

	public Map<String, ConstructedNetworks> getGroup2() {
		return group2;
	}

	public double getFDR() {
		return FDR;
	}

	public Map<String, Double> getG1Nodes() {
		return g1_nodes;
	}

	public Map<String, Double> getG2Nodes() {
		return g2_nodes;
	}

	public Map<String, Double> getG1Edges() {
		return g1_edges;
	}

	public Map<String, Double> getG2Edges() {
		return g2_edges;
	}

	public Map<String, Double> getP_rews() {
		return P_rews;
	}

	public Map<StrPair, Double> getDifferentialNetwork() {
		return differential_network;
	}

	public Map<StrPair, Map<String, List<String>>> getInteractionReasonsMap() {
		return interaction_reasons_map;
	}
	
	public Map<StrPair, Map<String, Integer>> getInteractionReasonsCountMap() {
		return interaction_reasons_count_map;
	}
	
	
	/**
	 * Write information about the differential network to a Cytoscape-usable file
	 * @param diffnet_out_path
	 */
	public void writeDiffnet(String diffnet_out_path) {
		
		List<String> to_write = new LinkedList<>();
		int groupwise_comparisons = this.group1.size() * this.group2.size();
		
		// header
		to_write.add("Protein1 Protein2 Type Count Probability p-val p-val_adj Reasons");
		for (StrPair pair:this.significantly_rewired_interactions) {
			
			String sign = "-";
			if (this.interaction_direction_map.get(pair))
				sign = "+";
			
			double v = this.differential_network.get(pair);
			double p = this.interaction_p_map.get(pair);
			List<String> sorted_reasons = this.interaction_sorted_reasons_map.get(pair);
			
			to_write.add(pair.getL() + " " + pair.getR() + " " + sign + " " 
			        + (int) Math.abs(v) + " " + Math.abs(v / groupwise_comparisons) + " " + p + " " + this.rawp2adjp.get(p) + " " + String.join(",", sorted_reasons));
		}
		
		Utilities.writeEntries(to_write, diffnet_out_path);
	}
	
	/**
	 * Outputs a map of proteins and affected interactions where alternative splicing events of those proteins are the reasons for a difference regarding specific interactions.
	 * A boolean switch defines if only the most appearing reason is counted or every instance
	 * @param only_major
	 * @return
	 */
	public Map<String, List<StrPair>> determineAltSplicingSwitches(boolean only_major) {
		Map<String, List<StrPair>> relevant_subset = new HashMap<>();
		
		for (StrPair pair: this.interaction_sorted_reasons_map.keySet()) {
			
			for (String unprocessed_reason:this.interaction_sorted_reasons_map.get(pair)) {
			
				// if very different reasons for this change: skip
				if (unprocessed_reason.startsWith("no_change"))
					continue;
				
				for (String s:unprocessed_reason.split("/")) {
					String[] split_temp = s.split(":")[0].split("\\(");
					String protein = split_temp[0];
					String reason = split_temp[1].substring(0, split_temp[1].length() - 1);
					
					if (reason.contains("->")) {
						if (!relevant_subset.containsKey(protein))
							relevant_subset.put(protein, new LinkedList<>());
						relevant_subset.get(protein).add(pair);
					}
				}
				
				if (only_major)
					break;
				
			}
		}
		
		return relevant_subset;
	}
}
