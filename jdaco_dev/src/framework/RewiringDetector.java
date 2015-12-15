package framework;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

public class RewiringDetector {

	private static BinomialTest binom_test = new BinomialTest();
	private int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
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
	private Map<StrPair, Double> interaction_alt_splicing_fraction_map = new HashMap<>();
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
		PPIN ppin1 = group1.get(sample1).getPPIN();
		PPIN ppin2 = group2.get(sample2).getPPIN();

		List<String> reasons = new LinkedList<>();

		// pre-check1: is a change found in the network-pair?
		if (ppin1.getWeights().containsKey(interaction) == ppin2.getWeights().containsKey(interaction))
			return "no_change";

		// pre-check2: is the specific change found in the network-pair?
		if (addition) { // interaction should not be in sample1 but in sample2; otherwise -> opposing change (no change already caught)
			if ( !(!ppin1.getWeights().containsKey(interaction) && ppin2.getWeights().containsKey(interaction)) )
				return "opposing_change";
		} else { // vice versa
			if ( !(ppin1.getWeights().containsKey(interaction) && !ppin2.getWeights().containsKey(interaction)) )
				return "opposing_change";
		}

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

		return String.join("/", reasons);
	}

	/**
	 * For all relevant reasons_counts determined for a protein pair,
	 * determines the fraction of relevant alternative splicing events. (1-this fraction=fraction of differential gene expression)
	 * If both interaction partners show changes, they contribute equally (a half count each) to the full count.
	 * @param reasons_count
	 * @return
	 */
	private double determineAltSplicingFraction(Map<String, Integer> reasons_count) {

		double all_events = 0.0;
		double alt_splicing_events = 0.0;

		for (String reason:reasons_count.keySet()) {
			double current_count = reasons_count.get(reason);

			// change could be in one or both proteins
			String[] contributions = reason.split("/");

			// full count if only one protein changed, half of counts for the reason of the individual changes if both changed
			for (String ind_reasons:contributions) {

				if (ind_reasons.contains("->"))
					alt_splicing_events += current_count / contributions.length;

				if (!ind_reasons.endsWith("_change")) // no_change / oppposing_change
					all_events += current_count / contributions.length;
			}
		}

		return alt_splicing_events / all_events;
	}
	
	public RewiringDetector(Map<String, ConstructedNetworks> group1, Map<String, ConstructedNetworks> group2, double FDR, int no_threads) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.no_threads = no_threads;
		
		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	public RewiringDetector(Map<String, ConstructedNetworks> group1, Map<String, ConstructedNetworks> group2, double FDR) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;

		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	public RewiringDetector(Map<String, ConstructedNetworks> group1, Map<String, ConstructedNetworks> group2, double FDR, int no_threads, boolean strict_denominator) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.no_threads = no_threads;
		this.strict_denominator = strict_denominator;

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

	@SuppressWarnings("unchecked")
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

		// parallel assessment of pairwise differences
		List<PPIComparatorTask> comparison_calculations = new LinkedList<>();
		for (String sample1:this.group1.keySet())
			for (String sample2:group2.keySet())
				comparison_calculations.add(new PPIComparatorTask(sample1, sample2));

		ExecutorService es = Executors.newFixedThreadPool(this.no_threads);

		try {
			for (Future<Object[]> f:es.invokeAll(comparison_calculations)) {
				Object[] obj = f.get();

				// added interactions
				for (StrPair pair: (Set<StrPair>) obj[0]) {
					if(!overall_added.containsKey(pair))
						overall_added.put(pair, 0);
					overall_added.put(pair, overall_added.get(pair) + 1 );
				}

				// lost interactions
				for (StrPair pair: (Set<StrPair>) obj[1]) {
					if(!overall_lost.containsKey(pair))
						overall_lost.put(pair, 0);
					overall_lost.put(pair, overall_lost.get(pair) + 1 );
				}

			}
		} catch (Exception e1) {
			System.err.println("Problem during assessment of groupwise differences.");
			e1.printStackTrace();
			System.exit(1);
		}
		es.shutdown();

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

				// compute ratio of differential gene/transcript expression
				this.interaction_alt_splicing_fraction_map.put(pair, determineAltSplicingFraction(reasons_count));

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
		to_write.add("Protein1 Protein2 Type Count Probability p-val p-val_adj Reasons AS_fraction");
		for (StrPair pair:this.significantly_rewired_interactions) {

			String sign = "-";
			if (this.interaction_direction_map.get(pair))
				sign = "+";

			double v = this.differential_network.get(pair);
			double p = this.interaction_p_map.get(pair);
			List<String> sorted_reasons = this.interaction_sorted_reasons_map.get(pair);
			double AS_fraction = this.interaction_alt_splicing_fraction_map.get(pair);

			to_write.add(pair.getL() + " " + pair.getR() + " " + sign + " "
			        + (int) Math.abs(v) + " " + Math.abs(v / groupwise_comparisons) + " " + p +
			        " " + this.rawp2adjp.get(p) + " " + String.join(",", sorted_reasons) + " " + AS_fraction);
		}

		Utilities.writeEntries(to_write, diffnet_out_path);
	}

	/**
	 * Outputs a map of proteins and affected interactions where alternative splicing events of those proteins are the reasons for a difference regarding specific interactions.
	 * The first boolean switch specifies if only the most appearing reason is counted or every instance,
	 * the second switch if the interaction should be only be counted if it is mainly (at least 50%) driven by AS events.
	 * @param only_major
	 * @return
	 */
	public Map<String, List<StrPair>> determineAltSplicingSwitches(boolean only_major, boolean only_mainly_AS) {
		Map<String, List<StrPair>> relevant_subset = new HashMap<>();

		for (StrPair pair: this.interaction_sorted_reasons_map.keySet()) {

			// if specified, kicks all diff. interactions that are not mainly (at least 50%) driven by AS events
			if (only_mainly_AS && this.interaction_alt_splicing_fraction_map.get(pair) < 0.5)
				continue;

			int count_last = 0;
			for (String unprocessed_reason:this.interaction_sorted_reasons_map.get(pair)) {
				String preprocessed_reason = unprocessed_reason.split(":")[0];

				if (only_major) {
					int count = this.interaction_reasons_count_map.get(pair).get(preprocessed_reason);

					// only major reason or any reason that appears equally often
					if (count < count_last)
						break;

					count_last = count;
				}


				// if no/wrong kind of change: skip
				if (preprocessed_reason.startsWith("no_") || preprocessed_reason.startsWith("opp"))
					continue;

				for (String s:preprocessed_reason.split("/")) {
					String[] split_temp = s.split("\\(");
					String protein = split_temp[0];
					String reason = split_temp[1].substring(0, split_temp[1].length() - 1);

					if (reason.contains("->")) {
						if (!relevant_subset.containsKey(protein))
							relevant_subset.put(protein, new LinkedList<>());
						relevant_subset.get(protein).add(pair);
					}
				} // within reason split

			} // reasons iteration

		} // per interaction iteration

		return relevant_subset;
	}

	// helper class for faster comparison
	class PPIComparatorTask implements Callable<Object[]> {
		private String sample1;
		private String sample2;

		public PPIComparatorTask(String sample1, String sample2) {
			this.sample1 = sample1;
			this.sample2 = sample2;
		}

		@Override
		public Object[] call() throws Exception {
			PPIN ppin1 = group1.get(sample1).getPPIN();
			PPIN ppin2 = group2.get(sample2).getPPIN();
			Set<StrPair> added_interactions = ppin2.removeAllIAs(ppin1);
			Set<StrPair> lost_interactions = ppin1.removeAllIAs(ppin2);
			double denominator = (double) ppin1.mergeAllIAs(ppin2).size();

			if (strict_denominator)
				denominator = Math.min(ppin1.getSizes()[1], ppin2.getSizes()[1]);

			P_rews.put(sample1 + "-" + sample2, ( (double) added_interactions.size() + lost_interactions.size() ) / denominator);

			return new Object[]{added_interactions, lost_interactions};
		}

	}
	
	/**
	 * Heuristically calculates the smallest set of reasons necessary to explain all changes,
	 * output as a list of strings PROTEIN(difference):count
	 * @return
	 */
	public List<String> getMinMostLikelyReasons() {
		
		Map<String, Set<StrPair>> reason_IA_map = new HashMap<>();
		Map<String, Integer> reason_count_map = new HashMap<>();
		
		// fills datastructures, afterwards we know which reasons affect which interactions and each reasons occurance is counted to give it an importance score
		for (StrPair IA:this.interaction_reasons_count_map.keySet()) {
			Map<String, Integer> count_map = this.interaction_reasons_count_map.get(IA);
			
			for (String reason:count_map.keySet()) {
				if (reason.startsWith("no_") || reason.startsWith("opp")) // kick non-helpful reasons
					continue;
				for (String ex_reason:reason.split("/")) {
					
					if (!reason_IA_map.containsKey(ex_reason))
						reason_IA_map.put(ex_reason, new HashSet<StrPair>());
					reason_IA_map.get(ex_reason).add(IA);
					
					if (!reason_count_map.containsKey(ex_reason))
						reason_count_map.put(ex_reason, 0);
					reason_count_map.put(ex_reason, reason_count_map.get(ex_reason) + count_map.get(reason));
				}
			}
		}
		
		List<String> reasons = new ArrayList<>(reason_IA_map.keySet());
		reasons.sort((String s1, String s2) -> reason_count_map.get(s2).compareTo(reason_count_map.get(s1))); // sorts from highest to lowest
		
		Set<StrPair> unsatisfied_IAs = new HashSet<>(this.interaction_reasons_map.keySet());
		List<String> result = new LinkedList<>();
		
		for (String reason:reasons) {
			int to_sat_pre = unsatisfied_IAs.size();
			unsatisfied_IAs.removeAll(reason_IA_map.get(reason));
			int to_sat_post = unsatisfied_IAs.size();
			
			if (to_sat_pre != to_sat_post)
				result.add(reason + ":" + reason_count_map.get(reason));
			
			if (to_sat_post == 0)
				break;
		}
		
		return result;
	}
}