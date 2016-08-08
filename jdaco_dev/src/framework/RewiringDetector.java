package framework;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

public class RewiringDetector {

	private final int max_in_thread_iterations = 500;
	private final int max_thread_cycles = 2;
	private final BinomialTest binom_test = new BinomialTest();
	private int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	private final Map<String, RewiringDetectorSample> group1;
	private final Map<String, RewiringDetectorSample> group2;
	private double P_rew;
	private double P_rew_std;
	private double FDR;
	private boolean strict_denominator = false;
	private boolean only_diffnet = false;
	private String organism_database;
	private PrintStream verbose;
	
	private List<Double> P_rews_temp = new LinkedList<>();
	
	// all detected changes
	private Map<StrPair, Integer> differential_network = new HashMap<>();
	private Map<Double, Double> rawp2adjp = new HashMap<>();

	// only stored for significant detected changes
	private List<StrPair> significantly_rewired_interactions = new LinkedList<>(); // sorted by importance, important ones first
	private Map<StrPair, Double> interaction_p_map = new HashMap<>();
	private Map<StrPair, Boolean> interaction_direction_map = new HashMap<>();
	private Map<StrPair, Map<String, Integer>> interaction_reasons_count_map = new HashMap<>();
	private Map<StrPair, Double> interaction_alt_splicing_fraction_map = new HashMap<>();
	private Map<StrPair, List<String>> interaction_sorted_reasons_map = new HashMap<>();
	
	/**
	 * Most powerful constructor
	 * @param group1
	 * @param group2
	 * @param FDR
	 * @param no_threads
	 * @param verbose
	 * @param strict_denominator
	 * @param only_diffnet -> many features not usable
	 */
	public RewiringDetector(Map<String, RewiringDetectorSample> group1, Map<String, RewiringDetectorSample> group2, double FDR, int no_threads, PrintStream verbose, boolean strict_denominator, boolean only_diffnet) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.no_threads = no_threads;
		this.verbose = verbose;
		this.strict_denominator = strict_denominator;
		this.only_diffnet = only_diffnet;
		
		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	/**
	 * Custom constructor
	 * @param group1
	 * @param group2
	 * @param FDR
	 * @param no_threads
	 * @param verbose
	 */
	public RewiringDetector(Map<String, RewiringDetectorSample> group1, Map<String, RewiringDetectorSample> group2, double FDR, int no_threads, PrintStream verbose) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		this.no_threads = no_threads;
		this.verbose = verbose;
		
		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	/**
	 * Simple constructor
	 * @param group1
	 * @param group2
	 * @param FDR
	 */
	public RewiringDetector(Map<String, RewiringDetectorSample> group1, Map<String, RewiringDetectorSample> group2, double FDR) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		
		this.determineGroupwiseDifferences();
		this.assessRewiring();
	}
	
	/*
	 * major init functions
	 */
	
	@SuppressWarnings("unchecked")
	private void determineGroupwiseDifferences() {
		
		/*
		 * preparations and distribution of work
		 */
		
		// build all comparison pairs
		int overall_comparisons = getNumberOfComparisons();
		List<String[]> comparisons = new ArrayList<>(overall_comparisons);
		for (String sample1:this.group1.keySet())
			for (String sample2:this.group2.keySet()) {
				String[] sample_pair = new String[]{sample1, sample2};
				comparisons.add(sample_pair);
			}
		
		// build worker objects from chunks that are at most max_in_thread_iterations big, but at least are well distributed
		List<PPIComparatorTask> comparison_calculations = new ArrayList<>( (overall_comparisons / Math.min(this.max_in_thread_iterations, Math.max(overall_comparisons / this.no_threads, 1))) );
		for ( List<String[]> temp:Utilities.partitionListIntoChunks(comparisons, Math.min(this.max_in_thread_iterations, Math.max(overall_comparisons / this.no_threads, 1))) )
			comparison_calculations.add(new PPIComparatorTask(temp));
		comparisons = null;
		
		/*
		 * calculations
		 */
		
		try {
			// after max_thread_cycles of no_threads workers, memory is cleared up and a report is given
			for (List<PPIComparatorTask> temp_calculations:Utilities.partitionListIntoChunks(comparison_calculations, this.max_thread_cycles * this.no_threads)) {
				
				ExecutorService es = Executors.newFixedThreadPool(this.no_threads);
				
				// compute and collect results
				double tasks = 0;
				Object[] obj = null;
				long start = System.currentTimeMillis();
				for (Future<Object[]> f:es.invokeAll(temp_calculations)) {
					obj = f.get();

					// collect diffnet results
					for (Entry<StrPair, Integer> entry:((Map<StrPair, Integer>) obj[0]).entrySet())
						this.differential_network.put(entry.getKey(), this.differential_network.getOrDefault(entry.getKey(), 0) + entry.getValue());

					// collect individual P_rews
					P_rews_temp.addAll((List<Double>) obj[1]);
					tasks += ((List<Double>) obj[1]).size();

					// some manual help for the GC
					obj[0] = null;
					obj[1] = null;
					obj = null;
				}
				long end = System.currentTimeMillis();
				
				// some occasionally encouraged heap cleaning
				es.shutdown();
				es = null;
				System.gc();
				
				// optional feedback
				if (this.verbose != null) {
					double per_task_duration = (end-start) / tasks;
					this.verbose.println(P_rews_temp.size() + " / " + getNumberOfComparisons() + " comparisons finished ("+ ((int) per_task_duration) + "ms per task, " + ( (int) ((per_task_duration*(getNumberOfComparisons()-P_rews_temp.size()))/1000/60) ) + "min left) ...");
					this.verbose.flush();
				}
			}
			
		} catch (Exception e) {
			System.err.println("Problem during assessment of groupwise differences.");
			e.printStackTrace();
			System.exit(1);
		}
		
		// clean from zeros (results of +1 and -1 ...)
		this.differential_network.entrySet().removeIf( e -> e.getValue().equals(0));
	}

	@SuppressWarnings("unchecked")
	private void assessRewiring() {
		
		if (this.verbose != null) {
			this.verbose.println("Assessing significance of rewiring ...");
			this.verbose.flush();
		}
		
		Map<StrPair, Double> test_map = new HashMap<>();
		Map<Double, LinkedList<StrPair>> p2pair = new HashMap<>();
		
		// compute P_rew statistics and delete temporary data
		this.P_rew = Utilities.getMean(this.P_rews_temp);
		this.P_rew_std = Utilities.getStd(this.P_rews_temp);
		this.P_rews_temp.clear();
		this.P_rews_temp = null;

		// catch this, encourage to use default methodology
		if (this.P_rew > 1.0) {
			System.err.println("P_rew too high, less strict denominator is advised.");
			return;
		}

		int groupwise_comparisons = this.getNumberOfComparisons();

		for (StrPair pair:this.differential_network.keySet()) {
			int v = Math.abs(this.differential_network.get(pair));

			double raw_p = binom_test.binomialTest(groupwise_comparisons, v, this.P_rew, AlternativeHypothesis.GREATER_THAN);
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

		if (this.verbose != null) {
			this.verbose.println("Evaluating significantly rewired interactions ...");
			this.verbose.flush();
		}
		
		k = 1;
		p_values = new ArrayList<>(new HashSet<>(p_values));
		Collections.sort(p_values);
		
		List<ReasonEvalTask> reason_calculations = new LinkedList<>();
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
				
				// reasoning part can be bypassed
				if (!this.only_diffnet)
					reason_calculations.add(new ReasonEvalTask(pair, addition));
				
				k++;
			}
		}
		
		// sort list of significantly rewired interactions
		this.significantly_rewired_interactions.sort( (e2,e1) -> Integer.compare(Math.abs(this.differential_network.get(e1)), Math.abs(this.differential_network.get(e2))) );
		
		// reasoning part can be bypassed
		if (this.only_diffnet)
			return;
		
		// calculate and store reasons
		ExecutorService es = Executors.newFixedThreadPool(this.no_threads);
		try {
			Object[] obj = null;
			for (Future<Object[]> f:es.invokeAll(reason_calculations)) {
				obj = f.get();
				StrPair pair = (StrPair) obj[0];
				this.interaction_reasons_count_map.put(pair, (Map<String,Integer>) obj[1]);
				this.interaction_sorted_reasons_map.put(pair, (List<String>) obj[2]);
				this.interaction_alt_splicing_fraction_map.put(pair, (double) obj[3]);
			}
		} catch (Exception e) {
			System.err.println("Problem during evaluation of rewiring reasons.");
			e.printStackTrace();
			System.exit(1);
		}
		
		// cleanup
		reason_calculations = null;
		es.shutdown();
		es = null;
		System.gc();
	}
	
	
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
		Set<StrPair> i1 = group1.get(sample1).getInteractions();
		Set<StrPair> i2 = group2.get(sample2).getInteractions();
		
		// detailed transcript info
		Map<String, Set<String>> pd1 = group1.get(sample1).getProteinToUsedDomains();
		Map<String, Set<String>> pd2 = group2.get(sample2).getProteinToUsedDomains();
		Map<String, Map<String, Set<String>>> ppd1 = group1.get(sample1).getProteinToProteinByDomains();
		Map<String, Map<String, Set<String>>> ppd2 = group2.get(sample2).getProteinToProteinByDomains();
		Set<String> decayed1 = group1.get(sample1).getDecayProteins();
		Set<String> decayed2 = group2.get(sample2).getDecayProteins();
		
		List<String> reasons = new LinkedList<>();
		
		
		// pre-check1: is a change found in the network-pair?
		if (i1.contains(interaction) == i2.contains(interaction))
			return "no_change";

		// pre-check2: is the specific change found in the network-pair?
		if (addition) { // interaction should not be in sample1 but in sample2; otherwise -> opposing change (no change already caught)
			if ( !(!i1.contains(interaction) && i2.contains(interaction)) )
				return "opposing_change";
		} else { // vice versa
			if ( !(i1.contains(interaction) && !i2.contains(interaction)) )
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

		
		/*
		 * check on transcript-level if indications are given, then go into details
		 */
		
		// first check for relevant splicing of p1: p1 in both samples and Ensembl transcript changed
		if (m1.containsKey(p1) && m2.containsKey(p1) && !m1.get(p1).equals(m2.get(p1))) {
			
			// check details: filter for change
			boolean has_domains1 = pd1.containsKey(p1);
			boolean has_domains2 = pd2.containsKey(p1);
			boolean is_decay1 = decayed1.contains(p1);
			boolean is_decay2 = decayed2.contains(p1);
			
			// if no non-FB domains in neither sample: switch cannot be relevant
			if (has_domains1 || has_domains2) {
				
				Set<String> relevant_changed_domains = null;
				
				if (addition) { // addition of an interaction -> relevant domain(s) in group 2 but not in 1
					if (has_domains1 && has_domains2 && ppd2.get(p1).containsKey(p2)) {
						relevant_changed_domains = new HashSet<>(pd2.get(p1));
						relevant_changed_domains.removeAll(pd1.get(p1));
					} else if (has_domains2 && ppd2.get(p1).containsKey(p2)) {// and consequently no in sample1, vice versa not helpful for an addition
						relevant_changed_domains = new HashSet<>(pd2.get(p1));
					}
				} else { // loss of an interaction -> relevant domain(s) in group 1 but not in 2
					if (has_domains1 && has_domains2 && ppd1.get(p1).containsKey(p2)) {
						relevant_changed_domains = new HashSet<>(pd1.get(p1));
						relevant_changed_domains.removeAll(pd2.get(p1));
					} else if (has_domains1 && ppd1.get(p1).containsKey(p2)) {// and consequently no in sample2, vice versa not helpful for a loss
						relevant_changed_domains = new HashSet<>(pd1.get(p1));
					}
				}
				
				
				// switch can only be relevant if a definitive change in that direction happened for p1
				if (relevant_changed_domains != null && relevant_changed_domains.size() > 0) {
					
					if (addition) { // addition of an interaction -> relevant domain(s) in group 2 but not in 1
						
						// if a domains that then supports the interaction is newly included, then it was a reason
						relevant_changed_domains.retainAll(ppd2.get(p1).get(p2));

						
					} else { // loss of an interaction -> relevant domain(s) in group 1 but not in 2
						
						// if a domain that maintained the interaction fell away, then it was a reason
						relevant_changed_domains.retainAll(ppd1.get(p1).get(p2));
						
					}
					
					if (relevant_changed_domains.size() > 0) {
						reasons.add( p1 + "(" + m1.get(p1) + "->" + m2.get(p1) + ")");
					}
					
				} // end of relevant domain check loop
			} // end of domain annotation check loop
			
			// switch could be due to decay in one sample
			if (is_decay1 || is_decay2) {
				reasons.add( p1 + "(" + m1.get(p1) + "->" + m2.get(p1) + ")");
			}
			
		} // end check for splicing of p1
		
		// analogous check for relevant splicing of p2: only check if p2 in both samples and Ensembl id changed
		if (m1.containsKey(p2) && m2.containsKey(p2) && !m1.get(p2).equals(m2.get(p2))) {
			
			// check details: filter for change
			boolean has_domains1 = pd1.containsKey(p2);
			boolean has_domains2 = pd2.containsKey(p2);
			boolean is_decay1 = decayed1.contains(p2);
			boolean is_decay2 = decayed2.contains(p2);
			
			// if no non-FB domains in neither sample: switch cannot be relevant
			if (has_domains1 || has_domains2) {
				
				Set<String> relevant_changed_domains = null;
				
				if (addition) { // addition of an interaction -> relevant domain(s) in group 2 but not in 1
					if (has_domains1 && has_domains2 && ppd2.get(p2).containsKey(p1)) {
						relevant_changed_domains = new HashSet<>(pd2.get(p2));
						relevant_changed_domains.removeAll(pd1.get(p2));
					} else if (has_domains2 && ppd2.get(p2).containsKey(p1)) {// and consequently no in sample1, vice versa not helpful for an addition
						relevant_changed_domains = new HashSet<>(pd2.get(p2));
					}
				} else { // loss of an interaction -> relevant domain(s) in group 1 but not in 2
					if (has_domains1 && has_domains2 && ppd1.get(p2).containsKey(p1)) {
						relevant_changed_domains = new HashSet<>(pd1.get(p2));
						relevant_changed_domains.removeAll(pd2.get(p2));
					} else if (has_domains1 && ppd1.get(p2).containsKey(p1)) {// and consequently no in sample2, vice versa not helpful for a loss
						relevant_changed_domains = new HashSet<>(pd1.get(p2));
					}
				}

				
				// switch can only be relevant if a definitive change in that direction happened for p2
				if (relevant_changed_domains != null && relevant_changed_domains.size() > 0) {
					
					if (addition) { // addition of an interaction -> relevant domain(s) in group 2 but not in 1
						
						// if a domains that then supports the interaction is newly included, then it was a reason
						relevant_changed_domains.retainAll(ppd2.get(p2).get(p1));

						
					} else { // loss of an interaction -> relevant domain(s) in group 1 but not in 2
						
						// if a domain that maintained the interaction fell away, then it was a reason
						relevant_changed_domains.retainAll(ppd1.get(p2).get(p1));
						
					}
					
					if (relevant_changed_domains.size() > 0) {
						reasons.add( p2 + "(" + m1.get(p2) + "->" + m2.get(p2) + ")");
					}
					
				} // end of relevant domain check loop
			} // end of domain annotation check loop
			
			// switch could be due to decay in one sample
			if (is_decay1 || is_decay2) {
				reasons.add( p2 + "(" + m1.get(p2) + "->" + m2.get(p2) + ")");
			}
			
		} // end check for splicing of p2

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

	public Map<String, RewiringDetectorSample> getGroup1() {
		return group1;
	}

	public Map<String, RewiringDetectorSample> getGroup2() {
		return group2;
	}

	public double getFDR() {
		return FDR;
	}

	public double getP_rew() {
		return P_rew;
	}
	
	public double getP_rew_std() {
		return P_rew_std;
	}
	
	public int getNumberOfComparisons() {
		return this.group1.size() * this.group2.size();
	}
	
	/**
	 * Returns the COMPLETE differential network, not only the significant rewiring events
	 * @return
	 */
	public Map<StrPair, Integer> getDifferentialNetwork() {
		return this.differential_network;
	}
	
	/**
	 * Returns a sorted listed of significantly rewired interactions (high to low)
	 * @return
	 */
	public List<StrPair> getSignificantlyRewiredInteractions() {
		return this.significantly_rewired_interactions;
	}

	public Map<StrPair, Map<String, Integer>> getInteractionReasonsCountMap() {
		return interaction_reasons_count_map;
	}

	public Map<StrPair, Boolean> getInteractionDiectionMap() {
		return this.interaction_direction_map;
	}

	/**
	 * Write information about the differential network to a Cytoscape-usable file
	 * @param diffnet_out_path
	 */
	public void writeDiffnet(String diffnet_out_path) {

		List<String> to_write = new LinkedList<>();
		int groupwise_comparisons = getNumberOfComparisons();

		// header
		to_write.add("Protein1 Protein2 Type Count Probability p-val p-val_adj Reasons AS_fraction");
		for (StrPair pair:this.significantly_rewired_interactions) {

			String sign = "-";
			if (this.interaction_direction_map.get(pair))
				sign = "+";

			double v = Math.abs(this.differential_network.get(pair));
			double p = this.interaction_p_map.get(pair);
			List<String> sorted_reasons = this.interaction_sorted_reasons_map.getOrDefault(pair, new ArrayList<>(1));
			
			if (only_diffnet)
				sorted_reasons.add("-");
			
			double AS_fraction = this.interaction_alt_splicing_fraction_map.getOrDefault(pair, 0.0);

			to_write.add(pair.getL() + " " + pair.getR() + " " + sign + " "
			        + (int) v + " " + v / groupwise_comparisons + " " + p +
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

	/**
	 * Determines matching Ensembl organism database from one of the networks
	 * @return
	 */
	public String getAppropriateOrganismDatabase() {
		
		// only determine once
		if (organism_database == null) {
			for (RewiringDetectorSample rds:this.group1.values()) {
				Set<String> proteins = new HashSet<>();
				for (StrPair pair:rds.getInteractions()) {
					proteins.add(pair.getL());
					proteins.add(pair.getR());
				}
				this.organism_database = DataQuery.getEnsemblOrganismDatabaseFromProteins(proteins);
				break;
			}
		}
		
		return this.organism_database;
	}
	
	public void writeProteinAttributes(String out_path) {
		
		// determine all proteins in diffnet
		Set<String> proteins = new HashSet<>();
		for (StrPair pair:this.significantly_rewired_interactions) {
			proteins.add(pair.getL());
			proteins.add(pair.getR());
		}
		
		// determine entries
		Map<String, Integer> reasons_count_map = determineReasonCountMap();
		List<String> min_mostlikely_reasons = getMinMostLikelyReasons();
		
		Map<String, Integer> expr_reasons = new HashMap<>();
		Map<String, Integer> AS_reasons = new HashMap<>();
		Map<String, Integer> all_reasons = new HashMap<>();
		
		// prefill datastructures
		for (String protein:proteins) {
			expr_reasons.put(protein, 0);
			AS_reasons.put(protein, 0);
			all_reasons.put(protein, 0);
		}
		
		// fill data
		for (String reason:reasons_count_map.keySet()) {
			
			String protein = reason.split("\\(")[0];
			all_reasons.put(protein, all_reasons.get(protein) + reasons_count_map.get(reason));
			
			// distinguish between just expression and alternative splicing
			if (reason.contains("->")) {
				AS_reasons.put(protein, AS_reasons.get(protein) + reasons_count_map.get(reason));
			} else {
				expr_reasons.put(protein, expr_reasons.get(protein) + reasons_count_map.get(reason));
			}
		}
		
		Set<String> min_mostl_reason_proteins = new HashSet<>();
		HashMap<String, String> min_most_reasons = new HashMap<>();
		for (String reason:min_mostlikely_reasons) {
			String protein = reason.split("\\(")[0];
			String tr_reason = reason.split("\\(")[1];
			min_mostl_reason_proteins.add(protein);
			min_most_reasons.put(protein, tr_reason.substring(0, tr_reason.length()-1));
		}
		
		// naming data
		Map<String, String> gene_to_name = DataQuery.getGenesCommonNames(getAppropriateOrganismDatabase());
		Map<String, String> up_to_name = new HashMap<>();
		for (String[] data:DataQuery.getGenesTranscriptsProteins(this.organism_database)) {
			String gene = data[0];
			String protein = data[2];
			up_to_name.put(protein, gene_to_name.get(gene));
		}
		
		/*
		 * prepare output
		 */
		
		List<String> to_write = new LinkedList<>();
		to_write.add("UniProt_ACC Gene_name Part_of_min_reasons Alteration_type Transcriptomic_alteration Score");
		
		for (String protein:proteins) {
			String min_reasons = "no";
			String cause = "/";
			String type = "/";
			if (min_mostl_reason_proteins.contains(protein)) {
				min_reasons = "yes";
				cause = min_most_reasons.get(protein);
				type = "DE";
				if (cause.contains("->"))
					type = "AS";
			}
			
			to_write.add(protein + " " + up_to_name.get(protein) + " " + min_reasons + " " + type + " " + cause + " " + all_reasons.get(protein));
		}
		
		Utilities.writeEntries(to_write, out_path);
	}
	
	private Map<String, Integer> determineReasonCountMap() {
		Map<String, Integer> reason_count_map = new HashMap<>();
		
		// basically the same code as in the heuristic
		for (StrPair IA:this.interaction_reasons_count_map.keySet()) {
			Map<String, Integer> count_map = this.interaction_reasons_count_map.get(IA);
			for (String reason:count_map.keySet()) {
				if (reason.startsWith("no_") || reason.startsWith("opp")) // kick non-helpful reasons
					continue;
				for (String ex_reason:reason.split("/")) {
					if (!reason_count_map.containsKey(ex_reason))
						reason_count_map.put(ex_reason, 0);
					reason_count_map.put(ex_reason, reason_count_map.get(ex_reason) + count_map.get(reason));
				}
			}
		}
		
		return reason_count_map;
	}
	
	/**
	 * Heuristically calculates the smallest set of reasons necessary to explain all changes,
	 * output as a list of strings PROTEIN(difference)
	 * @return
	 */
	public List<String> getMinMostLikelyReasons() {
		return getMinMostLikelyReasons(null);
	}
	
	/**
	 * Heuristically calculates the smallest set of reasons necessary to explain all changes of a subset of rewired interactions,
	 * output as a list of strings PROTEIN(difference)
	 * @return
	 */
	public List<String> getMinMostLikelyReasons(Set<StrPair> subset) {
		
		Map<String, Set<StrPair>> reason_IA_map = new HashMap<>();
		Map<String, Integer> reason_count_map = new HashMap<>();
		
		// fills datastructures, afterwards we know which reasons affect which interactions and each reasons occurrence is counted to give it an importance score
		for (StrPair IA:this.interaction_reasons_count_map.keySet()) {
			
			// if specified, only count rewired interactions of the specified subset
			if (subset != null && !subset.contains(IA))
				continue;
			
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
		
		// convert maximum to minimum weight problem
		int max = Collections.max(reason_count_map.values()) + 1;
		for (String reason:reason_count_map.keySet())
			reason_count_map.put(reason, max - reason_count_map.get(reason));
		
		Set<StrPair> unsatisfied_IAs = new HashSet<>(this.interaction_reasons_count_map.keySet());
		
		// if specified, only account for subset
		if (subset != null)
			unsatisfied_IAs.retainAll(subset);
		
		List<String> result = new LinkedList<>();
		
		// while not all rewiring events are explained
		while (unsatisfied_IAs.size() > 0) {
			String best_reason = null;
			double best_weight_per_IA = Double.MAX_VALUE;
			Set<StrPair> best_additionally_satisfied_IAs = null;
			
			// determine minimum weight per element choice
			for (String reason:reason_count_map.keySet()) {
				double weight = reason_count_map.get(reason);
				Set<StrPair> additionally_satisfied_IAs = new HashSet<>(reason_IA_map.get(reason));
				additionally_satisfied_IAs.retainAll(unsatisfied_IAs);
				
				if (additionally_satisfied_IAs.size() == 0)
					continue;
				double weight_per_IA = weight / additionally_satisfied_IAs.size();
				
				if (weight_per_IA < best_weight_per_IA) {
					best_weight_per_IA = weight_per_IA;
					best_additionally_satisfied_IAs = additionally_satisfied_IAs;
					best_reason = reason;
				}
			}
			
			// realize locally optimal choice
			result.add(best_reason);
			reason_count_map.remove(best_reason);
			unsatisfied_IAs.removeAll(best_additionally_satisfied_IAs);
		}
		
		return result;
	}
	
	
	/*
	 * Helper classes for parallel processing
	 */
	
	class PPIComparatorTask implements Callable<Object[]> {
		private List<String[]> comparisons;

		public PPIComparatorTask(List<String[]> comparisons) {
			this.comparisons = comparisons;
		}

		@Override
		public Object[] call() throws Exception {

			Map<StrPair, Integer> diff_temp = new HashMap<>(4096); // less number of size-change events
			List<Double> P_rews_temp = new ArrayList<>(comparisons.size());

			for (String[] samples:comparisons) {
				Set<StrPair> i1 = group1.get(samples[0]).getInteractions();
				Set<StrPair> i2 = group2.get(samples[1]).getInteractions();

				Set<StrPair> added_interactions = new HashSet<>(i2);
				added_interactions.removeAll(i1);

				for (StrPair pair:added_interactions)
					diff_temp.put(pair, diff_temp.getOrDefault(pair, 0) + 1 );

				Set<StrPair> lost_interactions = new HashSet<>(i1);
				lost_interactions.removeAll(i2);

				for (StrPair pair:lost_interactions)
					diff_temp.put(pair, diff_temp.getOrDefault(pair, 0) - 1 );

				Set<StrPair> union = new HashSet<>(i1);
				union.addAll(i2);
				double denominator = (double) union.size();

				if (strict_denominator)
					denominator = Math.min(i1.size(), i2.size());

				P_rews_temp.add(( (double) added_interactions.size() + lost_interactions.size() ) / denominator );
			}

			return new Object[]{diff_temp, P_rews_temp};
		}
	}
	
	class ReasonEvalTask implements Callable<Object[]> {
		private StrPair pair;
		private boolean addition;

		public ReasonEvalTask(StrPair pair, boolean addition) {
			this.pair = pair;
			this.addition = addition;
		}

		@Override
		public Object[] call() throws Exception {
			// check the exact reason for the difference between all samples and count them
			Map<String, Integer> reasons_count = new HashMap<>();
			for (String sample1:group1.keySet())
				for (String sample2:group2.keySet()) {
					String reason = checkReason(addition, pair, sample1, sample2);
					reasons_count.put(reason, reasons_count.getOrDefault(reason, 0) + 1);
				}
			
			// building sorted (descending) reasons-string for output
			List<String> sorted_reasons = new ArrayList<>(reasons_count.keySet());
			sorted_reasons.sort((e2,e1) -> reasons_count.get(e1).compareTo(reasons_count.get(e2)));
			sorted_reasons.replaceAll(e->e + ":" + reasons_count.get(e));
			
			// return everything
			return new Object[]{pair, reasons_count, sorted_reasons, determineAltSplicingFraction(reasons_count)};
		}
	}
}