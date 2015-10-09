package framework;

import java.io.PrintStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
/**
 * DACO algorithm execution
 * @author Thorsten Will
 */
public class DACO {
	
	// networks
	private final PPIN ppi;
	private final DDIN ddi;
	
	// collections for calculations need to be thread-safe
	private Set<HashSet<String>> temp_results = Collections.newSetFromMap(new ConcurrentHashMap<HashSet<String>, Boolean>());
	private int number_of_threads = Runtime.getRuntime().availableProcessors() / 2;
	private ThreadPoolExecutor pool = new ThreadPoolExecutor(number_of_threads, number_of_threads, 5, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
	private int max_depth_of_search = 10;
	private double pair_building_threshold = 0.75;
	private double prob_cutoff = 0.5;
	private PrintStream verbose;
	private int compute_timeout = 0; // minutes
	
	/**
	 * Some convenient constructors
	 */
	public DACO(DDIN ddi, PPIN ppi) {
		this.ddi = ddi;
		this.ppi = ppi;
	}
	
	public DACO(ConstructedNetworks constr) {
		this.ddi = constr.getDDIN();
		this.ppi = constr.getPPIN();
	}
	
	public DACO(DDIN ddi, PPIN ppi, int number_of_threads, int max_depth_of_search, double pair_building_threshold, double probability_cutoff, PrintStream verbose) {
		this.ddi = ddi;
		this.ppi = ppi;
		
		this.number_of_threads = number_of_threads;
		this.max_depth_of_search = max_depth_of_search;
		this.pair_building_threshold = pair_building_threshold;
		this.prob_cutoff = probability_cutoff;
		this.verbose = verbose;
	}
	
	public DACO(DDIN ddi, PPIN ppi, int number_of_threads, int max_depth_of_search, double pair_building_threshold, double probability_cutoff, PrintStream verbose, int compute_timeout) {
		this.ddi = ddi;
		this.ppi = ppi;
		
		this.number_of_threads = number_of_threads;
		this.max_depth_of_search = max_depth_of_search;
		this.pair_building_threshold = pair_building_threshold;
		this.prob_cutoff = probability_cutoff;
		this.verbose = verbose;
		this.compute_timeout = compute_timeout;
	}
	
	public DACO(ConstructedNetworks constr, int number_of_threads, int max_depth_of_search, double pair_building_threshold, double probability_cutoff, PrintStream verbose) {
		this.ddi = constr.getDDIN();
		this.ppi = constr.getPPIN();
		
		this.number_of_threads = number_of_threads;
		this.max_depth_of_search = max_depth_of_search;
		this.pair_building_threshold = pair_building_threshold;
		this.prob_cutoff = probability_cutoff;
		this.verbose = verbose;
	}
	
	/**
	 * Sets the general parameters
	 * @param number_of_threads
	 * @param max_depth_of_search
	 * @param pair_building_threshold
	 * @param probability_cutoff
	 * @param verbose
	 * @param compute_timeout
	 */
	public void setComputationParameters(int number_of_threads, int max_depth_of_search, double pair_building_threshold, double probability_cutoff, PrintStream verbose, int compute_timeout) {
		this.number_of_threads = number_of_threads;
		this.max_depth_of_search = max_depth_of_search;
		this.pair_building_threshold = pair_building_threshold;
		this.prob_cutoff = probability_cutoff;
		this.verbose = verbose;
		this.compute_timeout = compute_timeout;
	}
	
	/**
	 * StepFunction implements the parallel-part of the algorithm
	 */
	private final class StepFunction implements Runnable {
		
		final HashSet<String> internal_proteins;
		final HashSet<StrPair> domain_interactions; 
		final double c_w_in;
		final double c_w_out;
		final double P;
		
		public StepFunction(HashSet<String> internal_proteins, HashSet<StrPair> domain_interactions, double c_w_in, double c_w_out, double p) {
			this.internal_proteins = internal_proteins;
			this.domain_interactions = domain_interactions;
			this.c_w_in = c_w_in;
			this.c_w_out = c_w_out;
			P = p;
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime
					* result
					+ ((domain_interactions == null) ? 0 : domain_interactions
							.hashCode());
			result = prime
					* result
					+ ((internal_proteins == null) ? 0 : internal_proteins
							.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			StepFunction other = (StepFunction) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (domain_interactions == null) {
				if (other.domain_interactions != null)
					return false;
			} else if (!domain_interactions.equals(other.domain_interactions))
				return false;
			if (internal_proteins == null) {
				if (other.internal_proteins != null)
					return false;
			} else if (!internal_proteins.equals(other.internal_proteins))
				return false;
			return true;
		}
		
		/**
		 * actual step
		 */
		@Override
		public void run() {
			final double c_sum = c_w_in + c_w_out;
			double max_coh = c_w_in / c_sum;
			double n_w_in = 0.0;
			double n_w_out = 0.0;
			
			// determine occupied domains from edges
			HashSet<String> occupied_domains = new HashSet<String>();
			for (StrPair edge : domain_interactions) {
				occupied_domains.add(edge.getL());
				occupied_domains.add(edge.getR());
			}
			
			// search for additional protein that may max. cohesiveness
			HashMap<String, LinkedList<StrPair>> incident = getIncidentNodes(internal_proteins, occupied_domains);
			String add_max_protein = null;
			for (String protein : incident.keySet()) {
				final double[] coh = ppi.computeDeltaCohesiveness(protein, internal_proteins);
				final double temp_coh = (c_w_in + coh[0]) / (c_sum + coh[1]);
				if (temp_coh > max_coh) {
					max_coh = temp_coh;
					add_max_protein = protein;
					n_w_in = coh[0];
					n_w_out = coh[1];
				}
			}
			
			// search for removable protein that could max. cohesiveness
			HashMap<String, StrPair> boundary = getBoundaryNodes(internal_proteins, domain_interactions);
			String del_max_protein = null;
			for (String protein : boundary.keySet()) {
				final double[] coh = ppi.computeDeltaCohesiveness(protein, internal_proteins);
				final double temp_coh = (c_w_in - coh[0]) / (c_sum - coh[1]);
				if (temp_coh > max_coh) {
					max_coh = temp_coh;
					del_max_protein = protein;
					n_w_in = coh[0];
					n_w_out = coh[1];
				}
			}
			
			// choose coh-max. move
			if (del_max_protein == null) {
				if (add_max_protein == null) {
					// "return" result
					temp_results.add(internal_proteins);
					return;
				} else {
					// add case
					// filter alternatives
					final LinkedList<StrPair> domain_edge_alternatives = filterDomainInteractionAlternatives(incident.get(add_max_protein));
					
					if (domain_edge_alternatives.size() == 1) {
						// single addition
						final StrPair distinct_domain_edge = domain_edge_alternatives.get(0);
						final double P_n = P * ppi.getWeights().get(new StrPair(ddi.getDomain_to_protein().get(distinct_domain_edge.getL()), ddi.getDomain_to_protein().get(distinct_domain_edge.getR())));
						
						if (P_n >= prob_cutoff) {
							// add
							final HashSet<String> new_internal_proteins = new HashSet<String>(internal_proteins);
							final HashSet<StrPair> new_domain_interactions = new HashSet<StrPair>(domain_interactions);
							new_internal_proteins.add(add_max_protein);
							
							// cutoff, don't recurse
							if (new_internal_proteins.size() == max_depth_of_search) {
								temp_results.add(new_internal_proteins);
								return;
							}
							new_domain_interactions.add(distinct_domain_edge);
							
							// recurse,  necessary to make new objects
							try {
								pool.execute(new StepFunction(new_internal_proteins, new_domain_interactions, c_w_in + n_w_in, c_w_out + n_w_out - n_w_in, P_n));
							} catch (RejectedExecutionException e) {
								temp_results.add(new_internal_proteins);
							}
						} else {
							// "return" result
							temp_results.add(internal_proteins);
						}
						return;
						
					} else {
						// many possibilities to add the same protein
						final double d_in = c_w_in + n_w_in;
						final double d_out = c_w_out + n_w_out - n_w_in;
						boolean any_addition = false;
						
						final HashSet<String> new_internal_proteins = new HashSet<String>(internal_proteins);
						new_internal_proteins.add(add_max_protein);
						
						if (pool.isShutdown()) {
							temp_results.add(new_internal_proteins);
							return;
						}
						
						try {
							for (StrPair domain_edge : domain_edge_alternatives) {
								final double P_n = P * ppi.getWeights().get(new StrPair(ddi.getDomain_to_protein().get(domain_edge.getL()), ddi.getDomain_to_protein().get(domain_edge.getR())));
								if (P_n >= prob_cutoff) {
									// add
									any_addition = true;
									// cutoff, don't recurse
									if (new_internal_proteins.size() == max_depth_of_search) {
										temp_results.add(new_internal_proteins);
										return;// ok if only 1 choice is above threshold...
									}
									final HashSet<StrPair> new_domain_interactions = new HashSet<StrPair>(domain_interactions);
									new_domain_interactions.add(domain_edge);
									
									// recurse, only necessary to make new domain-IA objects
									pool.execute(new StepFunction(new_internal_proteins, new_domain_interactions, d_in, d_out, P_n));	
								} 
							}
						} catch (RejectedExecutionException e) {
							temp_results.add(new_internal_proteins);
							return;
						} 
						
						
						
						// if no extension possible, return current
						if (!any_addition)
							temp_results.add(internal_proteins);
						return;
					}
				}
			} else {
				// remove case
				final StrPair distinct_domain_edge = boundary.get(del_max_protein);
				final double P_n = P / ppi.getWeights().get(new StrPair(ddi.getDomain_to_protein().get(distinct_domain_edge.getL()), ddi.getDomain_to_protein().get(distinct_domain_edge.getR())));
				final HashSet<String> new_internal_proteins = new HashSet<String>(internal_proteins);
				final HashSet<StrPair> new_domain_interactions = new HashSet<StrPair>(domain_interactions);
				new_internal_proteins.remove(del_max_protein);
				new_domain_interactions.remove(distinct_domain_edge);
				
				// recurse, make new objects
				try {
					pool.execute(new StepFunction(new_internal_proteins, new_domain_interactions, c_w_in - n_w_in, c_w_out - n_w_out + n_w_in, P_n));
				} catch (RejectedExecutionException e) {
					temp_results.add(new_internal_proteins);
				}
					
				return;
			}
			
		}

		private DACO getOuterType() {
			return DACO.this;
		}
	}
	
	/**
	 * Implements a TimerTask that stops the execution of pathological subcases (if wanted) and returns the intermediate results instead
	 */
	private final class EarlyStopper extends TimerTask {

		@Override
		public void run() {
			
			if (verbose != null)
				verbose.println("Timeout of " + compute_timeout + " min exceeded: forced shutdown while collecting intermediate results.");
			
			// stops algorithm and retrieves remaining incomplete iterations
			List<Runnable> intermediates = pool.shutdownNow();
			for (Runnable r:intermediates) {
				StepFunction step = (StepFunction) r;
				temp_results.add(step.internal_proteins);
			}
		}
		
	}
	
	/**
	 * Starts the pair growing and DACO algorithm for a given seed node
	 * @param protein
	 * @param other_seed_proteins
	 * @return
	 */
	public HashSet<HashSet<String>> growPairs(String protein, Set<String> other_seed_proteins) {
		
		// re-init necessary stuff
		temp_results = Collections.newSetFromMap(new ConcurrentHashMap<HashSet<String>, Boolean>(Math.min(number_of_threads / 2, 2)));
		pool = new ThreadPoolExecutor(number_of_threads, number_of_threads, 5, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		
		// check if in PPI
		if (!ppi.contains(protein)) {
			if (verbose != null)
				verbose.println(protein+" not in PPIN.");
			return new HashSet<HashSet<String>>();
		}
		
		// determine start-pairs and DDIs
		HashSet<String> seed = new HashSet<String>();
		seed.add(protein);
		HashMap<String, LinkedList<StrPair>> incident = getIncidentNodes(seed, new HashSet<String>());
		LinkedList<StepFunction> jobs = new LinkedList<DACO.StepFunction>();
		
		// get all DDIs and proteins
		int no_pair_proteins = 0;
		for (String other_protein : incident.keySet()) {
			// if known that more will be processed, don't start symmetric
			if (other_seed_proteins.contains(other_protein) && protein.compareTo(other_protein) < 0)
				continue;
			double P = ppi.getWeights().get(new StrPair(protein, other_protein));
			if (P <= pair_building_threshold)
				continue;
			no_pair_proteins++;
			LinkedList<StrPair> relevant_ddis = filterDomainInteractionAlternatives(incident.get(other_protein));
			final HashSet<String> temp_set = new HashSet<String>(seed);
			temp_set.add(other_protein);
			double[] coh = ppi.computeClusterCohesiveness(temp_set);
			for (StrPair ddi_choice : relevant_ddis) {
				final HashSet<StrPair> new_domain_interactions = new HashSet<StrPair>();
				new_domain_interactions.add(ddi_choice);
				jobs.add(new StepFunction(temp_set, new_domain_interactions, coh[0], coh[1], P));
			}
		}
		// end job building
		
		// old Python implementation tried to have at least 2 pairs at this point, JDACO doesn't do that
		
		
		// start computing and output
		Timer timer = new Timer();
		
		if (verbose != null)
			verbose.println("Seeding from "+no_pair_proteins+" pairs with "+jobs.size()+ " unique domain interactions.");
		
		if (this.compute_timeout > 0)
			timer.schedule(new EarlyStopper(), this.compute_timeout * 60 * 1000);
		
		for (StepFunction job : jobs)
			pool.execute(job);
		do {
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		while (pool.getQueue().size() > 0 || pool.getActiveCount() > 0);
		
		timer.cancel();
		pool.shutdown();
		
		Runtime.getRuntime().gc();
		
		return new HashSet<HashSet<String>>(temp_results);
	}
	
	/**
	 * Starts the pair growing and DACO algorithm for all nodes automatically,
	 * merges all results and returns the final set of sets in a special results class.
	 * @param seed_proteins
	 * @return
	 */
	public DACOResultSet batchGrowPairs(Set<String> seed_proteins) {
		
		HashSet<HashSet<String>> results = new HashSet<HashSet<String>>();
		
		for (String tf : seed_proteins) {
			results.addAll(this.growPairs(tf, seed_proteins));
		}
		
		DACO.filterOutput(results);
		
		return new DACOResultSet(results, seed_proteins);
	}
	
	/**
	 * Returns boundary nodes and
	 * @param internal_proteins
	 * @param domain_interactions
	 * @return
	 */
	private HashMap<String, StrPair> getBoundaryNodes(HashSet<String> internal_proteins, HashSet<StrPair> domain_interactions) {
		HashMap<String, StrPair> boundary_nodes = new HashMap<String, StrPair>();
		HashMap<String, LinkedList<StrPair>> associated_dom_interactions = new HashMap<String, LinkedList<StrPair>>();
		
		for (StrPair domain_interaction : domain_interactions) {
			String d1 = domain_interaction.getL();
			String d2 = domain_interaction.getR();
			String p1 = ddi.getDomain_to_protein().get(d1);
			String p2 = ddi.getDomain_to_protein().get(d2);
			if (!associated_dom_interactions.containsKey(p1))
				associated_dom_interactions.put(p1, new LinkedList<StrPair>());
			associated_dom_interactions.get(p1).add(domain_interaction);
			if (!associated_dom_interactions.containsKey(p2))
				associated_dom_interactions.put(p2, new LinkedList<StrPair>());
			associated_dom_interactions.get(p2).add(domain_interaction);
		}
		for (String protein : associated_dom_interactions.keySet()) {
			List<StrPair> ddi_options = associated_dom_interactions.get(protein);
			if (ddi_options.size() == 1)
				boundary_nodes.put(protein, ddi_options.get(0));
		}
		associated_dom_interactions = null;
		
		return boundary_nodes;
	}
	
	/**
	 * Filters for types of interactions that occur several times
	 * @param usable_interactions
	 * @return
	 */
	private LinkedList<StrPair> filterDomainInteractionAlternatives(LinkedList<StrPair> usable_interactions) {
		HashSet<String> already_seen = new HashSet<String>();
		LinkedList<StrPair> filtered_list = new LinkedList<StrPair>();
		for (StrPair pair : usable_interactions) {
			// reverse ordering
			String[] sp1 = pair.getL().split("\\|");
			String[] sp2 = pair.getR().split("\\|");
			String s1 = sp1[2]+sp1[1];
			String s2 = sp2[2]+sp2[1];
			String hash;
			if (s1.compareTo(s2) < 0)
				hash = s1+s2;
			else
				hash = s2+s1;
			if (already_seen.contains(hash)) 
				continue;
			already_seen.add(hash);
			filtered_list.add(pair);
		}
		return filtered_list;
	}
	
	/**
	 * Returns incident proteins and how they can be added on the domain level
	 * @param internal_proteins
	 * @param occupied_domains
	 * @return
	 */
	private HashMap<String, LinkedList<StrPair>> getIncidentNodes(HashSet<String> internal_proteins, HashSet<String> occupied_domains) {
		HashMap<String, LinkedList<StrPair>> incident_nodes = new HashMap<String, LinkedList<StrPair>>();
		for (String protein : internal_proteins)
			for (String domain1: ddi.getProtein_to_domains().get(protein)) {
				if (occupied_domains.contains(domain1)) continue;
				for (String domain2: ddi.getDDIs().get(domain1)) {
					if (occupied_domains.contains(domain2)) continue;
					String protein2 = ddi.getDomain_to_protein().get(domain2);
					if (internal_proteins.contains(protein2)) continue;
					
					// add option
					if (!incident_nodes.containsKey(protein2))
						incident_nodes.put(protein2, new LinkedList<StrPair>());
					incident_nodes.get(protein2).add(new StrPair(domain1, domain2));
				}
			}
			
		return incident_nodes;
	}
	
	/**
	 * Merges and writes output to a csv-file
	 * @param out_file
	 * @param results
	 */
	public static void writeAndFilterOutput(String out_file, HashSet<HashSet<String>> results) {
		// clean from subsets
		DACO.filterOutput(results);
		
		List<String> to_write = new LinkedList<>();

		for (HashSet<String> cluster : results) 
			to_write.add( String.join(",", cluster) );
		
		Utilities.writeEntries(to_write, out_file);
	}
	
	/**
	 * Merges proper subsets in results (in-memory!)
	 * @param results
	 */
	public static void filterOutput(HashSet<HashSet<String>> results) {
		// clean from subsets
		HashSet<HashSet<String>> sub = new HashSet<HashSet<String>>();
		for (HashSet<String> i : results)
			for (HashSet<String> j : results) {
				if (i.equals(j))
					continue;
				if (j.containsAll(i))
					sub.add(i);
			}
		results.removeAll(sub);
	}
}
