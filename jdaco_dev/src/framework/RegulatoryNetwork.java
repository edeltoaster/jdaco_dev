package framework;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

/**
 * Class to handle regulatory networks
 * @author Thorsten Will
 */
public class RegulatoryNetwork {
	private final DACOResultSet TF_complexes;
	private final BindingDataHandler bdh;
	private int d_min = -50;
	private int d_max = 50;
	private int no_threads = Runtime.getRuntime().availableProcessors()/2;
	private int min_TFs = 1;
	
	// only stores TF in complex, not all partners
	private final Map<HashSet<String>, Set<String>> complex_to_targets = new HashMap<>();
	private final Map<String, List<HashSet<String>>> tf_to_complex = new HashMap<>();
	
	public RegulatoryNetwork(HashSet<HashSet<String>> TF_complexes, BindingDataHandler bdh, int min_TFs) {
		this.TF_complexes = new DACOResultSet(TF_complexes, bdh.getTFsWithBindingData());
		this.bdh = bdh;
		this.min_TFs = min_TFs;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(Set<HashSet<String>> TF_complexes, BindingDataHandler bdh, int d_min, int d_max, int no_threads, int min_TFs) {
		this.TF_complexes = new DACOResultSet(TF_complexes, bdh.getTFsWithBindingData());
		this.bdh = bdh;
		this.d_min = d_min;
		this.d_max = d_max;
		this.no_threads = no_threads;
		this.min_TFs = min_TFs;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(List<HashSet<String>> TF_complexes, BindingDataHandler bdh, int d_min, int d_max, int no_threads, int min_TFs) {
		this.TF_complexes = new DACOResultSet(new HashSet<HashSet<String>>(TF_complexes), bdh.getTFsWithBindingData());
		this.bdh = bdh;
		this.d_min = d_min;
		this.d_max = d_max;
		this.no_threads = no_threads;
		this.min_TFs = min_TFs;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(DACOResultSet complexes, BindingDataHandler bdh, int d_min, int d_max, int no_threads, int min_TFs) {
		this.TF_complexes = complexes;
		this.bdh = bdh;
		this.d_min = d_min;
		this.d_max = d_max;
		this.no_threads = no_threads;
		this.min_TFs = min_TFs;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(DACOResultSet complexes, BindingDataHandler bdh) {
		this.TF_complexes = complexes;
		this.bdh = bdh;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(Set<HashSet<String>> TF_complexes, BindingDataHandler bdh) {
		this.TF_complexes = new DACOResultSet(TF_complexes, bdh.getTFsWithBindingData());
		this.bdh = bdh;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(Set<HashSet<String>> TF_complexes, BindingDataHandler bdh, int d_min, int d_max, int no_threads) {
		this.TF_complexes = new DACOResultSet(TF_complexes, bdh.getTFsWithBindingData());
		this.bdh = bdh;
		this.d_min = d_min;
		this.d_max = d_max;
		this.no_threads = no_threads;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(DACOResultSet complexes, BindingDataHandler bdh, int d_min, int d_max, int no_threads) {
		this.TF_complexes = complexes;
		this.bdh = bdh;
		this.d_min = d_min;
		this.d_max = d_max;
		this.no_threads = no_threads;
		this.constructNetwork();
	}
	
	public RegulatoryNetwork(DACOResultSet complexes, BindingDataHandler bdh, int min_TFs) {
		this.TF_complexes = complexes;
		this.bdh = bdh;
		this.min_TFs = min_TFs;
		this.constructNetwork();
	}
	
	/**
	 * Constructs complex_to_targets and tf_to_complex from a written RegNet file, all other fields have default initialization.
	 * @param regnet_file
	 */
	public RegulatoryNetwork(String regnet_file) {
		this.TF_complexes = null;
		this.bdh = null;
		this.readRegulatoryNetwork(regnet_file);
	}
	
	/**
	 * Constructs actual TF complex to target map, only considers complexes with at least min_TFs TFs
	 * @param min_TFs
	 */
	private void constructNetwork() {
		
		// speed-up building process
		ExecutorService threadpool = Executors.newFixedThreadPool(this.no_threads);
		CompletionService<HashMap<HashSet<String>, Set<String>>> pool = new ExecutorCompletionService<>(threadpool);
		
		int i = 0;
		for (HashSet<String> tfs_in_complex:this.TF_complexes.getSeedToComplexMap().keySet()) {
			if (tfs_in_complex.size() < min_TFs)
				continue;
			pool.submit(new ConstructHelper(tfs_in_complex));
			i++;
		}
		int n = i;
		try {
			for (i = 0; i < n; i++) {
				HashMap<HashSet<String>, Set<String>> target_map;
				target_map = pool.take().get(); // always gets a finished job
				for (HashSet<String> tfs_in_complex:target_map.keySet()) {
					Set<String> targets = target_map.get(tfs_in_complex);
					if (targets.size() > 0) {
						this.complex_to_targets.put(tfs_in_complex, targets);
						for (String tf:tfs_in_complex) {
							if (!this.tf_to_complex.containsKey(tf))
								this.tf_to_complex.put(tf, new LinkedList<HashSet<String>>());
							this.tf_to_complex.get(tf).add(tfs_in_complex);
						}
					}
				}
			
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		threadpool.shutdown();
	}
	
	/**
	 * Constructs complex_to_targets and tf_to_complex from a written RegNet file
	 * @param regnet_file
	 */
	private void readRegulatoryNetwork(String regnet_file) {
		
		for (String line:Utilities.readFile(regnet_file)) {
			if (line.startsWith("TF(complex)")) // format: TF(complex) target edgetype
				continue;
			String[] split = line.trim().split("\\s+");
			String edgetype = split[2];
			String tf_compl = split[0];
			String target = split[1];
			
			if (edgetype.equals("TFC/target")) {
				HashSet<String> complex = new HashSet<String>(Arrays.asList(tf_compl.split("/")));
				if (complex.size() == 1 && !this.tf_to_complex.containsKey(tf_compl)) {
					this.tf_to_complex.put(tf_compl, new LinkedList<HashSet<String>>());
					this.tf_to_complex.get(tf_compl).add(complex); // for consistency
				}
				if (!this.complex_to_targets.containsKey(complex))
					this.complex_to_targets.put(complex, new HashSet<String>());
				this.complex_to_targets.get(complex).add(target);
			} else { // type TF/TFC
				HashSet<String> complex = new HashSet<String>(Arrays.asList(target.split("/")));
				if (!this.tf_to_complex.containsKey(tf_compl))
					this.tf_to_complex.put(tf_compl, new LinkedList<HashSet<String>>());
				this.tf_to_complex.get(tf_compl).add(complex);
			}
		}
		
	}
	
	/**
	 * Writes the TF complex -> target-association to a text file, only writes associations with at least write_min_TFs TFs in the complex.
	 * @param out_file
	 * @param write_min_TFs
	 */
	public void writeRegulatoryNetwork(String out_file, int write_min_TFs) {
		List<String> to_write = new LinkedList<>();
		to_write.add("TF(complex) target edgetype");

		// write complex -> target
		String tfc_target = " TFC/target";
		String tf_tfc = " TF/TFC";
		for (HashSet<String> tfs_in_complex:this.complex_to_targets.keySet()) {
			if (tfs_in_complex.size() < write_min_TFs)
				continue;
			
			String complex = String.join("/", tfs_in_complex);
			for (String target:this.complex_to_targets.get(tfs_in_complex)) {
				to_write.add(complex + " " + target + tfc_target);
			}
		}

		// write TF -> complex
		for (String tf:this.tf_to_complex.keySet())
			for (HashSet<String> complex:this.tf_to_complex.get(tf)) {
				// if complex involved only 1 TF this would result in a self-interaction
				if (complex.size() == 1 || complex.size() < write_min_TFs)
					continue;
				to_write.add(tf + " " + String.join("/", complex) + tf_tfc);
			}
		
		Utilities.writeEntries(to_write, out_file);

	}
	
	/**
	 * Writes the TF complex -> target-association to a text file
	 * @param out_file
	 */
	public void writeRegulatoryNetwork(String out_file) {
		writeRegulatoryNetwork(out_file, this.min_TFs);
	}
	
	/**
	 * Writes further information to all nodes in a text file
	 * @param out_file
	 */
	public void writeSimpleNodeTable(String out_file) {
		writeNodeTable(out_file, null, null);
	}
	
	/**
	 * Writes further information to all nodes in a text file, automatically retrieves gene names.
	 * @param out_file
	 */
	public void writeNodeTable(String out_file) {		
		writeNodeTable(out_file, DataQuery.getUniprotToGeneNameMap(tf_to_complex.keySet()), null);
	}
	
	/**
	 * Writes further information to all nodes in a text file, automatically retrieves gene names and can add additional custom information.
	 * @param out_file
	 * @param annotational_data
	 */
	public void writeNodeTable(String out_file, Map<String, Map<String, String>> annotational_data) {
		writeNodeTable(out_file, DataQuery.getUniprotToGeneNameMap(tf_to_complex.keySet()), annotational_data);
	}
	
	@Deprecated
	/**
	 * Writes further information to all nodes in a text file, adds HGNC identifiers and user-given boolean labels
	 * @param out_file
	 */
	public void writeHumanNodeTableOld(String out_file, Map<String, Set<String>> additional_data) {
		Set<String> seen_targets = new HashSet<>();
		List<String> to_write = new LinkedList<>();
		
		List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
		Map<String, String> up_hgnc = new HashMap<>();
		for (String[] s:HGNC_map) {
			up_hgnc.put(s[1], s[0]);
		}
		
		String header = "Node Nodetype HGNC";
		if (additional_data != null)
			for (String label:additional_data.keySet())
				header += " " + label;
		
		to_write.add(header);
		for (HashSet<String> tfs_in_complex:this.complex_to_targets.keySet()) {
			String complex = String.join("/", tfs_in_complex);
			String temp = "";
			boolean first = true;
			for (String s:tfs_in_complex) {
				if (first) {
					temp += up_hgnc.get(s);
					first = false;
				}
				else 
					temp += "/" + up_hgnc.get(s);
			}
			
			String data = complex + " complex " + temp;
			if (additional_data != null)
				for (String label:additional_data.keySet()) {
					Set<String> of_interest = new HashSet<String>(additional_data.get(label));
					of_interest.retainAll(tfs_in_complex);
					if (of_interest.size() > 0)
						data += " yes";
					else
						data += " no";
				}
			to_write.add(data);
			
			for (String target:this.complex_to_targets.get(tfs_in_complex)) {
				
				// don't do twice
				if (seen_targets.contains(target))
					continue;
				seen_targets.add(target);
				
				if (this.tf_to_complex.keySet().contains(target))
					data = target + " TF " + up_hgnc.get(target);
				else
					data = target + " protein " + up_hgnc.get(target);
				
				if (additional_data != null)
					for (String label:additional_data.keySet()) {
						if (additional_data.get(label).contains(target))
							data += " yes";
						else
							data += " no";
					}
				to_write.add(data);
			}
		}
		
		// ensure that all TFs were described, even if they are in no target set
		for (HashSet<String> tfs_in_complex:this.complex_to_targets.keySet())
			for (String tf:tfs_in_complex) {
				if (seen_targets.contains(tf))
					continue;
				seen_targets.add(tf);
				
				String data = tf + " TF " + up_hgnc.get(tf);
				
				if (additional_data != null)
					for (String label:additional_data.keySet()) {
						if (additional_data.get(label).contains(tf))
							data += " yes";
						else
							data += " no";
					}
				to_write.add(data);
			}
		
		Utilities.writeEntries(to_write, out_file);
	}
	
	@Deprecated
	/**
	 * Writes further information to all nodes in a text file, adds HGNC identifiers
	 * @param out_file
	 */
	public void writeHumanNodeTable(String out_file) {
		this.writeHumanNodeTableOld(out_file, null);
	}
	
	/**
	 * Writes further information to all nodes in a text file.
	 * @param out_file
	 * @param name_conversion
	 * @param annotational_data
	 */
	public void writeNodeTable(String out_file, Map<String, String> name_conversion, Map<String, Map<String, String>> annotational_data) {
		List<String> to_write = new LinkedList<>();
		
		// write header
		String header = "Node Nodetype";
		
		if (name_conversion != null)
			header += " Gene";
		
		if (annotational_data != null)
			for (String label:annotational_data.keySet())
				header += " " + label;
		
		to_write.add(header);
		String data = null;
		Set<String> seen_targets = new HashSet<>();
		for (HashSet<String> tfs_in_complex:this.complex_to_targets.keySet()) {
			String complex = String.join("/", tfs_in_complex);
			data = complex + " complex";
			
			// mark single TF complexes so that the are correctly annotated as complexes only rather than also as TFs
			if (tfs_in_complex.size() == 1)
				seen_targets.addAll(tfs_in_complex);
			
			// converts names, names protein by UniProt Acc if not in map
			if (name_conversion != null) {
				String genes = String.join("/", tfs_in_complex.stream().map(p->name_conversion.getOrDefault(p, p)).collect(Collectors.toList()));
				data += " " + genes;
			}
			
			// adds annotational data, annotates with / if there is no data for this TF combination
			if (annotational_data != null)
				for (String label:annotational_data.keySet()) {
					Map<String, String> annotation_map = annotational_data.get(label);
					data += " " + annotation_map.getOrDefault(tfs_in_complex.toString(), "/");
				}
			to_write.add(data);
		}
		
		// second iteration ensures that TFs that are also sole complexes are named correctly
		for (HashSet<String> tfs_in_complex:this.complex_to_targets.keySet())
			for (String target:this.complex_to_targets.get(tfs_in_complex)) {
				
				// don't do twice
				if (seen_targets.contains(target))
					continue;
				seen_targets.add(target);
				
				if (this.tf_to_complex.keySet().contains(target))
					data = target + " TF";
				else
					data = target + " protein";
				
				// converts name, names protein by UniProt Acc if not in map
				if (name_conversion != null) {
					String gene = name_conversion.getOrDefault(target, target);
					data += " " + gene;
				}
				
				// adds annotational data, annotates with / if there is no data for this target
				if (annotational_data != null)
					for (String label:annotational_data.keySet()) {
						Map<String, String> annotation_map = annotational_data.get(label);
						data += " " + annotation_map.getOrDefault(target, "/");
					}
				to_write.add(data);
			}
		
		// ensure that all TFs were described, even if they are in no target set
		for (HashSet<String> tfs_in_complex:this.complex_to_targets.keySet())
			for (String tf:tfs_in_complex) {
				if (seen_targets.contains(tf))
					continue;
				seen_targets.add(tf);
				
				data = tf + " TF";
				
				// converts name, names protein by UniProt Acc if not in map
				if (name_conversion != null) {
					String gene = name_conversion.getOrDefault(tf, tf);
					data += " " + gene;
				}
				
				// adds annotational data, annotates with / if there is no data for this TF
				if (annotational_data != null)
					for (String label:annotational_data.keySet()) {
						Map<String, String> annotation_map = annotational_data.get(label);
						data += " " + annotation_map.getOrDefault(tf, "/");
					}
				to_write.add(data);
			}
		
		Utilities.writeEntries(to_write, out_file);
	}
	
	/*
	 * graph-theoretic algorithms
	 */
	
	/**
	 * Prunes regulatory network to largest strongly connected component 
	 * and iteratively removes those parts of the network where not all TFs needed for a complex are within the network.
	 */
	public void pruneToLargestSCCs() {
		
		if (this.complex_to_targets.size() == 0)
			return;
		
		// convert to general graph format
		Map<String, List<String>> graph = new HashMap<>();
		for (HashSet<String> complex:this.complex_to_targets.keySet())
			graph.put("C:" + String.join("/", complex), new ArrayList<>(this.complex_to_targets.get(complex)));
		
		for (String tf:this.tf_to_complex.keySet())
			graph.put(tf, this.tf_to_complex.get(tf).stream().map( c -> "C:" + String.join("/", c) ).collect(Collectors.toList()));
		
		List<Set<String>> SCCs = Utilities.getSCCs(graph);
		
		Set<String> largest_SCC = null;
		int largest_size = 0;
		for (Set<String> SCC:SCCs) {
			int current_size = SCC.size();
			if (current_size > largest_size) {
				largest_size = current_size;
				largest_SCC = SCC;
			}
		}
		
		// convert back to TFs, targets, complexes
		Set<String> allowed_proteins = new HashSet<>();
		Set<HashSet<String>> allowed_complexes = new HashSet<>();
		for (String node:largest_SCC) {
			// complex
			if (node.startsWith("C:")) {
				allowed_complexes.add(new HashSet<>(Arrays.asList(node.substring(2).split("/"))));
			} else {
				allowed_proteins.add(node);
			}
		}
		
		this.pruneToConsistency(allowed_proteins, allowed_complexes);
	}
	
	/**
	 * Iteratively prune to ensure all TF -> complex links exist. Input sets are altered for consistency.
	 * @param allowed_proteins
	 * @param allowed_complexes
	 */
	private void pruneToConsistency(Set<String> allowed_proteins, Set<HashSet<String>> allowed_complexes) {
		
		boolean prune = false;
		int allowed_complexes_before = Integer.MAX_VALUE;
		int allowed_proteins_before = Integer.MAX_VALUE;
		
		do {
			prune = false;
			allowed_complexes_before = allowed_complexes.size();
			allowed_proteins_before = allowed_proteins.size();
			
			// remove complexes including TFs that are not allowed
			for (HashSet<String> complex:this.complex_to_targets.keySet())
				// if not all TFs are allowed, the complex cannot exist
				if (complex.stream().anyMatch(p -> !allowed_proteins.contains(p)))
					allowed_complexes.remove(complex);
			
			// prune complexes
			this.complex_to_targets.keySet().retainAll(allowed_complexes);
			// shrink to target sets that are allowed
			this.complex_to_targets.values().forEach(targets -> targets.retainAll(allowed_proteins));
			// remove empty target sets
			this.complex_to_targets.entrySet().removeIf(e -> e.getValue().size() == 0);
			
			// refine allowed complexes and proteins based on target set-pruning
			allowed_complexes.retainAll(this.complex_to_targets.keySet());
			allowed_proteins.retainAll(Utilities.getValueSetFromSetMultimap(this.complex_to_targets));
			
			// prune TFs
			this.tf_to_complex.keySet().retainAll(allowed_proteins);
			// prune complexes that are targeted
			this.tf_to_complex.values().forEach(targets -> targets.retainAll(allowed_complexes));
			// remove empty target sets
			this.tf_to_complex.entrySet().removeIf(e -> e.getValue().size() == 0);
			
			allowed_proteins.retainAll(this.tf_to_complex.keySet());
			allowed_complexes.retainAll(Utilities.getValueSetFromMultimap(this.tf_to_complex));
			
			if (allowed_complexes.size() < allowed_complexes_before || allowed_proteins.size() < allowed_proteins_before)
				prune = true;
			
			// if there was pruning, repeat once more
		} while (prune);
	}
	
	/**
	 * Remove given protein set from regulatory network
	 * @param proteins_to_prune
	 */
	public void removeProteinSet(Set<String> proteins_to_prune) {
		// initialize with current network
		Set<String> allowed_proteins = new HashSet<>(this.tf_to_complex.keySet());
		this.complex_to_targets.values().forEach(targets -> allowed_proteins.addAll(targets));
		Set<HashSet<String>> allowed_complexes = new HashSet<>(this.complex_to_targets.keySet());
		
		// remove stuff that should be pruned
		allowed_proteins.removeAll(proteins_to_prune);
		for (HashSet<String> complex:this.complex_to_targets.keySet())
			// if not all TFs are allowed, the complex cannot exist
			if (complex.stream().anyMatch(p -> !allowed_proteins.contains(p)))
				allowed_complexes.remove(complex);
		
		this.pruneToConsistency(allowed_proteins, allowed_complexes);
	}
	
	/*
	 * getters
	 */
	
	/**
	 * Returns [#complexes, #TFs, #compl_target_IAs, #TF_compl_IAs].
	 * Note that complexes consisting of single TFs are also counted, even if they do not appear in the network visualization.
	 * @return
	 */
	public int[] getSizes() {
		int[] sizes = new int[4];
		sizes[0] = this.complex_to_targets.size();
		sizes[1] = this.tf_to_complex.size();
		sizes[2] = this.complex_to_targets.values().stream().mapToInt(t -> t.size()).sum();
		sizes[3] = this.tf_to_complex.values().stream().mapToInt(c -> c.size()).sum();
		return sizes;
	}
	
	/**
	 * Returns "#complexes, #TFs, #compl_target_IAs, #TF_compl_IAs".
 	 * Note that complexes consisting of single TFs are also counted, even if they do not appear in the network visualization.
	 * @return
	 */
	public String getSizesStr() {
		int[] sizes = this.getSizes();
		return sizes[0] + " complexes / " + sizes[1] + " TFs / " + sizes[2] + " complex_target_IAs / " + sizes[3] + " TF_complex_IAs";
	}
	
	/**
	 * Returns all proteins that are targeted by regulation but are no TF themselves and are targeted by complex with at least sink_min_TFs
	 * @param sink_min_TFs
	 */
	public Set<String> getSinkProteins(int sink_min_TFs) {
		Set<String> tfs = this.bdh.getTFsWithBindingData();
		Set<String> sink_proteins = new HashSet<String>();
		for (HashSet<String> TF_complex:this.complex_to_targets.keySet()) {
			if (TF_complex.size() < sink_min_TFs)
				continue;
			for (String target:this.complex_to_targets.get(TF_complex)) {
				if (tfs.contains(target))
					continue;
				sink_proteins.add(target);
			}
		}
		return sink_proteins;
	}
	
	/**
	 * Returns all proteins that are targeted by regulation but are no TF themselves
	 * @return
	 */
	public Set<String> getSinkProteins() {
		return getSinkProteins(this.min_TFs);
	}
	
	public DACOResultSet getTFComplexes() {
		return TF_complexes;
	}

	/**
	 * Returns all proteins somehow included in the network
	 * @return
	 */
	public Set<String> getIncludedProteins() {
		Set<String> proteins = new HashSet<String>(this.tf_to_complex.keySet());
		for (Set<String> temp:this.complex_to_targets.values())
			proteins.addAll(temp);
		return proteins;
	}
	
	public BindingDataHandler getBindingDataHandler() {
		return bdh;
	}

	public int getDmin() {
		return d_min;
	}

	public int getDmax() {
		return d_max;
	}
	
	public int getMinTFs() {
		return min_TFs;
	}
	public Map<HashSet<String>, Set<String>> getComplexToTargets() {
		return complex_to_targets;
	}

	public Map<String, List<HashSet<String>>> getTFToComplex() {
		return tf_to_complex;
	}
	
	/**
	 * Returns the #TFs in the biggest TF complex with common targets
	 * @return
	 */
	public int getBiggestTFComplexWithTargets() {
		int max = 0;
		for (HashSet<String> tf_complex:this.complex_to_targets.keySet()) {
			if (tf_complex.size() > max)
				max = tf_complex.size();
		}
		return max;
	}
	
	/**
	 * Helper to speedup building
	 */
	private final class ConstructHelper implements Callable<HashMap<HashSet<String>, Set<String>>>{
		
		private final HashSet<String> tfs_in_complex;

		ConstructHelper(HashSet<String> tfs_in_complex) {
			this.tfs_in_complex = tfs_in_complex;
		}
		
		@Override
		public HashMap<HashSet<String>, Set<String>> call() throws Exception {
			HashMap<HashSet<String>, Set<String>> temp = new HashMap<HashSet<String>, Set<String>>();
			temp.put(tfs_in_complex, bdh.getAdjacencyPossibilities(tfs_in_complex, d_min, d_max, false));
			return temp;
		}
		
	}
}
