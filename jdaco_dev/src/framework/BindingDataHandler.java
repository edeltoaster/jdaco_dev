package framework;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

/**
 * Tools for TF->target-gene handling
 * @author Thorsten Will
 */
public class BindingDataHandler {
	private final HashMap<String, HashMap<String, LinkedList<BindingSite>>> TF_to_targets = new HashMap<>(1024);
	private final Map<String, Set<String>> target_to_TFs = new HashMap<>(8096);
	
	/**
	 * Unrestricted constructor:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there).
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile) {
		this.readFIMO(fimo_outputfile, null, 1.0, null, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	/**
	 * Constructor with p-value threshold:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there).
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile, double p_threshold) {
		this.readFIMO(fimo_outputfile, null, p_threshold, null, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	/**
	 * Restricted constructor:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there)
	 * but only care for the TFs in a given set.
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile, Collection<String> expressed_TFs) {
		this.readFIMO(fimo_outputfile, expressed_TFs, 1.0, null, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	/**
	 * Constructor with p-value threshold:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there)
	 * but only care for the TFs in a given set.
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile, double p_threshold, Collection<String> expressed_TFs) {
		this.readFIMO(fimo_outputfile, expressed_TFs, p_threshold, null, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	/**
	 * Restricted constructor:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there),
	 * but only care for the TFs in a given set and the acc. targets given in another.
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile, Collection<String> expressed_TFs, Collection<String> acc_targets) {
		this.readFIMO(fimo_outputfile, expressed_TFs, 1.0, acc_targets, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	/**
	 * Restricted constructor:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there),
	 * but only care for the TFs in a given set, BS below p_threshold and the acc. targets given in another.
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile, Collection<String> expressed_TFs, double p_threshold, Collection<String> acc_targets) {
		this.readFIMO(fimo_outputfile, expressed_TFs, p_threshold, acc_targets, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	/**
	 * Restricted constructor:
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there),
	 * but only care for the TFs in a given set, BS below p_threshold and the acc. targets given in another.
	 * Additionally, only allow binding sites in the range between min_start and max_end around the TSS, for example. 
	 * This strongly depends on the processing and annotation of the binding data.
	 * All proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	public BindingDataHandler(String fimo_outputfile, Collection<String> expressed_TFs, double p_threshold, Collection<String> acc_targets, int min_start, int max_end) {
		this.readFIMO(fimo_outputfile, expressed_TFs, p_threshold, acc_targets, min_start, max_end);
	}
	
	/**
	 * Read binding data from MEME-suite's FIMO txt-file (gzipped also fine, ending .gz assumed there)
	 * all proteins (and associated genes) are assumed to be in UniProt accs.
	 */
	private void readFIMO(String fimo_outputfile, Collection<String> expressed_TFs, double p_threshold, Collection<String> acc_targets, int min_start, int max_end) {
		
		BufferedReader in = null;
		try {
			if (fimo_outputfile.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fimo_outputfile))));
			else
				in = new BufferedReader(new FileReader(fimo_outputfile));
			while (in.ready()) {
				String line = in.readLine();
				
				if (line == null)
					break;
				
				// skip header and comments
				if (line.startsWith("#") || line.startsWith("motif_id"))
					continue;
				
				// parse
				String[] data = line.split("\\s+");
				String tf = data[0];
				
				// skip if TFs are restricted and TF is not in the list of relevant TFs
				if (expressed_TFs != null && !expressed_TFs.contains(tf))
					continue;
				
				String target = data[1];
				
				// skip if target not accessible
				if (acc_targets != null && !acc_targets.contains(target))
					continue;
				
				// p-value check
				double p = Double.parseDouble(data[6]);
				if (p > p_threshold)
					continue;
				
				int left = Integer.parseInt(data[2]);
				int right = Integer.parseInt(data[3]);
				
				if (left < min_start || right > max_end)
					continue;
				
				boolean plus_strand = true;
				if (data[4].equals("-"))
					plus_strand = false;
				BindingSite bs = new BindingSite(left, right, plus_strand);
				
				// TF to target
				if (!TF_to_targets.containsKey(tf))
					TF_to_targets.put(tf, new HashMap<>());
				if (!TF_to_targets.get(tf).containsKey(target))
					TF_to_targets.get(tf).put(target, new LinkedList<>());
				TF_to_targets.get(tf).get(target).add(bs);
				
				// target to TF
				if (!target_to_TFs.containsKey(target))
					target_to_TFs.put(target, new HashSet<>());
				target_to_TFs.get(target).add(tf);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
			
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		// sort binding data -> ensures every list is sorted
		for (HashMap<String, LinkedList<BindingSite>> second_map:this.TF_to_targets.values())
			for (LinkedList<BindingSite> individual_list:second_map.values())
				individual_list.sort(null);
	}
	
	/**
	 * Returns the set of TF with binding annotation
	 */
	public Set<String> getTFsWithBindingData() {
		return new HashSet<>(this.TF_to_targets.keySet());
	}
	
	/**
	 * Returns the set of targets the TF has (don't modify)
	 * @param TF
	 * @return
	 */
	public Set<String> getTargets(String TF) {
		if (!this.TF_to_targets.containsKey(TF))
			return new HashSet<>();
		
		return this.TF_to_targets.get(TF).keySet();
	}
	
	/**
	 * Returns the set of TF that regulate the given target (don't modify)
	 * @param target
	 * @return
	 */
	public Set<String> getRegulatingTFs(String target) {
		return this.target_to_TFs.getOrDefault(target, new HashSet<>());
	}
	
	public HashMap<String, HashMap<String, LinkedList<BindingSite>>> getTFToTargetsPositionalMap() {
		return TF_to_targets;
	}
	
	public Map<String, Set<String>> getTargetsToTFsMap() {
		return target_to_TFs;
	}
	
	/**
	 * Writes the TF->target-association to a textfile
	 * @param out_file
	 */
	public void writeRegulatoryNetwork(String out_file) {
		List<String> to_write = new LinkedList<>();

		to_write.add("TF target");
		for (String TF:this.TF_to_targets.keySet())
			for (String target:this.getTargets(TF) )
				to_write.add(TF+" "+target);
			
		Utilities.writeEntries(to_write, out_file);
	}
	
	/**
	 * Returns the set of targets that is shared by several TFs;
	 * TFs that have no annotated BSs lead to complexes without common targets.
	 * @param TFs
	 * @return
	 */
	public Set<String> getCommonTargets(Collection<String> TFs) {
		
		Set<String> common_targets = new HashSet<>();
		
		if (TFs.size() == 0)
			return common_targets;
		
		boolean initialized = false;
		
		for (String TF:TFs) {
			
			// add everything from the first one
			if (!initialized) {
				common_targets.addAll(this.TF_to_targets.getOrDefault(TF, new HashMap<>()).keySet());
				initialized = true;
			}
			// retain only the overlap
			else {
				common_targets.retainAll(this.TF_to_targets.getOrDefault(TF, new HashMap<>()).keySet());
			}
			
			// stop early if necessary
			if (common_targets.size() == 0)
				return common_targets;
		}
		
		return common_targets;
	}
	
	/**
	 * Checks for given TFs and pairwise distance constraints if there are target genes with BSs that are compatible (does NOT return BSs)
	 * @param TFs
	 * @param d_min
	 * @param d_max
	 * @return
	 */
	public Set<String> getAdjacencyPossibilities(Collection<String> TFs, int d_min, int d_max, boolean allow_self_IAs) {
		
		List<String> tf_list = new ArrayList<>(TFs);
		Set<String> result = new HashSet<>();
		
		// catch special cases
		if (TFs.size() == 0)
			return result;
		else if (TFs.size() == 1) {
			for (String target:this.TF_to_targets.getOrDefault(tf_list.get(0), new HashMap<>()).keySet()) {
				result.add(target);
			}
			return result;
		} else if (TFs.size() == 2) {
			return this.getPairwiseAdjacent(tf_list.get(0), tf_list.get(1), d_min, d_max).keySet();
		}
		
		// precompute all pairwise, memorize targets than can really be found by pairwise start
		Set<String> common_targets = new HashSet<>();
		Map<StrPair, HashMap<String, LinkedList<BindingSite[]>>> prec_pw_adj = new HashMap<>();
		for (String tf1:TFs)
			for (String tf2:TFs) {
				if ( (tf1.compareTo(tf2) < 0) || (tf1.compareTo(tf2) == 0 && !allow_self_IAs) )
					continue;
				StrPair pair = new StrPair(tf1, tf2);
				HashMap<String, LinkedList<BindingSite[]>> temp = getPairwiseAdjacent(tf1, tf2, d_min, d_max);
				prec_pw_adj.put(pair, temp);
				common_targets.addAll(temp.keySet());
			}
		
		// check for all common/shared targets
		common_targets.retainAll(getCommonTargets(TFs));
		for (String target:common_targets) {
			// simply stores if there was a valid combination found
			boolean found = false;
			// for all pairwise starting positions
			for (String tf1:TFs) {
				for (String tf2:TFs) {
					if ( (tf1.compareTo(tf2) < 0) || (tf1.compareTo(tf2) == 0 && !allow_self_IAs) )
						continue;
					HashMap<String, LinkedList<BindingSite[]>> calc_start = prec_pw_adj.get(new StrPair(tf1, tf2));
					
					// pairwise comparison already doesn't contain the target
					if (!calc_start.containsKey(target)) 
						continue;
					
					// add next TFs, but try every order of permutations
					List<String> remaining_TFs = new ArrayList<>(TFs);
					remaining_TFs.remove(tf1);
					remaining_TFs.remove(tf2);
					for (List<Integer> permutation:Utilities.getAllIntPermutations(remaining_TFs.size())) {
						LinkedList<BindingSite[]> temp_result = calc_start.get(target);
						for (int index_tf3:permutation) {
							String tf3 = remaining_TFs.get(index_tf3);
							temp_result = getAdjacentToComplex(temp_result, tf3, target, d_min, d_max);
							
							if (temp_result.size() == 0) // no need to try further in this tf-ordering
								break;
						}
						
						if (temp_result.size() == 0) // try next ordering
							continue;
						
						found = true;
					}
					
					// shortcut
					if (found)
						break;
				} // end TF2
				
				// shortcut
				if (found)
					break;
			}// end TF1 iteration
			
			if (found) 
				result.add(target);
			
		}
		
		return result;
	}
	
	/**
	 * Checks for given TFs and pairwise distance constraints if there are target genes and BSs that are compatible, 
	 * also returns BS-combinations (may have duplicates); self-interactions may or may not be allowed.
	 * @param TFs
	 * @param d_min
	 * @param d_max
	 * @return
	 */
	public HashMap<String, LinkedList<BindingSite[]>> getAdjacencyPossibilitiesExact(Collection<String> TFs, int d_min, int d_max, boolean allow_self_IAs) {
		
		List<String> tf_list = new ArrayList<String>(TFs);
		HashMap<String, LinkedList<BindingSite[]>> result = new HashMap<>();
		
		// catch special cases
		if (TFs.size() == 0)
			return result;
		else if (TFs.size() == 1) {
			for (String target:this.TF_to_targets.getOrDefault(tf_list.get(0), new HashMap<>()).keySet()) {
				result.put(target, new LinkedList<>());
				for (BindingSite bs:this.TF_to_targets.getOrDefault(tf_list.get(0), new HashMap<>()).get(target))
					result.get(target).add(new BindingSite[]{bs});
			}
			return result;
		} else if (TFs.size() == 2) {
			return this.getPairwiseAdjacent(tf_list.get(0), tf_list.get(1), d_min, d_max);
		}
		
		// precompute all pairwise
		Set<String> common_targets = new HashSet<String>();
		Map<StrPair, HashMap<String, LinkedList<BindingSite[]>>> prec_pw_adj = new HashMap<>();
		for (String tf1:TFs)
			for (String tf2:TFs) {
				if ( (tf1.compareTo(tf2) < 0) || (tf1.compareTo(tf2) == 0 && !allow_self_IAs) )
					continue;
				StrPair pair = new StrPair(tf1, tf2);
				HashMap<String, LinkedList<BindingSite[]>> temp = getPairwiseAdjacent(tf1, tf2, d_min, d_max);
				prec_pw_adj.put(pair, temp);
				common_targets.addAll(temp.keySet());
			}
		
		// check for all common/shared targets
		common_targets.retainAll(getCommonTargets(TFs));
		for (String target:common_targets) {
			LinkedList<BindingSite[]> all_pos = new LinkedList<>();
			
			// for all pairwise starting positions
			for (String tf1:TFs) 
				for (String tf2:TFs) {
					if ( (tf1.compareTo(tf2) < 0) || (tf1.compareTo(tf2) == 0 && !allow_self_IAs) )
						continue;
					HashMap<String, LinkedList<BindingSite[]>> calc_start = prec_pw_adj.get(new StrPair(tf1, tf2));
					
					// pairwise comparison already doesn't contain the target
					if (!calc_start.containsKey(target)) 
						continue;
					
					// add next TFs, but try every order of permutations
					List<String> remaining_TFs = new ArrayList<>(TFs);
					remaining_TFs.remove(tf1);
					remaining_TFs.remove(tf2);
					for (List<Integer> permutation:Utilities.getAllIntPermutations(remaining_TFs.size())) {
						LinkedList<BindingSite[]> temp_result = calc_start.get(target);
						for (int index_tf3:permutation) {
							String tf3 = remaining_TFs.get(index_tf3);
							temp_result = getAdjacentToComplex(temp_result, tf3, target, d_min, d_max);
							
							if (temp_result.size() == 0) // no need to try further in this tf-ordering
								break;
						}
						
						if (temp_result.size() == 0) // try next ordering
							continue;
						
						// add to found alternative positions
						all_pos.addAll(temp_result);
					}
				}
			// end TF1 iteration
			
			if (all_pos.size() > 0) {
				result.put(target, all_pos);
			}
		}
		
		return result;
	}
	
	/**
	 * Compute special case of pairwise adjacency, more specific
	 * @param tf1
	 * @param tf2
	 * @param d_min
	 * @param d_max
	 * @return
	 */
	private HashMap<String, LinkedList<BindingSite[]>> getPairwiseAdjacent(String tf1, String tf2, String regulated_target, int d_min, int d_max) {
		HashMap<String, LinkedList<BindingSite[]>> pw_results = new HashMap<>();
		Set<String> targets_to_check;
		
		// check if a directed query or a general query
		if (regulated_target == null) {
			Set<String> tf_set = new HashSet<>();
			tf_set.add(tf1);
			tf_set.add(tf2);
			targets_to_check = this.getCommonTargets(tf_set);
		} else {
			targets_to_check = new HashSet<>();
			targets_to_check.add(regulated_target);
		}
		
		for (String target:targets_to_check) {
			LinkedList<BindingSite[]> adj_intervals = new LinkedList<>();
			
			for (BindingSite bs1:this.TF_to_targets.get(tf1).get(target))
				for (BindingSite bs2:this.TF_to_targets.get(tf2).get(target)) {
					int adj_state = bs1.adjacentToRelation(bs2, d_min, d_max);
					// add valid interval as an possibility
					if (adj_state == 1) {
						adj_intervals.add(new BindingSite[]{bs1, bs2});
					} else if (adj_state == -1) {
						adj_intervals.add(new BindingSite[]{bs2, bs1});
					}
				}
			
			// if compatible intervals found
			if (adj_intervals.size() > 0) {
				if (adj_intervals.size() > 1)
					adj_intervals.sort(new Comparator<BindingSite[]>() {

					@Override
					public int compare(BindingSite[] o1, BindingSite[] o2) {
						return o1[0].compareTo(o2[0]);
					}
				});
				pw_results.put(target, adj_intervals);
			}
		}
		return pw_results;
	}
	/**
	 * Compute special case of pairwise adjacency
	 * @param tf1
	 * @param tf2
	 * @param d_min
	 * @param d_max
	 * @return
	 */
	private HashMap<String, LinkedList<BindingSite[]>> getPairwiseAdjacent(String tf1, String tf2, int d_min, int d_max) {
		return this.getPairwiseAdjacent(tf1, tf2, null, d_min, d_max);
	}
	
	private LinkedList<BindingSite[]> getAdjacentToComplex(LinkedList<BindingSite[]> prev_possibilities, String tf, String target, int d_min, int d_max) {
		
		LinkedList<BindingSite[]> adj_intervals = new LinkedList<>();
		LinkedList<BindingSite> pos_to_add = this.TF_to_targets.get(tf).get(target);
		
		for (BindingSite[] fixed_positions:prev_possibilities)  {
			
			// precompute ranges of interest
			int limit_left = fixed_positions[0].getLeft() - d_max;
			int limit_right = fixed_positions[fixed_positions.length-1].getRight() + d_max;
			int crit_left = fixed_positions[0].getLeft() - d_min;
			int crit_right = fixed_positions[fixed_positions.length-1].getRight() + d_min;
			
			for (BindingSite add_bs:pos_to_add) {
				
				// if below -> try next
				if (add_bs.getRight() < limit_left)
					continue;
				
				// if above -> all the further ones will be (since sorted)
				if (add_bs.getLeft() > limit_right)
					break;
				
				// if overlap -> try next
				if ( (add_bs.getLeft() > crit_left && add_bs.getLeft() < crit_right)
						|| (add_bs.getRight() > crit_left && add_bs.getRight() < crit_right))
					continue;
				
				// if possible: add usable position
				BindingSite[] new_fixed = new BindingSite[fixed_positions.length+1];
				for (int i = 0; i < fixed_positions.length; i++)
					new_fixed[i] = fixed_positions[i];
				new_fixed[fixed_positions.length] = add_bs;
				Arrays.sort(new_fixed, null);
				adj_intervals.add(new_fixed);
			}
		}
		
		if (adj_intervals.size() > 1)
			adj_intervals.sort(new Comparator<BindingSite[]>() {

			@Override
			public int compare(BindingSite[] o1, BindingSite[] o2) {
				return o1[0].compareTo(o2[0]);
			}
		});
		
		return adj_intervals;
	}
	
	public Map<Set<String>, LinkedList<BindingSite[]>> getPossibleComplexPartnersAtTarget(String target, int d_min, int d_max, boolean allow_self_IAs) {
		Set<String> targeting_tfs = this.target_to_TFs.get(target);
		Map<Set<String>, LinkedList<BindingSite[]>> found_TF_tuples = new HashMap<>();
		
		for (String tf1:targeting_tfs)
			for (String tf2:targeting_tfs) {
				if ( (tf1.compareTo(tf2) < 0) || (tf1.compareTo(tf2) == 0 && !allow_self_IAs) )
					continue;
				HashMap<String, LinkedList<BindingSite[]>> result = getPairwiseAdjacent(tf1, tf2, target, d_min, d_max);
				if (!result.containsKey(target))
					continue;
				Set<String> involved_tfs = new HashSet<>();
				involved_tfs.add(tf1);
				involved_tfs.add(tf2);
				found_TF_tuples.put(involved_tfs, result.get(target));
			}
		
		return found_TF_tuples;
	}
	
	/**
	 * Writes data to a BED file
	 * @param trackname
	 * @param chromosome
	 * @param output_file
	 * @param binding_data
	 */
	public static void writeBED(String trackname, String chromosome, String output_file, Map<Collection<String>, LinkedList<BindingSite[]>> binding_data) {
		List<String> output_content = new LinkedList<>();
		
		// write header
		String header = "track name=" + trackname;
		output_content.add(header);
		
		// write content
		Map<Integer, LinkedList<String>> sorted_data = new TreeMap<>();
		for (Collection<String> tfs:binding_data.keySet()) {
			String tf_complex = String.join("/", tfs);
			for (BindingSite[] bss:binding_data.get(tfs)) {
				int left = bss[0].getLeft();
				int right = bss[bss.length-1].getRight();
				String strand = "+";
				int score = 0;
				String temp = chromosome + " " + left + " " + right + " " + tf_complex + " " + score + " " + strand;
				if (!sorted_data.containsKey(left))
					sorted_data.put(left, new LinkedList<>());
				sorted_data.get(left).add(temp);
			}
		}
		
		for (int left:sorted_data.keySet()) {
			List<String> ensured_sorting = new ArrayList<>(new HashSet<>(sorted_data.get(left)));
			// sort by right
			ensured_sorting.sort(new Comparator<String>() {

				@Override
				public int compare(String s1, String s2) {
					int r1 = Integer.parseInt(s1.split(" ")[2]);
					int r2 = Integer.parseInt(s2.split(" ")[2]);
					return Integer.compare(r1, r2);
				}
			});
			for (String temp:ensured_sorting)
				output_content.add(temp);
		}
		// write to file
		Utilities.writeEntries(output_content, output_file);
	}
}
