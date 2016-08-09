package framework;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class DACOResultSet {
	
	private final HashSet<HashSet<String>> result;
	private final Set<String> abundant_seed_poteins = new HashSet<>();
	private final Map<HashSet<String>, LinkedList<HashSet<String>>> seed_to_complex_map = new HashMap<>();
	
	public DACOResultSet(String daco_out_file, String seed_file) {
		
		this.result = new HashSet<>();
		this.readResultCSV(daco_out_file);
		
		this.buildData(Utilities.readEntryFile(seed_file));
	}
	
	public DACOResultSet(String daco_out_file, Set<String> seed) {
		
		this.result = new HashSet<>();
		this.readResultCSV(daco_out_file);
		
		this.buildData(seed);
	}

	public DACOResultSet(HashSet<HashSet<String>> results, Set<String> seed) {
		this.result = results;
		this.buildData(seed);
	}

	/**
	 * Actually reads the file
	 * @param daco_out_file
	 */
	private void readResultCSV(String daco_out_file) {
		
		BufferedReader in = null;
		try {
			if (daco_out_file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(daco_out_file))));
			else
				in = new BufferedReader(new FileReader(daco_out_file));
			
			while (in.ready()) {
				String line = in.readLine();
				this.result.add(new HashSet<>(Arrays.asList(line.trim().split(","))));
			}
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
			
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
	}
	
	/**
	 * Filters to the subset of seed proteins
	 * @param to_test
	 * @param seed
	 * @return
	 */
	private static HashSet<String> filterSeedProteins(Set<String> to_test, Set<String> seed) {
		HashSet<String> from_seed = new HashSet<>(to_test);
		from_seed.retainAll(seed);
		return from_seed;
	}
	
	/**
	 * Buildup some data structures that may be useful
	 * @param seed
	 */
	private void buildData(Set<String> seed) {
		
		// build seed-set to complex map and note seed proteins that are abundant in the result
		for (HashSet<String> set:this.result) {
			HashSet<String> from_seed = filterSeedProteins(set, seed);
			if (from_seed.size() > 0) {
				abundant_seed_poteins.addAll(from_seed);
				if (!this.seed_to_complex_map.containsKey(from_seed))
					this.seed_to_complex_map.put(from_seed, new LinkedList<>());
				this.seed_to_complex_map.get(from_seed).add(set);
			}
		}
		
	}
	
	/**
	 * Rebuilds useful data-structures on the basis of relevant seed
	 * @param seed
	 */
	public void rebuildData(Set<String> refined_seed) {
		abundant_seed_poteins.clear();
		this.seed_to_complex_map.clear();
		
		// build seed-set to complex map and note seed proteins that are abundant in the result
		for (HashSet<String> set:this.result) {
			HashSet<String> from_seed = filterSeedProteins(set, refined_seed);
			if (from_seed.size() > 0) {
				abundant_seed_poteins.addAll(from_seed);
				if (!this.seed_to_complex_map.containsKey(from_seed))
					this.seed_to_complex_map.put(from_seed, new LinkedList<HashSet<String>>());
				this.seed_to_complex_map.get(from_seed).add(set);
			}
		}
		
	}
	
	/**
	 * Writes output to a csv-file
	 * @param out_file
	 * @param results
	 */
	public void writeCSV(String out_file) {
		List<String> to_write = new LinkedList<>();
		
		for (HashSet<String> cluster : this.result)
			to_write.add( String.join(",", cluster) );
		
		Utilities.writeEntries(to_write, out_file);
	}
	
	
	public double getComplexSetsSimilarity(DACOResultSet result_set2) {
		Set<HashSet<String>> result2 = result_set2.getResult();
		
		double sum = 0.0;
		for (HashSet<String> set1:this.result) {
			double max_Jsim = 0.0;
			for (HashSet<String> set2:result2) {
				double Jsim = Utilities.getJaccardSimilarity(set1, set2);
				
				if (Jsim > max_Jsim)
					max_Jsim = Jsim;
				
				if (max_Jsim == 1.0)
					break;
			}
			sum += max_Jsim;
		}
		
		for (HashSet<String> set1:result2) {
			double max_Jsim = 0.0;
			for (HashSet<String> set2:this.result) {
				double Jsim = Utilities.getJaccardSimilarity(set1, set2);
				
				if (Jsim > max_Jsim)
					max_Jsim = Jsim;
				
				if (max_Jsim == 1.0)
					break;
			}
			sum += max_Jsim;
		}
		
		return sum / (this.result.size() + result2.size()); // similar to dice-similarity
	}
	
	public double getSeedVariantSetsSimilarity(DACOResultSet result_set2) {
		Set<HashSet<String>> result2 = result_set2.getSeedToComplexMap().keySet();
		
		double sum = 0.0;
		for (HashSet<String> set1:this.seed_to_complex_map.keySet()) {
			double max_Jsim = 0.0;
			for (HashSet<String> set2:result2) {
				double Jsim = Utilities.getJaccardSimilarity(set1, set2);
				
				if (Jsim > max_Jsim)
					max_Jsim = Jsim;
				
				if (max_Jsim == 1.0)
					break;
			}
			sum += max_Jsim;
		}
		
		for (HashSet<String> set1:result2) {
			double max_Jsim = 0.0;
			for (HashSet<String> set2:this.seed_to_complex_map.keySet()) {
				double Jsim = Utilities.getJaccardSimilarity(set1, set2);
				
				if (Jsim > max_Jsim)
					max_Jsim = Jsim;
				
				if (max_Jsim == 1.0)
					break;
			}
			sum += max_Jsim;
		}
		
		return sum / (this.seed_to_complex_map.size() + result2.size()); // similar to dice-similarity
	}
	
	public HashSet<HashSet<String>> getResult() {
		return this.result;
	}
	
	public Set<String> getAbundantSeedProteins() {
		return this.abundant_seed_poteins;
	}

	public Map<HashSet<String>, LinkedList<HashSet<String>>> getSeedToComplexMap() {
		return this.seed_to_complex_map;
	}
	
	public Map<Set<String>, List<Set<String>>> getGeneralSeedToComplexMap() {
		Map<Set<String>, List<Set<String>>> temp = new HashMap<>();
		for (Set<String> seed_proteins:this.seed_to_complex_map.keySet()) {
			temp.put(seed_proteins, new LinkedList<Set<String>>());
			temp.get(seed_proteins).addAll(this.seed_to_complex_map.get(seed_proteins));
		}
		return temp;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((seed_to_complex_map == null) ? 0 : seed_to_complex_map.hashCode());
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
		DACOResultSet other = (DACOResultSet) obj;
		if (seed_to_complex_map == null) {
			if (other.seed_to_complex_map != null)
				return false;
		} else if (!seed_to_complex_map.equals(other.seed_to_complex_map))
			return false;
		return true;
	}
}
