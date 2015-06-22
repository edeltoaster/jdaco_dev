package framework;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class DACOResultSet {
	
	private final HashSet<HashSet<String>> result;
	private final Set<String> abundant_seed_poteins = new HashSet<String>();
	private final Map<HashSet<String>, LinkedList<HashSet<String>>> seed_to_complex_map = new HashMap<HashSet<String>, LinkedList<HashSet<String>>>();
	
	public DACOResultSet(String daco_out_file, String seed_file) {
		
		this.result = new HashSet<HashSet<String>>();
		this.readResultCSV(daco_out_file);
		
		this.buildData(Utilities.readEntryFile(seed_file));
	}
	
	public DACOResultSet(String daco_out_file, Set<String> seed) {
		
		this.result = new HashSet<HashSet<String>>();
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
		try {
			BufferedReader in = null;
			if (daco_out_file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(daco_out_file))));
			else
				in = new BufferedReader(new FileReader(daco_out_file));
			
			while (in.ready()) {
				String line = in.readLine();
				this.result.add(new HashSet<String>(Arrays.asList(line.trim().split(","))));
			}
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	/**
	 * Filters to the subset of seed proteins
	 * @param to_test
	 * @param seed
	 * @return
	 */
	private static HashSet<String> filterSeedProteins(Set<String> to_test, Set<String> seed) {
		HashSet<String> from_seed = new HashSet<String>(to_test);
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
					this.seed_to_complex_map.put(from_seed, new LinkedList<HashSet<String>>());
				this.seed_to_complex_map.get(from_seed).add(set);
			}
		}
		
	}
	
	/**
	 * Rebuilds useful data-structures on the basis of relevant seed (proteins with BS-annotated)
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

	public HashSet<HashSet<String>> getResult() {
		return result;
	}
	
	/**
	 * Writes output to a csv-file
	 * @param out_file
	 * @param results
	 */
	public void writeCSV(String out_file) {
		
		// write
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(out_file));
			for (HashSet<String> cluster : this.result) {
				String temp = String.join(",", cluster);
				bw.write(temp);
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public Set<String> getAbundantSeedProteins() {
		return abundant_seed_poteins;
	}

	public Map<HashSet<String>, LinkedList<HashSet<String>>> getSeedToComplexMap() {
		return seed_to_complex_map;
	}
	
}
