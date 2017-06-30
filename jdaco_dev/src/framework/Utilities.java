package framework;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Mixed helper functions
 * @author Thorsten Will
 */
public class Utilities {
	
	
	/*
	 * IO helpers
	 */
	
	/**
	 * Reads text file line by line and stores as a set (also works for .gz)
	 * lines are trimmed and comment lines marked with # are discarded
	 * @param in_file
	 * @return entry per line
	 */
	public static Set<String> readEntryFile(String in_file) {
		HashSet<String> input = new HashSet<>();
		BufferedReader in = null;
		try {
			if (in_file.endsWith(".gz") || in_file.endsWith(".gzip"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(in_file))));
			else
				in = new BufferedReader(new FileReader(in_file));
			
			while (in.ready()) {
				String line = in.readLine();
				if (line.startsWith("#"))
					continue;
				input.add(line.trim());
			}
			
		} catch (Exception e) {
			if (e instanceof FileNotFoundException)
				System.err.println("Problem while opening file " + in_file + ".");
			else
				System.err.println("Problem while parsing file " + in_file + ".");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
				// not helpful to have any output here
			}
		}
		
		return input;
	}
	
	/**
	 * Reads text file line by line and retains order (also works for .gz),
	 * lines are trimmed and comment lines marked with # are discarded
	 * @param in_file
	 * @return entry per line
	 */
	public static List<String> readFile(String in_file) {
		List<String> input = new LinkedList<String>();
		BufferedReader in = null;
		try {
			if (in_file.endsWith(".gz") || in_file.endsWith(".gzip"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(in_file))));
			else
				in = new BufferedReader(new FileReader(in_file));
			
			while (in.ready()) {
				String line = in.readLine();
				if (line.startsWith("#"))
					continue;
				input.add(line.trim());
			}
			
		} catch (Exception e) {
			if (e instanceof FileNotFoundException)
				System.err.println("Problem while opening file " + in_file + ".");
			else
				System.err.println("Problem while parsing file " + in_file + ".");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
				// not helpful to have any output here
			}
		}
		
		return input;
	}
	
	/**
	 * Given a directory [path] and a filename, 
	 * returns the absolute paths of all files in any directory below the given one with a certain filename.
	 * Skips hidden dirs/files and the ones starting with a dot.
	 * @param path
	 * @param filename
	 * @return
	 */
	public static List<File> getAllMatchingFilesInSubfolders(String path, String filename) {
		List<File> paths = new LinkedList<>();
		
		for (File f:listAllFilesWithinFolder(new File(path))) {
			if (f.getName().equals(filename))
				paths.add(f);
		}
		
		paths.sort(null);
		
		return paths;
	}
	
	/**
	 * Given a directory [path] and a filename, 
	 * returns the absolute paths of all files in any directory below the given one starting with a certain prefix.
	 * Skips hidden dirs/files and the ones starting with a dot.
	 * @param path
	 * @param filename
	 * @return
	 */
	public static List<File> getAllPrefixMatchingFilesInSubfolders(String path, String prefix) {
		List<File> paths = new LinkedList<>();
		
		for (File f:listAllFilesWithinFolder(new File(path))) {
			if (f.getName().startsWith(prefix))
				paths.add(f);
		}
		
		paths.sort(null);
		
		return paths;
	}
	
	/**
	 * Given a directory [path] and a filename, 
	 * returns the absolute paths of all files in any directory below the given one ending with a certain suffix or the suffix + ".gz".
	 * Skips hidden dirs/files and the ones starting with a dot.
	 * @param path
	 * @param filename
	 * @return
	 */
	public static List<File> getAllSuffixMatchingFilesInSubfolders(String path, String suffix) {
		List<File> paths = new LinkedList<>();
		String gz_suffix = suffix + ".gz";
		for (File f:listAllFilesWithinFolder(new File(path))) {
			if (f.getName().endsWith(suffix) || f.getName().endsWith(gz_suffix))
				paths.add(f);
		}
		
		paths.sort(null);
		
		return paths;
	}
	
	/**
	 * Recursively lists all files in any folder within the given directory, excludes hidden files/dirs and those starting with a dot
	 * @param path
	 * @return
	 */
	public static List<File> listAllFilesWithinFolder(File path) {
	
		List<File> paths = new LinkedList<>();
		
		if (path == null || !path.exists())
			return paths;
		
		for (File f:path.listFiles()) {
			if (f.isHidden() || f.getName().startsWith("."))
				continue;
			if (f.isDirectory())
				paths.addAll(listAllFilesWithinFolder(f));
			else if (f.isFile()) {
				paths.add(f);
			}
		}
		return paths;
	}
	
	/**
	 * Lists files and folders within the given directory, excludes hidden files/dirs and those starting with a dot
	 * @param path
	 * @return
	 */
	public static List<File> listDirectoriesAndFilesWithinFolder(File path) {
		List<File> paths = new LinkedList<>();

		for (File f:path.listFiles()) {
			if (f.isHidden() || f.getName().startsWith("."))
				continue;
			paths.add(f);
		}
		return paths;
	}
	
	/**
	 * Writes entries to a file using the entries' toString() function, uses GZIP if file to write to ends with .gz/.gzip
	 * @param <T>
	 * @param file
	 */
	public static <T> void writeEntries(Collection<T> data, String file) {
		BufferedWriter bw = null;
		try {
			if (file.endsWith(".gz") || file.endsWith(".gzip"))
				bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file))));
			else
				bw = new BufferedWriter(new FileWriter(file));
			
			for (T entry:data) {
				bw.write(entry.toString());
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			System.err.println("Problem while trying to write " + file);
		} finally {
			try {
				bw.close();
			} catch (Exception e) {
				// not helpful to add something here
			}
		}
	}
	
	
	/*
	 * Statistic helpers
	 */
	
	/**
	 * Computes the mean
	 * @param data
	 * @return
	 */
	public static double getMean(Collection<Double> data) {
		double sum = 0;
		for (Double d:data)
			sum += d;
		
		return sum / data.size();
	}
	
	/**
	 * Computes the median
	 * @param data
	 * @return
	 */
	public static double getMedian(Collection<Double> data) {
		
		// if the list only includes one element
		if (data.size() == 1)
			for (Double d:data)
				return d;
		
		List<Double> listed_data = new ArrayList<>(data);
		
		// sorts in ascending order
		Collections.sort(listed_data);
		
		int middle = listed_data.size() / 2;
		
		// if there is a middle element, take that; otherwise return mean of the two elements in the middle
		if (listed_data.size() % 2 == 1) {
			return listed_data.get(middle);
		} else {
			return (listed_data.get(middle-1) + listed_data.get(middle)) / 2.0;
		}
	}
	
	/**
	 * Computes the variance
	 * @param data
	 * @return
	 */
	public static double getVariance(Collection<Double> data) {
		double ssum = 0;
		double mean = getMean(data);
		
		for (Double d:data)
			ssum += (mean-d) * (mean-d);
		
		return ssum / data.size();
	}
	
	/**
	 * Computes the standard deviation
	 * @param data
	 * @return
	 */
	public static double getStd(Collection<Double> data) {
		return Math.sqrt(getVariance(data));
	}
	
	/**
	 * Computes the RMSD between two List<Doubles> of equal length
	 * @param data1
	 * @param data2
	 * @return
	 */
	public static double getRMSD(List<Double> data1, List<Double> data2) {
		double rmsd = 0;
		List<Double> data1_array = new ArrayList<>(data1);
		List<Double> data2_array = new ArrayList<>(data2);
		for (int i = 0; i < data1_array.size(); i++) {
			rmsd += Math.pow((data1_array.get(i) - data2_array.get(i)), 2);
		}
		rmsd /= data1_array.size();
		return Math.sqrt(rmsd);
	}
	/**
	 * Converts a list of Doubles to a double[]
	 * @param list
	 * @return
	 */
	public static double[] getDoubleArray(List<Double> list) {
		if (list.size() == 0)
			return new double[]{0.0};
		return list.stream().mapToDouble(d->d).toArray();
	}
	
	/**
	 * Converts raw pvalues to Benjamini-Hochberg adjusted pvalues according to given FDR,
	 * only returns significant objects
	 * @param raw_pvalues
	 * @param FDR
	 * @return
	 */
	public static <T> Map<T, Double> convertRawPValuesToBHFDR(Map<T, Double> raw_pvalues, double FDR) {
		// initialize reverse map
		Map<Double, List<T>> p_to_obj = new HashMap<>();
		for (T obj:raw_pvalues.keySet()) {
			double p = raw_pvalues.get(obj);
			if (!p_to_obj.containsKey(p))
				p_to_obj.put(p, new LinkedList<T>());
			p_to_obj.get(p).add(obj);
		}
		
		// find largest k and convert
		List<Double> p_values = new ArrayList<>(raw_pvalues.values());
		int m = p_values.size();
		Collections.sort(p_values);
		int k = 1;
		Map<Double, Double> raw_to_adj_p = new HashMap<>();
		for (double p:p_values) {
		    raw_to_adj_p.put(p, (p* m) / k); // if multiple have the same rank, take largest k
		    k++;
		}
		
		p_values = new ArrayList<>(new HashSet<>(p_values));
		Collections.sort(p_values);
		Collections.reverse(p_values);
		
		Map<T, Double> sign_obj_to_adjusted_p = new HashMap<>();
		double adj_p_before = 1.0;
		for (double p:p_values) {
			double adj_p = Math.min(raw_to_adj_p.get(p), adj_p_before);
			
			raw_to_adj_p.put(p, adj_p);
			adj_p_before = adj_p;
			
			if (adj_p >= FDR) //skip calculations for non-significant
				continue;
			
			for (T obj:p_to_obj.get(p))
				sign_obj_to_adjusted_p.put(obj, adj_p);
		}
		
		return sign_obj_to_adjusted_p;
	}
	
	/*
	 * Permutation helpers
	 */
	
	
	/**
	 * Builds all permutations of the integers 0...n-1 (of size n) according to the Steinhaus–Johnson–Trotter algorithm
	 * @param size
	 * @return
	 */
	public static List<List<Integer>> getAllIntPermutations(int n) {
		List<List<Integer>> output = new LinkedList<>();
		List<Integer> start = new ArrayList<>(1);
		start.add(0);
		output.add(start);
		for (int i = 1; i < n; i++)
			output = putInEveryPosition(i, output);
		
		return output;
	}
	
	private static List<List<Integer>> putInEveryPosition(int x, List<List<Integer>> list_of_lists) {
		List<List<Integer>> new_list_of_lists = new LinkedList<>();
		int n = list_of_lists.get(0).size();
		// for every permutation so far
		for (List<Integer> sublist:list_of_lists) 
			for (int i = 0; i <= n; i++) { // for every position so far
				List<Integer> temp_list = new ArrayList<>(sublist);
				temp_list.add(i, x);
				new_list_of_lists.add(temp_list);
			}
		
		return new_list_of_lists;
	}
	
	/*
	 * Collection helpers
	 */
	
	/**
	 * Jaccard similarity between sets of arbitrary classes
	 * @param a
	 * @param b
	 * @return
	 */
	public static <S> double getJaccardSimilarity(Set<S> a, Set<S>b) {
		HashSet<S> intersection = new HashSet<S>(a);
		intersection.retainAll(b);
		HashSet<S> union = new HashSet<S>(a);
		union.addAll(b);
		
		return ( (double) intersection.size() ) / union.size();
	} 
	
	/**
	 * Builds a set of all values in a multimap
	 * @param map
	 * @return
	 */
	public static <K, V> Set<V> getValueSetFromMultimap(Map<K, List<V>> map) {
		Set<V> value_set = new HashSet<>();
		
		for (List<V> collection:map.values())
			value_set.addAll(collection);
		
		return value_set;
	}
	
	/**
	 * Breaks a list into chunks; for best performance, input should be an ArrayList
	 * @param input
	 * @param chunkSize
	 * @return
	 */
    public static <T> List<List<T>> partitionListIntoChunks(List<T> input, int chunkSize) {
    	// implementation from http://stackoverflow.com/questions/12026885/common-util-to-break-a-list-into-batch
        int inputSize = input.size();
        int chunkCount = (int) Math.ceil(inputSize / (double) chunkSize);

        Map<Integer, List<T>> map = new HashMap<>(chunkCount);
        List<List<T>> chunks = new ArrayList<>(chunkCount);

        for (int i = 0; i < inputSize; i++) {

            map.computeIfAbsent(i / chunkSize, (ignore) -> {

                List<T> chunk = new ArrayList<>();
                chunks.add(chunk);
                return chunk;

            }).add(input.get(i));
        }

        return chunks;
    }
    
    /**
     * Constructs powerset
     * @param data
     * @param min
     * @param max
     * @return
     */
    public static <T> Set<Set<T>> getPowerSet(Set<T> data, int min, int max) {
    	// implementation adjusted from http://stackoverflow.com/questions/1670862/obtaining-a-powerset-of-a-set-in-java
    	List<T> list = new ArrayList<T>(data);
    	int n = list.size();

    	Set<Set<T>> powerSet = new HashSet<Set<T>>();

    	for( long i = 0; i < (1 << n); i++) {
    	    Set<T> element = new HashSet<T>();
    	    for( int j = 0; j < n; j++ )
    	    	if( (i >> j) % 2 == 1 )
    	    		element.add(list.get(j));
    	    
    	    if (element.size() >= min && element.size() <= max)
    	    	powerSet.add(element); 
    	}

    	return powerSet;
    }
    
    /*
     * Graph algorithms
     */
    
    private static int tarjan_index = -1;
    
    /**
     * Tarjan's strongly connected component algorithm, see https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
     * @param graph
     */
    public static List<Set<String>> getSCCs(Map<String, List<String>> graph) {
    	
    	// determine all nodes
    	Set<String> nodes = new HashSet<>();
    	nodes.addAll(graph.keySet());
    	graph.values().stream().forEach(s -> nodes.addAll(s));
    	
    	// necessary data structures
    	Map<String, Integer> node_index = new HashMap<>();
    	Map<String, Integer> node_lowlink = new HashMap<>();
    	Deque<String> stack = new LinkedList<>();
    	Set<String> on_stack = new HashSet<>();
    	
    	List<Set<String>> SCCs = new LinkedList<>();
    	tarjan_index = 0;
    	for (String node:nodes) {
    		if (node_index.containsKey(node))
    			continue;
    		tarjan_recurse(node, graph, SCCs, node_index, node_lowlink, stack, on_stack);
    	}
    	
    	return SCCs;
    }
    
    private static void tarjan_recurse(String node, Map<String, List<String>> graph, List<Set<String>> SCCs, Map<String, Integer> node_index, Map<String, Integer> node_lowlink, Deque<String> stack, Set<String> on_stack) {
    	// sets depth to smallest unused index
    	node_index.put(node, tarjan_index);
		node_lowlink.put(node, tarjan_index);
		tarjan_index++;
		stack.push(node);
		on_stack.add(node);
		
		// depth first search
		if (graph.containsKey(node))
			for (String node2:graph.get(node)) {
				if (!node_index.containsKey(node2)) {
					tarjan_recurse(node2, graph, SCCs, node_index, node_lowlink, stack, on_stack);
					node_lowlink.put(node, Math.min(node_lowlink.get(node), node_lowlink.get(node2)));
				}
				else if (on_stack.contains(node2)) {
					node_lowlink.put(node, Math.min(node_lowlink.get(node), node_index.get(node2)));
				}
			}
		
		// if node is a root node, pop from stack and generate SCC
		if (node_lowlink.get(node).equals(node_index.get(node))) {
			
			Set<String> SCC = new HashSet<>();
			String node2 = null;
			do {
				node2 = stack.pop();
				on_stack.remove(node2);
				SCC.add(node2);
			} while (!node.equals(node2));
			SCCs.add(SCC);
		}
    }
}
