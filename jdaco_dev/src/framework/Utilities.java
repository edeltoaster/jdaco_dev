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
	 * Reads text file line by line with, for example, proteins (also works for .gz)
	 * @param in_file
	 * @return entry per line
	 */
	public static HashSet<String> readEntryFile(String in_file) {
		HashSet<String> input = new HashSet<>();
		try {
			BufferedReader in = null;
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
			in.close();
		} catch (Exception e) {
			if (e instanceof FileNotFoundException)
				System.err.println("Problem while opening file " + in_file + ".");
			else
				System.err.println("Problem while parsing file " + in_file + ".");
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
		
		for (File f:listDirectoriesAndFiles(new File(path))) {
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
		
		for (File f:listDirectoriesAndFiles(new File(path))) {
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
		for (File f:listDirectoriesAndFiles(new File(path))) {
			if (f.getName().endsWith(suffix) || f.getName().endsWith(gz_suffix))
				paths.add(f);
		}
		
		paths.sort(null);
		
		return paths;
	}
	
	/**
	 * Lists all files in any folder within the given directory
	 * @param path
	 * @return
	 */
	private static List<File> listDirectoriesAndFiles(File path) {
		List<File> paths = new LinkedList<>();

		for (File f:path.listFiles()) {
			if (f.isHidden() || f.getName().startsWith("."))
				continue;
			if (f.isDirectory())
				paths.addAll(listDirectoriesAndFiles(f));
			else if (f.isFile()) {
				paths.add(f);
			}
		}
		return paths;
	}
	
	/**
	 * Writes entries to a file
	 * @param file
	 */
	public static void writeEntries(Collection<String> data, String file) {
		try {
			
			BufferedWriter bw = null;
			if (file.endsWith(".gz") || file.endsWith(".gzip"))
				bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file))));
			else
				bw = new BufferedWriter(new FileWriter(file));
			
			for (String entry:data) {
				bw.write(entry);
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			System.err.println("Problem while trying to write " + file);
		}
	}
	
	
	/*
	 * Statistic helpers
	 */
	
	
	public static double getMean(Collection<Double> data) {
		double sum = 0;
		for (Double d:data)
			sum += d;
		
		return sum / data.size();
	}
	
	public static double getMedian(Collection<Double> data) {
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
	
	public static double getVariance(Collection<Double> data) {
		double ssum = 0;
		double mean = getMean(data);
		
		for (Double d:data)
			ssum += (mean-d) * (mean-d);
		
		return ssum / data.size();
	}
	
	public static double getStd(Collection<Double> data) {
		return Math.sqrt(getVariance(data));
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
	
	public static List<List<Integer>> putInEveryPosition(int x, List<List<Integer>> list_of_lists) {
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
    
}
