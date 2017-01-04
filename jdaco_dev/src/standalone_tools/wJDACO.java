package standalone_tools;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import framework.ConstructedNetworks;
import framework.DACO;
import framework.Utilities;
import framework.wDACO;

/**
 * JDACO cmd-tool
 * @author Thorsten Will
 */
public class wJDACO {
	private static String ppin_file;
	private static String ddin_file;
	private static String abundance_file;
	private static String seed_file;
	
	private static ConstructedNetworks cn;
	private static int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1);
	private static String output_file;
	private static int max_depth = 0;
	private static double percentile = 5;
	private static Set<String> seed;
	private static double pair_threshold = -1.0;
	private static double prob_threshold = -1.0;
	private static PrintStream out = System.out;
	private static int compute_timeout = 60;
	
	public static void printHelp() {
		System.out.println("usage: java -jar JDACO.jar ([OPTIONS]) [PPI-NETWORK] [DDI-NETWORK] [MAJ-TRANSCRIPT-FILE] [SEED-FILE] [MAX-DEPTH] [OUT-FILE]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-t=[#threads] : number of threads to use (default: #cores/2)");
		System.out.println("	-p=[percentile] : percentile of distribution used to automatically construct start pairs and complex probability cutoff (default: 5 as 5th highest percentile)");
		System.out.println("	-pb=[threshold] : overwrite pair-building cutoff with [threshold] (default: based on distribution)");
		System.out.println("	-cp=[probability] : overwrite internal complex probability cutoff with [probability], set to zero to turn off (default: determined pair threshold**([MAX-DEPTH]-1))");
		System.out.println("	-ct=[minutes] : limits compute time per seed-protein to [minutes] (default: 60min)");
		System.out.println("	-s : silent mode, no output during actual computation");
		
		System.out.println();
		
		System.out.println("[PPI-NETWORK] / [DDI-NETWORK] / [MAJ-TRANSCRIPT-FILE]:");
		System.out.println("	Any protein-protein interaction network in SIF-format: Protein1 Protein2 (weight).");
		System.out.println("	Proteins are assumed to be given as UniProt or HGNC accessions.");
		System.out.println("	Any matching domain-domain interaction network and major transcript file as generated by PPIXpress.");
		
		System.out.println();
		
		System.out.println("[SEED-FILE] :");
		System.out.println("	Linewise list of seed proteins given as UniProt accessions.");
		
		System.out.println();
		
		System.out.println("[MAX-DEPTH] :");
		System.out.println("	Max. depth of search.");
		
		System.out.println();
		
		System.out.println("[OUT-FILE] :");
		System.out.println("	The outcome is written as a csv.");
		
		System.out.println();
		
		System.exit(0);
	}
	
	/**
	 * Parse arguments
	 * @param args
	 */
	public static void parseInput(String[] args) {
		
		for (String arg:args) {
			
			// help needed?
			if (arg.equals("-h") || arg.equals("-help"))
				printHelp();
			
			// threads
			else if (arg.startsWith("-t")) 
				no_threads = Math.max(Integer.parseInt(arg.split("=")[1]), 1);
			
			// pair building threshold
			else if (arg.startsWith("-pb"))
				pair_threshold = Double.parseDouble(arg.split("=")[1]);
			
			// percentile
			else if (arg.startsWith("-p"))
				percentile = Double.parseDouble(arg.split("=")[1]);
			
			// probability threshold
			else if (arg.startsWith("-cp"))
				prob_threshold = Double.parseDouble(arg.split("=")[1]);
			
			// compute time
			else if (arg.startsWith("-ct"))
				compute_timeout = Integer.parseInt(arg.split("=")[1]);
			
			// silent mode
			else if (arg.equals("-s")) 
				out = null;
			
			// other parameters
			else {
				if (ppin_file == null) {
					ppin_file = arg;
				}
				else if (ddin_file == null) {
					ddin_file = arg;
				}
				else if (abundance_file == null) {
					abundance_file = arg;
				}
				else if (seed_file == null) {
					seed_file = arg;
					seed = Utilities.readEntryFile(seed_file);
				}
				else if (max_depth == 0)
					max_depth = Integer.parseInt(arg);
				else if (output_file == null)
					output_file = arg;
			}
		}
		
		// check if input is set
		if (ppin_file == null) {
			System.err.println("No protein interaction network given.");
			System.exit(1);
		}
		
		if (ddin_file == null) {
			System.err.println("No domain interaction network given.");
			System.exit(1);
		}
		
		if (abundance_file == null) {
			System.err.println("No major transcript and abundance data given.");
			System.exit(1);
		}
		
		if (seed_file == null) {
			System.err.println("No seed protein file given.");
			System.exit(1);
		}
		if (max_depth < 2) {
			System.err.println("Max. depth should be at least two.");
			System.exit(1);
		}
		
		if (output_file == null) {
			System.err.println("No output file set.");
			System.exit(1);
		}
		
		cn = new ConstructedNetworks(ppin_file, ddin_file, abundance_file, null, false);
	}
	
	public static void main(String[] args) {

		if (args.length < 6 || args.length > 11) {
			printHelp();
		}
		
		// parse parameters
		parseInput(args);
		
		// output parameters
		System.out.println("Running JDACO using:");
		
		System.out.println("PPIN: " + ppin_file + " ("+ cn.getPPIN().getSizesStr() + ")");
		System.out.println("DDIN: " + ddin_file + " ("+ cn.getDDIN().getSizesStr() + ")");
		
		System.out.println("Seed protein list: " + seed_file + " (" + seed.size() + " proteins)");
		
		// determine thresholds
		if (pair_threshold == -1.0) {
			System.out.println("Percentile: " + percentile);
			pair_threshold = cn.getPPIN().getPercentile(percentile);
		}
		
		if (prob_threshold == -1.0)
			prob_threshold = Math.pow(pair_threshold, max_depth-1);
		
		System.out.println("Pair-building threshold: " + pair_threshold);
		System.out.println("Probability threshold: " + prob_threshold);
		System.out.println("Max. depth of search: " + max_depth);
		
		if (compute_timeout > 0) 
			System.out.println("Compute timeout: " + compute_timeout + " min");
		
		System.out.println("Output file: " + output_file);

		System.out.println("#threads: " + no_threads);
		System.out.println("");
		System.out.flush();

		// initialize calculation and set parameters
		wDACO wdaco = new wDACO(cn, no_threads, max_depth, pair_threshold, prob_threshold, out);
		
		
		// carry out computation
		long start = System.currentTimeMillis();
		System.out.println("Computing ...");
		HashSet<HashSet<String>> results = new HashSet<>();
		int n = seed.size();
		int i = 1;
		for (String tf : seed) {
			if (out != null)
				out.println(i + "/" + n + ": " + tf);
			results.addAll(wdaco.growPairs(tf, seed));
			++i;
		}
		long duration = System.currentTimeMillis() - start;
		
		
		// aftermath
		wdaco.ensurePoolShutdown();
		
		// filter and write output
		DACO.writeAndFilterOutput(output_file, results, seed);
		System.out.println(results.size() + " candidates written to output.");
		
		// below 2 min show results in seconds
		if (duration < 120000)
			System.out.println("Overall time: " + TimeUnit.MILLISECONDS.toSeconds(duration) + " sec");
		else
			System.out.println("Overall time: " + TimeUnit.MILLISECONDS.toMinutes(duration) + " min");
	}
}
