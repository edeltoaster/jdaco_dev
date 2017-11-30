package standalone_tools;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DiffComplexDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;

/**
 * DiffComp cmd-tool
 * @author Thorsten Will
 */
public class DiffComp {
	static String version_string = "DiffComp 1.0pre";
	
	private static String path_group1 = "[GROUP1-FOLDER]";
	private static String path_group2 = "[GROUP2-FOLDER]";
	private static Map<String, QuantDACOResultSet> group1;
	private static Map<String, QuantDACOResultSet> group2;
	private static String path_seed = "no seed file defined";
	private static Set<String> seed = null;
	private static String output_folder;
	
	private static double FDR = 0.05;
	private static boolean parametric = false;
	private static boolean paired = false;
	private static boolean incorporate_supersets = false;
	private static double min_variant_fraction = 0.75;
	private static boolean human_readable = false;
	private static int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	
	/**
	 * Reads all matching output pairs of JDACO results/major transcripts from a certain folder: 
	 * [folder]/[sample-string]_complexes.txt(.gz) and [folder]/[sample-string]_major-transcripts.txt(.gz)
	 * @param folder
	 * @return
	 */
	public static Map<String, QuantDACOResultSet> readQuantDACOComplexes(String folder) {
		Map<String, QuantDACOResultSet> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(folder, "_major-transcripts.txt")) {
			String gz = "";
			if (f.getName().endsWith(".gz"))
				gz = ".gz";
			String pre = f.getAbsolutePath().split("_major-transcripts")[0];
			String sample = f.getName().split("_major-transcripts")[0];
			QuantDACOResultSet qdr = new QuantDACOResultSet(pre + "_complexes.txt" + gz, seed, pre + "_major-transcripts.txt" + gz);
			data.put(sample, qdr);
		}

		return data;
	}
	
	/**
	 * Prints the help message
	 */
	public static void printHelp() {
		System.out.println("usage: java -jar DiffComp.jar ([OPTIONS]) [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-fdr=[FDR] : false discovery rate (default: 0.05)");
		System.out.println("	-mf=[MIN_VAR_FRACTION] : fraction of group a complex must be part of to be considered in the analysis (default: 0.75)");
		System.out.println("	-s=[SEED-FILE] : seed file used for (J)DACO complex predictions (default: none)");
		System.out.println("	-t=[#threads] : number of threads to use (default: #cores/2)");
		System.out.println("	-nd : assume normal distribution when testing (default: nonparametric)");
		System.out.println("	-p : assume paired/dependent data (default: unpaired/independent)");
		System.out.println("	-ss : also associate supersets (default: no supersets)");
		System.out.println("	-hr : additionally output human readable files with gene names and rounded numbers (default: no output)");
		
		System.out.println();
		
		System.out.println("[GROUPx-FOLDER] :");
		System.out.println("	Standard PPIXpress/JDACO-output is read from those folders. JDACO results and major transcripts are needed.");
		
		System.out.println();
		
		System.out.println("[OUTPUT-FOLDER]");
		System.out.println("	The outcome is written to this folder. If it does not exist, it is created.");
		
		System.out.println();
		
		System.exit(0);
	}
	
	
	/**
	 * Prints the version of the program
	 */
	public static void printVersion() {
		System.out.println(version_string);
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
			
			// output version
			else if (arg.equals("-version"))
				printVersion();
			
			// parse FDR
			else if (arg.startsWith("-fdr="))
				FDR = Double.parseDouble(arg.split("=")[1]);
			
			// parse min_variant_fraction
			else if (arg.startsWith("-mf="))
				min_variant_fraction = Double.parseDouble(arg.split("=")[1]);
			
			// add seed file
			else if (arg.startsWith("-s=")) {
				path_seed = arg.split("=")[1];
				seed = Utilities.readEntryFile(path_seed);
				
				if (seed == null)
					System.exit(1); // error will be thrown by the utility-function
				
				if (seed.isEmpty()) {
					System.err.println("Seed file " + path_seed + " is empty.");
					System.exit(1);
				}
			}
			
			// parse #threads
			else if (arg.startsWith("-t="))
				no_threads = Math.max(Integer.parseInt(arg.split("=")[1]), 1);
			
			// parametric?
			else if (arg.equals("-nd"))
				parametric = true;
			
			// paired?
			else if (arg.equals("-p"))
				paired = true;
			
			// supersets?
			else if (arg.equals("-ss"))
				incorporate_supersets = true;
			
			// output human readable files?
			else if (arg.equals("-hr"))
				human_readable = true;
			
			// read groupwise input
			else if (group1 == null) {
				path_group1 = arg;
				if (!path_group1.endsWith("/"))
					path_group1 += "/";
				group1 = readQuantDACOComplexes(path_group1);
			}
			else if (group2 == null) {
				path_group2 = arg;
				if (!path_group2.endsWith("/"))
					path_group2 += "/";
				group2 = readQuantDACOComplexes(path_group2);
			}
			
			// set output-folder
			else if (output_folder == null) {
				output_folder = arg;
				if (!output_folder.endsWith("/"))
					output_folder += "/";
			}
			
		}
		
		// some final checks
		if (group1 == null || group1.isEmpty() || group2 == null || group2.isEmpty()) {
			
			if (group1 == null || group1.isEmpty())
				System.err.println(path_group1 + " does not exist or contains no useable data.");
			if (group2 == null || group2.isEmpty())
				System.err.println(path_group2 + " does not exist or contains no useable data.");
			
			System.exit(1);
		}
		
		if (group1.size() < 2 || group2.size() < 2) {
			
			if (group1.size() < 2)
				System.err.println(path_group1 + " does not contain enough data. At least two samples per group are needed for a statistical evaluation.");
			if (group2.size() < 2)
				System.err.println(path_group2 + " does not contain enough data. At least two samples per group are needed for a statistical evaluation.");
			
			System.exit(1);
		}
		
		if (output_folder == null) {
			System.err.println("Please add an output folder.");
			System.exit(1);
		}
	}
	
	
	public static void main(String[] args) {
		
		if (args.length == 1 && args[0].equals("-version"))
			printVersion();
		
		if (args.length < 3) {
			printHelp();
		}
		
		// parse cmd-line and set all parameters
		try {
			parseInput(args);
		} catch (Exception e) {
			System.out.println("Something went wrong while reading the command-line, please check your input!");
			System.out.println();
			printHelp();
		}
		
		// preface
		if (group1.size() < 5 || group2.size() < 5) {
			System.out.println("Computations will run as intended, but be aware that at least five samples per group are recommended to gather meaningful results.");
		}
		
		System.out.println("Seed: " + path_seed);
		System.out.println("FDR: " + FDR);
		
		String test = "Wilcoxon rank sum test (unpaired, nonparametric)";
		if (paired && parametric)
			test = "paired t-test (paired, parametric)";
		else if (paired)
			test = "Wilcoxon signed-rank test (paired, nonparametric)";
		else if (parametric)
			test = "Welch's unequal variances t-test (unpaired, parametric)";
		
		System.out.println("Statistical test: " + test);
		System.out.println("Min. fraction: " + min_variant_fraction);
		System.out.println("Incorporate supersets: " + (incorporate_supersets ? "yes" : "no"));
		System.out.println();
		
		// computations
		System.out.println("Processing " + path_group1 + " (" + group1.keySet().size() + ") vs " + path_group2 + " (" + group2.keySet().size() + ") ...");
		System.out.flush();
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, FDR, parametric, paired, incorporate_supersets, min_variant_fraction, no_threads);
		
		if (seed != null)
			System.out.println(dcd.getRawPValues().size() + " complexes tested, " + dcd.getSignificanceSortedComplexes().size() +" (" + dcd.getSignSortedVariants(false, false).size() + " seed combination variants) significant.");
		else
			System.out.println(dcd.getRawPValues().size() + " complexes tested, " + dcd.getSignificanceSortedComplexes().size() + " significant.");
		System.out.flush();
		
		/*
		 * write file-output
		 */
		
		if (dcd.getSignificanceSortedComplexes().size() == 0) {
			System.out.println("Since there is nothing to report, no output is written.");
			System.exit(0);
		}
		
		// check if output-folder exists, otherwise create it
		File f = new File(output_folder);
		if (!f.exists())
			f.mkdir();
		
		// write output files
		System.out.println("Writing results ...");
		dcd.writeSignSortedComplexes(output_folder + "diff_complexes.txt", false);
		
		if (seed != null)
			dcd.writeSignSortedVariants(output_folder + "diff_seed_variants.txt", false);
		
		if (human_readable) {
			System.out.println("Retrieving data on gene names ...");
			System.out.flush();
			System.out.println("Output human readable results ...");
			dcd.writeSignSortedComplexes(output_folder + "diff_complexes_hr.txt", true);
			if (seed != null)
				dcd.writeSignSortedVariants(output_folder + "diff_seed_variants_hr.txt", true);
		}
	}
}