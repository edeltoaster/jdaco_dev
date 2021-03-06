package standalone_tools;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPEnrichment;
import framework.DiffSeedCombDetector;
import framework.QuantDACOResultSet;
import framework.Utilities;

/**
 * CompleXChange cmd-tool
 * @author Thorsten Will
 */
public class CompleXChange {
	static String version_string = "CompleXChange 1.01";
	
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
	private static boolean also_seedcomb_calcs = false;
	private static boolean also_seed_enrich = false;
	private static boolean filter_allosome = false;
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
		System.out.println("usage: java -jar CompleXChange.jar ([OPTIONS]) [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-fdr=[FDR] : false discovery rate (default: 0.05)");
		System.out.println("	-mf=[MIN_VAR_FRACTION] : fraction of either group a complex must be part of to be considered in the analysis (default: 0.75)");
		System.out.println("	-s=[SEED-FILE] : file listing proteins around which the complexes are centered, e.g. seed file used for JDACO complex predictions (default: none)");
		System.out.println("	-t=[#threads] : number of threads to be used (default: #cores/2)");
		System.out.println("	-nd : assume normal distribution when testing (default: no assumption, non-parametric tests applied)");
		System.out.println("	-p : assume paired/dependent data (default: unpaired/independent tests applied)");
		System.out.println("	-ss : also associate supersets in the analyses (default: no supersets)");
		System.out.println("	-enr : determine seed proteins enriched in up/down-regulated complexes (as in GSEA)");
		System.out.println("	-sc : additional analysis based on seed combinations rather than sole complexes (as in GSEA)");
		System.out.println("	-hr : additionally output human readable files with gene names and rounded numbers (default: no output)");
		System.out.println("	-f : filter complexes with allosome proteins (default: no filtering)");
		
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
			
			// seed combinations?
			else if (arg.equals("-sc"))
				also_seedcomb_calcs = true;
			
			// seed enrichment?
			else if (arg.equals("-enr"))
				also_seed_enrich = true;
			
			// output human readable files?
			else if (arg.equals("-hr"))
				human_readable = true;
			
			// filter allosome proteins?
			else if (arg.equals("-f"))
				filter_allosome = true;
			
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
		
		// check for invalid usage of seedcomb mode
		if (also_seedcomb_calcs && seed == null) {
			also_seedcomb_calcs = false;
			System.out.println("Analysis of seed combination variants in complexes requires a seed file, this additional analysis is skipped.");
		}
		
		if (also_seed_enrich && seed == null) {
			also_seed_enrich = false;
			System.out.println("Analysis of seed protein enrichment in complexes requires a seed file, this additional analysis is skipped.");
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
		
		String test = "Wilcoxon rank sum test (unpaired, non-parametric)";
		if (paired && parametric)
			test = "paired t-test (paired, parametric)";
		else if (paired)
			test = "Wilcoxon signed-rank test (paired, non-parametric)";
		else if (parametric)
			test = "Welch's unequal variances t-test (unpaired, parametric)";
		
		System.out.println("Statistical test: " + test);
		System.out.println("Min. fraction: " + min_variant_fraction);
		System.out.println("Incorporate supersets: " + (incorporate_supersets ? "yes" : "no"));
		
		if (filter_allosome) {
			System.out.print("Filtering of complexes with allosome proteins enabled. Downloading data ... ");
			String db = DataQuery.getEnsemblOrganismDatabaseFromProteins(group1.values().iterator().next().getAbundantSeedProteins());
			DataQuery.getAllosomeProteins(db);
			System.out.println("done.");
		}	
		System.out.println();
		
		// computations
		boolean output_to_write = false;
		System.out.println("Processing " + path_group1 + " (" + group1.keySet().size() + ") vs " + path_group2 + " (" + group2.keySet().size() + ") ...");
		System.out.flush();
		
		DiffComplexDetector dcd = new DiffComplexDetector(group1, group2, FDR, parametric, paired, incorporate_supersets, min_variant_fraction, no_threads, filter_allosome);
		
		if (seed != null)
			System.out.println(dcd.getRawPValues().size() + " complexes tested, " + dcd.getSignificanceSortedComplexes().size() +" (" + dcd.getSignSortedVariants(false, false).size() + " seed combination variants) significant.");
		else
			System.out.println(dcd.getRawPValues().size() + " complexes tested, " + dcd.getSignificanceSortedComplexes().size() + " significant.");
		System.out.flush();
		
		output_to_write = !dcd.getSignificanceSortedComplexes().isEmpty();
		
		DiffSeedCombDetector dscd = null;
		if (also_seedcomb_calcs) {
			dscd = new DiffSeedCombDetector(group1, group2, FDR, parametric, paired, incorporate_supersets, min_variant_fraction, no_threads, filter_allosome);
			System.out.println(dscd.getVariantsRawPValues().size() + " seed combinations tested, " + dscd.getSignificanceSortedVariants().size() + " significant.");
			
			if (!dscd.getSignificanceSortedVariants().isEmpty())
				output_to_write = true;
			
			System.out.flush();
		}
		
		SPEnrichment spe = null;
		if (seed == null && also_seed_enrich) {
			System.out.println("Please specify a seed protein file if seed protein enrichment should be determined.");
			System.out.flush();
		}
		else if (seed != null && also_seed_enrich) {
			System.out.println("Calculating seed protein enrichment with 10000 permutations to approximate null distribution and only considering seed proteins participating in at least 10 complexes.");
			System.out.flush();
			spe = dcd.calculateSPEnrichment(FDR, 10000, 10);
			
			if (!spe.getSignificanceSortedSeedProteins().isEmpty())
				output_to_write = true;
		}
		
		/*
		 * write file-output
		 */
		
		if (!output_to_write) {
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
		
		if (seed != null) {
			dcd.writeSignSortedVariants(output_folder + "diff_seedcomb_in_complexes.txt", false);
			
			if (also_seedcomb_calcs)
				dscd.writeSignSortedVariants(output_folder + "diff_seed_combinations.txt", false);
			
			if (also_seed_enrich)
				spe.writeSignificantSeedProteins(output_folder + "enriched_seed_proteins.txt");
		}
		
		if (human_readable) {
			System.out.println("Retrieving data on gene names ...");
			System.out.flush();
			System.out.println("Output human readable results ...");
			dcd.writeSignSortedComplexes(output_folder + "diff_complexes_hr.txt", true);
			if (seed != null) {
				dcd.writeSignSortedVariants(output_folder + "diff_seedcomb_in_complexes_hr.txt", true);
				
				if (also_seedcomb_calcs)
					dscd.writeSignSortedVariants(output_folder + "diff_seed_combinations_hr.txt", true);
			}
		}
	}
}