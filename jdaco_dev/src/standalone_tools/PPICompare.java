package standalone_tools;

import java.io.File;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

import framework.DataQuery;
import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.Utilities;

/**
 * PPICompare cmd-tool
 * @author Thorsten Will
 */
public class PPICompare {
	static String version_string = "PPICompare 1.04";
	
	private static String path_group1 = "[GROUP1-FOLDER]";
	private static String path_group2 = "[GROUP2-FOLDER]";
	private static String output_folder;
	private static double FDR = 0.05;
	private static boolean output_protein_attributes = true;
	
	private static int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1);
	private static PrintStream verbose = null;
	
	private static Map<String, RewiringDetectorSample> group1;
	private static Map<String, RewiringDetectorSample> group2;

	
	/**
	 * Prints the help message
	 */
	public static void printHelp() {
		System.out.println("usage: java -jar PPICompare.jar ([OPTIONS]) [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-fdr=[FDR] : false discovery rate (default: 0.05)");
		System.out.println("	-n : no output of protein attribute table (no data retrieval necessary)");
		System.out.println("	-r=[release] : try to use data from a certain Ensembl release, uses newest if specific release not found");		
		System.out.println("	-s=[US, UK, AS, 'specific URL'] : change initial server, note that US and asian mirrors only store the last two releases (default: US)");
		System.out.println("	-t=[#threads] : number of threads to use (default: #cores/2)");
		System.out.println("	-v : be more verbose");
		
		System.out.println();
		
		System.out.println("[GROUPx-FOLDER] :");
		System.out.println("	Standard PPIXpress-output is read from those folders. PPIN, DDIN and major transcript files are needed.");
		
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
			else if (arg.startsWith("-fdr"))
				FDR = Double.parseDouble(arg.split("=")[1]);
			
			// toggle writing of protein output table
			else if (arg.equals("-n"))
				output_protein_attributes = false;
			
			// try to enforce a release manually
			else if (arg.startsWith("-r="))
				DataQuery.enforceSpecificEnsemblRelease(arg.split("=")[1]);
			
			// server switching
			else if (arg.startsWith("-s=")){
				String option = arg.split("=")[1];
				if (option.equals("UK"))
					DataQuery.switchServer("ensembldb.ensembl.org:3306");
				else if (option.equals("US"))
					DataQuery.switchServer("useastdb.ensembl.org:3306");
				else if (option.equals("AS"))
					DataQuery.switchServer("asiadb.ensembl.org:3306");
				else
					DataQuery.switchServer(option);
			}
			
			// parse #threads
			else if (arg.startsWith("-t="))
				no_threads = Math.max(Integer.parseInt(arg.split("=")[1]), 1);
			
			else if (arg.equals("-v"))
				verbose = System.out;
			
			// read groupwise input
			else if (group1 == null) {
				path_group1 = arg;
				if (!path_group1.endsWith("/"))
					path_group1 += "/";
				group1 = RewiringDetectorSample.readNetworks(path_group1);	
			}
			else if (group2 == null) {
				path_group2 = arg;
				if (!path_group2.endsWith("/"))
					path_group2 += "/";
				group2 = RewiringDetectorSample.readNetworks(path_group2);
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
		
		// some preface
		System.out.println("Processing " + path_group1 + " (" + group1.keySet().size() + ") vs " + path_group2 + " (" + group2.keySet().size() + ") : ");
		
		if (group1.size() < 3 || group2.size() < 3) {
			System.out.println("PPICompare will run as intended, but be aware that at least three samples per group are recommended to gather meaningful results.");
		}
		
		System.out.flush();
		
		// actual calculations
		RewiringDetector rd = new RewiringDetector(group1, group2, FDR, no_threads, verbose);
		double P_rew_rounded = Math.round(rd.getP_rew() * 1000d) / 1000d;
		
		// stdout-output
		System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getSignificantlyRewiredInteractions().size() + " significant rewiring events." );
		
		List<String> minReasons = null;
		if (rd.getSignificantlyRewiredInteractions().size() > 0) {
			minReasons = rd.getMinMostLikelyReasons();
			System.out.println(minReasons.size() + " transcriptomic alterations can explain all significant changes.");
		}
		
		System.out.flush();
		
		
		/*
		 * write file-output
		 */
		
		// check if output-folder exists, otherwise create it
		File f = new File(output_folder);
		if (!f.exists())
			f.mkdir();
		
		// write differential interaction attribute table
		rd.writeDiffnet(output_folder + "differential_network.txt");
		
		// write output of optimization approach
		if (minReasons != null)
			Utilities.writeEntries(minReasons, output_folder + "min_reasons.txt");
		
		// if necessary, write protein attribute table
		if (output_protein_attributes) {
			
			// version output and data retrieval
			String organism_database = rd.getAppropriateOrganismDatabase();
			String ensembl_version = organism_database.split("_")[organism_database.split("_").length-2];
			
			System.out.print("Retrieving ENSEMBL " + ensembl_version + " data from database " + organism_database + " (may take some minutes) ... ");
			DataQuery.getGenesCommonNames(organism_database);
			System.out.print("33% ... ");
			System.out.flush();
			
			DataQuery.getGenesTranscriptsProteins(organism_database);
			System.out.println("100%.");
			System.out.flush();
			
			// write node table
			rd.writeProteinAttributes(output_folder + "protein_attributes.txt");
		}
	}
}