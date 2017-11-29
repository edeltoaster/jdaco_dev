package standalone_tools;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import framework.QuantDACOResultSet;
import framework.Utilities;

/**
 * DiffComp cmd-tool
 * @author Thorsten Will
 */
public class DiffComp {
	static String version_string = "DiffComp dev version";
	
	private static String path_group1 = "[GROUP1-FOLDER]";
	private static String path_group2 = "[GROUP2-FOLDER]";
	private static Map<String, QuantDACOResultSet> group1;
	private static Map<String, QuantDACOResultSet> group2;
	private static String output_folder;
	private static double FDR = 0.05;
	private static boolean parametric = false;
	private static boolean paired = false;
	private static boolean incorporate_supersets = false;
	private static int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1); // assuming HT/SMT systems
	// TODO: think about seed
	
	// TODO: edit doc
	/**
	 * Reads all matching PPIXpress output pairs of PPINs/major transcripts from a certain folder: 
	 * [folder]/[sample-string]_ppin.txt(.gz), [folder]/[sample-string]_ddin.txt(.gz) and [folder]/[sample-string]_major-transcripts.txt(.gz)
	 * Assumes the standard file endings "_ppin.txt(.gz), "_ddin.txt(.gz) and _major-transcripts.txt(.gz)".
	 * @param folder
	 * @return
	 */
	public static Map<String, QuantDACOResultSet> readQuantDACOComplexes(String folder) {
		Map<String, QuantDACOResultSet> data = new HashMap<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(folder, "_ppin.txt")) {
			String gz = "";
			if (f.getName().endsWith(".gz"))
				gz = ".gz";
			String pre = f.getAbsolutePath().split("_ppin")[0];
			String sample = f.getName().split("_ppin")[0];
			//TODO: fill (pre + "_ppin.txt" + gz, pre + "_ddin.txt" + gz, pre + "_major-transcripts.txt" + gz);
			
			//data.put(sample, rds);
		}

		return data;
	}
	
	/**
	 * Prints the help message
	 */
	public static void printHelp() {
		System.out.println("usage: java -jar CompAre.jar ([OPTIONS]) [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-fdr=[FDR] : false discovery rate (default: 0.05)");
		System.out.println("	-t=[#threads] : number of threads to use (default: #cores/2)");
		System.out.println("	-nd : assume normal distribution when testing (default: nonparametric)");
		System.out.println("	-p : assume paired/dependent data (default: unpaired/independent)");
		System.out.println("	-s : also associate supersets (default: no supersets)");
		
		System.out.println();
		
		System.out.println("[GROUPx-FOLDER] :");
		System.out.println("	Standard PPIXpress/JDACO-output is read from those folders. PPIN, major transcript and JDACO results are needed.");
		
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
			else if (arg.equals("-s"))
				incorporate_supersets = true;
			
			// read groupwise input
			else if (group1 == null) {
				path_group1 = arg;
				if (!path_group1.endsWith("/"))
					path_group1 += "/";
				group1 = null;	// TODO: write method, see also below
			}
			else if (group2 == null) {
				path_group2 = arg;
				if (!path_group2.endsWith("/"))
					path_group2 += "/";
				group2 = null;// TODO: write method, see also above
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
		
		if (group1.size() < 3 || group2.size() < 3) {
			
			if (group1.size() < 3)
				System.err.println(path_group1 + " does not contain enough data. At least three samples per group are needed for a statistical evaluation.");
			if (group2.size() < 3)
				System.err.println(path_group2 + " does not contain enough data. At least three samples per group are needed for a statistical evaluation.");
			
			System.exit(1);
		}
		
		if (output_folder == null) {
			System.err.println("Please add an output folder.");
			
			System.exit(1);
		}
	}
	
	
	public static void main(String[] args) {
		
		if (args.length < 3) { // TODO: check length when all parameters are added
			printHelp();
		}
		
		// parse cmd-line and set all parameters
		try {
			parseInput(args);
		} catch (Exception e){
			printHelp();
		}
		
		// TODO: implement cmd-tool stub for DiffComp
	}
}