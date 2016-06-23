package standalone_tools;

import java.io.File;
import java.util.List;
import java.util.Map;

import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.StrPair;
import framework.Utilities;

/**
 * PPICompare cmd-tool
 * @author Thorsten Will
 */
public class PPICompare {
	
	private static String path_group1 = "[GROUP1-FOLDER]";
	private static String path_group2 = "[GROUP2-FOLDER]";
	private static Map<String, RewiringDetectorSample> group1;
	private static Map<String, RewiringDetectorSample> group2;
	private static String output_folder;
	private static double FDR = 0.05;
	private static boolean output_protein_attributes = true;
	private static int no_threads = Math.max(Runtime.getRuntime().availableProcessors() / 2, 1);
	
	public static void printHelp() {
		System.out.println("usage: java -jar PPICompare.jar ([OPTIONS]) [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-fdr=[FDR] : false discovery rate (default: 0.05)");
		System.out.println("	ENSEMBL VERSION?");
		System.out.println("	-t=[#threads] : number of threads to use (default: #cores/2)");
		
		System.out.println();
		
		System.out.println("[GROUPx-FOLDER] :");
		System.out.println("	Standard PPIXpress-output is read from those folders. PPIN and major transcript files are needed.");
		
		System.out.println();
		
		System.out.println("[OUTPUT-FOLDER]");
		System.out.println("	The outcome is written to this folder. If it does not exist, it is created.");
		
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
			
			// parse FDR
			else if (arg.startsWith("-fdr"))
				FDR = Double.parseDouble(arg.split("=")[1]);
			
			// parse #threads
			else if (arg.startsWith("-t"))
				no_threads = Math.max(Integer.parseInt(arg.split("=")[1]), 1);
			
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
		
		if (output_folder == null) {
			System.err.println("Please add an output folder.");
			
			System.exit(1);
		}
		
			
	}
	
	public static void main(String[] args) {
		
		if (args.length < 3) {
			printHelp();
		}
		
		// parse cmd-line and set all parameters
		try {
			parseInput(args);
		} catch (Exception e){
			e.printStackTrace();
			printHelp();
		}
		
		// retrieval with output
		// TODO: retrieval with output PPIXpress style
		
		// some preface
		System.out.print("Processing " + path_group1 + " (" + group1.keySet().size() + ") vs " + path_group2 + " (" + group2.keySet().size() + ") : ");
		
		// actual calculations
		RewiringDetector rd = new RewiringDetector(group1, group2, FDR, no_threads, null);
		Map<String, List<StrPair>> major_alt_splice_switches = rd.determineAltSplicingSwitches(true, false);
		Map<String, List<StrPair>> all_alt_splice_switches = rd.determineAltSplicingSwitches(false, true);
		double P_rew_rounded = (double) Math.round(rd.getP_rew() * 1000d) / 1000d;
		List<String> minReasons = rd.getMinMostLikelyReasons();
		
		// stdout-output
		System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getSignificantlyRewiredInteractions().size() + " dIAs" );
		System.out.println(major_alt_splice_switches.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches).size() + " diff. interactions.");
		System.out.println(all_alt_splice_switches.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches).size() + " diff. interactions that are mainly driven by AS events.");
		System.out.println(minReasons.size() + " alterations can explain all significant changes.");
		
		/*
		 * write file-output
		 */
		
		// check if output-folder exists, otherwise create it
		File f = new File(output_folder);
		if (!f.exists())
			f.mkdir();
		
		rd.writeDiffnet(output_folder + "differential_net.txt");
		Utilities.writeEntries(major_alt_splice_switches.keySet(), output_folder + "major_AS_proteins.txt");
		Utilities.writeEntries(all_alt_splice_switches.keySet(), output_folder + "contributing_AS_proteins.txt");
		Utilities.writeEntries(minReasons, output_folder + "min_reasons.txt");
		
		if (output_protein_attributes)
			rd.writeProteinAttributes(output_folder + "protein_attributes.txt");
		
	}
}
