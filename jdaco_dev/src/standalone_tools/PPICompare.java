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
	
	private static String path_group1;
	private static String path_group2;
	private static Map<String, RewiringDetectorSample> group1;
	private static Map<String, RewiringDetectorSample> group2;
	private static String output_folder;
	private static double FDR = 0.05;
	
	public static void printHelp() {
		System.out.println("usage: java -jar PPICompare.jar [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER] ([FDR])");
		
		System.out.println("[GROUPx-FOLDER] :");
		System.out.println("	Standard PPIXpress-output is read from those folders. PPIN and major transcript files are needed.");
		
		System.out.println("");
		
		System.out.println("[OUTPUT-FOLDER]");
		System.out.println("	The outcome is written to this folder. If it does not exist, it is created.");
		
		System.out.println("");
		
		System.out.println("optional: [FDR]");
		System.out.println("	False discovery rate (default: 0.05)");
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
			
			// parse FDR
			else {
				FDR = Double.parseDouble(arg);
			}
			
		}
	}
	
	public static void main(String[] args) {
		if (args.length != 3 && args.length != 4) {
			printHelp();
		}
		
		// parse cmd-line and set all parameters
		try {
			parseInput(args);
		} catch (Exception e){
			printHelp();
		}
		
		// do
		System.out.print("Processing " + path_group1 + " (" + group1.keySet().size() + ") vs " + path_group2 + " (" + group2.keySet().size() + ") : ");
		RewiringDetector rd = new RewiringDetector(group1, group2, FDR);
		Map<String, List<StrPair>> major_alt_splice_switches = rd.determineAltSplicingSwitches(true, false);
		Map<String, List<StrPair>> all_alt_splice_switches = rd.determineAltSplicingSwitches(false, true);
		double P_rew_rounded = (double) Math.round(rd.getP_rew() * 1000d) / 1000d;
		
		// stdout-output
		System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getInteractionReasonsCountMap().size() + " dIAs" );
		System.out.println(major_alt_splice_switches.keySet().size() + " alt. spliced proteins are the major reason that affect " + Utilities.getValueSetFromMultimap(major_alt_splice_switches).size() + " diff. interactions.");
		System.out.println(all_alt_splice_switches.keySet().size() + " alt. spliced proteins contribute to a change in the " + Utilities.getValueSetFromMultimap(all_alt_splice_switches).size() + " diff. interactions that are mainly driven by AS events.");
		
		List<String> minReasons = rd.getMinMostLikelyReasons();
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
		rd.writeProteinAttributes(output_folder + "protein_attributes.txt");
		
	}
}
