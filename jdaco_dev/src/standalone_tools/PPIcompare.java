package standalone_tools;

import java.io.File;

import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;

/**
 * PPIcompare cmd-tool
 * @author Thorsten Will
 */
public class PPIcompare {
	
	private static String original_ppin_path;
	private static String group1_folder;
	private static String group2_folder;
	private static double expr_threshold = 1.0;
	private static String output_folder;
	
	public static void printHelp() {
		System.out.println("usage: java -jar PPIXpress.jar [reference PPIN] [folder1] [folder2] [cutoff] [output-folder]");

		System.exit(0);
	}
	
	/**
	 * Parse arguments
	 * @param args
	 */
	public static void parseInput(String[] args) {
		
		for (String arg:args) {
			// help needed?
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("--h") || arg.equals("--help"))
				printHelp();
		}
		
		// check if all there
		original_ppin_path = args[0];
		group1_folder = args[1];
		group2_folder = args[2];
		
		try {
			expr_threshold = Double.parseDouble(args[3]);
		} catch (Exception e) {
			printHelp();
		}
		
		output_folder = args[4];
	}
	
	public static void main(String[] args) {
		
		if (args.length < 4) {
			printHelp();
		}
		
		// parse commandline and set all parameters
		parseInput(args);
		
		// check if output-folder exists, otherwise create it
		if (!output_folder.endsWith("/"))
			output_folder += "/";
		File f = new File(output_folder);
		if (!f.exists())
			f.mkdir();
		
		PPIN original_ppin = new PPIN(original_ppin_path);
		NetworkBuilder builder = new NetworkBuilder(original_ppin);
		
		System.out.println("BLUEPRINT net-builder, TPM threshold: " + expr_threshold);
		System.out.println("Original PPIN: " + original_ppin_path);
		System.out.println("Ensembl version: " + DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		System.out.println("3did:" + DataQuery.get3didVersion());
		System.out.println("iPfam:" + DataQuery.getIPfamVersion());
		
		
	}
}
