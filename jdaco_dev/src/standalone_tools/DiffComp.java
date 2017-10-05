package standalone_tools;

/**
 * DiffComp cmd-tool
 * @author Thorsten Will
 */
public class DiffComp {
	static String version_string = "DiffComp dev version";
	
	
	/**
	 * Prints the help message
	 */
	public static void printHelp() {
		System.out.println("usage: java -jar CompAre.jar ([OPTIONS]) [GROUP1-FOLDER] [GROUP2-FOLDER] [OUTPUT-FOLDER]");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-fdr=[FDR] : false discovery rate (default: 0.05)");
		System.out.println("	-t=[#threads] : number of threads to use (default: #cores/2)");

		
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
			printHelp();
		}
		
		// TODO: implement cmd-tool stub for DiffComp
	}
}