package standalone_tools;

/**
 * PPIcompare cmd-tool
 * @author Thorsten Will
 */
public class PPIcompare {
	
	public static void printHelp() {
		System.out.println("usage: java -jar PPIcompare.jar [folder1] [folder2]");
		System.exit(0);
	}
	
	public static void main(String[] args) {
		
		if (args.length < 2) {
			printHelp();
		}
		
		// TODO: all
		
	}
}
