package standalone_tools;

import framework.PPIN;

/**
 * PPIXpress helper tool for Evans Kataka :-)
 * @author Thorsten Will
 */
public class PPIXpress_preout {

	
	public static void printHelp() {
		System.out.println("usage: java -jar PPIXpress_preout.jar [INPUT-NETWORK]");
		System.exit(0);
	}
	
	public static void main(String[] args) {
		
		if (args.length < 1)
			printHelp();
		
		PPIN original_ppin = new PPIN(args[0]);
		
		// write to file what is used in PPIXpress
		original_ppin.writePPIN("ppixpress_ppin.txt");
		
		// update accessions and write to file
		PPIN original_ppin_updated = original_ppin.updateUniprotAccessions();
		original_ppin_updated.writePPIN("ppixpress_ppin_updated.txt");
	}
}
