package standalone_tools;

import java.util.HashSet;
import java.util.concurrent.TimeUnit;

import framework.ConstructedNetworks;
import framework.DACO;
import framework.DDIN;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.Utilities;
/**
 * JDACO cmd-tool
 * @author Thorsten Will
 */
public class JDACO_old_style {
	
	public static void main(String[] args) {
		
		if (args.length < 5 || args.length > 7) {
			System.out.println("usage: java -jar JDACO.jar PPIN SEED-FILE PAIR-THRESHOLD DEPTH OUT-FILE [PROBABILITY-THRESHOLD=0.5] [THREADS]");
			System.exit(0);
		}
		
		// get parameters
		String ppi_file = args[0];
		String seed_protein_file = args[1];
		double pair_threshold = Double.parseDouble(args[2]);
		int max_depth = Integer.parseInt(args[3]);
		String output_file = args[4];
		double prob_threshold = 0.5;
		int threads = 2;
		
		if (args.length >= 6)
			prob_threshold = Double.parseDouble(args[5]);
		
		if (args.length == 7)
			threads = Integer.parseInt(args[6]);
		
		// output parameters
		System.out.println("Running JDACO using:");
		System.out.println("PPIN: "+ppi_file);
		System.out.println("Seed protein list: "+seed_protein_file);
		System.out.println("Pair threshold: "+pair_threshold);
		System.out.println("Max. depth of search: "+max_depth);
		System.out.println("Output file: "+output_file);
		System.out.println("Probability threshold: "+prob_threshold);
		System.out.println("#threads: "+threads);
		System.out.println("");
		
		// read input
		PPIN ppi = new PPIN(ppi_file);
		System.out.println("Constructing domain-domain interaction network with ...");
		ConstructedNetworks networks = NetworkBuilder.constructAssociatedIsoformNetworks(ppi);
		DDIN ddi = networks.getDDIN();
		System.out.println("... "+ddi.getNumberOfProteins()+" proteins,");
		System.out.println("... "+ddi.getNumberOfDomains()+" domains,");
		System.out.println("... "+ddi.getNumberOfDDIs()+" domain-domain interactions");
		System.out.println("using Ensembl's '"+networks.getDB()+"' data.");
		System.out.println("");
		System.out.flush();
		HashSet<String> seed = Utilities.readEntryFile(seed_protein_file);
		
		// initialize calculation and set parameters
		DACO daco = new DACO(networks, threads, max_depth, pair_threshold, prob_threshold, System.out);
		
		// some enforced cleaning, just in case :-)
		Runtime.getRuntime().gc();
		
		// carry out computation
		long start = System.currentTimeMillis();
		System.out.println("Computing ...");
		HashSet<HashSet<String>> results = new HashSet<HashSet<String>>();
		int n = seed.size();
		int i = 1;
		for (String tf : seed) {
			System.out.println(i + "/" + n + ": " + tf);
			results.addAll(daco.growPairs(tf, seed));
			++i;
		}
		
		long duration = System.currentTimeMillis() - start;
		System.out.println("Comp. time: " + TimeUnit.MILLISECONDS.toMinutes(duration) + " min");
		
		// filter and write output
		DACO.writeAndFilterOutput(output_file, results);
		System.out.println(results.size() + " candidates written to output.");
	}
}
