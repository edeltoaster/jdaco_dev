package diff_compl_mono_geu;

import java.io.File;
import java.util.List;
import java.util.Map;

import framework.RewiringDetector;
import framework.RewiringDetectorSample;
import framework.Utilities;

public class check_monocyte_rewiring {
	static String output_folder = "PPICompare_out/";
	static String group1_folder = "classical/";
	static String group2_folder = "nonclassical/";
	static double q_threshold = 0.05;
	
	public static void main(String[] args) {
		Map<String, RewiringDetectorSample> group1 = RewiringDetectorSample.readNetworks(group1_folder);
		Map<String, RewiringDetectorSample> group2 = RewiringDetectorSample.readNetworks(group2_folder);
		
		// some preface
		System.out.println("Processing " + group1_folder + " (" + group1.keySet().size() + ") vs " + group2_folder + " (" + group2.keySet().size() + ") : ");
		System.out.flush();
		RewiringDetector rd = new RewiringDetector(group1, group2, q_threshold, 64, System.out);
		double P_rew_rounded = Math.round(rd.getP_rew() * 1000d) / 1000d;
		
		// stdout-output
		System.out.println(rd.getNumberOfComparisons()  + " comparisons, " + "P_rew: " + P_rew_rounded + ", " + rd.getSignificantlyRewiredInteractions().size() + " significant rewiring events." );
		
		List<String> minReasons = null;
		if (rd.getSignificantlyRewiredInteractions().size() > 0) {
			minReasons = rd.getMinMostLikelyReasons();
			System.out.println(minReasons.size() + " transcriptomic alterations can explain all significant changes.");
		}
		
		System.out.flush();
		
		// check if output-folder exists, otherwise create it
		File f = new File(output_folder);
		if (!f.exists())
			f.mkdir();
		
		// write differential interaction attribute table
		rd.writeDiffnet(output_folder + "differential_network.txt");
		
		// write output of optimization approach
		if (minReasons != null)
			Utilities.writeEntries(minReasons, output_folder + "min_reasons.txt");
		
	}
}
