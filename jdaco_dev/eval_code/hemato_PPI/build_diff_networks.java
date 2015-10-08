package hemato_PPI;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import framework.PPIN;
import framework.RewiringDetector;
import framework.Utilities;

public class build_diff_networks {
	
	static double FDR = 0.05;
	static String network_folder = "/Users/tho/Desktop/BLUEPRINT_networks_0.03125/";
	static String results_root = "/Users/tho/Desktop/test/";
	
	public static void main(String[] args) {
		
		// read all data
		
		List<PPIN> HSCs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "HSC/", ".tsv.gz"))
			HSCs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> MPPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "MPP/", ".tsv.gz"))
			MPPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> CLPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "CLP/", ".tsv.gz"))
			CLPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> CMPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "CMP/", ".tsv.gz"))
			CMPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> GMPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "GMP/", ".tsv.gz"))
			GMPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> MEPs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "MEP/", ".tsv.gz"))
			MEPs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> EBs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "EB/", ".tsv.gz"))
			EBs.add(new PPIN(f.getAbsolutePath()));
		
		List<PPIN> MKs = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(network_folder + "MK/", ".tsv.gz"))
			MKs.add(new PPIN(f.getAbsolutePath()));
		
		// process data
		
		RewiringDetector.comparePPINs(HSCs, MPPs, FDR, results_root + "HSC_MPP/");
		RewiringDetector.comparePPINs(MPPs, CMPs, FDR, results_root + "MPP_CMP/");
		RewiringDetector.comparePPINs(MPPs, CLPs, FDR, results_root + "MPP_CLP/");
		RewiringDetector.comparePPINs(CMPs, MEPs, FDR, results_root + "CMP_MEP/");
		RewiringDetector.comparePPINs(CMPs, GMPs, FDR, results_root + "CMP_GMP/");
		RewiringDetector.comparePPINs(MEPs, MKs, FDR, results_root + "MEP_MK/");
		RewiringDetector.comparePPINs(MEPs, EBs, FDR, results_root + "MEP_EB/");
	}
}
