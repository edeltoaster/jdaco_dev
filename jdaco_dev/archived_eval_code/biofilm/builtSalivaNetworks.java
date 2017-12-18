package biofilm;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import framework.ConstructedNetworks;
import framework.DDIN;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.Utilities;

public class builtSalivaNetworks {
	static String project_folder = "/Users/tho/GDrive/Work/biofilm/";
	
	public static void main(String[] args) {
		System.out.println("run at 03.07.15, PrePPI");
		Map<String, Set<String>> data_prot_map = new HashMap<String, Set<String>>();
		
		for (File file:Utilities.getAllSuffixMatchingFilesInSubfolders(project_folder + "merged_proteins/", "_merged.txt")) {
			String filename = file.getName();
			String sample = filename.split("_")[0];
			data_prot_map.put(sample, Utilities.readEntryFile(file.getAbsolutePath()));
		}
		
		// initialize stuff
		PPIN org_ppi = new PPIN("mixed_data/human_ppi.tsv.gz");
		NetworkBuilder builder = new NetworkBuilder(org_ppi);
		
	
		// versioning output
		System.out.println("versions:");
		System.out.println("Ensembl: " + builder.getDB());
		System.out.println("3did: " + DataQuery.get3didVersion());
		System.out.println("iPfam: " + DataQuery.getIPfamVersion());
		// process
		for (String sample:data_prot_map.keySet()) {
			System.out.println("Processing " + sample + " ...");
			ConstructedNetworks constr = builder.constructAssociatedNetworksFromProteinSet(data_prot_map.get(sample), false);
			PPIN ppi = constr.getPPIN();
			DDIN ddi = constr.getDDIN();
			
			System.out.println("PPI size: "+ ppi.getSizesStr());
			
			ppi.writePPIN(project_folder + "networks/" + sample +"_ppi.tsv");
			ddi.writeDDIN(project_folder + "networks/" + sample +"_ddi.tsv");
		}
		
	}

}
