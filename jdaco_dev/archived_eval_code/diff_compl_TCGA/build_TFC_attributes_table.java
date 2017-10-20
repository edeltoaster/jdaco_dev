package diff_compl_TCGA;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import framework.BindingDataHandler;
import framework.Utilities;

public class build_TFC_attributes_table {

	static String local_results = "/Users/tho/Dropbox/Work/projects/cancer_complexome/data/diffcompl_results_95_5_-25-25/";
	static String local_output_folder = "/Users/tho/Desktop/";
	
	public static Set<HashSet<String>> readVariants(String results_path) {
		Set<HashSet<String>> TFCs = new HashSet<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(results_path, "sign.txt.gz")) {
			
			// only read ALL TFCs
			if (!f.getName().endsWith("_compl_vsign.txt.gz") && !f.getName().endsWith("_tfc_sign.txt.gz"))
				continue;
			
			for (String s:Utilities.readFile(f.getAbsolutePath())) {
				String[] spl = s.split("\\s+");
				
				// split header
				if (spl[0].endsWith("seed_comb"))
					continue;
				
				HashSet<String> tfc = new HashSet<>(Arrays.asList(spl[0].split("/")));
				TFCs.add(tfc);
			}
		}
		
		return TFCs;
	}
	
	
	public static void main(String[] args) {
		// determine relevant TFCs and TFs
		Set<HashSet<String>> TFCs = readVariants(local_results);
		Set<String> TFs = new HashSet<>();
		TFCs.stream().forEach(tfc -> TFs.addAll(tfc));
		
		System.out.println("DiffTFC/DiffCompl Results contained " + TFCs.size() + " TFCs of " + TFs.size() + " TFs.");
		
		definitions.printBindingDataParameters();
		
		// determine all tfc-targets
		BindingDataHandler bdh = new BindingDataHandler(definitions.binding_data);
		List<String> output_to_write = new LinkedList<>();
		for (HashSet<String> tfc:TFCs) {
			Set<String> targets = bdh.getAdjacencyPossibilities(tfc, definitions.d_min, definitions.d_max, false);
			if (targets.size() > 0) {
				String tfc_string = String.join("/", tfc);
				String target_string = String.join(",", targets);
				output_to_write.add(tfc_string + " " + target_string); 
			}
		}
		
		// write all targets
		Utilities.writeEntries(bdh.getTargetsToTFsMap().keySet(), local_output_folder + "bdh_targets.txt.gz");
		Utilities.writeEntries(output_to_write, local_output_folder + "tfc_targets.txt.gz");
		
	}
}
