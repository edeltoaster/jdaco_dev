package PPIComp_hemato;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class parse_TOP2A_expression {
	
	
	static String BLUEPRINT_expr_folder = "/Users/tho/Desktop/BLUEPRINT_expr/";
	static Map<String, String> folder_type_map = new HashMap<>();
	static String top2a_tr = "ENST00000423485";
	
	public static void main(String[] args) {
		
		for (String s:Utilities.readEntryFile("/Users/tho/git/jdaco_dev/jdaco_dev/eval_code/PPIComp_hemato/cell_types.txt")) {
			if (s.startsWith("#"))
				continue;
			String[] spl = s.trim().split(" ");
			folder_type_map.put(spl[0], spl[1]);
		}
		
		List<String> to_write = new LinkedList<String>();
		to_write.add("cell type\tTOP2A expression");
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders(BLUEPRINT_expr_folder, ".rsem.tsv.gz")) {
			String path = f.getAbsolutePath();
			String[] path_split = path.split("/");
			String file_name = path_split[path_split.length-1].split("\\.")[0];
			
			String cell_type = path_split[path_split.length-2];
			
			if (!folder_type_map.containsKey(cell_type))
				continue;
			
			cell_type = folder_type_map.get(cell_type);
			
			// filter N, M, NK, CD4 to venous blood samples
			if (cell_type.equals("M") || cell_type.equals("N") || cell_type.equals("NK") || cell_type.equals("CD4"))
				if (file_name.startsWith("Cord"))
					continue;
			
			Map<String, Float> transcr_expr = TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, 0.0);
			
			to_write.add(cell_type + "\t" + transcr_expr.get(top2a_tr));
		}
		
		Utilities.writeEntries(to_write, "/Users/tho/GDrive/Work/projects/hemato_rewiring/diffnet_analysis/scripts/consequence/TOP2A_expr.txt.gz");
	}
}
