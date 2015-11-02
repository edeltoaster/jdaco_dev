package mixed;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class TFmir_adds {
	
	public static void main(String[] args) {
		
		List<String> output = new LinkedList<>();
		for (File f:Utilities.getAllSuffixMatchingFilesInSubfolders("/Users/tho/Dropbox/Work/projects/TFmir_adds/", "gtf.gz")) {
			String sample = f.getName().split("RnaSeq")[1].split("Cell")[0].split("GeneEns")[0];
			System.out.println("Processing " + sample);
			
			Map<String, Float> abundances = TranscriptAbundanceReader.readCSHLData(f.getAbsolutePath());
			
			if (abundances.get("no_samples") < 2) {
				System.out.println("only 1 sample -> withdraw");
				continue;
			}
			abundances.remove("no_samples");
			
			output.add(sample + ":" + String.join(",", abundances.keySet()));
		}
		
		Utilities.writeEntries(output, "/Users/tho/Dropbox/Work/projects/TFmir_adds/mouse_tissues_ENS.txt.gz");
	}
}
