package old_mixed;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

public class compareTissuePPINs {

	public static void computeForNetwork(String network, String results_out_folder) {
		/**
		 * Compute core-network and builder
		 */
		System.out.println("Reading original network ...");
		PPIN core_abundance_network = new PPIN("mixed_data/"+network);
		PPIN core_splice_network = new PPIN("mixed_data/"+network);
		NetworkBuilder builder = new NetworkBuilder(core_abundance_network);
		int[] sizes = core_abundance_network.getSizes();
		System.out.println("Sizes: "+ sizes[0]+ "/" + sizes[1]);
		
		Set<String> core_abundance_proteins = new HashSet<String>(core_abundance_network.getProteins());
		Set<String> core_splice_proteins = new HashSet<String>(core_abundance_network.getProteins());
		
		// make output folders
		new File(results_out_folder).mkdir();
		new File(results_out_folder+"networks").mkdir();
		new File(results_out_folder+"proteins").mkdir();
		
		// write unpruned
		core_abundance_network.writePPIN(results_out_folder + "unpruned.tsv");
		
		Utilities.writeEntries(core_abundance_proteins, results_out_folder+network.split("\\.")[0]+"_proteins.txt");
		Map<String, PPIN[]> sample_mapping = new HashMap<String, PPIN[]>();
		// read and map all files
		System.out.println("Reading and processing samples ...");
		for (File file:Utilities.getAllMatchingFilesInSubfolders("bodymap", "isoforms.fpkm_tracking.gz")) {
			
			String sample = file.getParentFile().getName();
			String sample_path = file.getParentFile().getAbsolutePath();
			
			System.out.println("Reading "+sample);
			PPIN abundance_network = builder.constructAssociatedNetworksFromGeneAbundance(TranscriptAbundanceReader.getGeneAbundanceFromCufflinksFPKM(sample_path+"/isoforms.fpkm_tracking.gz", 1.0), true).getPPIN();
			PPIN splice_network = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readCufflinksTranscriptsFPKM(sample_path+"/isoforms.fpkm_tracking.gz", 1.0), true).getPPIN();
			
			sample_mapping.put(sample, new PPIN[]{abundance_network, splice_network});
		}
		
		System.out.println("Computing core-network ...");
		for (String sample:sample_mapping.keySet()) {
			System.out.println("Processing: "+sample);
			PPIN abundance_network = sample_mapping.get(sample)[0];
			PPIN splice_network = sample_mapping.get(sample)[1];
			
			core_abundance_network = core_abundance_network.retainAll(abundance_network);
			core_abundance_proteins.retainAll(abundance_network.getProteins());
			
			core_splice_network = core_splice_network.retainAll(splice_network);
			core_splice_proteins.retainAll(splice_network.getProteins());
			
		}
		
		// write output
		core_abundance_network.writePPIN(results_out_folder+"networks/core_abundance_ppin.tsv");
		Utilities.writeEntries(core_abundance_proteins, results_out_folder+"proteins/core_abundance_proteins.txt");
		
		core_splice_network.writePPIN(results_out_folder+"networks/core_splice_ppin.tsv");
		Utilities.writeEntries(core_splice_proteins, results_out_folder+"proteins/core_splice_proteins.txt");
		
		/**
		 * Compute tissue-specific subnetworks
		 */
		System.out.println("Computing tissue-specific networks ...");
		for (String sample:sample_mapping.keySet()) {
			System.out.println("Processing: "+sample);
			PPIN abundance_network = sample_mapping.get(sample)[0];
			PPIN splice_network = sample_mapping.get(sample)[1];
			Set<String> abundance_proteins = new HashSet<String>(abundance_network.getProteins());
			Set<String> splice_proteins = new HashSet<String>(splice_network.getProteins());
			abundance_network.writePPIN(results_out_folder+"networks/"+sample+"_abundance_ppin.tsv");
			splice_network.writePPIN(results_out_folder+"networks/"+sample+"_splice_ppin.tsv");
			
			// remove all others
			for (String other_sample:sample_mapping.keySet()) {
				if (sample.equals(other_sample))
					continue;
				
				PPIN other_abundance_network = sample_mapping.get(other_sample)[0];
				PPIN other_splice_network = sample_mapping.get(other_sample)[1];
				
				abundance_network = abundance_network.removeAll(other_abundance_network);
				abundance_proteins.removeAll(other_abundance_network.getProteins());
				
				splice_network = splice_network.removeAll(other_splice_network);
				splice_proteins.removeAll(other_splice_network.getProteins());
			}
			
			// write output
			abundance_network.writePPIN(results_out_folder+"networks/"+sample+"_abundance_ppin_restricted.tsv");
			Utilities.writeEntries(abundance_proteins, results_out_folder+"proteins/"+sample+"_abundance_proteins.txt");
			Utilities.writeEntries(abundance_network.getProteins(), results_out_folder+"proteins/"+sample+"_abundance_IA_proteins.txt");
			
			splice_network.writePPIN(results_out_folder+"networks/"+sample+"_splice_ppin_restricted.tsv");
			Utilities.writeEntries(splice_proteins, results_out_folder+"proteins/"+sample+"_splice_proteins.txt");
			Utilities.writeEntries(splice_network.getProteins(), results_out_folder+"proteins/"+sample+"_splice_IA_proteins.txt");
		}
	}
	
	public static void main(String[] args) {
		System.out.println("human merged");
		computeForNetwork("human_merged.tsv.gz", "/Users/tho/Desktop/merged/");
	}
}
