package eval_old_mixed;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;
// robustness test:
// define TS-proteins as those in AT MOST 2 tissues / vs only in one
public class compareTissuePPINs2 {

	public static void computeForNetwork(String network, String results_out_folder) {
		/**
		 * Compute core-network and builder
		 */
		System.out.println("Reading original network ...");
		PPIN core_abundance_network = new PPIN("mixed_data/"+network);
		NetworkBuilder builder = new NetworkBuilder(core_abundance_network);
		int[] sizes = core_abundance_network.getSizes();
		System.out.println("Sizes: "+ sizes[0]+ "/" + sizes[1]);
		
		Set<String> core_abundance_proteins = new HashSet<String>(core_abundance_network.getProteins());
		
		// make output folders
		new File(results_out_folder).mkdir();
		new File(results_out_folder+"proteins").mkdir();
				
		Utilities.writeEntries(core_abundance_proteins, results_out_folder+network.split("\\.")[0]+"_proteins.txt");
		Map<String, PPIN[]> sample_mapping = new HashMap<String, PPIN[]>();
		
		Map<String, Set<String>> ab_protein_to_tissue = new HashMap<String, Set<String>>();
		Map<String, Set<String>> sp_protein_to_tissue = new HashMap<String, Set<String>>();
		
		// read and map all files
		System.out.println("Reading and processing samples ...");
		for (File file:Utilities.getAllMatchingFilesInSubfolders("bodymap", "isoforms.fpkm_tracking.gz")) {
			
			String sample = file.getParentFile().getName();
			String sample_path = file.getParentFile().getAbsolutePath();
			
			System.out.println("Reading "+sample);
			PPIN abundance_network = builder.constructAssociatedNetworksFromGeneAbundance(TranscriptAbundanceReader.getGeneAbundanceFromCufflinksFPKM(sample_path+"/isoforms.fpkm_tracking.gz", 1.0)).getPPIN();
			PPIN splice_network = builder.constructAssociatedNetworksFromTranscriptAbundance(TranscriptAbundanceReader.readCufflinksTranscriptsFPKM(sample_path+"/isoforms.fpkm_tracking.gz", 1.0)).getPPIN();
			
			sample_mapping.put(sample, new PPIN[]{abundance_network, splice_network});
		}
		
		System.out.println("Count and map ...");
		for (String sample:sample_mapping.keySet()) {
			System.out.println("Processing: "+sample);
			PPIN abundance_network = sample_mapping.get(sample)[0];
			PPIN splice_network = sample_mapping.get(sample)[1];
			
			// count protein occurances
			// abundance
			for (String protein:abundance_network.getProteins()) {
				if (!ab_protein_to_tissue.containsKey(protein))
					ab_protein_to_tissue.put(protein, new HashSet<String>());
				ab_protein_to_tissue.get(protein).add(sample);
			}
			// splice
			for (String protein:splice_network.getProteins()) {
				if (!sp_protein_to_tissue.containsKey(protein))
					sp_protein_to_tissue.put(protein, new HashSet<String>());
				sp_protein_to_tissue.get(protein).add(sample);
			}
		}
		
		/**
		 * Compute tissue-specific subnetworks
		 */
		System.out.println("Computing tissue-specific networks ...");
		for (String sample:sample_mapping.keySet()) {
			System.out.println("Processing: "+sample);
			PPIN abundance_network = sample_mapping.get(sample)[0];
			PPIN splice_network = sample_mapping.get(sample)[1];
			
			// check for all proteins if there and at most in two
			Set<String> TS_proteins = new HashSet<String>();
			for (String protein:abundance_network.getProteins()) {
				if (ab_protein_to_tissue.get(protein).size() > 2)
					continue;
				if (ab_protein_to_tissue.get(protein).contains(sample))
					TS_proteins.add(protein);
			}
			Utilities.writeEntries(TS_proteins, results_out_folder+"proteins/"+sample+"_abundance_proteins.txt");
			
			TS_proteins.clear();
			for (String protein:splice_network.getProteins()) {
				if (sp_protein_to_tissue.get(protein).size() > 2)
					continue;
				if (sp_protein_to_tissue.get(protein).contains(sample))
					TS_proteins.add(protein);
			}
			Utilities.writeEntries(TS_proteins, results_out_folder+"proteins/"+sample+"_splice_proteins.txt");
		}
	}
	
	public static void main(String[] args) {

	}
}
