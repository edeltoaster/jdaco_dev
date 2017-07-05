package stem_cell_complexeomes;

import java.util.HashSet;
import java.util.Set;

import framework.BindingDataHandler;
import framework.QuantDACOResultSet;
import framework.RegulatoryNetwork;

public class check_example_regnets {

	public static String data_path = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/regnet_examples/input_data/";
	public static String binding_data_path = "/Users/tho/Dropbox/Work/data_general/binding_sites/hocomoco_v10_EPD_v4_5k.txt.gz";
	public static String results_out = "/Users/tho/Dropbox/Work/projects/stem_cell_complexeome/regnet_examples/regnets/";
	
	public static void main(String[] args) {
		
		QuantDACOResultSet qdr_H1 = new QuantDACOResultSet(data_path + "H1-hESC_1_1_ENCSR000COU.csv.gz", definitions.seed_file, data_path + "H1-hESC_1_1_ENCSR000COU_major-transcripts.txt.gz");
		QuantDACOResultSet qdr_CD14 = new QuantDACOResultSet(data_path + "CD14-positive-monocyte_1_1_ENCSR000CUC.csv.gz", definitions.seed_file, data_path + "CD14-positive-monocyte_1_1_ENCSR000CUC_major-transcripts.txt.gz");
		
		Set<String> involved_tfs = new HashSet<>();
		qdr_H1.getSeedToComplexMap().keySet().stream().forEach(c -> involved_tfs.addAll(c));
		qdr_CD14.getSeedToComplexMap().keySet().stream().forEach(c -> involved_tfs.addAll(c));
		
		System.out.println("Reading binding data for " + involved_tfs.size() + " TFs.");
		BindingDataHandler bdh = new BindingDataHandler(binding_data_path, involved_tfs, 0.0001, involved_tfs);
		
		
		System.out.println("Building regnets ...");
		RegulatoryNetwork regnet_H1 = new RegulatoryNetwork(qdr_H1, bdh, definitions.d_min, definitions.d_max, 3);
		regnet_H1.writeRegulatoryNetwork(results_out + "H1_regnet.txt");
		regnet_H1.writeNodeTable(results_out + "H1_nodetable.txt");
		System.out.println("H1: " + regnet_H1.getSizesStr());
		
		RegulatoryNetwork regnet_CD14 = new RegulatoryNetwork(qdr_CD14, bdh, definitions.d_min, definitions.d_max, 3);
		regnet_CD14.writeRegulatoryNetwork(results_out + "CD14_regnet.txt");
		regnet_CD14.writeNodeTable(results_out + "CD14_nodetable.txt");
		System.out.println("CD14: " + regnet_CD14.getSizesStr());
	}
}
