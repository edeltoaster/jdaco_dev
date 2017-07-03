package stem_cell_complexeomes;

import framework.RegulatoryNetwork;

public class check_regnets {
	
	// testing and checking some stuff first
	public static void main(String[] args) {
		RegulatoryNetwork regnet = new RegulatoryNetwork("/Users/tho/Desktop/95/plurisub_regnet_only.txt");
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/1.txt");
		regnet.pruneToLargestSCCs();
		System.out.println(regnet.getSizesStr());
		regnet.writeRegulatoryNetwork("/Users/tho/Desktop/2.txt");
	}
}
