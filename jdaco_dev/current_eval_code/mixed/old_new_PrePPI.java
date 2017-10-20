package mixed;

import framework.PPIN;
import framework.StrPair;

public class old_new_PrePPI {
	
	public static void main(String[] args) {
		PPIN old_preppi = new PPIN("mixed_data/human_preppi_initial_release.tsv.gz");
		System.out.println(old_preppi.getSizesStr());
		PPIN new_preppi = new PPIN("mixed_data/human_PrePPI_17_01_17.txt.gz");
		System.out.println(new_preppi.getSizesStr());
		PPIN new_preppi_hc = new PPIN("mixed_data/human_PrePPI_17_01_17_hc.txt.gz");
		System.out.println(new_preppi_hc.getSizesStr());
		
		StrPair oct4_sox2 = new StrPair("Q01860", "P48431");
		StrPair nanog_oct4 = new StrPair("Q9H9S0", "Q01860");
		StrPair nanog_sox2 = new StrPair("Q9H9S0", "P48431");
		
		System.out.println("oct4/sox2 old:" + " " + old_preppi.getWeights().get(oct4_sox2));
		System.out.println("oct4/sox2 new:" + " " + new_preppi.getWeights().get(oct4_sox2));
		System.out.println("oct4/sox2 new hc:" + " " + new_preppi_hc.getWeights().get(oct4_sox2));
		
		System.out.println("nanog_oct4 old:" + " " + old_preppi.getWeights().get(nanog_oct4));
		System.out.println("nanog_oct4 new:" + " " + new_preppi.getWeights().get(nanog_oct4));
		System.out.println("nanog_oct4 new hc:" + " " + new_preppi_hc.getWeights().get(nanog_oct4));
		
		System.out.println("nanog_sox2 old:" + " " + old_preppi.getWeights().get(nanog_sox2));
		System.out.println("nanog_sox2 new:" + " " + new_preppi.getWeights().get(nanog_sox2));
		System.out.println("nanog_sox2 new hc:" + " " + new_preppi_hc.getWeights().get(nanog_sox2));
	}
}
