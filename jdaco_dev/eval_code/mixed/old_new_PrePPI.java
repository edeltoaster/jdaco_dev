package mixed;

import framework.PPIN;

public class old_new_PrePPI {
	
	public static void main(String[] args) {
		PPIN old_preppi = new PPIN("mixed_data/human_preppi_initial_release.tsv.gz");
		System.out.println(old_preppi.getSizesStr());
		PPIN new_preppi = new PPIN("mixed_data/human_PrePPI_17_01_17.txt.gz");
		System.out.println(new_preppi.getSizesStr());
		PPIN new_preppi_hc = new PPIN("mixed_data/human_PrePPI_17_01_17_hc.txt.gz");
		System.out.println(new_preppi_hc.getSizesStr());
	}
}
