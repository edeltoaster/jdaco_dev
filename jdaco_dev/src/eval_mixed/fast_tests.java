package eval_mixed;

import framework.DataQuery;
import framework.PPIN;

public class fast_tests {
	
	public static void main(String[] args) {
		PPIN bg = new PPIN("mixed_data/human_BG.tsv");
		System.out.println("bg:" + bg.getSizesStr());
		PPIN hprd = new PPIN("mixed_data/human_hprd_r9.tsv");
		System.out.println("hprd:" + hprd.getSizesStr());
		PPIN intact = DataQuery.getIntActNetwork("9606");
		System.out.println("intact:" + intact.getSizesStr());
		
		intact.mergeAllIAs(bg).mergeAllIAs(hprd).writePPIN("mixed_data/human_merged.tsv");
		
		intact = new PPIN("mixed_data/human_merged.tsv");
		System.out.println("merged:" + intact.getSizesStr());
	}
}
