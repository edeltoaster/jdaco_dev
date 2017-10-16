package mixed;

import framework.DataQuery;

public class fast_tests {
	
	public static void main(String[] args) {
		System.out.println(DataQuery.getProteinsWithGO("GO:0006915", "9606").size());
	}
}
