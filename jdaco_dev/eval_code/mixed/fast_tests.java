package mixed;

import framework.DataQuery;
import framework.DiffComplexDetector;
import framework.DiffComplexDetector.SPEnrichment;

public class fast_tests {
	
	public static void main(String[] args) {
		DiffComplexDetector dcd = new DiffComplexDetector("raw_pvalues_medians.txt.gz", 0.05, 64);
		SPEnrichment enrich = dcd.calculateTFEnrichment(0.05, 5000, 10);
		for (String tf:enrich.getSignificanceSortedSeedProteins()) {
			System.out.println(enrich.getSignificantSeedProteinDirections().get(tf) + " " + tf + " " + DataQuery.getHGNCNameFromProtein(tf) + " " + enrich.getSignificantSeedProteinQvalues().get(tf));
		}
	}
}
