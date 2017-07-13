package mixed;

import framework.DiffComplexDetector;
import framework.DiffComplexDetector.TFEnrichment;

public class fast_tests {
	
	public static void main(String[] args) {
		DiffComplexDetector dcd = new DiffComplexDetector("raw_pvalues_medians.txt.gz", 0.05, 64);
		TFEnrichment tf_enrich = dcd.calculateTFEnrichment(0.05, 2000);
	}
}
