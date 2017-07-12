package mixed;

import java.util.Map;

import framework.DataQuery;
import framework.DiffComplexDetector;

public class fast_tests {
	
	public static void main(String[] args) {
		DiffComplexDetector dcd = new DiffComplexDetector("raw_pvalues_medians.txt.gz", 0.05, 64);
		Map<String, Double> sign_TFs = dcd.determineGSEAStyle(0.05, 5000);
		System.out.println(sign_TFs.size());
		for (String TF:sign_TFs.keySet()) {
			System.out.println(DataQuery.getHGNCNameFromProtein(TF) + " " + TF + " " + sign_TFs.get(TF));
		}
	}
}
