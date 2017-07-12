package mixed;

import java.util.Map;

import framework.DataQuery;
import framework.DiffComplexDetector;

public class fast_tests {
	
	public static void main(String[] args) {
		DiffComplexDetector dcd = new DiffComplexDetector("raw_pvalues_medians.txt.gz", 0.05, 64);
		Map<String, Double> sign_TFs = dcd.determineGSEAStyle(0.05, 5000);
		System.out.println(sign_TFs.size());
		for (String TFdir:sign_TFs.keySet()) {
			String tf = TFdir.substring(1);
			String dir = TFdir.substring(0, 1);
			System.out.println(dir + " " + DataQuery.getHGNCNameFromProtein(tf) + " " + tf + " " + sign_TFs.get(TFdir));
		}
	}
}
