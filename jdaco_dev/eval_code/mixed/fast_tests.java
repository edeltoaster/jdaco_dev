package mixed;

import framework.TranscriptAbundanceReader;

public class fast_tests {
	
	public static void main(String[] args) {
		String file = "/Users/tho/Downloads/results/SRR493366/kallisto/abundance.h5";
		System.out.println(TranscriptAbundanceReader.readKallistoH5(file, 0.0, false));
	}
}
