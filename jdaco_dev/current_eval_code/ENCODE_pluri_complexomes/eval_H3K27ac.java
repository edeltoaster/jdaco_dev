package ENCODE_pluri_complexomes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

import framework.BindingDataHandler;
import framework.Utilities;

public class eval_H3K27ac {
	static String diff_compl_file = "/Users/tho/GDrive/Work/projects/stem_cell_complexome/diff_compl_0.75/preppi_sign.txt";
	static String fimo_data = "/Users/tho/GDrive/Work/data_general/ENCODE_and_ROADMAP/H3K27ac/hocomoco_genehancer.tsv.gz";
	static Set<String> TFs = Utilities.readEntryFile("mixed_data/hocomoco_human_core_TFs_v11.txt.gz");
	static Map<String, Double> pluri = readBED("/Users/tho/GDrive/Work/data_general/ENCODE_and_ROADMAP/H3K27ac/pluri_avg.bed.gz");
	static Map<String, Double> nonpluri = readBED("/Users/tho/GDrive/Work/data_general/ENCODE_and_ROADMAP/H3K27ac/nonpluri_max.bed.gz");
	
	/**
	 * Reads BED file, returns map of regions to score.
	 * @param bed_file
	 * @return
	 */
	public static Map<String, Double> readBED(String bed_file) {
		Map<String, Double> data = new HashMap<>();
		for (String line : Utilities.readFile(bed_file)) {
			String[] spl = line.trim().split("\\t");
			String region = spl[0] + ":" + spl[1] + "-" + spl[2];
			Double score = Double.parseDouble(spl[4]);
			data.put(region, score);
		}
		return data;
	}
	
	public static String regionCheck(Set<String> regions) {
		
		double[] pluri_scores = new double[regions.size()];
		double[] nonpluri_scores = new double[regions.size()];
		int i = 0;
		double diff = 0;
		for (String region : regions) {
			pluri_scores[i] = pluri.get(region);
			nonpluri_scores[i] = nonpluri.get(region);
			diff += pluri_scores[i];
			diff -= nonpluri_scores[i];
			i++;
		}
		
		String diff_sign = "+";
		if (diff < 0)
			diff_sign = "-";
		
		// above 20 samples, test statistic W is approx. normally distributed
		boolean compute_exact_p = true;
		if (regions.size() > 20)
			compute_exact_p = false;
		
		WilcoxonSignedRankTest wsrt = new WilcoxonSignedRankTest();
		double pval = wsrt.wilcoxonSignedRankTest(pluri_scores, nonpluri_scores, compute_exact_p);
		
		return diff_sign + " " + pval;
	}
	
	public static void main(String[] args) {
		
		// read upregulated diff complexes
		List<Set<String>> pluri_upCs = new LinkedList<>();
		List<Set<String>> pluri_downCs = new LinkedList<>();
		Set<String> pluri_TFs = new HashSet<>();
		for (String s : Utilities.readFile(diff_compl_file)) {
			if (s.startsWith("(sub)"))
					continue;
			String[] spl = s.split(" ");
			Set<String> tfc = new HashSet<String>(Arrays.asList(spl[0].split("/")));
			
			if (spl[1].equals("+")) {
				pluri_upCs.add(tfc);
			} else {
				pluri_downCs.add(tfc);
			}
			pluri_TFs.addAll(tfc);
		}
		pluri_TFs.retainAll(TFs);
		
		System.out.println(pluri_upCs.size() + " sign upregulated pluri complexes");
		System.out.println(pluri_downCs.size() + " sign downregulated pluri complexes");
		System.out.println(pluri_TFs.size() + " TFs therein");
		
//		Set<String> Hac_proteins = DataQuery.getProteinsWithGO("GO:001657", "9606");
//		Set<String> Hdac_proteins = DataQuery.getProteinsWithGO("GO:0016575", "9606");
		Set<String> Hac_proteins = new HashSet<>();
		for (String s:Utilities.readFile("mixed_data/H_ac.tsv")) {
			if (s.startsWith("GENE"))
					continue;
			Hac_proteins.add(s.split("\t")[1]);
		}
		System.out.println("GO data retrieved.");
		
		List<Set<String>> pluri_Hac_upCs = new LinkedList<>();
		for (Set<String> complex : pluri_upCs) {
			Set<String> ov = new HashSet<>(complex);
			ov.retainAll(Hac_proteins);
			if (ov.size() > 0) {
				pluri_Hac_upCs.add(complex);
			}
		}
		System.out.println(pluri_Hac_upCs.size() + " Hac complexes in up");
		
		List<Set<String>> pluri_Hac_downCs = new LinkedList<>();
		for (Set<String> complex : pluri_downCs) {
			Set<String> ov = new HashSet<>(complex);
			ov.retainAll(Hac_proteins);
			if (ov.size() > 0) {
				pluri_Hac_downCs.add(complex);
			}
		}
		System.out.println(pluri_Hac_downCs.size() + " Hac complexes in down");
		
		
		// gather "pluri enhancers"
		Set<String> pluri_regions = new HashSet<>();
		for (String region : pluri.keySet()) {
			if (nonpluri.get(region) == 0.0 && pluri.get(region) > 0)
				pluri_regions.add(region);
		}
		System.out.println(pluri_regions.size() + " pluri enhancer regions");
		BindingDataHandler bdh = new BindingDataHandler(fimo_data, TFs, pluri_regions);
		System.out.println("FIMO data read.");
		System.out.println();
		
		System.out.println("++++");
		double n = 0;
		double score = 0;
		double hitscore = 0;
		for (Set<String> complex : pluri_Hac_upCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			
			if (tfc.size() < 2)
				continue;
			
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, -20, +10, false);
			//System.out.println("Checking " + complex + "/" + tfc + " targetting " + target_regions.size() + " enhancers");
			n += target_regions.size();
			if (target_regions.size() > 0) {
				double sub_score = 0;
				for (String reg : target_regions) {
					sub_score += pluri.get(reg);
					hitscore += pluri.get(reg);
				}
				sub_score /= target_regions.size();
				score += sub_score;
			}
		}
		hitscore /= n;
		n /= pluri_Hac_upCs.size();
		score /= pluri_Hac_upCs.size();
		System.out.println(n + "  " + score + "  " + hitscore);
		
		
		System.out.println();
		System.out.println("----");
		n = 0;
		score = 0;
		hitscore = 0;
		for (Set<String> complex : pluri_Hac_downCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			
			if (tfc.size() < 2)
				continue;
			
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, -20, +10, false);
			//System.out.println("Checking " + complex + "/" + tfc + " targetting " + target_regions.size() + " enhancers");
			n += target_regions.size();
			if (target_regions.size() > 0) {
				double sub_score = 0;
				for (String reg : target_regions) {
					sub_score += pluri.get(reg);
					hitscore += pluri.get(reg);
				}
				sub_score /= target_regions.size();
				score += sub_score;
			}
		}
		hitscore /= n;
		n /= pluri_Hac_downCs.size();
		score /= pluri_Hac_downCs.size();
		System.out.println(n + "  " + score + "  " + hitscore);
		
		
		System.out.println();
		System.out.println("all+-+-+-+-+-+-+-+-+-");
		n = 0;
		score = 0;
		hitscore = 0;
		List<String> tfl = new ArrayList<>(TFs);
		int sample_size = 1000;
		for (int i = 0; i < sample_size; i++) {
			Set<String> tfc = new HashSet<>();
			Collections.shuffle(tfl);
			tfc.add(tfl.get(0));
			tfc.add(tfl.get(1));
			
			if (i % 2 == 0)
				tfc.add(tfl.get(2));
			
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, -20, +10, false);
			//System.out.println("Checking " + tfc + " targetting " + target_regions.size() + " enhancers");
			n += target_regions.size();
			if (target_regions.size() > 0) {
				double sub_score = 0;
				for (String reg : target_regions) {
					sub_score += pluri.get(reg);
					hitscore += pluri.get(reg);
				}
				sub_score /= target_regions.size();
				score += sub_score;
			}
		}
		hitscore /= n;
		n /= sample_size;
		score /= sample_size;
		System.out.println(n + "  " + score + "  " + hitscore);
		
		System.out.println();
		System.out.println("TFCTFsonly+-+-+-+-+-+-+-+-+-");
		n = 0;
		score = 0;
		hitscore = 0;
		tfl = new ArrayList<>(pluri_TFs);
		for (int i = 0; i < sample_size; i++) {
			Set<String> tfc = new HashSet<>();
			Collections.shuffle(tfl);
			tfc.add(tfl.get(0));
			tfc.add(tfl.get(1));
			
			if (i % 2 == 0)
				tfc.add(tfl.get(2));
			
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, -20, +10, false);
			//System.out.println("Checking " + tfc + " targetting " + target_regions.size() + " enhancers");
			n += target_regions.size();
			if (target_regions.size() > 0) {
				double sub_score = 0;
				for (String reg : target_regions) {
					sub_score += pluri.get(reg);
					hitscore += pluri.get(reg);
				}
				sub_score /= target_regions.size();
				score += sub_score;
			}
		}
		hitscore /= n;
		n /= sample_size;
		score /= sample_size;
		System.out.println(n + "  " + score + "  " + hitscore);
	}

}
