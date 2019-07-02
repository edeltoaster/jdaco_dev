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
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

import framework.BindingDataHandler;
import framework.DataQuery;
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
		Set<String> pluri_upTFs = new HashSet<>();
		for (String s : Utilities.readFile(diff_compl_file)) {
			if (s.startsWith("(sub)"))
					continue;
			String[] spl = s.split(" ");
			Set<String> tfc = new HashSet<String>(Arrays.asList(spl[0].split("/")));
			
			if (spl[1].equals("+")) {
				pluri_upCs.add(tfc);
				pluri_upTFs.addAll(tfc);
			} else {
				pluri_downCs.add(tfc);
			}
			
		}
		pluri_upTFs.retainAll(TFs);
		
		System.out.println(pluri_upCs.size() + " sign upregulated pluri complexes");
		System.out.println(pluri_downCs.size() + " sign downregulated pluri complexes");
		System.out.println(pluri_upTFs.size() + " TFs in upregulated");

		Set<String> Hac_proteins = new HashSet<>();
		System.out.println("P300/CBP only");
		Hac_proteins.add("Q09472"); // EP300
		Hac_proteins.add("Q92793"); // CBP
		String P300 = "Q09472"; // EP300
		String CBP ="Q92793"; // CBP
		
		boolean at_least_two_TFs = false;
		System.out.println("at least two TFs: "+ at_least_two_TFs);
		int d_min = -20;
		int d_max = 10;
		System.out.println("d_min / d_max: " + d_min + "-" + d_max);
		int repetitions = 10000;
		System.out.println("#repetitions: " + repetitions);
		
		List<Set<String>> pluri_Hac_upCs = new LinkedList<>();
		List<Set<String>> pluri_P300_upCs = new LinkedList<>();
		List<Set<String>> pluri_CBP_upCs = new LinkedList<>();
		List<Set<String>> pluri_both_upCs = new LinkedList<>();
		for (Set<String> complex : pluri_upCs) {
			Set<String> ov = new HashSet<>(complex);
			ov.retainAll(Hac_proteins);
			if (ov.size() > 0) {
				pluri_Hac_upCs.add(complex);
				if (ov.contains(P300) && ov.contains(CBP))
					pluri_both_upCs.add(complex);
				if (ov.contains(P300))
					pluri_P300_upCs.add(complex);
				else
					pluri_CBP_upCs.add(complex);
			}
		}
		System.out.println(pluri_Hac_upCs.size() + " Hac complexes in up");
		System.out.println(pluri_P300_upCs.size() + " P300 complexes in up");
		System.out.println(pluri_CBP_upCs.size() + " CBP complexes in up");
		System.out.println(pluri_both_upCs.size() + " P300+CBP complexes in up");
		
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
		double samples = 0;
		double score = 0;
		double hitscore = 0;
		for (Set<String> complex : pluri_Hac_upCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			
			if (at_least_two_TFs && tfc.size() < 2)
				continue;
			samples++;
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, d_min, d_max, false);
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
		n /= samples;
		score /= samples;
		System.out.println(samples + "  " + n + "  " + score + "  " + hitscore);
		double ref_n = n;
		double ref_score = score;
		double ref_hitscore = hitscore;
		
		
		System.out.println();
		System.out.println("----");
		n = 0;
		samples = 0;
		score = 0;
		hitscore = 0;
		for (Set<String> complex : pluri_Hac_downCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			
			if (at_least_two_TFs && tfc.size() < 2)
				continue;
			samples++;
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, d_min, d_max, false);
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
		n /= samples;
		score /= samples;
		System.out.println(samples + "  " + n + "  " + score + "  " + hitscore);

		
		// prepare sampling distribution
		List<Integer> distr = new ArrayList<>();
		for (Set<String> complex : pluri_Hac_upCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			if (at_least_two_TFs && tfc.size() < 2)
				continue;
			distr.add(tfc.size());
		}
		
		System.out.println();
		System.out.println("all+-+-+-+-+-+-+-+-+-");
		double tot_n = 0;
		double tot_score = 0;
		double tot_hitscore = 0;
		double p_n = 0;
		double p_score = 0;
		double p_hitscore = 0;
		List<String> tfl = new ArrayList<>(TFs);
		for (int rep = 0; rep < repetitions; rep++) {
			n = 0;
			score = 0;
			hitscore = 0;
			for (int s : distr) {
				Set<String> tfc = new HashSet<>();
				Collections.shuffle(tfl);
				for (int i = 0; i < s ; i++)
					tfc.add(tfl.get(i));
				
				Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, d_min, d_max, false);
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
			n /= distr.size();
			score /= distr.size();
			tot_n += n;
			tot_score += score;
			tot_hitscore += hitscore;
			
			// check if better than ref
			if (n > ref_n)
				p_n++;
			if (score > ref_score)
				p_score++;
			if(hitscore > ref_hitscore)
				p_hitscore++;
		}
		tot_n /= repetitions;
		tot_score /= repetitions;
		tot_hitscore /= repetitions;
		System.out.println(distr.size() + "  " + tot_n + "  " + tot_score + "  " + tot_hitscore);
		p_n /= repetitions;
		p_score /= repetitions;
		p_hitscore /= repetitions;
		System.out.println("     " + "  " + p_n + "  " + p_score + "  " + p_hitscore);
		
		System.out.println();
		System.out.println("upTFs only+-+-+-+-+-+-+-+-+-");
		tot_n = 0;
		tot_score = 0;
		tot_hitscore = 0;
		p_n = 0;
		p_score = 0;
		p_hitscore = 0;
		tfl = new ArrayList<>(pluri_upTFs);
		for (int rep = 0; rep < repetitions; rep++) {
			n = 0;
			score = 0;
			hitscore = 0;
			for (int s : distr) {
				Set<String> tfc = new HashSet<>();
				Collections.shuffle(tfl);
				for (int i = 0; i < s ; i++)
					tfc.add(tfl.get(i));
				
				Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, d_min, d_max, false);
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
			n /= distr.size();
			score /= distr.size();
			tot_n += n;
			tot_score += score;
			tot_hitscore += hitscore;
			
			// check if better than ref
			if (n > ref_n)
				p_n++;
			if (score > ref_score)
				p_score++;
			if(hitscore > ref_hitscore)
				p_hitscore++;
		}
		tot_n /= repetitions;
		tot_score /= repetitions;
		tot_hitscore /= repetitions;
		System.out.println(distr.size() + "  " + tot_n + "  " + tot_score + "  " + tot_hitscore);
		p_n /= repetitions;
		p_score /= repetitions;
		p_hitscore /= repetitions;
		System.out.println("     " + "  " + p_n + "  " + p_score + "  " + p_hitscore);
		
		
		/*
		 *  detailed output
		 */
		
		Map<String, String> up_to_gene = DataQuery.getUniprotToGeneNameMap(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		Set<String> P300_targets = new HashSet<>();
		Set<String> CBP_targets = new HashSet<>();
		Set<String> all_targets = new HashSet<>();
		Map<String, List<String>> region_P300_map = new HashMap<>();
		Map<String, List<String>> region_CBP_map = new HashMap<>();
		
		System.out.println();
		System.out.println("P300++++");
		n = 0;
		samples = 0;
		score = 0;
		hitscore = 0;
		for (Set<String> complex : pluri_P300_upCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			
			if (at_least_two_TFs && tfc.size() < 2)
				continue;
			samples++;
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, d_min, d_max, false);
			P300_targets.addAll(target_regions);
			String complex_string = String.join("/", complex);
			String complex_genes = String.join("/", complex.stream().map(p -> up_to_gene.getOrDefault(p, p)).collect(Collectors.toSet()));
			String targets = String.join(",", target_regions);
			
			System.out.println(complex_genes + " " + complex_string + " " +  target_regions.size() + " " + targets);
			for (String target : target_regions) {
				if (!region_P300_map.containsKey(target))
					region_P300_map.put(target, new LinkedList<>());
				region_P300_map.get(target).add(complex_genes);
			}
			
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
		n /= samples;
		score /= samples;
		System.out.println(samples + "  " + n + "  " + score + "  " + hitscore);
		
		System.out.println();
		System.out.println("CBP++++");
		n = 0;
		samples = 0;
		score = 0;
		hitscore = 0;
		for (Set<String> complex : pluri_CBP_upCs) {
			Set<String> tfc = new HashSet<>(complex);
			tfc.retainAll(TFs);
			
			if (at_least_two_TFs && tfc.size() < 2)
				continue;
			samples++;
			Set<String> target_regions = bdh.getAdjacencyPossibilities(tfc, d_min, d_max, false);
			CBP_targets.addAll(target_regions);
			String complex_string = String.join("/", complex);
			String complex_genes = String.join("/", complex.stream().map(p -> up_to_gene.getOrDefault(p, p)).collect(Collectors.toSet()));
			String targets = String.join(",", target_regions);
			
			System.out.println(complex_genes + " " + complex_string + " " +  target_regions.size() + " " + targets);
			for (String target : target_regions) {
				if (!region_CBP_map.containsKey(target))
					region_CBP_map.put(target, new LinkedList<>());
				region_CBP_map.get(target).add(complex_genes);
			}
			
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
		n /= samples;
		score /= samples;
		System.out.println(samples + "  " + n + "  " + score + "  " + hitscore);
		
		all_targets.addAll(P300_targets);
		all_targets.addAll(CBP_targets);
		double pluri_region_size = pluri_regions.size();
		System.out.println();
		System.out.println("Enhancer targets:");
		System.out.println("P300 targets: " + P300_targets.size() + " covering " + P300_targets.size()/pluri_region_size);
		System.out.println("CBP targets: " + CBP_targets.size() + " covering " + CBP_targets.size()/pluri_region_size);
		System.out.println("joined targets: " + all_targets.size() + " covering " + all_targets.size()/pluri_region_size);
		
		List<String> to_write = new LinkedList<>();
		for (String target : all_targets) {
			String P300_size = "0";
			String P300_complexes = "/";
			if (region_P300_map.containsKey(target)) {
				P300_size = Integer.toString(region_P300_map.get(target).size());
				P300_complexes = String.join(",", region_P300_map.get(target));
			}
			String CBP_size = "0";
			String CBP_complexes = "/";
			if (region_CBP_map.containsKey(target)) {
				CBP_size = Integer.toString(region_CBP_map.get(target).size());
				CBP_complexes = String.join(",", region_CBP_map.get(target));
			}
			
			to_write.add(target + " " + P300_size + " " + CBP_size + " " + P300_complexes + " " + CBP_complexes);
		}
		Utilities.writeEntries(to_write, "/Users/tho/Desktop/region_complex_mapping.txt");
	}

}
