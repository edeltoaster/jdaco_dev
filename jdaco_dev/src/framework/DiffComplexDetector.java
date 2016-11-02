package framework;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

/**
 * Class for differential DACO complexes
 * @author Thorsten Will
 */
public class DiffComplexDetector {
	
	// given data
	private final Map<String, QuantDACOResultSet> group1;
	private final Map<String, QuantDACOResultSet> group2;
	private final double FDR;
	
	// processed / determined data
	private final Set<HashSet<String>> seed_combination_variant = new HashSet<>();
	private final Map<HashSet<String>, LinkedList<Double>> group1_abundances;
	private final Map<HashSet<String>, LinkedList<Double>> group2_abundances;
	
	// differential analysis results
	private final Map<HashSet<String>, Double> significance_variants_pvalues;
	private final List<HashSet<String>> significance_sorted_variants;
	
	public DiffComplexDetector(Map<String, QuantDACOResultSet> group1, Map<String, QuantDACOResultSet> group2, double FDR) {
		this.group1 = group1;
		this.group2 = group2;
		this.FDR = FDR;
		
		// determine potentially relevent seed variant combinations
		for (QuantDACOResultSet qdr:group1.values())
			this.seed_combination_variant.addAll(qdr.getSeedToComplexMap().keySet());
		for (QuantDACOResultSet qdr:group2.values())
			this.seed_combination_variant.addAll(qdr.getSeedToComplexMap().keySet());
		
		// determine abundance values
		this.group1_abundances = this.determineAbundanceOfSeedVariantsComplexes(group1);
		this.group2_abundances = this.determineAbundanceOfSeedVariantsComplexes(group2);
		
		// determine differential abundance, apply multiple hypothesis correction and filter to significant seed variants
		this.significance_variants_pvalues = this.determinePValues();
		// sort from most to least significant
		this.significance_sorted_variants = new ArrayList<>(this.significance_variants_pvalues.keySet());
		this.significance_sorted_variants.sort( (v1, v2) -> this.significance_variants_pvalues.get(v1).compareTo(this.significance_variants_pvalues.get(v2)));
	}
	
	/**
	 * Determines the abundance values of each seed variant found with the overall abundance of complexes containing it in the samples of the groups
	 * @param group
	 * @return map of seed variant and abundance in each sample
	 */
	private Map<HashSet<String>, LinkedList<Double>> determineAbundanceOfSeedVariantsComplexes(Map<String, QuantDACOResultSet> group) {
		Map<HashSet<String>, LinkedList<Double>> group_abundances = new HashMap<>();
		
		for (String sample:group.keySet()) {
			Map<HashSet<String>, Double> sample_abundances = group.get(sample).getAbundanceOfSeedVariantsComplexes();
			for (HashSet<String> variant:this.seed_combination_variant) {
				double abundance = sample_abundances.getOrDefault(variant, 0.0);
				if (!group_abundances.containsKey(variant))
					group_abundances.put(variant, new LinkedList<Double>());
				group_abundances.get(variant).add(abundance);
			}
		}
		
		return group_abundances;
	}
	
	private Map<HashSet<String>, Double> determinePValues() {
		MannWhitneyUTest mwu = new MannWhitneyUTest();
		Map<HashSet<String>, Double> test_results = new HashMap<>();
		for (HashSet<String> variant:this.seed_combination_variant) {
			double pm = mwu.mannWhitneyUTest(Utilities.getDoubleArray(this.group1_abundances.get(variant)), Utilities.getDoubleArray(this.group2_abundances.get(variant)));
			test_results.put(variant, pm);
		}
		
		return Utilities.convertRawPValuesToBHFDR(test_results, this.FDR);
	}

	public void printResults() {
		for (HashSet<String> variant: this.significance_sorted_variants) {
			System.out.println(variant + " " + this.significance_variants_pvalues.get(variant));
		}
	}

	// TODO: think about what is helpful, probably inherit classes that are more specific, like TF complexes
	
	public Map<String, QuantDACOResultSet> getGroup1() {
		return group1;
	}

	public Map<String, QuantDACOResultSet> getGroup2() {
		return group2;
	}

	public double getFDR() {
		return FDR;
	}

	public Set<HashSet<String>> getSeedCombinationVariant() {
		return seed_combination_variant;
	}

	public Map<HashSet<String>, LinkedList<Double>> getGroup1Abundances() {
		return group1_abundances;
	}

	public Map<HashSet<String>, LinkedList<Double>> getGroup2Abundances() {
		return group2_abundances;
	}

	public Map<HashSet<String>, Double> getSignificanceVariantsPValues() {
		return significance_variants_pvalues;
	}

	public List<HashSet<String>> getSignificanceSortedVariants() {
		return significance_sorted_variants;
	}
	
}
