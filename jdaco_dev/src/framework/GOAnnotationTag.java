package framework;

import java.util.HashSet;
import java.util.Set;

public class GOAnnotationTag {
	private final String taxon;
	private final String tag_name;
	private final Set<String> positive_GO_terms;
	private final Set<String> negative_GO_terms;
	private final boolean include_IEA;
	private final Set<String> positive_proteins = new HashSet<>();
	private final Set<String> negative_proteins = new HashSet<>();
	private final Set<String> mixed_proteins = new HashSet<>();
	
	public GOAnnotationTag(String taxon, String tag_name, Set<String> positive_GO_terms, Set<String> negative_GO_terms, boolean include_IEA) {
		this.taxon = taxon;
		this.tag_name = tag_name;
		this.positive_GO_terms = positive_GO_terms;
		this.negative_GO_terms = negative_GO_terms;
		this.include_IEA = include_IEA;
		
		this.retrieveAndProcessData();
	}
	
	/**
	 * Retrieve annotation data from QuickGO
	 */
	private void retrieveAndProcessData() {
		for (String GO_term:this.positive_GO_terms)
			this.positive_proteins.addAll(DataQuery.getProteinsWithGO(GO_term, this.taxon, this.include_IEA, false));
		for (String GO_term:this.negative_GO_terms)
			this.negative_proteins.addAll(DataQuery.getProteinsWithGO(GO_term, this.taxon, this.include_IEA, false));
		
		this.mixed_proteins.addAll(positive_proteins);
		this.mixed_proteins.retainAll(negative_proteins);
		
		this.positive_proteins.removeAll(mixed_proteins);
		this.negative_proteins.removeAll(mixed_proteins);
	}

	/**
	 * Return taxon of organism associated with retrieved proteins
	 * @return
	 */
	public String getTaxon() {
		return taxon;
	}

	/**
	 * Return the description string of the annotation tag
	 * @return
	 */
	public String getTagName() {
		return tag_name;
	}

	/**
	 * Returns associated positive GO terms (+)
	 * @return
	 */
	public Set<String> getPositiveGOTerms() {
		return positive_GO_terms;
	}

	/**
	 * Returns associated negative GO terms (-)
	 * @return
	 */
	public Set<String> getNegativeGOTerms() {
		return negative_GO_terms;
	}

	/**
	 * Returns if proteins annotated by purely computationally inferred evidence were included  
	 * @return
	 */
	public boolean IncludesIEA() {
		return include_IEA;
	}

	/**
	 * Returns all strictly positive associated proteins
	 * @return
	 */
	public Set<String> getPositiveProteins() {
		return positive_proteins;
	}

	/**
	 * Returns all strictly negative associated proteins
	 * @return
	 */
	public Set<String> getNegativeProteins() {
		return negative_proteins;
	}

	/**
	 * Returns all proteins with both positive and negative annotation regarding this tag
	 * @return
	 */
	public Set<String> getMixedProteins() {
		return mixed_proteins;
	}

}
