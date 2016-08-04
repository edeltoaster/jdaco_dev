package framework;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
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
	
	/*
	 * Constructor if data is given within framework
	 */
	public GOAnnotationTag(String taxon, String tag_name, Set<String> positive_GO_terms, Set<String> negative_GO_terms, boolean include_IEA) {
		this.taxon = taxon;
		this.tag_name = tag_name;
		this.positive_GO_terms = positive_GO_terms;
		this.negative_GO_terms = negative_GO_terms;
		this.include_IEA = include_IEA;
		
		this.retrieveAndProcessData();
	}
	
	/*
	 * Constructor if data is given from a string (definition file)
	 */
	public GOAnnotationTag(String taxon, boolean include_IEA, String data_string) {
		this.taxon = taxon;
		this.include_IEA = include_IEA;
		
		String[] data = data_string.trim().split("\\s+");
		this.tag_name = data[0];
		
		this.positive_GO_terms = new HashSet<>();
		for (String GO_term:data[1].split(","))
			this.positive_GO_terms.add(GO_term);
		
		this.negative_GO_terms = new HashSet<>();
		for (String GO_term:data[2].split(","))
			this.negative_GO_terms.add(GO_term);
		
		this.retrieveAndProcessData();
	}
	
	/*
	 * Constructor if data is given from a string (GOAnnotator file)
	 */
	public GOAnnotationTag(String data_string) {
		String[] data = data_string.trim().split("\\s+");
		
		this.taxon = data[0];
		this.tag_name = data[1];
		
		this.positive_GO_terms = new HashSet<>();
		for (String GO_term:data[2].split(","))
			this.positive_GO_terms.add(GO_term);
		
		this.negative_GO_terms = new HashSet<>();
		for (String GO_term:data[3].split(","))
			this.negative_GO_terms.add(GO_term);
		
		if (data[4].equals("IEA"))
			this.include_IEA = true;
		else
			this.include_IEA = false;
		
		for (String protein:data[5].split(","))
			this.positive_proteins.add(protein);
		
		for (String protein:data[6].split(","))
			this.negative_proteins.add(protein);
		
		for (String protein:data[7].split(","))
			this.mixed_proteins.add(protein);
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
	 * Returns space-separated string representation of all data
	 * @return
	 */
	public String getDataString() {
		List<String> to_string = new LinkedList<>();
		to_string.add(this.taxon);
		to_string.add(this.tag_name);
		to_string.add(String.join(",", this.positive_GO_terms));
		to_string.add(String.join(",", this.negative_GO_terms));
		
		if (this.include_IEA)
			to_string.add("IEA");
		else
			to_string.add("no IEA");
		
		to_string.add(String.join(",", this.positive_proteins));
		to_string.add(String.join(",", this.negative_proteins));
		to_string.add(String.join(",", this.mixed_proteins));
		
		return String.join(" ", to_string);
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
