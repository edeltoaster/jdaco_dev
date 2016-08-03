package framework;

import java.util.LinkedList;
import java.util.List;

public class GOAnnotator {
	private final String taxon;
	List<GOAnnotationTag> tags = new LinkedList<>();
	
	/*
	 * Constructor to load definition data from file
	 */
	public GOAnnotator(String taxon, String file) {
		this.taxon = taxon;
		Utilities.readFile(file);
	}
	
	/*
	 * Constructor to load old GOAnnotator object from file
	 */
	public GOAnnotator(String file) {
		Utilities.readFile(file);
		this.taxon = "";
	}
	
	/**
	 * Return taxon of organism associated with retrieved data
	 * @return
	 */
	public String getTaxon() {
		return taxon;
	}
	
}
