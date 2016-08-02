package framework;

import java.util.LinkedList;
import java.util.List;

public class GOAnnotator {
	private final String taxon;
	List<GOAnnotationTag> tags = new LinkedList<>();
	
	public GOAnnotator(String taxon) {
		this.taxon = taxon;
	}
	
	
	/**
	 * Return taxon of organism associated with retrieved data
	 * @return
	 */
	public String getTaxon() {
		return taxon;
	}
	
}
