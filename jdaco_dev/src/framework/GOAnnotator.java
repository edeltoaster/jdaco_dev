package framework;

import java.util.LinkedList;
import java.util.List;

public class GOAnnotator {
	private String taxon = "not set";
	private boolean include_IEA;
	List<GOAnnotationTag> tags = new LinkedList<>();
	
	/*
	 * Constructor to load definition data from definition file
	 */
	public GOAnnotator(String taxon, boolean include_IEA, String file) {
		this.taxon = taxon;
		this.include_IEA = include_IEA;
		
		this.readInputFile(file);
	}
	
	/*
	 * Constructor to load definition data from definition file
	 */
	public GOAnnotator(String taxon, String file) {
		this.taxon = taxon;
		this.include_IEA = true;
		
		this.readInputFile(file);
	}
	
	/*
	 * Constructor to load old GOAnnotator object from file, taxon and IEA settings are overwritten
	 */
	public GOAnnotator(String file) {
		this.readInputFile(file);
	}
	
	private void readInputFile(String file) {
		List<String> data = Utilities.readFile(file);
		
		// empty file not helpful
		if (data.size() == 0)
			return;
		
		int columns = data.get(0).split("\\s+").length;
		boolean is_definition = false;
		
		// number of columns will help to distinguish different types of data files
		if (columns == 3)
			is_definition = true;
		
		for (String line:data) {
			// thread as definition file
			if (is_definition)
				tags.add(new GOAnnotationTag(this.taxon, this.include_IEA, line));
			// thread as output file with retrieved data
			else
				tags.add(new GOAnnotationTag(line));
		}
		
		// taxon and include_IEA still need to be set
		if (!is_definition) {
			// tags is non-empty
			this.taxon = tags.get(0).getTaxon();
			this.include_IEA = tags.get(0).IncludesIEA();
		}
	}
	
	/**
	 * Write retrieved data to file for later usage
	 * @param output_file
	 */
	public void writeRetrievedData(String output_file) {
		List<String> to_write = new LinkedList<>();
		for (GOAnnotationTag tag:this.tags)
			to_write.add(tag.getDataString());
		Utilities.writeEntries(to_write, output_file);
	}
	
	/**
	 * Return taxon of organism associated with retrieved data
	 * @return
	 */
	public String getTaxon() {
		return taxon;
	}
	
	/**
	 * Returns if proteins annotated by purely computationally inferred evidence were included  
	 * @return
	 */
	public boolean IncludesIEA() {
		return include_IEA;
	}
	
	/**
	 * Returns included GOAnnotationTags
	 * @return
	 */
	public List<GOAnnotationTag> getGOAnnotationTags() {
		return this.tags;
	}
	
}
