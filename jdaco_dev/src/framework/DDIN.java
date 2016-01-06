package framework;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Domain-domain interaction network suitable for performant DACO algorithm execution
 * @author Thorsten Will
 */
public class DDIN {
	
	// okay with many threads since only written once
	private HashMap<String, String[]> ddis;
	private HashMap<String, String[]> protein_to_domains;
	private HashMap<String, String> domain_to_protein;
	
	/**
	 * Some constructors
	 */
	public DDIN(HashMap<String, String[]> ddis, HashMap<String, String[]> protein_to_domains, HashMap<String, String> domain_to_protein) {
		this.ddis = ddis;
		this.protein_to_domains = protein_to_domains;
		this.domain_to_protein = domain_to_protein;
	}
	
	public DDIN(Map<String, List<String>> ddis, Map<String, List<String>> protein_to_domains, Map<String, String> domain_to_protein) {
		// transform to right data structures
		this.ddis = new HashMap<>();
		for (String domain:ddis.keySet())
			this.ddis.put(domain, ddis.get(domain).toArray(new String[ddis.get(domain).size()]));
		
		this.protein_to_domains = new HashMap<String, String[]>();
		for (String protein:protein_to_domains.keySet())
			this.protein_to_domains.put(protein, protein_to_domains.get(protein).toArray(new String[protein_to_domains.get(protein).size()]));
		
		this.domain_to_protein = new HashMap<>(domain_to_protein);
	}
	
	// format [Protein/Domain1 Domain2 IAType Weight] assumed
	public DDIN(String file) {
		this.domain_to_protein = new HashMap<>();
		Map<String, List<String>> protein_to_domains = new HashMap<>();
		Map<String, List<String>> ddis = new HashMap<>();
		
		// read file
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				String line = in.readLine();
				String[] split = line.split("\\s+");
				String partner1 = split[0];
				String partner2 = split[1];
				String type = split[2];
				
				// skip first line
				if (partner1.equals("Protein/Domain1"))
					continue;
				
				if (type.equals("pd")) { // protein-domain
					
					// protein_to_domain
					if (!protein_to_domains.containsKey(partner1))
						protein_to_domains.put(partner1, new LinkedList<String>());
					protein_to_domains.get(partner1).add(partner2);
					
					// domain_to_protein
					this.domain_to_protein.put(partner2, partner1);
					
				} else { // domain-domain interaction
					// add to DDI
					if (!ddis.containsKey(partner1))
						ddis.put(partner1, new LinkedList<String>());
					ddis.get(partner1).add(partner2);
					if (!ddis.containsKey(partner2))
						ddis.put(partner2, new LinkedList<String>());
					ddis.get(partner2).add(partner1);
				}
			}
			
		} catch (Exception e) {
			if (e instanceof FileNotFoundException)
				System.err.println("Problem while opening domain interaction network " + file + ".");
			else
				System.err.println("Problem while parsing domain interaction network " + file + ".");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
				// no output necessary at this point
			}
		}
		
		this.ddis = new HashMap<>();
		for (String domain:ddis.keySet())
			this.ddis.put(domain, ddis.get(domain).toArray(new String[ddis.get(domain).size()]));
		
		this.protein_to_domains = new HashMap<>();
		for (String protein:protein_to_domains.keySet())
			this.protein_to_domains.put(protein, protein_to_domains.get(protein).toArray(new String[protein_to_domains.get(protein).size()]));
		
		this.domain_to_protein = new HashMap<>(domain_to_protein);
	}
	
	public HashMap<String, String[]> getDDIs() {
		return ddis;
	}

	public HashMap<String, String[]> getProtein_to_domains() {
		return protein_to_domains;
	}

	public HashMap<String, String> getDomain_to_protein() {
		return domain_to_protein;
	}
	
	/**
	 * Writes the network to a file in SIF format
	 * @param file
	 */
	public void writeDDIN(String file) {
		List<String> to_write = new LinkedList<>();
		
		to_write.add("Protein/Domain1 Domain2 IAType Weight");


		// write proteins -> domains
		for (String protein:this.protein_to_domains.keySet())
			for (String domain:this.protein_to_domains.get(protein))
				to_write.add(protein + " " + domain + " pd 10");

		//write DDIs
		for (String domain1:this.ddis.keySet())
			for (String domain2:this.ddis.get(domain1)) {
				// no double occurance
				if (domain1.compareTo(domain2) <= 0)
					continue;
				to_write.add(domain1 + " " + domain2 + " dd 1");
			}
		
		Utilities.writeEntries(to_write, file);
	}
	
	/**
	 * @return number of proteins in DDIN
	 */
	public int getNumberOfProteins() {
		return this.protein_to_domains.keySet().size();
	}
	
	/**
	 * @return number of domains in DDIN
	 */
	public int getNumberOfDomains() {
		return this.domain_to_protein.keySet().size();
	}
	
	/**
	 * @return number of domain-domain interactions in DDIN
	 */
	public int getNumberOfDDIs() {
		int n = 0;
		for (String domain1:this.ddis.keySet())
			for (String domain2:this.ddis.get(domain1))
				if (domain1.compareTo(domain2) < 0) // count only once
					n++;
		return n;
	}
	
	/**
	 * Returns [#proteins,#domains,#interactions]
	 * @return
	 */
	public int[] getSizes() {
		int[] sizes = new int[3];
		sizes[0] = this.getNumberOfProteins();
		sizes[1] = this.getNumberOfDomains();
		sizes[2] = this.getNumberOfDDIs();
		return sizes;
	}
	
	/**
	 * Returns "#proteins / #domains / #interactions"
	 * @return
	 */
	public String getSizesStr() {
		int[] sizes = this.getSizes();
		return sizes[0] + " proteins / " + sizes[1] + " domains / " + sizes[2] + " interactions";
	}
}
