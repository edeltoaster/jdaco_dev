package PPIXpress_BRCA;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.StrPair;

public class eval_coverage {
	
	public static void main(String[] args) {
		System.out.println("11.5.2015");
		System.out.println("IntAct r189");
		Map<String, String> taxon_map = new HashMap<String, String>();
		taxon_map.put("human (homo sapiens)", "9606");
		taxon_map.put("mouse (mus musculus)", "10090");
		taxon_map.put("fruit fly (Drosophila melanogaster)", "7227");
		taxon_map.put("yeast (S.cerevisiae 288c)", "559292");
		taxon_map.put("human (homo sapiens) BG", "9606");
		
		Map<String, PPIN> ppis = new HashMap<String, PPIN>();
		
		Map<String, Double> all = new HashMap<String, Double>();
		Map<String, Double> local_only = new HashMap<String, Double>();
		Map<String, Double> stricter = new HashMap<String, Double>();
		Map<String, Double> retrieved = new HashMap<String, Double>();
		
		Map<String, Double> d_all = new HashMap<String, Double>();
		Map<String, Double> d_local_only = new HashMap<String, Double>();
		Map<String, Double> d_stricter = new HashMap<String, Double>();
		Map<String, Double> d_retrieved = new HashMap<String, Double>();
		
		System.out.println("ALL");
		Map<String, List<String>> ddis = DataQuery.getKnownDDIs();
		System.out.println("#domains: "+ ddis.size());
		HashSet<StrPair> test = new HashSet<StrPair>();
		for (String d1:ddis.keySet())
			for (String d2:ddis.get(d1))
				test.add(new StrPair(d1, d2));
		System.out.println("#IAs: "+ test.size());
		
		for (String organism:taxon_map.keySet()) {
			String taxon = taxon_map.get(organism);
			
			PPIN ppi = null;
			if (organism.equals("human (homo sapiens) BG"))
				ppi = new PPIN("mixed_data/human_biogrid.txt.gz");
			else
				ppi = DataQuery.getIntActNetwork(taxon);
			ppis.put(organism, ppi);
			NetworkBuilder netb = new NetworkBuilder(ppi);
			double map = netb.getMappingPercentage();
			all.put(organism, map);
			map = netb.getMappingDomainPercentage();
			d_all.put(organism, map);
		}
		
		DataQuery.localDDIsOnly();
		System.out.println("LOCAL-ONLY");
		ddis = DataQuery.getKnownDDIs();
		System.out.println("#domains: "+ ddis.size());
		test = new HashSet<StrPair>();
		for (String d1:ddis.keySet())
			for (String d2:ddis.get(d1))
				test.add(new StrPair(d1, d2));
		System.out.println("#IAs: "+ test.size());
		
		for (String organism:taxon_map.keySet()) {
			PPIN ppi = ppis.get(organism);
			NetworkBuilder netb = new NetworkBuilder(ppi);
			double map = netb.getMappingPercentage();
			local_only.put(organism, map);
			map = netb.getMappingDomainPercentage();
			d_local_only.put(organism, map);
		}
		
		DataQuery.stricterLocalDDIs();
		System.out.println("STRICTER LOCAL-ONLY");
		ddis = DataQuery.getKnownDDIs();
		System.out.println("#domains: "+ ddis.size());
		test = new HashSet<StrPair>();
		for (String d1:ddis.keySet())
			for (String d2:ddis.get(d1))
				test.add(new StrPair(d1, d2));
		System.out.println("#IAs: "+ test.size());
		
		for (String organism:taxon_map.keySet()) {
			PPIN ppi = ppis.get(organism);
			NetworkBuilder netb = new NetworkBuilder(ppi);
			double map = netb.getMappingPercentage();
			stricter.put(organism, map);
			map = netb.getMappingDomainPercentage();
			d_stricter.put(organism, map);
		}
		
		DataQuery.onlyRetrievedDDIs();
		System.out.println("ONLY RETRIEVED");
		ddis = DataQuery.getKnownDDIs();
		System.out.println("#domains: "+ ddis.size());
		test = new HashSet<StrPair>();
		for (String d1:ddis.keySet())
			for (String d2:ddis.get(d1))
				test.add(new StrPair(d1, d2));
		System.out.println("#IAs: "+ test.size());
		for (String organism:taxon_map.keySet()) {
			PPIN ppi = ppis.get(organism);
			NetworkBuilder netb = new NetworkBuilder(ppi);
			double map = netb.getMappingPercentage();
			retrieved.put(organism, map);
			map = netb.getMappingDomainPercentage();
			d_retrieved.put(organism, map);
		}
		
		System.out.println();
		System.out.println();
		System.out.println("SUMMARY: frac. of proteins with usable domains / frac. of interactions");
		for (String organism:taxon_map.keySet()) {
			System.out.println(organism);
			System.out.println("size: " + ppis.get(organism).getSizesStr());
			System.out.println("all: " + String.format("%.3g", d_all.get(organism)) + " / " + String.format("%.3g",all.get(organism)));
			System.out.println("local: " + String.format("%.3g",d_local_only.get(organism)) + " / " + String.format("%.3g",local_only.get(organism)));
			System.out.println("retrieved: " + String.format("%.3g",d_retrieved.get(organism)) + " / " + String.format("%.3g",retrieved.get(organism)));
			System.out.println("stricter: " + String.format("%.3g",d_stricter.get(organism)) + " / " + String.format("%.3g",stricter.get(organism)));
			System.out.println();
		}
		
	}
}
