package framework;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * Protein-protein interaction network implementation suitable for DACO implementation
 * @author Thorsten Will
 */
public class PPIN {
	
	// okay with many threads since only written once
	private HashMap<String, Set<String>> partners = new HashMap<String, Set<String>>();
	private HashMap<StrPair, Double> weights = new HashMap<StrPair, Double>();
	private HashMap<String, Double> w_whole = new HashMap<String, Double>();
	private boolean is_weighted = false;
	
	/**
	 * Reads PPIN from file, assumes SIF-file of [PROTEIN1 PROTEIN2 WEIGHT] or [PROTEIN1 PROTEIN2], may also be gzipped.
	 * @param file
	 */
	public PPIN(String file) {
		readPPINFile(file, 0.0);
	}
	
	/**
	 * Reads PPIN from file, assumes SIF-file of [PROTEIN1 PROTEIN2 WEIGHT] or [PROTEIN1 PROTEIN2],
	 * only adds interactions with weight >= cutoff. Assumes UniProt Accs.
	 * @param file
	 * @param cutoff
	 */
	public PPIN(String file, double cutoff) {
		readPPINFile(file, cutoff);
	}
	
	
	/**
	 * Reads PPIN from stream, assumes layout of [PROTEIN1 PROTEIN2 WEIGHT] or [PROTEIN1 PROTEIN2],
	 * only adds interactions with weight >= cutoff.
	 * @param file
	 * @param cutoff
	 * @param protein_format
	 */
	public PPIN(BufferedReader in, double cutoff) {
		processPPINInput(in, cutoff);
	}
	
	/**
	 * Internal method that maps file to generic reader
	 * @param file
	 * @param cutoff
	 */
	private void readPPINFile(String file, double cutoff) {
		try {
			BufferedReader in = null;
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			// parsing
			processPPINInput(in, cutoff);
			
			in.close();
		} catch (Exception e) {
			System.err.println("Problem while opening/parsing " + file + ".");
		}
	}
	
	/**
	 * Internal method to avoid multiplicity of code
	 * @param file
	 * @param cutoff
	 */
	private void processPPINInput(BufferedReader in, double cutoff) {
		boolean check = true;
		Map<String, String> conversion_map = null;
		try {
			while (in.ready()) {
				String line = in.readLine();
				
				// skip first line
				if (line == null || line.startsWith("Protein1") || line.isEmpty() || line.startsWith("#"))
					continue;
				String[] split = line.split("\\s+");
				String p1 = split[0];
				String p2 = split[1];
				
				// check if Uniprot or other, implemented: HGNC
				if (check) {
					
					// check if Uniprot
					if (!p1.matches("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")) {
						
						// check for ensembl GENES: general, yeast, fruitfly, c.elegans
						if ( (p1.length() > 6 && p1.startsWith("ENS") || p1.startsWith("YDL") || p1.startsWith("FBgn") || p1.startsWith("WBGene")) ) {
							System.out.println("Retrieving ENSEMBL conversion data ... ");
							Set<String> test_set = new HashSet<String>();
							test_set.add(p1);
							test_set.add(p2);
							
							String org_db = DataQuery.getEnsemblOrganismDatabaseFromProteins(DataQuery.getUniprotFromEnsemblGenes(test_set).values());
							if (org_db.equals("")) {
								System.err.println("None of the ENSEMBL genes in the first row of the PPIN data could be mapped to an organism.");
								return;
							}
							
							conversion_map = new HashMap<String, String>();
							for (String[] temp:DataQuery.getGenesTranscriptsProteins(org_db)) {
								conversion_map.put(temp[0], temp[2]);
							}
						} else {// otherwise assume HGNC, build map
							System.out.println("Retrieving HGNC conversion data ... ");
							conversion_map = new HashMap<String, String>();
							for (String[] temp:DataQuery.getHGNCProteinsGenes()) {
								conversion_map.put(temp[0], temp[1]);
							}
						}
					}
					check = false;
				} 
				
				if (conversion_map != null) {
					p1 = conversion_map.get(p1);
					p2 = conversion_map.get(p2);
					if (p1.equals("") || p2.equals(""))
						continue;
				}
				
				double w = 1.0;
				
				if (this.is_weighted || split.length == 3) {
					this.is_weighted = true;
					w = Double.parseDouble(split[2]);
				}
				
				// no self-interactions and a certain cutoff
				if (p1.equals(p2) || w < cutoff) 
					continue;
				
				// every interaction only once
				if (!this.partners.containsKey(p1)) {
					this.w_whole.put(p1, 0.0);
					this.partners.put(p1, new HashSet<String>());
				}
				
				this.partners.get(p1).add(p2);
				this.w_whole.put(p1, this.w_whole.get(p1) + w);
				
				if (!this.partners.containsKey(p2)) {
					this.w_whole.put(p2, 0.0);
					this.partners.put(p2, new HashSet<String>());
				}
				this.partners.get(p2).add(p1);
				this.w_whole.put(p2, this.w_whole.get(p2) + w);
				
				this.weights.put(new StrPair(p1, p2), w);
			}
			
			in.close();
		} catch (Exception e) {
			System.err.println("Problem while parsing protein interaction network.");
		}
	}
	
	/**
	 * Alternative constructor that builds subnetworks
	 * @param unpruned_ppi
	 * @param abundant_proteins
	 */
	public PPIN(PPIN unpruned_ppi, Set<String> abundant_proteins) {
		
		this.is_weighted = unpruned_ppi.is_weighted;
		
		for (String p1:unpruned_ppi.getPartners().keySet()) {
			// only if there
			if (!abundant_proteins.contains(p1))
				continue;
			for (String p2:unpruned_ppi.getPartners().get(p1)) {
				// only if 2 there and every pair only one time
				if (!abundant_proteins.contains(p2) || p1.compareTo(p2) <= 0)
					continue;
				StrPair pair = new StrPair(p1, p2);
				double w = unpruned_ppi.getWeights().get(pair);
				
				// every interaction only once -> add both directions
				if (!this.partners.containsKey(p1)) {
					this.w_whole.put(p1, 0.0);
					this.partners.put(p1, new HashSet<String>());
				}
				this.partners.get(p1).add(p2);
				this.w_whole.put(p1, this.w_whole.get(p1) + w);
				
				if (!this.partners.containsKey(p2)) {
					this.w_whole.put(p2, 0.0);
					this.partners.put(p2, new HashSet<String>());
				}
				this.partners.get(p2).add(p1);
				this.w_whole.put(p2, this.w_whole.get(p2) + w);
				
				this.weights.put(pair, w);
			}
		}
	}
	
	/**
	 * Alternative constructor that builds subnetworks of isoforms
	 * @param unpruned_ppi
	 * @param abundant_proteins
	 */
	public PPIN(PPIN unpruned_ppi, Map<String, Set<String>> supported_interactions) {
		
		this.is_weighted = unpruned_ppi.isWeighted();
		
		for (String p1:supported_interactions.keySet()) 
			for (String p2:supported_interactions.get(p1)) {
				// every pair only one time
				if (p1.compareTo(p2) <= 0)
					continue;
				StrPair pair = new StrPair(p1, p2);
				double w = unpruned_ppi.getWeights().get(pair);
				
				// every interaction only once -> add both directions
				if (!this.partners.containsKey(p1)) {
					this.w_whole.put(p1, 0.0);
					this.partners.put(p1, new HashSet<String>());
				}
				this.partners.get(p1).add(p2);
				this.w_whole.put(p1, this.w_whole.get(p1) + w);
				
				if (!this.partners.containsKey(p2)) {
					this.w_whole.put(p2, 0.0);
					this.partners.put(p2, new HashSet<String>());
				}
				this.partners.get(p2).add(p1);
				this.w_whole.put(p2, this.w_whole.get(p2) + w);
				
				this.weights.put(pair, w);
			}
	}
	
	/**
	 * Alternative constructor
	 * @param abundant_proteins
	 */
	public PPIN(Map<String, Set<String>> interactions) {
		
		this.is_weighted = false;
		
		for (String p1:interactions.keySet()) 
			for (String p2:interactions.get(p1)) {
				// every pair only one time
				if (p1.compareTo(p2) <= 0)
					continue;
				StrPair pair = new StrPair(p1, p2);
				double w = 1.0;
				
				// every interaction only once -> add both directions
				if (!this.partners.containsKey(p1)) {
					this.w_whole.put(p1, 0.0);
					this.partners.put(p1, new HashSet<String>());
				}
				this.partners.get(p1).add(p2);
				this.w_whole.put(p1, this.w_whole.get(p1) + w);
				
				if (!this.partners.containsKey(p2)) {
					this.w_whole.put(p2, 0.0);
					this.partners.put(p2, new HashSet<String>());
				}
				this.partners.get(p2).add(p1);
				this.w_whole.put(p2, this.w_whole.get(p2) + w);
				
				this.weights.put(pair, w);
			}
	}
	
	/**
	 * Alternative constructor that combines two networks where one defines the interactions and one the weights,
	 * if wanted, edges for which no weights can be assigned from the other network are simply added with a weight of 1.
	 */
	public PPIN(PPIN interaction_network, PPIN weight_network, boolean add_unweighted) {
		this.is_weighted = true;
		
		for (String p1:interaction_network.getPartners().keySet()) {
			// only if there
			if (!add_unweighted && !weight_network.contains(p1))
				continue;
			for (String p2:interaction_network.getPartners().get(p1)) {
				// only if p2 also partner and every pair only one time
				if ((!add_unweighted && !weight_network.getPartners().get(p1).contains(p2)) || p1.compareTo(p2) <= 0)
					continue;
				StrPair pair = new StrPair(p1, p2);
				
				// weight from other network
				double w = 1.0;
				if (weight_network.getWeights().containsKey(pair))
					w = weight_network.getWeights().get(pair);
				
				// every interaction only once -> add both directions
				if (!this.partners.containsKey(p1)) {
					this.w_whole.put(p1, 0.0);
					this.partners.put(p1, new HashSet<String>());
				}
				this.partners.get(p1).add(p2);
				this.w_whole.put(p1, this.w_whole.get(p1) + w);
				
				if (!this.partners.containsKey(p2)) {
					this.w_whole.put(p2, 0.0);
					this.partners.put(p2, new HashSet<String>());
				}
				this.partners.get(p2).add(p1);
				this.w_whole.put(p2, this.w_whole.get(p2) + w);
				
				this.weights.put(pair, w);
			}
		}
	}
	
	/**
	 * Alternative merging constructor
	 * @param PPINs_to_merge
	 */
	public PPIN(List<PPIN> PPINs_to_merge) {
		
		for (PPIN p:PPINs_to_merge) {
			for (String p1:p.partners.keySet()) 
				for (String p2:p.partners.get(p1)) {
					// every pair only one time
					if (p1.compareTo(p2) <= 0)
						continue;
					StrPair pair = new StrPair(p1, p2);
					double w = p.weights.get(pair);
					
					if (w != 1.0)
						this.is_weighted = true;
					
					// every interaction only once -> add both directions
					if (!this.partners.containsKey(p1)) {
						this.w_whole.put(p1, 0.0);
						this.partners.put(p1, new HashSet<String>());
					}
					
					// interaction already added, don't add weights twice
					if (this.partners.get(p1).contains(p2))
						continue;
					
					this.partners.get(p1).add(p2);
					this.w_whole.put(p1, this.w_whole.get(p1) + w);
					
					if (!this.partners.containsKey(p2)) {
						this.w_whole.put(p2, 0.0);
						this.partners.put(p2, new HashSet<String>());
					}
					this.partners.get(p2).add(p1);
					this.w_whole.put(p2, this.w_whole.get(p2) + w);
					
					this.weights.put(pair, w);
				}
		}
	}
	
	public boolean contains(String protein) {
		return this.w_whole.keySet().contains(protein);
	}
	
	public HashMap<String, Set<String>> getPartners() {
		return partners;
	}

	/**
	 * Returns all interactions as a set of StrPair
	 * @return
	 */
	public Set<StrPair> getInteractions() {
		return new HashSet<StrPair>(this.weights.keySet());
	}
	
	public HashMap<StrPair, Double> getWeights() {
		return weights;
	}
	
	public HashMap<String, Double> getWholeWeights() {
		return w_whole;
	}
	
	public boolean isWeighted() {
		return this.is_weighted;
	}
	
	/**
	 * Retrieves edge-weights range [0..1] from the functional association network STRING,
	 * if no weights can be assigned, a weight of 1.0 is given for add_unweighted = true or the edge isn't added for false.
	 * @param add_unweighted
	 * @return
	 */
	public PPIN getAsSTRINGWeighted(boolean add_unweighted) {
		String taxon_id = DataQuery.getTaxonFromProteins(this.getProteins());
		PPIN STRING = DataQuery.getSTRINGNetwork(taxon_id);
		
		return new PPIN(this, STRING, add_unweighted);
	}
	
	/**
	 * Computes cohesiveness for a set of proteins
	 * @param proteins
	 * @return [w_in, w_out]
	 */
	public double[] computeClusterCohesiveness(HashSet<String> proteins) {
		double[] coh = new double[2];
		for (String p1 : proteins)
			for (String p2 : partners.get(p1)) {
				// unequality checked during construction
				if ((proteins.contains(p2)) && (p1.compareTo(p2) < 0) ) // only count once
					coh[0] += weights.get(new StrPair(p1, p2)); // internal edge
				else
					coh[1] += weights.get(new StrPair(p1, p2)); // external edge
			}
		return coh;
	}
	
	/**
	 * Computes the change in cohesiveness
	 * @param protein
	 * @param proteins
	 * @return [w_in, w_out]
	 */
	public double[] computeDeltaCohesiveness(String protein, HashSet<String> proteins) {
		double[] coh = new double[2];
		Set<String> temp = partners.get(protein);
		for (String p2 : proteins)
			if (temp.contains(p2))
				coh[0] += weights.get(new StrPair(protein, p2));
		coh[1] = w_whole.get(protein) - coh[0];
		return coh;
	}
	
	/**
	 * @return All proteins in the network
	 */
	public Set<String> getProteins() {
		return new HashSet<String>(this.w_whole.keySet());
	}
	
	/**
	 * Returns [#proteins,#interactions]
	 * @return
	 */
	public int[] getSizes() {
		int[] sizes = new int[2];
		sizes[0] = this.getProteins().size();
		sizes[1] = this.getAllWeights().size();
		return sizes;
	}
	
	/**
	 * Returns "#proteins / #interactions"
	 * @return
	 */
	public String getSizesStr() {
		int[] sizes = this.getSizes();
		return sizes[0] + " proteins / " + sizes[1] + " interactions";
	}
	/**
	 * Writes the network to a file in SIF format
	 * @param file
	 */
	public void writePPIN(String file) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			bw.write("Protein1 Protein2 weight");
			bw.newLine();
			for (StrPair pair:this.weights.keySet()) {
				String temp = pair.getL() + " " + pair.getR() + " " + this.weights.get(pair);
				bw.write(temp);
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			System.err.println("Problem while trying to write PPIN.");
		}
	}
	
	/**
	 * Permutes edges,
	 * weights not relevant anymore afterwards.
	 * @param iterations
	 */
	public void permutateEdges(int iterations) {
		// algorithm from "Specificity and Stability in Topology of Protein Networks", Science 2002
		Random rnd = new Random(System.currentTimeMillis());
		List<StrPair> edges = new ArrayList<StrPair>(this.weights.size());
		StrPair AB = null;
		StrPair CD = null;
		StrPair AD = null;
		StrPair BC = null;
		String A = null;
		String B = null;
		String C = null;
		String D = null;
		int e = this.weights.size();
		for (int i = 0; i< iterations; i++) {
			// init
			edges.clear();
			edges.addAll(this.weights.keySet());
			
			// search for 2 edges A<->B, C<->D
			Collections.shuffle(edges, rnd);
			while (AB == null) {
				
				AB = edges.get(rnd.nextInt(e));
				CD = edges.get(rnd.nextInt(e));
				
				if (AB.equals(CD))
					continue;
				
				// check
				AD = new StrPair(AB.getL(), CD.getR());
				BC = new StrPair(AB.getR(), CD.getL());
				
				// new ones not in current network, but in principle realizable
				if (this.weights.containsKey(AD) || this.weights.containsKey(BC)) {
					AB = null;
				}
			}
			
			// rewire to A<->D, B<->C
			A = AB.getL();
			B = AB.getR();
			C = CD.getL();
			D = CD.getR();
			
			// rewire weights
			this.weights.remove(AB);
			this.weights.remove(CD);
			this.weights.put(AD, 1.0);
			this.weights.put(BC, 1.0);
			
			// rewire partners
			this.partners.get(A).remove(B);
			this.partners.get(B).remove(A);
			this.partners.get(C).remove(D);
			this.partners.get(D).remove(C);
			
			this.partners.get(A).add(D);
			this.partners.get(D).add(A);
			this.partners.get(B).add(C);
			this.partners.get(C).add(B);
		}
	}
	
	private List<Double> getAllWeights() {
		List<Double> all_weights = new ArrayList<Double>(this.weights.values());
		// sorts in ascending order
		Collections.sort(all_weights);
		return all_weights;
	}
	
	/**
	 * Suggests a probability-cutoff for the pair building, best 5% interactions -> 5
	 * @param percentile
	 * @return
	 */
	public double getPercentile(double percentile) {
		List<Double> all_weights = this.getAllWeights();
		int n = all_weights.size();
		double p = (100-percentile) / 100.0 * n;
		
		return all_weights.get((int) p);
	}
	
	/**
	 * Suggests a probability-cutoff for the growth of complexes
	 * @param percentile
	 * @param depth
	 * @return
	 */
	public double getSampledComplexCutoff(double percentile, int depth) {
		List<Double> all_weights = this.getAllWeights();
		List<Double> high_weights = new LinkedList<Double>();
		int no_samples = 10000;
		double cutoff = this.getPercentile(percentile);
		
		for (Double P:all_weights)
			if (P > cutoff)
				high_weights.add(P);
		double P_sampled = 0.0;
		for (int i = 0; i< no_samples; i++) {
			double P_sample = 1.0;
			Collections.shuffle(high_weights);
			for (Double P:high_weights.subList(0, depth))
				P_sample *= P;
			P_sampled += P_sample;
		}
		return P_sampled / no_samples;
	}
	
	// set-like operations
	/**
	 * Returns a new PPIN where all interactions from the given PPIN are removed
	 * @param other_ppin
	 * @return
	 */
	public PPIN removeAllIAs(PPIN other_ppin) {
		Map< String, Set<String>> supported_interactions = new HashMap<String, Set<String>>();
		
		// do a deep copy
		for (String protein:this.partners.keySet())
			supported_interactions.put(protein, new HashSet<String>(this.partners.get(protein)));
				
		// remove all IAs from given set
		for (String protein1:other_ppin.partners.keySet())
			for (String protein2:other_ppin.partners.get(protein1)) {
				if (protein1.compareTo(protein2) <= 0) // no IA twice
					continue;
				
				if (supported_interactions.containsKey(protein1)) {
					supported_interactions.get(protein1).remove(protein2);
					
					if (supported_interactions.get(protein1).size() == 0)
						supported_interactions.remove(protein1);
				}
				
				if (supported_interactions.containsKey(protein2)) {
					supported_interactions.get(protein2).remove(protein1);
					if (supported_interactions.get(protein2).size() == 0)
						supported_interactions.remove(protein2);
				}
			}
		
		return new PPIN(this, supported_interactions);
	}
	
	/**
	 * Returns a new PPIN where only interactions from the given PPIN are retain
	 * @param other_ppin
	 * @return
	 */
	public PPIN retainAllIAs(PPIN other_ppin) {
		Map< String, Set<String>> supported_interactions = new HashMap<String, Set<String>>();
		
		// do a deep copy
		for (String protein:this.partners.keySet())
			supported_interactions.put(protein, new HashSet<String>(this.partners.get(protein)));
		
		// remove IAs that are not in given set
		for (String protein1:this.partners.keySet())
			for (String protein2:this.partners.get(protein1)) {
				if (protein1.compareTo(protein2) <= 0) // no IA twice
					continue;
				
				// if protein1 in other PPIN
				if (!other_ppin.partners.containsKey(protein1)) {
					supported_interactions.remove(protein1);
				} else {
					if (!other_ppin.partners.get(protein1).contains(protein2)) {
						supported_interactions.get(protein1).remove(protein2);
						if (supported_interactions.get(protein1).size() == 0)
							supported_interactions.remove(protein1);
					}
				}
				
				// if protein2 in other PPIN
				if (!other_ppin.partners.containsKey(protein2)) {
					supported_interactions.remove(protein2);
				} else {
					if (!other_ppin.partners.get(protein2).contains(protein1)) {
						supported_interactions.get(protein2).remove(protein1);
						if (supported_interactions.get(protein2).size() == 0)
							supported_interactions.remove(protein2);
					}
				}
				
			}
		
		return new PPIN(this, supported_interactions);
	}
	
	/**
	 * Returns a new PPIN where all interactions from this and the given PPIN are merged; weights are lost and set to 1
	 * @param other_ppin
	 * @return
	 */
	public PPIN mergeAllIAs(PPIN other_ppin) {
		Map< String, Set<String>> all_interactions = new HashMap<String, Set<String>>();
		
		// do a deep copy
		for (String protein:this.partners.keySet())
			all_interactions.put(protein, new HashSet<String>(this.partners.get(protein)));
		
		// add more IAs
		for (String protein1:other_ppin.partners.keySet())
			for (String protein2:other_ppin.partners.get(protein1)) {
				if (protein1.compareTo(protein2) <= 0) // no IA twice
					continue;
				
				//protein1
				if (all_interactions.containsKey(protein1)) {
					all_interactions.get(protein1).add(protein2);
				} else {
					all_interactions.put(protein1, new HashSet<String>());
					all_interactions.get(protein1).add(protein2);
				}
				
				// protein2
				if (all_interactions.containsKey(protein2)) {
					all_interactions.get(protein2).add(protein1);
				} else {
					all_interactions.put(protein2, new HashSet<String>());
					all_interactions.get(protein2).add(protein1);
				}
			}
		
		return new PPIN(all_interactions);
	}
	
	public double getPairWeight(String protein1, String protein2) {
		if (this.getPartners().get(protein1).contains(protein2))
			return this.getWeights().get(new StrPair(protein1, protein2));
		
		return 0.0;
	}
}
