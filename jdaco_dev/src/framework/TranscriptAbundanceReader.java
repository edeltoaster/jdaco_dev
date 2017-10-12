package framework;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * Helper class with parsers for plenty of file types describing expression data
 * @author Thorsten Will
 */
public class TranscriptAbundanceReader {
	
	/**
	 * Reads Cufflinks FPKM_tracking files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only transcripts/genes with transcripts above threshold or if a gene file is parsed genes above threshold.
	 * @param file
	 * @param threshold
	 * @param return_gene_level
	 * @return map of gene/transcript -> FPKM (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readCufflinksFPKMFile(String file, double threshold, boolean return_gene_level) {
		Map<String, Float> transcript_abundance = new HashMap<>(1024);
		Map<String, Float> gene_abundance = new HashMap<>(1024);
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				if (line.startsWith("tracking_id"))
					continue;
				String[] split = line.split("\\s+");
				String transcript = split[0].split("\\.")[0];
				String gene = split[3].split("\\.")[0];
				float fpkm = Float.parseFloat(split[9]);
				if (fpkm <= threshold) // common threshold = 1
					continue;
				
				transcript_abundance.put(transcript, fpkm);
				
				// gene abundance is the sum of transcript abundances
				if (!gene_abundance.containsKey(gene))
					gene_abundance.put(gene, 0f);
				gene_abundance.put(gene, gene_abundance.get(gene) + fpkm);
			}
			
		} catch (Exception e) {
			System.err.println("Problem while trying to parse Cufflinks FPKM file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		if (return_gene_level)
			return gene_abundance;
		
		return transcript_abundance;
	}
	
	/**
	 * Reads RSEM output files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only transcripts/genes with transcripts above threshold or if a gene file is parsed genes above threshold.
	 * @param file
	 * @param threshold
	 * @param return_gene_level
	 * @return map of gene/transcript -> TPM (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readRSEMFile(String file, double threshold, boolean return_gene_level) {
		Map<String, Float> transcript_abundance = new HashMap<>(1024);
		Map<String, Float> gene_abundance = new HashMap<>(1024);
		int transcript_index = 0;
		int gene_index = 1;
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				
				// check actual type of file
				if (line.startsWith("transcript_id")) // leave indices as already set
					continue;
				if (line.startsWith("gene_id")) { // change indices
					transcript_index = 1;
					gene_index = 0;
					continue;
				}
				
				String[] split = line.split("\\s+");
				String transcript = split[transcript_index].split("\\.")[0]; // shear off version
				String gene = split[gene_index].split("\\.")[0];
				
				float tpm = Float.parseFloat(split[5]);
				if (tpm <= threshold) // common threshold = 1
					continue;
				
				transcript_abundance.put(transcript, tpm);
				
				// gene abundance is the sum of transcript abundances
				if (!gene_abundance.containsKey(gene))
					gene_abundance.put(gene, 0f);
				gene_abundance.put(gene, gene_abundance.get(gene) + tpm);
			}
			
		} catch (Exception e) {
			System.err.println("Problem while trying to parse RSEM file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		if (return_gene_level)
			return gene_abundance;
		
		return transcript_abundance;
	}
	
	/**
	 * Reads Kallisto output files of genes/transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only genes/transcripts with genes/transcripts above TPM threshold.
	 * @param file
	 * @param threshold
	 * @return map of genes/transcript -> TPM (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readKallistoFile(String file, double threshold) {
		return readKallistoFile(file, threshold, false);
	}
	
	/**
	 * Reads Kallisto output files of genes/transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only genes/transcripts with genes/transcripts above TPM or count threshold.
	 * @param file
	 * @param threshold
	 * @param get_counts
	 * @return map of genes/transcript -> TPM or estimated count (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readKallistoFile(String file, double threshold, boolean get_counts) {
		Map<String, Float> transcript_abundance = new HashMap<>(1024);
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz") || file.endsWith(".gzip"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				
				if (line.startsWith("target_id")) //
					continue;
				
				String[] split = line.split("\\s+");
				String transcript = split[0].split("\\.")[0]; // shear off version
				float tpm = Float.parseFloat(split[4]);
				
				// tpm stores estinated counts instead
				if (get_counts)
					tpm = Float.parseFloat(split[3]);
				
				if (tpm <= threshold)
					continue;
				
				transcript_abundance.put(transcript, tpm);
			}

		} catch (Exception e) {
			System.err.println("Problem while trying to parse Kallisto file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		return transcript_abundance;
	}
	
	/**
	 * Reads Kallisto output files of transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only transcripts above TPM threshold, convert to gene abundances.
	 * @param file
	 * @param threshold
	 * @return map of genes -> TPM (for all transcripts above > threshold)
	 */
	public static Map<String, Float> readKallistoTranscriptsAsGenes(String file, double threshold, String organism_database) {
		
		Map<String, Float> abundance = readKallistoFile(file, threshold);
		Map<String, String> transcript_to_gene = new HashMap<String, String>(1024);
			
		// get association by query
		for (String[] data:DataQuery.getGenesTranscriptsProteins(organism_database)) {
			String gene = data[0];
			String transcript = data[1];
			transcript_to_gene.put(transcript, gene);
		}
			
		Map<String, Float> gene_abundance = new HashMap<String, Float>(1024);
		for (String transcript:abundance.keySet()) {
			if (!transcript_to_gene.containsKey(transcript))
				continue;
			String gene = transcript_to_gene.get(transcript);
				
			// gene = sum of all its transcripts
			if (!gene_abundance.containsKey(gene))
				gene_abundance.put(gene, 0f);
			gene_abundance.put(gene, gene_abundance.get(gene) + abundance.get(transcript));
				
		}
		
		return gene_abundance;
	}
	
	/**
	 * Reads Kallisto h5 output files. Creates temporary results that are cleaned afterwards.
	 * A UNIX-system and a working installation of kallisto in /usr/local/bin/kallisto are assumed.
	 * @param h5_file
	 * @param threshold
	 * @param get_counts
	 * @return
	 */
	public static Map<String, Float> readKallistoH5(String h5_file, double threshold, boolean get_counts) {
		String temp_folder = h5_file+"_temp";
		String[] call1 = {"/usr/local/bin/kallisto", "h5dump", "-o", temp_folder, h5_file};
		String[] call2 = {"rm", "-rf", temp_folder};
		Map<String, Float> expr_data = null;
		try {
			Process prc1 = Runtime.getRuntime().exec(call1);
			BufferedReader br = new BufferedReader(new InputStreamReader(prc1.getErrorStream()));
			String line = "";
			String abundance_file = "";
			while ((line = br.readLine()) != null) {
				if (line.startsWith("[h5dump] writing abundance file:")) {
					abundance_file = line.split("\\s+")[4];
				}
			}
			prc1.waitFor();
			
			// read data
			expr_data = TranscriptAbundanceReader.readKallistoFile(abundance_file, threshold, get_counts);
			
			// remove generated data
			Process prc2 = Runtime.getRuntime().exec(call2);
			prc2.waitFor();
		} catch (Exception e) {
			System.err.println("Problem while trying to convert Kallisto h5 file. "
					+ "Be aware that this method only works on UNIX systems for which 'kallisto' can be called in de terminal without further specification.");
		}
		
		return expr_data;
	}
	
	/**
	 * Reads GTF files of genes or transcripts (gzipped also fine, ending .gz assumed there) as from Cufflinks/ENCODE
	 * and allow all only transcripts/genes with transcripts above threshold or if a gene file is parsed genes above threshold.
	 * @param file
	 * @param threshold
	 * @param return_gene_level
	 * @return map of gene/transcript -> FPKM (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readGencodeGTFFile(String file, double threshold, boolean return_gene_level) {
		Map<String, Float> transcript_abundance = new HashMap<>(1024);
		Map<String, Float> gene_abundance = new HashMap<>(1024);
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				String[] split = line.split("\\t");
				
				String type = split[2];
				
				if (!type.equals("transcript") && !type.equals("gene"))
					continue;
				
				// parse attributes
				String transcript = "";
				String gene = "";
				float fpkm = 0;
				float sample_count = 1;
				for (String s:split[8].split(";")) {
					String[] temp = s.trim().split("\\s+");
					
					// special case: FPKM "       nan"
					if (temp.length > 2)
						continue;
					
					// cut out "" if there are any
					if (temp[1].startsWith("\""))
						temp[1] = temp[1].substring(1, temp[1].length()-1);
					
					if (temp[0].equals("gene_id"))
						gene = temp[1].split("\\.")[0];
					else if (temp[0].equals("transcript_id"))
						transcript = temp[1].split("\\.")[0];
					else if (temp[0].equals("FPKM") || temp[0].equals("RPKM") || temp[0].equals("TPM")) { // if there is one specifically elaborated datapoint -> take it and don't think about further ones
						fpkm = Float.parseFloat(temp[1]);
						break;
					}
					else if (temp[0].startsWith("FPKM") || temp[0].startsWith("RPKM") || temp[0].startsWith("TPM")) { // startswith but NOT equals
						fpkm += Float.parseFloat(temp[1]);
						sample_count++;
					}
				}
				// average if there are several in the same file (the case for some ENCODE data)
				fpkm /= sample_count;
				
				if (fpkm <= threshold) // common threshold = 1
					continue;
				
				transcript_abundance.put(transcript, fpkm);
				
				// gene abundance is the sum of transcript abundances
				if (!gene_abundance.containsKey(gene))
					gene_abundance.put(gene, 0f);
				gene_abundance.put(gene, gene_abundance.get(gene) + fpkm);
			}

		} catch (Exception e) {
			System.err.println("Problem while trying to parse Gencode GTF file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		if (return_gene_level)
			return gene_abundance;
		
		return transcript_abundance;
	}
	
	/**
	 * Reads GTF files of genes/transcripts as specified by CSHL ENCODE paper: 2 expressed replicates (>0) and iIDR < 0.1,
	 * also adds a transcript "no_samples" that specifies number of samples.
	 * @param file
	 * @param threshold
	 * @param return_gene_level
	 * @return map of transcript/no_samples -> RPKM/FPKM (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readCSHLData(String file) {
		Map<String, Float> transcript_abundance = new HashMap<>(1024);
		int max_sample_count = 1;
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				String[] split = line.split("\\t");
				
				String type = split[2];
				
				if (!type.equals("transcript") && !type.equals("gene"))
					continue;
				
				// parse attributes
				String transcript = "";
				float fpkm = 0;
				int sample_count = 0;
				for (String s:split[8].split(";")) {
					String[] temp = s.trim().split("\\s+");
					
					// special case: FPKM "       nan"
					if (temp.length > 2 || !temp[1].startsWith("\""))
						continue;
					
					temp[1] = temp[1].substring(1, temp[1].length()-1);
					
					if (temp[0].equals("transcript_id"))
						transcript = temp[1].split("\\.")[0];
					else if (temp[0].equals("gene_id"))
						transcript = temp[1].split("\\.")[0];
					else if (temp[0].equals("FPKM") || temp[0].equals("RPKM")) { // if there is one specifically elaborated datapoint -> take it and don't think about further ones
						fpkm = Float.parseFloat(temp[1]);
						sample_count = 1;
						break;
					}
					else if (temp[0].startsWith("FPKM") || temp[0].startsWith("RPKM")) { // startswith but NOT equals
						if (Float.parseFloat(temp[1]) == 0) {
							fpkm = 0;
							break;
						}
						fpkm += Float.parseFloat(temp[1]);
						sample_count++;
					}
					else if (temp[0].equals("iIDR")) {
						if (temp[1].equals("NA") || Float.parseFloat(temp[1]) >= 0.1)
							fpkm = 0;
						break;
					}
				}
				max_sample_count = Math.max(max_sample_count, sample_count);
				
				if (fpkm == 0)
					continue;
				
				fpkm /= sample_count;
				
				transcript_abundance.put(transcript, fpkm);
			}

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		transcript_abundance.put("no_samples", (float) max_sample_count);
		
		return transcript_abundance;
	}
	
	/**
	 * Reads simple expression files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only transcripts/genes with transcripts/genes above threshold. 
	 * @param file
	 * @param threshold
	 * @return map of gene/transcript -> expression (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readExpressionFile(String file, double threshold, boolean header) {
		return readExpressionFile(file, threshold, null, header);
	}
	
	/**
	 * Reads simple expression files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and allow all only transcripts/genes with transcripts/genes above threshold. 
	 * @param file
	 * @param threshold
	 * @param organism_database -> to convert to genes if necessary
	 * @return map of gene/transcript -> expression (if a transcript/gene above > threshold)
	 */
	public static Map<String, Float> readExpressionFile(String file, double threshold, String organism_database, boolean header) {
		Map<String, Float> abundance = new HashMap<>(1024);
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				
				if (line.startsWith("#"))
					continue;
				if (header) {
					header = false;
					continue;
				}
				String[] split = line.split("\\s+");
				String transcript = split[0].split("\\.")[0];
				float expr = Float.parseFloat(split[1]);
				if (expr <= threshold) // common threshold = 1
					continue;
				
				abundance.put(transcript, expr);
			}

		} catch (Exception e) {
			System.err.println("Problem while trying to parse expression file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		if (organism_database != null) {
			Map<String, String> transcript_to_gene = new HashMap<String, String>(1024);
			
			// get association by query
			for (String[] data:DataQuery.getGenesTranscriptsProteins(organism_database)) {
				String gene = data[0];
				String transcript = data[1];
				transcript_to_gene.put(transcript, gene);
			}
			
			Map<String, Float> gene_abundance = new HashMap<String, Float>();
			for (String transcript:abundance.keySet()) {
				if (!transcript_to_gene.containsKey(transcript))
					continue;
				String gene = transcript_to_gene.get(transcript);
				
				// gene = sum of all its transcripts
				if (!gene_abundance.containsKey(gene))
					gene_abundance.put(gene, 0f);
				gene_abundance.put(gene, gene_abundance.get(gene) + abundance.get(transcript));
				
			}
			abundance = gene_abundance;
		}
		
		return abundance;
	}
	
	/**
	 * Reads Cufflinks FPKM_tracking files of transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns all that are above a certain threshold.
	 * @param file
	 * @return map transcript -> FPKM (if > threshold)
	 */
	public static Map<String, Float> readCufflinksTranscriptsFPKM(String file, double threshold) {
		return readCufflinksFPKMFile(file, threshold, false);
	}
	
	/**
	 * Reads Cufflinks FPKM_tracking files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns all genes that are above a certain threshold. If a transcript file is used a gene is expressed if at least one transcript is above the threshold
	 * @param file
	 * @return map gene -> FPKM (if > threshold)
	 */
	public static Map<String, Float> readCufflinksGenesFPKM(String file, double threshold) {
		return readCufflinksFPKMFile(file, threshold, true);
	}
	
	/**
	 * Reads RSEM files of transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns all that are above a certain threshold.
	 * @param file
	 * @return map transcript -> TPM (if > threshold)
	 */
	public static Map<String, Float> readRSEMTranscriptsTPM(String file, double threshold) {
		return readRSEMFile(file, threshold, false);
	}
	
	/**
	 * Reads RSEM files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns all genes that are above a certain threshold. If a transcript file is used a gene is expressed if at least one transcript is above the threshold
	 * @param file
	 * @return map gene -> TPM (if > threshold)
	 */
	public static Map<String, Float> readRSEMGenesTPM(String file, double threshold) {
		return readRSEMFile(file, threshold, true);
	}
	
	/**
	 * Reads GTF files of transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns all that are above a certain threshold.
	 * @param file
	 * @return map transcript -> FPKM (if > threshold)
	 */
	public static Map<String, Float> readGencodeGTFTranscripts(String file, double threshold) {
		return readGencodeGTFFile(file, threshold, false);
	}
	
	/**
	 * Reads GTF files of genes or transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns all genes that are above a certain threshold. If a transcript file is used a gene is expressed if at least one transcript is above the threshold
	 * @param file
	 * @return map gene -> FPKM (if > threshold)
	 */
	public static Map<String, Float> readGencodeGTFGenes(String file, double threshold) {
		return readGencodeGTFFile(file, threshold, true);
	}
	
	/**
	 * Reads TCGA quantile-normalized isoform RSEM files of UCSC transcripts (gzipped also fine, ending .gz assumed there)
	 * and returns either transcripts or genes and their expression level; uses version of UCSC id.
	 * @param file
	 * @param transcript_threshold
	 * @param boolean return_gene_level
	 * @return map Ensembl transcript or gene -> RSEM (if transcript > threshold)
	 */
	public static Map<String, Float> readTCGAIsoformRSEMFile(String file, double transcript_threshold, boolean return_gene_level) {
		Map<String, Float> abundance = new HashMap<>(1024);
		
		Map<String, String> ucsc_to_ensembl = DataQuery.getUCSChg19toTranscriptMap();
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				String line = in.readLine();
				if (line.startsWith("isoform_id"))
					continue;
				String[] split = line.split("\\s+");
				String transcript = split[0];
				float rsem = Float.parseFloat(split[1]);
				
				// ensures consistency across transcr./genes.
				if (!ucsc_to_ensembl.containsKey(transcript))
					continue;
				
				// UCSC -> Ensembl
				transcript = ucsc_to_ensembl.get(transcript);
				
				// may map several times, sum is taken
				if (!abundance.containsKey(transcript))
					abundance.put(transcript, 0f);
				
				abundance.put(transcript, abundance.get(transcript) + rsem );
				
			}

		} catch (Exception e) {
			System.err.println("Problem while trying to parse TCGA file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		// check for threshold
		abundance.entrySet().removeIf(e -> e.getValue().floatValue() <= transcript_threshold);
		
		// if gene-level wanted: convert
		if (return_gene_level) {
			String db = DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens");
			// get association by query
			Map<String, String> transcript_to_gene = new HashMap<>(1024);
			for (String[] data:DataQuery.getGenesTranscriptsProteins(db)) {
				String gene = data[0];
				String transcript = data[1];
				transcript_to_gene.put(transcript, gene);
			}
			Map<String, Float> gene_abundance = new HashMap<>();
			for (String transcript:abundance.keySet()) {
				// every transcript should be mappable
				String gene = transcript_to_gene.get(transcript);
				
				// gene = sum of all its transcripts (not important for the context but generally correct)
				if (!gene_abundance.containsKey(gene))
					gene_abundance.put(gene, 0f);
				gene_abundance.put(gene, gene_abundance.get(gene) + abundance.get(transcript));
				
			}
			
			abundance = gene_abundance;
		}
		
		return abundance;
	}
	
	/**
	 * Reads TCGA normalized isoform RSEM files of UCSC transcripts (gzipped also fine, ending .gz assumed there)
	 * @param file
	 * @param threshold
	 * @return map Ensembl transcript -> RSEM (if > threshold)
	 */
	public static Map<String, Float> readTCGAIsoformRSEM(String file, double threshold) {
		return readTCGAIsoformRSEMFile(file, threshold, false);
	}
	
	/**
	 * Reads TCGA normalized isoform RSEM files of UCSC transcripts (gzipped also fine, ending .gz assumed there)
	 * @param file
	 * @param threshold
	 * @return map Ensembl gene -> RSEM (if any of its transcripts above threshold)
	 */
	public static Map<String, Float> readTCGAIsoformRSEMAsGenes(String file, double transcript_threshold) {
		return readTCGAIsoformRSEMFile(file, transcript_threshold, true);
	}
	
	/**
	 * Reads TCGA normalized gene RSEM files (gzipped also fine, ending .gz assumed there)
	 * @param file
	 * @return Ensembl gene -> RSEM (if > threshold)
	 */
	public static Map<String, Float> readTCGAGeneRSEM(String file, double threshold) {
		Map<String, Float> gene_abundance = new HashMap<>(1024);
		Map<String, String> HGCN_to_ensembl = new HashMap<>(1024);
		for (String[] entry:DataQuery.getHGNCProteinsGenes()) {
			if (entry[2].equals(""))
				continue;
			HGCN_to_ensembl.put(entry[0], entry[2]);
		}
		
		BufferedReader in = null;
		try {
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				String line = in.readLine();
				if (line.startsWith("gene_id"))
					continue;
				String[] split = line.split("\\s+");
				String HGCN = split[0].split("\\|")[0];
				if (HGCN.equals("?"))
					continue;
				Float rsem = Float.parseFloat(split[1]);
				
				if (rsem <= threshold) // common cutoff
					continue;
				
				if (HGCN_to_ensembl.containsKey(HGCN))
					gene_abundance.put(HGCN_to_ensembl.get(HGCN), rsem);
			}

		} catch (Exception e) {
			System.err.println("Problem while trying to parse TCGA file.");
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		return gene_abundance;
	}
	
	/**
	 * Reads Cufflinks FPKM_tracking file of genes or transcripts and returns the expressed ones (gzipped also fine, ending .gz assumed there).
	 * If a genes files is used the cutoff is used on the gene's levels, if a transcript file is given as the input, a gene is said to be expressed 
	 * if at least one of its transcripts is expressed.
	 * @param file
	 * @return set of transcribed genes (FPKM > threshold)
	 */
	public static Set<String> getGeneAbundanceFromCufflinksFPKM(String file, double threshold) {
		return readCufflinksGenesFPKM(file, threshold).keySet();
	}
	
	/**
	 * Reads TCGA normalized transcript file and returns the expressed genes (gzipped also fine, ending .gz assumed there)
	 * @param file
	 * @return set of transcribed Ensembl genes with transcripts above threshold
	 */
	public static Set<String> getGeneAbundanceFromTCGAIsoformRSEM(String file, double threshold) {
		return readTCGAIsoformRSEMFile(file, threshold, true).keySet();
	}
	
	/**
	 * Reads TCGA normalized gene file and returns the expressed genes (gzipped also fine, ending .gz assumed there)
	 * @param file
	 * @return set of transcribed Ensembl genes
	 */
	public static Set<String> getGeneAbundanceFromTCGAGeneRSEM(String file, double threshold) {
		return readTCGAGeneRSEM(file, threshold).keySet();
	}
	
	/**
	 * Reads TCGA normalized gene file and returns the expressed proteins (gzipped also fine, ending .gz assumed there)
	 * @param file
	 * @return set of expressed Uniprot proteins
	 */
	public static Set<String> getProteinAbundanceFromTCGAGeneRSEM(String file, double threshold) {
		Set<String> expressed_HGCN = readTCGAGeneRSEM(file, threshold).keySet();
		Set<String> expressed_proteins = new HashSet<>(1024);
		
		Map<String, String> HGCN_to_up = new HashMap<>();
		for (String[] entry:DataQuery.getHGNCProteinsGenes()) {
			if (entry[1].equals(""))
				continue;
			HGCN_to_up.put(entry[0], entry[1]);
		}
		
		for (String hgcn:expressed_HGCN)
			if (HGCN_to_up.containsKey(hgcn))
				expressed_proteins.add(HGCN_to_up.get(hgcn));
		 
		return expressed_proteins;
	}
	
	
	/*
	 * some little helpers
	 */
	
	
	/**
	 * Infers from the Ensembl identifier if it is a (T)ranscript, (G)ene or (P)rotein.
	 * Special cases: 
	 * - yeast genes starts with Y. SGD ORF.
	 * - fruitfly from flybase -> FBgn0032956 gene, FBtr0085899. FlyBase
	 * - C. elegans -> WBGene00004893, F53H8.4. WormBase
	 * - return 'o' for identifiers used in TCGA
	 * @param identifier
	 * @return
	 */
	public static char getEnsemblIDType(String identifier) {
		
		// is TCGA identifier, UCSC transcript or HGNC|Entrez
		if (identifier.startsWith("uc") || identifier.contains("|"))
			return 'o';
		
		// yeast gene
		if (identifier.startsWith("Y"))
			return 'G';
		else if (identifier.startsWith("F")) // C.elegans transcript
			return 'T';
		
		int i = 0;
		while (!Character.isDigit(identifier.charAt(i))) {
			i++;
			if (identifier.length() == i) // for non-digit identifiers
				return 'o';
		}
		
		char t = identifier.charAt(--i);
		
		if (t == 'n')
			return 'G';
		else if (t == 'r')
			return 'T';
		
		return t;
	}
	
	/**
	 * Tries to infer the type of input file. Supported are Cufflinks, TCGA and linewise files.
	 * At the moment Ensembl-identifiers are assumed.
	 * @param path
	 * @return type of file
	 */
	public static String inferTranscriptAbundanceFileType(String file) {
		
		BufferedReader in = null;
		try {
			boolean second_row = false;
			String may_be = "";
			
			if (file.endsWith(".gz"))
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			else
				in = new BufferedReader(new FileReader(file));
			
			while (in.ready()) {
				// parse
				String line = in.readLine();
				if (line.startsWith("#"))
					continue;
				String[] split = line.split("\\t");
				
				// should be Kallisto
				if (split.length == 5) {
					if (!second_row) {
						second_row = true;
						continue;
					}
					char type = getEnsemblIDType(split[0]);
					return "Kallisto " + type;
					
				} else if (split.length == 9) { // assume GTF output due to length
					String temp = split[2];
					if (temp.equals("gene"))
						return "GTF G";
					else
						return "GTF T";
					
				} else if (split.length == 13) {// assume Cufflinks due to length
					
					// get next line, no information here
					if (line.startsWith("tracking_id"))
						continue;
					
					char type = getEnsemblIDType(split[0]);
					return "Cufflinks " + type;
					
				} else if (split.length > 13) {// assume RSEM due to length
					
					// assume type from header
					if (line.startsWith("gene_id"))
						return "RSEM G";
					else if (line.startsWith("transcript_id"))
						return "RSEM T";
					
				} else {
					
					// assume TCGA from header
					if (line.startsWith("gene_id")) {
						may_be = "TCGA G";
						second_row = true;
						continue;
					}
					else if (line.startsWith("isoform_id")) {
						may_be = "TCGA T";
						second_row = true;
						continue;
					}
					
					if (line.startsWith("Protein")) // PPIN / DDIN ?
						return "other";
					
					// assume simple linewise file
					char type = getEnsemblIDType(split[0]);
					
					if (type == 'o') {
						if (second_row)
							if (may_be.equals(""))
								return "unknown";
							else
								return may_be; // some TCGA determined earlier
						else { // header for linewise?
							second_row = true;
							continue;
						}
					}
					
					if (second_row)
						return "linewise (header) " + type;
					else
						return "linewise " + type;
				}
				
			}
			
		} catch (Exception e) {
			System.err.println("Problem while trying to infer type of expression file.");
			return "unknown";
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		return "unknown";
	}
	
	/**
	 * Averages expression data from several files,
	 * note that the data should have been retrieved without any threshold to not introduce any artefacts.
	 * @param data
	 * @param threshold
	 * @return
	 */
	public static Map<String, Float> averageAbundances(Collection<Map<String, Float>> data, double threshold) {
		
		int datasets = data.size();
		
		// find all genes/transcripts that were noted
		Set<String> transcripts = new HashSet<>(1024);
		for (Map<String, Float> d:data) {
			transcripts.addAll(d.keySet());
		}
		
		Map<String, Float> temp_abundances = new HashMap<>(1024);
		for (String s:transcripts) {
			float expr = 0;
			for (Map<String, Float> d:data)
				expr += d.getOrDefault(s, 0f);
			temp_abundances.put(s, expr / datasets);
		}
		
		// filter for genes/transcripts that meet the threshold after averaging
		Map<String, Float> final_abundances = new HashMap<>(1024);
		for (String s:temp_abundances.keySet()) {
			float avg_expr = temp_abundances.get(s);
			if (avg_expr <= threshold)
				continue;
			final_abundances.put(s, avg_expr);
		}
		
		return final_abundances;
	}
	
	/**
	 * Convert expression data to counts/transcripts per million (better comparability across samples) 
	 * @param expr_data (non-pruned expression values)
	 * @param threshold
	 * @return normalized "per million" values of transcripts/genes with expression above threshold
	 */
	public static Map<String, Float> convertToPMMeasure(Map<String, Float> expr_data, double threshold) {
		float norm_sum = expr_data.values().stream().reduce(0f, Float::sum);
		Map<String, Float> tpms = new HashMap<>();
		
		for (String t:expr_data.keySet()) {
			float f = (float) (1000000 * expr_data.get(t) / norm_sum);
			
			if (f > threshold)
				tpms.put(t, f);
		}
		
		return tpms;
	}
	
	/**
	 * Normalizes transcript count data by the length of the mapped transcript
	 * @param expr_data
	 * @return length-normalized data
	 */
	public static Map<String, Float> normalizeByTranscriptLength(Map<String, Float> expr_data, String ensembl_db) {
		Map<String, Integer> transcript_length = DataQuery.getTranscriptsCDNALength(ensembl_db);
		Map<String, Float> norm_data = new HashMap<>();
		
		for (String transcript:expr_data.keySet()) {
			if (!transcript_length.containsKey(transcript))
				continue;
			norm_data.put(transcript, expr_data.get(transcript) / transcript_length.get(transcript));
		}
		
		return norm_data;
	}
	
	/**
	 * Read arbitrary but supported sample file, returns null if not supported
	 * @param path
	 * @return
	 */
	public static Map<String, Float> readSample(String path, double threshold, boolean gene_level_only, String organism_database) {
		String type = TranscriptAbundanceReader.inferTranscriptAbundanceFileType(path);
		Map<String, Float> abundance = null;
		
		if (type.equals("Cufflinks G")) {
			abundance = TranscriptAbundanceReader.readCufflinksGenesFPKM(path, threshold);
		} else if (type.equals("Cufflinks T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readCufflinksGenesFPKM(path, threshold);
			else
				abundance = TranscriptAbundanceReader.readCufflinksTranscriptsFPKM(path, threshold);
			
		} else if (type.equals("RSEM G")) {
			abundance = TranscriptAbundanceReader.readRSEMGenesTPM(path, threshold);
		} else if (type.equals("RSEM T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readRSEMGenesTPM(path, threshold);
			else
				abundance = TranscriptAbundanceReader.readRSEMTranscriptsTPM(path, threshold);
			
		} else if (type.equals("GTF G")) {
			abundance = TranscriptAbundanceReader.readGencodeGTFGenes(path, threshold);
		} else if (type.equals("GTF T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readGencodeGTFGenes(path, threshold);
			else
				abundance = TranscriptAbundanceReader.readGencodeGTFTranscripts(path, threshold);
			
		} else if (type.equals("linewise G")) {
			abundance = TranscriptAbundanceReader.readExpressionFile(path, threshold, false);
		} else if (type.equals("linewise T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readExpressionFile(path, threshold, organism_database, false);
			else
				abundance = TranscriptAbundanceReader.readExpressionFile(path, threshold, false);
			
		} else if (type.equals("linewise (header) G")) {
			abundance = TranscriptAbundanceReader.readExpressionFile(path, threshold, true);
		} else if (type.equals("linewise (header) T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readExpressionFile(path, threshold, organism_database, true);
			else
				abundance = TranscriptAbundanceReader.readExpressionFile(path, threshold, true);
			
		} else if (type.equals("TCGA G")) {
			abundance = TranscriptAbundanceReader.readTCGAGeneRSEM(path, threshold);
		} else if (type.equals("TCGA T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readTCGAIsoformRSEMAsGenes(path, threshold);
			else
				abundance = TranscriptAbundanceReader.readTCGAIsoformRSEM(path, threshold);
			
		} else if (type.equals("Kallisto G")) {
			abundance = TranscriptAbundanceReader.readKallistoFile(path, threshold);
		} else if (type.equals("Kallisto T")) {
			if (gene_level_only)
				abundance = TranscriptAbundanceReader.readKallistoTranscriptsAsGenes(path, threshold, organism_database);
			else
				abundance = TranscriptAbundanceReader.readKallistoFile(path, threshold);
		}
		
		return abundance;
	}
	
	/**
	 * Determines an expression threshold from the percentile of the expression value distribution, 5% as threshold -> input 5
	 * @param file
	 * @param percentile
	 * @return
	 */
	public static double getExprThresholdFromPercentile(String file, double percentile, boolean gene_level_only, String organism_database) {
		
		// read input
		List<Float> values = new LinkedList<Float>(readSample(file, -1, gene_level_only, organism_database).values());
		
		// get percentile in list
		Collections.sort(values);
		double p = (percentile) / 100.0 * values.size();
		
		return values.get((int) p);
	}
}