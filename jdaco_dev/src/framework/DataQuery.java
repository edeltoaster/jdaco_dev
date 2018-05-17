package framework;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

/**
 * Data retrieval helper
 * @author Thorsten Will
 */
public class DataQuery {

	// Servers
	// ensembldb.ensembl.org:3306 from Ensembl 48 onwards, 3337: GRCh37 current/previous
	// useastdb.ensembl.org:3306 current und previous version only
	private static String ensembl_mysql = "ensembldb.ensembl.org:3306";
	private static Map<String, List<String>> known_DDIs;
	private static boolean up2date_DDIs = true;
	private static String STRING_version = "10";
	private static String specific_ensembl_release = ""; // global option
	private static String specific_3did_release = "current";
	private static int retries = 0;
	private static int timeout = 3*60; // 3min
	private static long last_server_change = System.currentTimeMillis(); // holds track of server switches, important for parallel retrieval
	private static PrintStream err_out = System.err;
	private static boolean stricter_local_DDIs = false;
	private static boolean no_local_DDIs = false;
	
	// caches
	private static Map<String, List<String[]>> cache_genestransprots = new HashMap<>();
	private static Map<String, Map<String, String>> cache_isoformtranscr = new HashMap<>();
	private static Map<String,Map<String, List<String>>> cache_transcrdom = new HashMap<>();
	private static Map<String, Set<String>> cache_decay_transcripts = new HashMap<>();
	private static Map<String, Map<String, List<String>>> cache_ensembl_proteins = new HashMap<>();
	private static Map<String, Map<String, String>> cache_ensembl_names = new HashMap<>();
	private static Map<String, String> cache_db = new HashMap<>();
	private static Map<String, PPIN> cache_STRING = new HashMap<>();
	private static Map<String, String> cache_ucsc = new HashMap<>();
	private static Map<String, Map<String, String>> cache_refseq = new HashMap<>();
	private static Map<String, Map<String, Integer>> cache_transcript_length = new HashMap<>();
	private static Map<String, String> cache_ensembl_db = new HashMap<>();
	private static Map<String, Map<String, String>> cache_uniprot_gene_map = new HashMap<>();
	private static List<String[]> cache_HGNC;
	private static Map<String, String> cache_uniprot_HGNC;
	private static Map<String, String> uniprot_sec_accs;
	private static String uniprot_release;
	
	
	/**
	 * Queries UniProt to find organism name from a given protein
	 * @param uniprot_acc
	 * @return organism string
	 */
	public static String getOrganismFromUniprot(String uniprot_acc) {
		
		String organism = "";
		BufferedReader datastream = null;
		try {
			URL server = new URL("http://www.uniprot.org/uniprot/"+uniprot_acc+".txt");
			URLConnection connection = server.openConnection();
			// read
			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty())
					continue;
				if (line.startsWith("OS")) {
					line = line.split("   ")[1];
					organism = line.substring(0, line.length()-1);
					break;
				}
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("UniProt");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to infer organism from UniProt in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
			return getOrganismFromUniprot(uniprot_acc);
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
				// no output worthwhile
			}
		}
		
		return organism;
	}
	
	/**
	 * Queries UniProt to find NCBI taxon from a given protein
	 * @param uniprot_acc
	 * @return organism string
	 */
	public static String getTaxonFromUniprot(String uniprot_acc) {
		
		String taxon = "";
		BufferedReader datastream = null;
		try {
			URL server = new URL("http://www.uniprot.org/uniprot/"+uniprot_acc+".txt");
			URLConnection connection = server.openConnection();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty())
					continue;
				if (line.startsWith("OX")) {
					line = line.split("\\s+")[1];
					taxon = line.replace("NCBI_TaxID=", "");
					taxon = taxon.replace(";", "");
					break;
				}
			}

		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("UniProt");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to infer taxon from UniProt in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			return getTaxonFromUniprot(uniprot_acc);
			
		} finally {
			try {
				datastream.close();
			} catch (IOException e) {
				// no output necessary
			}
		}
		
		return taxon;
	}
	
	/***
	 * Detects taxon id for organism from given proteins
	 * @param collection of proteins
	 * @return taxon
	 */
	public static String getTaxonFromProteins(Collection<String> proteins) {
		String taxon = "";
		for (String protein:proteins) {
			taxon = DataQuery.getTaxonFromUniprot(protein);
			if (!taxon.equals(""))
				break;
		}
		
		// fix for s. cerevisiae
		if (taxon.equals("559292"))
			taxon = "4932";
		
		return taxon;
	}
	
	/**
	 * Queries UniProt to find organism name from NCBI taxon
	 * @param taxon_id
	 * @return organism string
	 */
	public static String getOrganismFromTaxon(String taxon_id) {
		
		String organism = "";
		BufferedReader datastream = null;
		try {
			URL server = new URL("http://www.uniprot.org/taxonomy/"+taxon_id+".rdf");
			URLConnection connection = server.openConnection();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty() || !line.startsWith("<scientificName>"))
					continue;
				organism = line.trim();
				organism = organism.replace("<scientificName>", "");
				organism = organism.replace("</scientificName>", "");
				break;
			}

		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("UniProt");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to infer organism from UniProt in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			return getOrganismFromTaxon(taxon_id);
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
				// no output necessary
			}
		}
		
		return organism;
	}
	
	/***
	 * Detects most recent database version for certain organism or tries to get a specific version if it can be found
	 * @param organism
	 * @return Ensembl database name as string
	 */
	public static String getEnsemblOrganismDatabaseFromName(String organism) {
		
		if (DataQuery.cache_ensembl_db.containsKey(organism))
			return DataQuery.cache_ensembl_db.get(organism);
		
		String query = "'"+organism.split(" ")[0].toLowerCase()+"\\_"+organism.split(" ")[1].toLowerCase()+"\\_core\\_%\\_%'";
		String db = "";
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("show databases like "+query);
			
			// take newest version
			if (DataQuery.specific_ensembl_release.equals("")) {
				while (rs.next()) {
					db = rs.getString(1); // sorted by name
				}
			} else { // take a specific version
				while (rs.next()) {
					db = rs.getString(1); // sorted by name
					String version = db.split("_")[db.split("_").length-2];
					if (version.equals(DataQuery.specific_ensembl_release))
						break; // db will be the right one or the newest
				}
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to infer organism database from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getEnsemblOrganismDatabaseFromName(organism);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
				// no output necessary
			}
		}
		
		DataQuery.cache_ensembl_db.put(organism, db);
		
		return db;
	}
	
	/***
	 * Detects most recent database version for organism from given proteins
	 * @param collection of proteins
	 * @return Ensembl database name as string
	 */
	public static String getEnsemblOrganismDatabaseFromProteins(Collection<String> proteins) {
		String db = "";
		String organism = "";
		String query_protein = "";
		
		for (String protein:proteins) {
			query_protein = protein;
			if (DataQuery.cache_db.containsKey(protein))
				return DataQuery.cache_db.get(protein);
			
			organism = DataQuery.getOrganismFromUniprot(protein);
			if (!organism.equals(""))
				break;
		}
		
		if (!organism.equals(""))
			db = DataQuery.getEnsemblOrganismDatabaseFromName(organism);
		
		if (!db.equals(""))
			DataQuery.cache_db.put(query_protein, db);
		
		return db;
	}
	
	/**
	 * Queries Ensembl for association of UniProt Accs. to transcripts and genes
	 * @param organism_core_database
	 * @returnreturn list of arrays {GENE, TRANSCRIPT, UNIPROT-ACC}
	 */
	public static List<String[]> getGenesTranscriptsProteins(String organism_core_database) {
		return DataQuery.getGenesTranscriptsProteinsSingleQuery(organism_core_database); // defines which implementation is used
	}
	
	/**
	 * Queries Ensembl for association of UniProt Accs. to transcripts and genes
	 * @param organism_core_database
	 * @returnreturn list of arrays {GENE, TRANSCRIPT, UNIPROT-ACC}
	 */
	public static List<String[]> getGenesTranscriptsProteinsSingleQuery(String organism_core_database) {
		
		if (DataQuery.cache_genestransprots.containsKey(organism_core_database))
			return DataQuery.cache_genestransprots.get(organism_core_database);
		List<String[]> associations = new LinkedList<>();
		
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT gene.stable_id, transcript.stable_id, trans_table.dbprimary_acc FROM gene, transcript, (SELECT translation.translation_id, xref.dbprimary_acc FROM translation, object_xref, xref WHERE translation.translation_id=object_xref.ensembl_id AND object_xref.xref_id=xref.xref_id AND xref.external_db_id='2200') AS trans_table WHERE gene.gene_id=transcript.gene_id AND transcript.canonical_translation_id = trans_table.translation_id");
			while (rs.next()) {
				String[] row = {rs.getString(1), rs.getString(2), rs.getString(3)};
				associations.add(row);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			//e.printStackTrace();
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get gene/transcript/protein data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getGenesTranscriptsProteinsSingleQuery(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		DataQuery.cache_genestransprots.put(organism_core_database, associations);
		return associations;
	}
	
	/**
	 * Queries Ensembl for association of UniProt Accs. to transcripts and genes
	 * @param organism_core_database
	 * @returnreturn list of arrays {GENE, TRANSCRIPT, UNIPROT-ACC}
	 */
	public static List<String[]> getGenesTranscriptsProteinsParallelQuery(String organism_core_database) {
		
		if (DataQuery.cache_genestransprots.containsKey(organism_core_database))
			return DataQuery.cache_genestransprots.get(organism_core_database);
		
		List<String[]> associations = new LinkedList<>();
		
		ExecutorService es = Executors.newFixedThreadPool(2);
		
		try {
			// compute
			Future<Map<String,List<String>>> t_p = es.submit(new getTranscriptProteinsTask(organism_core_database));
			Future<Map<String,List<String>>> g_t = es.submit(new getGenesTranscriptsTask(organism_core_database));
			
			// get results
			Map<String,List<String>> transcript_proteins = t_p.get();
			Map<String,List<String>> gene_transcripts = g_t.get();
			es.shutdown();
			
			// write into array
			for (String gene:gene_transcripts.keySet())
				for (String transcript:gene_transcripts.get(gene)) {
					if (!transcript_proteins.containsKey(transcript))
						continue;
					for (String protein:transcript_proteins.get(transcript)) {
						associations.add(new String[]{gene, transcript, protein});
					}
				}
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			
			//e.printStackTrace();
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get gene/transcript/protein association from ENSEMBL in 10 seconds ..." );
			es.shutdownNow();
			
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getGenesTranscriptsProteinsParallelQuery(organism_core_database);
			
		}
		
		DataQuery.cache_genestransprots.put(organism_core_database, associations);
		return associations;
	}
	
	/**
	 * Queries Ensembl for association of Enseml genes with their common names as given in Entrez (for compatibility)
	 * @param organism_core_database
	 */
	public static Map<String, String> getGenesCommonNames(String organism_core_database) {
		
		// return directly if in cache
		if (DataQuery.cache_ensembl_names.containsKey(organism_core_database))
			return DataQuery.cache_ensembl_names.get(organism_core_database);
		
		Map<String,String> gene_to_names = new HashMap<>();
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT gene.stable_id, xref.display_label "
					+ "FROM gene, object_xref, xref "
					+ "WHERE gene.gene_id=object_xref.ensembl_id AND object_xref.xref_id=xref.xref_id AND xref.external_db_id='1300'");
			
			while (rs.next()) {
				gene_to_names.put(rs.getString(1), rs.getString(2));
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			e.printStackTrace();
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get common names from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getGenesCommonNames(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		// store in cache
		DataQuery.cache_ensembl_names.put(organism_core_database, gene_to_names);
		
		return gene_to_names;
	}
	
	/**
	 * Queries Ensembl for association of Ensembl proteins to UniProt Accs., returns all accs per case
	 * @param organism_core_database
	 * @returnreturn map Ensembl -> UniProt accs
	 */
	public static Map<String, List<String>> getEnsemblToAllUniprotProteins(String organism_core_database) {
		
		// return directly if in cache
		if (DataQuery.cache_ensembl_proteins.containsKey(organism_core_database))
			return DataQuery.cache_ensembl_proteins.get(organism_core_database);
		
		Map<String, List<String>> ensembl_to_uniprot = new HashMap<>();
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT translation.stable_id, xref.dbprimary_acc "
					+ "FROM translation, object_xref, xref "
					+ "WHERE translation.translation_id=object_xref.ensembl_id AND object_xref.xref_id=xref.xref_id AND xref.external_db_id='2200'");
			
			while (rs.next()) {
				if (!ensembl_to_uniprot.containsKey(rs.getString(1)))
					ensembl_to_uniprot.put(rs.getString(1), new LinkedList<String>());
				ensembl_to_uniprot.get(rs.getString(1)).add(rs.getString(2));
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get protein data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getEnsemblToAllUniprotProteins(organism_core_database);
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		// store in cache
		DataQuery.cache_ensembl_proteins.put(organism_core_database, ensembl_to_uniprot);
		
		return ensembl_to_uniprot;
	}
	
	/**
	 * Given the organism_core_database, retrieves all Uniprot Accs of Allosome genes (X, Y, Z or W).
	 * @param organism_core_database
	 * @return
	 */
	public static Set<String> getAllosomeProteins(String organism_core_database) {
		
		Set<String> allosome_proteins = new HashSet<String>();
		
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			
			// first, determine the right coord_system_id for the organism database
			ResultSet rs1 = st.executeQuery("SELECT coord_system.coord_system_id FROM coord_system WHERE coord_system.name = 'chromosome' AND coord_system.attrib='default_version'");
			String chromosome_attrib = "";
			while (rs1.next()) {
				chromosome_attrib = rs1.getString(1);
			}
			
			// second, query proteins and where the genes are located
			ResultSet rs2 = st.executeQuery("SELECT trans_table.dbprimary_acc, seq_region.name FROM transcript, seq_region, (SELECT translation.translation_id, xref.dbprimary_acc FROM translation, object_xref, xref WHERE translation.translation_id=object_xref.ensembl_id AND object_xref.xref_id=xref.xref_id AND xref.external_db_id='2200') AS trans_table WHERE transcript.canonical_translation_id = trans_table.translation_id AND transcript.seq_region_id = seq_region.seq_region_id AND seq_region.coord_system_id='" + chromosome_attrib + "'");
			while (rs2.next()) {
				String protein = rs2.getString(1);
				String chromosome = rs2.getString(2);
				if (chromosome.equals("X") || chromosome.equals("Y") || chromosome.equals("Z") || chromosome.equals("W"))
					allosome_proteins.add(protein);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			//e.printStackTrace();
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get allosome protein data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getAllosomeProteins(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		return allosome_proteins;
	}
	
	/**
	 * Queries Ensembl for transcripts with RefSeq mRNA id
	 * @param organism_core_database
	 * @returm ENST->RefSeq
	 */
	public static Map<String, String> getTranscriptsWithRefSeq(String organism_core_database) {
		Map<String, String> refseq_annotated_transcripts = new HashMap<>();
		
		// get cached result
		if (DataQuery.cache_refseq.containsKey(organism_core_database))
			return DataQuery.cache_refseq.get(organism_core_database);
		
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			// query using database name from https://www.biostars.org/p/106470/
			ResultSet rs = st.executeQuery("SELECT transcript.stable_id, xref.display_label FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Transcript' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = 'RefSeq_mRNA'");
			while (rs.next()) {
				refseq_annotated_transcripts.put(rs.getString(1), rs.getString(2));
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			//e.printStackTrace();
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get RefSeq-annotated transcript data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getTranscriptsWithRefSeq(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}

		// store in cache
		DataQuery.cache_refseq.put(organism_core_database, refseq_annotated_transcripts);
		
		return refseq_annotated_transcripts;
	}
	
	/**
	 * Queries HGNC gene symbols to UniProt/Ensembl gene, may have ""
	 * @param organism_core_database
	 * @return
	 */
	public static List<String[]> getHGNCProteinsGenes() {
		
		// get from cache
		if (DataQuery.cache_HGNC != null)
			return DataQuery.cache_HGNC;
		
		// queries HGNC
		List<String[]> entries = new LinkedList<>();
		BufferedReader datastream = null;
		try {
			URL server = new URL("http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=md_prot_id&col=md_ensembl_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit");
			URLConnection connection = server.openConnection();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty() || line.startsWith("Approved"))
					continue;
				String[] temp = line.trim().split("\t");
				if (temp.length <= 2)
					temp = new String[]{temp[0],"",""};
				entries.add(temp);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("HGNC");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get naming data from HGNC in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getHGNCProteinsGenes();
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}
		
		// store
		DataQuery.cache_HGNC = entries;
		
		return entries;
	}
	
	/**
	 * Translates collection of human UniProt accessions to a list of HGNC gene symbols.
	 * @param input
	 * @return
	 */
	public static List<String> batchHGNCNamesFromProteins(Collection<String> input) {
		
		// if not cached, get data
		if (cache_uniprot_HGNC == null) {
			List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
			cache_uniprot_HGNC = new HashMap<String, String>();
			for (String[] s:HGNC_map) {
				cache_uniprot_HGNC.put(s[1], s[0]);
			}
		}
		
		List<String> output = new LinkedList<String>();
		for (String in:input)
			output.add(cache_uniprot_HGNC.getOrDefault(in, in));
		
		return output;
	}
	
	/**
	 * Translates a human UniProt accessions to a HGNC gene symbol,
	 * returns the UniProt accession if not mappable.
	 * @param Uniprot_acc
	 * @return
	 */
	public static String getHGNCNameFromProtein(String Uniprot_acc) {
		
		// if not cached, get data
		if (cache_uniprot_HGNC == null) {
			List<String[]> HGNC_map = DataQuery.getHGNCProteinsGenes();
			cache_uniprot_HGNC = new HashMap<String, String>();
			for (String[] s:HGNC_map) {
				cache_uniprot_HGNC.put(s[1], s[0]);
			}
		}
		
		return cache_uniprot_HGNC.getOrDefault(Uniprot_acc, Uniprot_acc);
	}
	
	/**
	 * Retrieves naming data to provide Uniprot -> gene name map.
	 * @return
	 */
	public static Map<String, String> getUniprotToGeneNameMap(Collection<String> sample_proteins) {
		
		if (sample_proteins.size() == 0)
			return new HashMap<String, String>();
		
		String organism_db = DataQuery.getEnsemblOrganismDatabaseFromProteins(sample_proteins);
		
		// get from cache
		if (DataQuery.cache_uniprot_gene_map.containsKey(organism_db))
			return DataQuery.cache_uniprot_gene_map.get(organism_db);
		
		Map<String, String> gene_to_name = DataQuery.getGenesCommonNames(organism_db);
		Map<String, String> up_to_name = new HashMap<>();
		for (String[] data:DataQuery.getGenesTranscriptsProteins(organism_db)) {
			String gene = data[0];
			String protein = data[2];
			up_to_name.put(protein, gene_to_name.get(gene));
		}
		
		// remove broken naming information
		up_to_name.entrySet().removeIf(e -> e.getValue() == null);
		
		// save in cache
		DataQuery.cache_uniprot_gene_map.put(organism_db, up_to_name);
		
		return up_to_name;
	}
	
	/**
	 * Retrieves naming data to provide Uniprot -> gene name map.
	 * @return
	 */
	public static Map<String, String> getUniprotToGeneNameMap(String organism_db) {
		
		// get from cache
		if (DataQuery.cache_uniprot_gene_map.containsKey(organism_db))
			return DataQuery.cache_uniprot_gene_map.get(organism_db);
		
		Map<String, String> gene_to_name = DataQuery.getGenesCommonNames(organism_db);
		Map<String, String> up_to_name = new HashMap<>();
		for (String[] data:DataQuery.getGenesTranscriptsProteins(organism_db)) {
			String gene = data[0];
			String protein = data[2];
			up_to_name.put(protein, gene_to_name.get(gene));
		}
		
		// remove broken naming information
		up_to_name.entrySet().removeIf(e -> e.getValue() == null);
		
		// save in cache
		DataQuery.cache_uniprot_gene_map.put(organism_db, up_to_name);
		
		return up_to_name;
	}
	
	/**
	 * Translates a human UniProt accessions to a gene symbol,
	 * returns the UniProt accession if not mappable.
	 * @param Uniprot_acc
	 * @return
	 */
	public static String getGeneNameFromProtein(String Uniprot_acc) {
		String organism_db = DataQuery.getEnsemblOrganismDatabaseFromName(DataQuery.getOrganismFromUniprot(Uniprot_acc));
		
		// if not cached, get data
		if (!DataQuery.cache_uniprot_gene_map.containsKey(organism_db)) {
			DataQuery.getUniprotToGeneNameMap(organism_db);
		}
		
		return DataQuery.cache_uniprot_gene_map.get(organism_db).getOrDefault(Uniprot_acc, Uniprot_acc);
	}
	
	/**
	 * Translates collection of human UniProt accessions to a list of gene symbols.
	 * @param input
	 * @return
	 */
	public static List<String> batchGeneNamesFromProteins(Collection<String> input) {
		
		String organism_db = DataQuery.getEnsemblOrganismDatabaseFromProteins(input);
		
		// if not cached, get data
		if (!DataQuery.cache_uniprot_gene_map.containsKey(organism_db)) {
			DataQuery.getUniprotToGeneNameMap(organism_db);
		}
		
		List<String> output = new LinkedList<String>();
		for (String in:input)
			output.add(cache_uniprot_gene_map.get(organism_db).getOrDefault(in, in));
		
		return output;
	}
	
	/**
	 * Queries Uniprot biomart to determine organism, then get mapping to ensembl
	 * @return
	 */
	public static Map<String,String> getUniprotFromEnsemblGenes(Collection<String> ensembl_genes) {
		
		// Java string join in 1.8+, not even in 1.7
		String genes = "";
		boolean first = true;
		for (String gene:ensembl_genes) {
			if (first) {
				genes += gene;
				first = false;
				continue;
			} 
			genes += "," + gene;
		}
		
		String query_xml = "<!DOCTYPE Query><Query client=\"true\" processor=\"TSV\" limit=\"-1\" header=\"1\"><Dataset name=\"uniprot\" config=\"uniprot_config\"><Filter name=\"ensembl_id\" value=\"" + genes + "\" filter_list=\"\"/><Filter name=\"entry_type\" value=\"Swiss-Prot\" filter_list=\"\"/><Attribute name=\"accession\"/><Attribute name=\"ensembl_id\"/></Dataset></Query>";
		Map<String, String> up_ens_map = new HashMap<>();
		
		BufferedReader datastream = null;
		try {
			URL server = new URL("http://www.biomart.org/biomart/martservice?query="+ URLEncoder.encode(query_xml, "UTF-8"));
			URLConnection connection = server.openConnection();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty() || line.startsWith("Accession"))
					continue;
				String[] temp = line.trim().split("\t");
				
				if (temp.length < 2)
					break;
				
				String up = temp[0];
				String ens = temp[1];
				up_ens_map.put(ens, up);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("BioMART");
			e.printStackTrace();
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get naming data from BioMART in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getUniprotFromEnsemblGenes(ensembl_genes);
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}
		
		return up_ens_map;
	}
	
	/**
	 * Queries UCSC hg19 for association of UCSC ids to Ensembl transcripts
	 * @return map USCS_id.version -> ENST
	 */
	public static Map<String, String> getUCSChg19toTranscriptMap() {
		
		if (!DataQuery.cache_ucsc.isEmpty())
			return DataQuery.cache_ucsc;

		Connection connection = null;
		try {
			// ucsc -> ens
			connection = DriverManager.getConnection("jdbc:mysql://"+ "genome-euro-mysql.soe.ucsc.edu" +"/" + "hg19", "genomep", "password");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT name, value FROM knownToEnsembl");
			
			while (rs.next()) {
				String ucsc = rs.getString(1); 
				String ens = rs.getString(2);
				DataQuery.cache_ucsc.put(ucsc, ens);
			}

			// store UCSC internal updates
			st = connection.createStatement();
			st.setQueryTimeout(timeout);
			rs = st.executeQuery("SELECT oldId, newId FROM kg6ToKg7");
			
			while (rs.next()) {
				String ucsc_old = rs.getString(1); 
				String ucsc_new = rs.getString(2);
				if (DataQuery.cache_ucsc.containsKey(ucsc_new))
					DataQuery.cache_ucsc.put(ucsc_old, DataQuery.cache_ucsc.get(ucsc_new));
			}
			
			st = connection.createStatement();
			st.setQueryTimeout(timeout);
			rs = st.executeQuery("SELECT oldId, newId FROM kg5ToKg6");
			
			while (rs.next()) {
				String ucsc_old = rs.getString(1); 
				String ucsc_new = rs.getString(2);
				if (DataQuery.cache_ucsc.containsKey(ucsc_new))
					DataQuery.cache_ucsc.put(ucsc_old, DataQuery.cache_ucsc.get(ucsc_new));
			}
			
			st = connection.createStatement();
			st.setQueryTimeout(timeout);
			rs = st.executeQuery("SELECT oldId, newId FROM kg4ToKg5");
			
			while (rs.next()) {
				String ucsc_old = rs.getString(1); 
				String ucsc_new = rs.getString(2);
				if (DataQuery.cache_ucsc.containsKey(ucsc_new))
					DataQuery.cache_ucsc.put(ucsc_old, DataQuery.cache_ucsc.get(ucsc_new));
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("UCSC");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get USCS data in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getUCSChg19toTranscriptMap();
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}

		return cache_ucsc;
	}
	
	/**
	 * Queries Ensembl for associated isoform transcript (longest CDS) of UniProt Accs.
	 * @param organism_core_database
	 * @returnreturn map of UniProt -> longest Ensembl transcript
	 */
	public static Map<String, String> getIsoformTranscriptsOfProteins(String organism_core_database) {
		
		if (DataQuery.cache_isoformtranscr.containsKey(organism_core_database))
			return DataQuery.cache_isoformtranscr.get(organism_core_database);
		
		Map<String, String> longest_transcript = new HashMap<>();
		Map<String, Integer> longest_length = new HashMap<>();
		
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT transcript.stable_id, trans_table.dbprimary_acc, translation_attrib.value FROM transcript, (SELECT translation.translation_id, xref.dbprimary_acc FROM translation, object_xref, xref WHERE translation.translation_id=object_xref.ensembl_id AND object_xref.xref_id=xref.xref_id AND xref.external_db_id='2200') AS trans_table, translation_attrib WHERE transcript.canonical_translation_id = trans_table.translation_id AND trans_table.translation_id = translation_attrib.translation_id AND translation_attrib.attrib_type_id = '167'");
			
			while (rs.next()) {
				String transcript = rs.getString(1);
				String uniprot = rs.getString(2);
				int length = Integer.parseInt(rs.getString(3));
				
				if (longest_length.containsKey(uniprot)) {
					if (length > longest_length.get(uniprot)) {
						longest_length.put(uniprot, length);
						longest_transcript.put(uniprot, transcript);
					}
				} else {
					longest_length.put(uniprot, length);
					longest_transcript.put(uniprot, transcript);
				}
				
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get transcript data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getIsoformTranscriptsOfProteins(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		DataQuery.cache_isoformtranscr.put(organism_core_database, longest_transcript);
		return longest_transcript;
	}
	
	/**
	 * Queries Ensembl for domain information of transcripts
	 * @param organism_core_database
	 * @return map of Ensembl transcripts to Pfam domains
	 */
	public static Map<String, List<String>> getTranscriptsDomains(String organism_core_database) {
		
		if (DataQuery.cache_transcrdom.containsKey(organism_core_database))
			return DataQuery.cache_transcrdom.get(organism_core_database);
		
		Map <String, List<String>> transcript_domain_mapping = new HashMap<>();
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT transcript.stable_id, protein_feature.hit_name, protein_feature.hit_start, protein_feature.hit_end, transcript.biotype "
					+ "FROM transcript, translation, protein_feature, analysis "
					+ "WHERE transcript.canonical_translation_id = translation.translation_id AND translation.translation_id=protein_feature.translation_id AND protein_feature.analysis_id = analysis.analysis_id AND analysis.logic_name = 'Pfam'");
			
			while (rs.next()) {
				if (!transcript_domain_mapping.containsKey(rs.getString(1)))
					transcript_domain_mapping.put(rs.getString(1), new LinkedList<String>());
				transcript_domain_mapping.get(rs.getString(1)).add(rs.getString(2));
				
				String biotype = rs.getString(5);
				if (biotype.equals("nonsense_mediated_decay") || biotype.equals("non_stop_decay")) { // see http://www.gencodegenes.org/gencode_biotypes.html
					if (!DataQuery.cache_decay_transcripts.containsKey(organism_core_database))
						DataQuery.cache_decay_transcripts.put(organism_core_database, new HashSet<String>());
					DataQuery.cache_decay_transcripts.get(organism_core_database).add(rs.getString(1));
				}
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get transcript domain data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getTranscriptsDomains(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		// sort for faster checks later
		for (List<String> domain_annotations:transcript_domain_mapping.values())
			Collections.sort(domain_annotations);
		
		DataQuery.cache_transcrdom.put(organism_core_database, transcript_domain_mapping);
		return transcript_domain_mapping;
	}
	
	/**
	 * Queries Ensembl for length information of transcripts
	 * @param organism_core_database
	 * @return map of Ensembl transcripts to length of exons in bp
	 */
	public static Map<String, Integer> getTranscriptsCDNALength(String organism_core_database) {
		
		// check cache
		if (DataQuery.cache_transcript_length.containsKey(organism_core_database))
			return DataQuery.cache_transcript_length.get(organism_core_database);
		
		Map <String, Integer> transcript_bp_count = new HashMap<>();
		Connection connection = null;
		try {
			connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/"+organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			ResultSet rs = st.executeQuery("SELECT transcript.stable_id, exon.seq_region_start, exon.seq_region_end FROM transcript, exon_transcript, exon "
					+ "WHERE transcript.transcript_id = exon_transcript.transcript_id AND exon.exon_id = exon_transcript.exon_id");
			
			while (rs.next()) {
				String transcript = rs.getString(1);
				int start = rs.getInt(2);
				int end = rs.getInt(3);
				int length = 1 + end - start;
				transcript_bp_count.put(transcript, transcript_bp_count.getOrDefault(transcript, 0) + length);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("ENSEMBL");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get transcript length data from ENSEMBL in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			DataQuery.switchServer();
			
			return getTranscriptsCDNALength(organism_core_database);
			
		} finally {
			try {
				connection.close();
			} catch (Exception e) {
			}
		}
		
		// caching
		DataQuery.cache_transcript_length.put(organism_core_database, transcript_bp_count);
		
		return transcript_bp_count;
	}
	
	/**
	 * Returns those translated Ensembl transcripts that are associated with either "nonsense_mediated_decay" or "non_stop_decay".
	 * @param organism_core_database
	 * @return
	 */
	public static Set<String> getDecayTranscripts(String organism_core_database) {
		if (!DataQuery.cache_transcrdom.containsKey(organism_core_database))
			getTranscriptsDomains(organism_core_database);
		
		return DataQuery.cache_decay_transcripts.get(organism_core_database);
	}
	
	/**
	 * Builds a map of UniProt -> all domain types found in transcripts
	 * @param organism_core_database
	 * @returnreturn map of UniProt -> all domain types found in transcripts
	 */
	public static Map<String, Set<String>> getHolisticProteinDomainMap(String organism_core_database) {
		Map<String, List<String>> transcript_domain_mapping = getTranscriptsDomains(organism_core_database);
		List<String[]> associations = getGenesTranscriptsProteins(organism_core_database);
		Map<String, Set<String>> holistic_protein_domain_mapping = new HashMap<>();
		for (String[] array:associations) {
			String protein = array[2];
			String transcript = array[1];
			
			if (!holistic_protein_domain_mapping.containsKey(protein))
				holistic_protein_domain_mapping.put(protein, new HashSet<String>());
			
			if (transcript_domain_mapping.containsKey(transcript)) {
				holistic_protein_domain_mapping.get(protein).addAll(transcript_domain_mapping.get(transcript));
			} else {
				holistic_protein_domain_mapping.get(protein).addAll(new LinkedList<String>());
			}
			
		}
		
		return holistic_protein_domain_mapping;
	}
	
	/**
	 * Builds a map of UniProt -> all domain types found in isoform transcripts
	 * @param organism_core_database
	 * @returnreturn map of UniProt -> all domain types found in isoform transcripts
	 */
	public static Map<String, Set<String>> getIsoformProteinDomainMap(String organism_core_database) {
		Map<String, List<String>> transcript_domain_mapping = getTranscriptsDomains(organism_core_database);
		Map<String, String> isoform_map = getIsoformTranscriptsOfProteins(organism_core_database);
		Map<String, Set<String>> isoform_protein_domain_mapping = new HashMap<>();
		
		for (String protein:isoform_map.keySet()) {
			String isoform = isoform_map.get(protein);
			
			if (transcript_domain_mapping.containsKey(isoform)) {
				isoform_protein_domain_mapping.put(protein, new HashSet<>(transcript_domain_mapping.get(isoform)));
			} else {
				isoform_protein_domain_mapping.put(protein, new HashSet<String>());
			}
			
		}
		
		return isoform_protein_domain_mapping;
	}
	
	/**
	 * Queries EBI QuickGO to get all UniProt proteins that are annotated with or descendants of a certain GO_id of an organism with certain taxon
	 * @param GO_id
	 * @param taxon
	 * @param include_IEA
	 * @return
	 */
	public static Set<String> getProteinsWithGO(String GO_id, String taxon, boolean include_IEA, boolean only_experimental, boolean report_genes) {
		
		Set<String> entries = new HashSet<>();
		BufferedReader datastream = null;
		try {
			URL server = new URL("https://www.ebi.ac.uk/QuickGO-Old/GAnnotation?format=tsv&gz=true&limit=-1&db=UniProtKB&tax="+ taxon +"&goid=" + GO_id);
			HttpURLConnection connection = (HttpURLConnection) server.openConnection();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(new GZIPInputStream(connection.getInputStream())));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty() || line.startsWith("DB"))
					continue;
				String[] temp = line.trim().split("\t");
				String inf = temp[9];
				
				// not IEA
				if (!include_IEA && inf.equals("IEA"))
					continue;
				
				// manual experimental: IMP,IGI,IPI,IDA,IEP,EXP
				if (only_experimental && !inf.equals("IMP") && !inf.equals("IGI") && !inf.equals("IPI") && !inf.equals("IDA") && !inf.equals("IEP") && !inf.equals("EXP"))
					continue;
				
				if (report_genes)
					entries.add(temp[3]);
				else
					entries.add(temp[1]);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			
			if (DataQuery.retries == 10)
				terminateRetrieval("QuickGO");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get GO data from QuickGO in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getProteinsWithGO(GO_id, taxon, include_IEA, only_experimental, report_genes);
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}
		
		return entries;
	}
	
	/**
	 * Queries EBI QuickGO to get all UniProt proteins that are annotated with or descendants of a certain GO_id of an organism with certain taxon,
	 * automatically includes also the ones with evidence code IEA
	 * @param GO_id
	 * @param taxon
	 * @return
	 */
	public static Set<String> getProteinsWithGO(String GO_id, String taxon) {
		return getProteinsWithGO(GO_id, taxon, true, false, false);
	}
	
	
	/**
	 * Retrieves STRING v10 network
	 * @param taxon_id
	 * @return
	 */
	public static PPIN getSTRINGNetwork(String taxon_id) {
		
		// try to get from cache
		if (DataQuery.cache_STRING.containsKey(taxon_id))
			return DataQuery.cache_STRING.get(taxon_id);
		
		StringBuilder output = new StringBuilder();
		String organism_core_database = getEnsemblOrganismDatabaseFromName(getOrganismFromTaxon(taxon_id));
		Map<String, List<String>> ensembl_to_uniprots = DataQuery.getEnsemblToAllUniprotProteins(organism_core_database);
		
		BufferedReader datastream = null;
		try {
			// URL is version-dependent
			String version = DataQuery.STRING_version.split("\\.")[0];
			if (DataQuery.STRING_version.split("\\.").length > 1)
				version += DataQuery.STRING_version.split("\\.")[1];
			URL server = new URL("http://string" + version +".embl.de/newstring_download/protein.links.v"+DataQuery.STRING_version+"/"+taxon_id+".protein.links.v"+DataQuery.STRING_version+".txt.gz");
			URLConnection connection = server.openConnection();
			
			// read and build pipe
			datastream = new BufferedReader(new InputStreamReader(new GZIPInputStream(connection.getInputStream())));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty() || line.startsWith("protein1"))
					continue;
				String[] temp = line.trim().split("\\s+");
				String protein1 = temp[0].split("\\.")[1];
				String protein2 = temp[1].split("\\.")[1];
				float P = Float.parseFloat(temp[2]) / 1000f;
				
				if (!ensembl_to_uniprots.containsKey(protein1) || !ensembl_to_uniprots.containsKey(protein2))
					continue;
				for (String up1:ensembl_to_uniprots.get(protein1))
					for (String up2:ensembl_to_uniprots.get(protein2))
						output.append(up1 + " " + up2 + " " + P + "\n");
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("STRING");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get interaction data from STRING in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getSTRINGNetwork(taxon_id);
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}
		
		// convert, stream is closed in PPIN constructor
		InputStream inputStream = new ByteArrayInputStream(output.toString().getBytes(Charset.forName("UTF-8")));
		
		// build
		PPIN STRING_network = new PPIN(new BufferedReader(new InputStreamReader(inputStream)), 0.0);
		
		// store in cache
		DataQuery.cache_STRING.put(taxon_id, STRING_network);
		
		return STRING_network;
	}
	
	/**
	 * Retrieves current IntAct network for a given organism; if there is no data for the taxon, an empty network is returned
	 * @param taxon_id
	 * @return
	 */
	public static PPIN getIntActNetwork(String taxon_id) {
		return getIntActNetwork(taxon_id, null);
	}
	
	/**
	 * Retrieves current IntAct network for a given organism, prints feedback to ps;  if there is no data for the taxon, an empty network is returned
	 * @param taxon_id
	 * @param ps
	 * @return
	 */
	public static PPIN getIntActNetwork(String taxon_id, PrintStream ps) {
		
		StringBuilder output = new StringBuilder();
		String line = null;
		ZipInputStream zis = null;
		try {
			URL server = new URL("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip");
			URLConnection connection = server.openConnection();
			
			// read and build pipe
			zis = new ZipInputStream(connection.getInputStream());
			ZipEntry file_in_zip = null;
			
			// for percentage output
			long size = (long) (connection.getContentLengthLong() * 9.97); // approximation, gzip etc ...
			long read = 0;
			int line_sep_len = System.getProperty("line.separator").getBytes().length;
			byte state = 0;
			while ( (file_in_zip = zis.getNextEntry()) != null)
				if (file_in_zip.getName().equals("intact.txt")) {
					BufferedReader datastream = new BufferedReader(new InputStreamReader(zis));
					while ( (line = datastream.readLine()) != null ) {
						read += line.getBytes().length + line_sep_len;
						
						// for output
						if (ps != null) {
							if (state == 0 && read > 0.24 * size) {
								state = 1;
								ps.print("25% ... ");
							} else if (state == 1 && read > 0.49 * size) {
								state = 2;
								ps.print("50% ... ");
							} else if (state == 2 && read > 0.74 * size) {
								state = 3;
								ps.print("75% ... ");
							}
						}
						
						if (line.isEmpty() || line.startsWith("#"))
							continue;
						
						String[] temp = line.split("\\t");
						
						// filter small molecules?
						if (temp[9].equals("-") || temp[10].equals("-"))
							continue;
						
						String tax1 = temp[9].split("\\|")[0].split(":")[1].split("\\(")[0];
						String tax2 = temp[10].split("\\|")[0].split(":")[1].split("\\(")[0];
						
						// both proteins must be in organism
						if (!tax1.equals(taxon_id) || !tax2.equals(taxon_id))
							continue;
						
						String[] p1 = temp[0].split(":");
						String[] p2 = temp[1].split(":");
						
						// both proteins must be given as uniprot accs
						if (!p1[0].equals("uniprotkb") || !p2[0].equals("uniprotkb"))
							continue;
						
						output.append(p1[1].split("-")[0] + " " + p2[1].split("-")[0] + "\n");
					}
					datastream.close();
					break;
				}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("IntAct");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get interaction data from IntAct in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getIntActNetwork(taxon_id, ps);
			
		} finally {
			try {
				zis.close();
			} catch (Exception e) {
			}
		}
		
		// convert, stream closed in PPIN constructor
		InputStream inputStream = new ByteArrayInputStream(output.toString().getBytes(Charset.forName("UTF-8")));
		
		// build
		
		return new PPIN(new BufferedReader(new InputStreamReader(inputStream)), 0.0);
	}
	
	/**
	 * Retrieves current iRefIndex network of physical interactions for a given organism;  if there is no data for the taxon, an empty network is returned
	 * @param taxon_id
	 * @return
	 */
	public static PPIN getIRefIndexNetwork(String taxon_id) {
		return getIRefIndexNetwork(taxon_id, null);
	}
	
	/**
	 * Retrieves current iRefIndex network of physical interactions for a given organism, prints feedback to ps;  if there is no data for the taxon, an empty network is returned
	 * @param taxon_id
	 * @param ps
	 * @return
	 */
	public static PPIN getIRefIndexNetwork(String taxon_id, PrintStream ps) {
		
		StringBuilder output = new StringBuilder();
		String line = null;
		
		ZipInputStream zis = null;
		try {
			// get actual filename
			URL server = new URL("http://irefindex.org/download/irefindex/data/current/psi_mitab/MITAB2.6/");
			URLConnection connection = server.openConnection();
			
			BufferedReader datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			
			while ( (line = datastream.readLine()) != null ) {
				output.append(line);
			}
			
			// parse html data
			String net_file = "";
			for (String s:output.toString().split("\"")) {
				if (s.startsWith(taxon_id) && s.contains("mitab") && s.endsWith(".txt.zip"))
					net_file = s;
			}
			
			// use the knowledge of the filename to get the data
			server = new URL("http://irefindex.org/download/irefindex/data/current/psi_mitab/MITAB2.6/" + net_file);
			connection = server.openConnection();
			output = new StringBuilder();
			
			// read and build pipe
			zis = new ZipInputStream(connection.getInputStream());
			ZipEntry file_in_zip = null;
			
			
			// for percentage output
			long size = (long) (connection.getContentLengthLong() * 9.97); // approximation, gzip etc ...
			long read = 0;
			int line_sep_len = System.getProperty("line.separator").getBytes().length;
			byte state = 0;
			
			while ( (file_in_zip = zis.getNextEntry()) != null)
				if (file_in_zip.getName().endsWith(".txt")) {
					datastream = new BufferedReader(new InputStreamReader(zis));
					while ( (line = datastream.readLine()) != null ) {
						read += line.getBytes().length + line_sep_len;
						
						// for output
						if (ps != null) {
							if (state == 0 && read > 0.24 * size) {
								state = 1;
								ps.print("25% ... ");
							} else if (state == 1 && read > 0.49 * size) {
								state = 2;
								ps.print("50% ... ");
							} else if (state == 2 && read > 0.74 * size) {
								state = 3;
								ps.print("75% ... ");
							}
						}
						
						// actual parsing
						if (line.isEmpty() || line.startsWith("#"))
							continue;
						
						String[] temp = line.split("\\t");
						
						String[] p1 = temp[0].split(":");
						String[] p2 = temp[1].split(":");
						
						// both proteins must be given as uniprot accs
						if (!p1[0].equals("uniprotkb") || !p2[0].equals("uniprotkb"))
							continue;
						// both proteins must be from organism
						if (!temp[9].contains(taxon_id) || !temp[10].contains(taxon_id))
							continue;
						
						// no unclear and unwanted interaction types
						String type = temp[11];
						if (!type.equals("MI:0914(association)") && !type.equals("MI:0407(direct interaction)") && !type.equals("MI:0914(association)") && !type.equals("MI:0195(covalent binding)") && !type.equals("MI:0915(physical association)") && !type.equals("MI:0000(psi-mi:\"MI:0407\")"))
							continue;

						output.append(p1[1].split("-")[0] + " " + p2[1].split("-")[0] + "\n");
					}
					break;
				}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("iRefIndex");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get interaction data from iRefIndex in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getIRefIndexNetwork(taxon_id, ps);
			
		} finally {
			try {
				zis.close();
			} catch (Exception e) {
			}
		}
		
		// convert, stream closed in PPIN constructor
		InputStream inputStream = new ByteArrayInputStream(output.toString().getBytes(Charset.forName("UTF-8")));
		
		// build
		
		return new PPIN(new BufferedReader(new InputStreamReader(inputStream)), 0.0);
	}

	/**
	 * Retrieves current weighted mentha network for a given organism; if there is no data for the taxon, an empty network is returned
	 * @param taxon_id
	 * @return
	 */
	public static PPIN getMenthaNetwork(String taxon_id) {
		return getMenthaNetwork(taxon_id, null);
	}
	
	/**
	 * Retrieves current weighted mentha network for a given organism, prints feedback to ps;  if there is no data for the taxon, an empty network is returned
	 * @param taxon_id
	 * @param ps
	 * @return
	 */
	public static PPIN getMenthaNetwork(String taxon_id, PrintStream ps) {
		
		StringBuilder output = new StringBuilder();
		String line = null;
		ZipInputStream zis = null;
		try {
			URL server = new URL("http://mentha.uniroma2.it/dumps/organisms/" + taxon_id + ".zip");
			URLConnection connection = server.openConnection();
			
			// read and build pipe
			zis = new ZipInputStream(connection.getInputStream());
			ZipEntry file_in_zip = null;
			
			// for percentage output
			long size = (long) (connection.getContentLengthLong() * 3.6); // approximation
			long read = 0;
			int line_sep_len = System.getProperty("line.separator").getBytes().length;
			byte state = 0;
			while ( (file_in_zip = zis.getNextEntry()) != null)
				if (file_in_zip.getName().equals(taxon_id)) {
					BufferedReader datastream = new BufferedReader(new InputStreamReader(zis));
					while ( (line = datastream.readLine()) != null ) {
						read += line.getBytes().length + line_sep_len;
						
						// for output
						if (ps != null) {
							if (state == 0 && read > 0.24 * size) {
								state = 1;
								ps.print("25% ... ");
							} else if (state == 1 && read > 0.49 * size) {
								state = 2;
								ps.print("50% ... ");
							} else if (state == 2 && read > 0.74 * size) {
								state = 3;
								ps.print("75% ... ");
							}
						}
						
						if (line.startsWith("Protein A") || line.isEmpty() || line.startsWith("#"))
							continue;
						
						String[] temp = line.split(";");
						String p1 = temp[0];
						String p2 = temp[2];
						String score = temp[4];
						
						output.append(p1 + " " + p2 + " " + score + "\n");
					}
					datastream.close();
					break;
				}
			
		} catch (FileNotFoundException e) {
			// data not available -> return empty PPIN
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("Mentha");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get interaction data from mentha in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getMenthaNetwork(taxon_id, ps);
			
		} finally {
			try {
				zis.close();
			} catch (Exception e) {
			}
		}
		
		// convert, stream closed in PPIN constructor
		InputStream inputStream = new ByteArrayInputStream(output.toString().getBytes(Charset.forName("UTF-8")));
		
		// build
		
		return new PPIN(new BufferedReader(new InputStreamReader(inputStream)), 0.0);
	}
	
	/**
	 * Retrieves a mapping of secondary to primary accessions from Uniprot,
	 * enables to filter/update data with old annotation
	 * @return
	 */
	public static Map<String, String> getUniprotSecondaryAccMap() {
		
		// cache
		if (DataQuery.uniprot_sec_accs != null)
			return DataQuery.uniprot_sec_accs;
		
		Map<String, String> sec_to_primary = new HashMap<>();
		BufferedReader datastream = null;
		try {
			URL server = new URL("ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/sec_ac.txt");
			URLConnection connection = server.openConnection();
			
			// read and build pipe
			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			String line;
			
			boolean no_data_yet = true;
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty())
					continue;
				
				if (no_data_yet) {
					if (line.startsWith("Release:")) {
						// parse version
						String[] temp = line.trim().split("\\s+");
						DataQuery.uniprot_release = temp[1];
					} else if (line.startsWith("_")) {
						no_data_yet = false;
					}
					
					continue;
				}
				
				// parse data
				String[] temp = line.trim().split("\\s+");
				sec_to_primary.put(temp[0], temp[1]);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("UniProt");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get accession data from Uniprot in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
			}
			
			return getUniprotSecondaryAccMap();
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}
			
		// write to cache
		DataQuery.uniprot_sec_accs = sec_to_primary;
		
		return sec_to_primary;
	}
	
	/**
	 * Returns release of Uniprot accession mapping data
	 */
	public static String getUniprotRelease() {
		if (DataQuery.uniprot_release == null)
			DataQuery.getUniprotSecondaryAccMap();
		
		return DataQuery.uniprot_release;
	}
	
	/**
	 * Generic stuff
	 */
	
	/**
	 * Enforces to use databases of a specific Ensembl release, if possible.
	 * @param release
	 */
	public static void enforceSpecificEnsemblRelease(String release) {
		DataQuery.specific_ensembl_release = release;
	}
	
	/**
	 * Enforces to use a specific 3did release, default: current -> "current"
	 * @param release
	 */
	public static void enforceSpecific3didRelease(String release) {
		DataQuery.specific_3did_release = release;
	}
	
	/**
	 * Clears cached results
	 */
	public static void resetCaches() {
		DataQuery.cache_genestransprots.clear();
		DataQuery.cache_isoformtranscr.clear();
		DataQuery.cache_transcrdom.clear();
		DataQuery.cache_decay_transcripts.clear();
		DataQuery.cache_ensembl_proteins.clear();
		DataQuery.cache_ensembl_names.clear();
		DataQuery.cache_db.clear();
		DataQuery.cache_STRING.clear();
		DataQuery.cache_ucsc.clear();
		DataQuery.cache_ensembl_db.clear();
	}
	
	/**
	 * Switches to another Ensembl server
	 */
	public static void switchServer() {
		
		long now = System.currentTimeMillis();
		
		// ensures that last change is at least 10 seconds ago
		if ((now - DataQuery.last_server_change) < 10000)
			return;
		
		DataQuery.last_server_change = now;
		
		if (DataQuery.ensembl_mysql.equals("ensembldb.ensembl.org:3337")) {
			err_out.println("Special GRCh37 ENSEMBL server is used and cannot be switched.");
			return;
		}
		
		String server1 = "ensembldb.ensembl.org:3306";
		String server2 = "useastdb.ensembl.org:3306";
		String server3 = "asiadb.ensembl.org:3306";
		
		if (DataQuery.ensembl_mysql.equals(server1))
			DataQuery.ensembl_mysql = server2;
		else if (DataQuery.ensembl_mysql.equals(server2))
			DataQuery.ensembl_mysql = server3;
		else // server 3
			DataQuery.ensembl_mysql = server1;
		
		err_out.println("ENSEMBL server changed to " + DataQuery.ensembl_mysql + ".");
	}
	
	/**
	 * Switches to another Ensembl server
	 */
	public static void switchServer(String server) {
		
		if (!DataQuery.ensembl_mysql.equals(server)) {
			DataQuery.ensembl_mysql = server;
			System.out.println("ENSEMBL server manually changed to " + DataQuery.ensembl_mysql + ".");
		}
		
	}
	
	/**
	 * Switches to GRCh37 Ensembl server
	 */
	public static void switchServerGRCh37() {
		
		if (!DataQuery.ensembl_mysql.equals("ensembldb.ensembl.org:3337")) {
			DataQuery.ensembl_mysql = "ensembldb.ensembl.org:3337";
			System.out.println("ENSEMBL server automatically changed to " + DataQuery.ensembl_mysql + " to retrieve GRCh37 data.");
		}
		
	}
	
	/**
	 * Is called after too many connection attempts.
	 * @param service
	 */
	public static void terminateRetrieval(String service) {
		err_out.println("Too many connection attempts while querying " + service + ".");
		err_out.println("Please check your internet connection or try again later since service might undergo update or maintenance.");
		System.exit(1);
	}
	
	/**
	 * Extends MySQL timeout
	 */
	public static void extendTimeout(int timeout) {
		DataQuery.timeout = timeout;
	}
	
	/**
	 * Set error out manually
	 */
	public static void setErrOut(PrintStream err_out) {
		DataQuery.err_out = err_out;
	}
	
	/**
	 * Retrieve 3did/iPfam and use normal local DDI data
	 */
	public static void defaultDDIs() {
		DataQuery.up2date_DDIs = true;
		DataQuery.no_local_DDIs = false;
		if (DataQuery.known_DDIs != null)
			DataQuery.known_DDIs = null;
	}
	
	/**
	 * Turn off iPfam retrieval, only local DOMINE/IDDI data used
	 */
	public static void localDDIsOnly() {
		DataQuery.up2date_DDIs = false;
		if (DataQuery.known_DDIs != null)
			DataQuery.known_DDIs = null;
	}
	
	/**
	 * No local data, only iPfam/3did
	 */
	public static void onlyRetrievedDDIs() {
		DataQuery.up2date_DDIs = true;
		DataQuery.no_local_DDIs = true;
		if (DataQuery.known_DDIs != null)
			DataQuery.known_DDIs = null;
	}
	
	/**
	 * Sets stricter threshold to precompiled IDDI/DOMINE data
	 */
	public static void stricterLocalDDIs() {
		DataQuery.stricter_local_DDIs = true;
		if (DataQuery.known_DDIs != null)
			DataQuery.known_DDIs = null;
	}
	
	/*
	 *  Reading DDIs from here
	 */
	
	/**
	 * Uses IDDI 1.0 and DOMINE 2.0, additionally most current iPfam if wanted (default: yes)
	 * @return compiled map of known interactions among Pfam domains
	 */
	public static Map<String, List<String>> getKnownDDIs() {
		
		// never read twice
		if (known_DDIs != null)
			return known_DDIs;
		
		known_DDIs = new HashMap<String, List<String>>();
		
		HashMap<String, HashSet<String>> domine_ddis = getDOMINE();
		HashMap<String, HashSet<String>> iddi_ddis = getIDDI();
		HashMap<String, HashSet<String>> ipfam_ddis = getIPfam();
		
		// merge in iddi_ddis
		for (String domain:domine_ddis.keySet()) {
			// merge if there
			if (iddi_ddis.containsKey(domain)) {
				iddi_ddis.get(domain).addAll(domine_ddis.get(domain));
			} else { // just add if not there
				iddi_ddis.put(domain, domine_ddis.get(domain));
			}
		}
		
		// merge iPfam DDIs
		for (String domain:ipfam_ddis.keySet()) {
			// merge if there
			if (iddi_ddis.containsKey(domain)) {
				iddi_ddis.get(domain).addAll(ipfam_ddis.get(domain));
			} else { // just add if not there
				iddi_ddis.put(domain, ipfam_ddis.get(domain));
			}
		}
		
		if (no_local_DDIs)
			iddi_ddis.clear();
		
		// add current 3did data
		if (DataQuery.up2date_DDIs) {
			
			HashMap<String, HashSet<String>> did_ddis = get3did();
			
			for (String domain:did_ddis.keySet()) {
				// merge if there
				if (iddi_ddis.containsKey(domain)) {
					iddi_ddis.get(domain).addAll(did_ddis.get(domain));
				} else { // just add if not there
					iddi_ddis.put(domain, did_ddis.get(domain));
				}
			}
			
		}
		
		// to list
		for (String domain:iddi_ddis.keySet()) {
			known_DDIs.put(domain, new ArrayList<>(iddi_ddis.get(domain)));
		}
		
		return known_DDIs;
	}
	
	/**
	 * Reads DOMINE data
	 */
	public static HashMap<String, HashSet<String>> getDOMINE() {
		HashMap<String, HashSet<String>> ddis = new HashMap<>();
		
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(DataQuery.class.getResourceAsStream("/data/DOMINE_2.0_HC.txt.gz"))));
			while (in.ready()) {
				String line = in.readLine();
				if (line.startsWith("PFAM1"))
					continue;
				String[] temp = line.split("\\|");
				String confidence = temp[temp.length-2];
				
				// high confidence and known interactions only
				if (!confidence.equals("HC") && !confidence.equals("NA"))
					continue;
				if (stricter_local_DDIs && !confidence.equals("NA"))
					continue;
				
				String d1 = temp[0];
				String d2 = temp[1];
				if (!ddis.containsKey(d1))
					ddis.put(d1, new HashSet<String>());
				if (!ddis.containsKey(d2))
					ddis.put(d2, new HashSet<String>());
				ddis.get(d1).add(d2);
				ddis.get(d2).add(d1);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		return ddis;
	}
	
	/**
	 * Reads IDDI data
	 */
	public static HashMap<String, HashSet<String>> getIDDI() {
		HashMap<String, HashSet<String>> ddis = new HashMap<>();
		
		// from original paper: 0.329 for 98% accuracy, 0.102 for 90% ; one MET4-partner in yeast: 0.20580978 needed, for example
		double threshold = 0.102;
		if (stricter_local_DDIs)
			threshold = 0.329;
		
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(DataQuery.class.getResourceAsStream("/data/IDDI_1.0_90.txt.gz"))));
			while (in.ready()) {
				String line = in.readLine();
				if (line.startsWith("PFAM1"))
					continue;
				String[] temp = line.split("\\s+");
				float score = Float.parseFloat(temp[temp.length-1]);
				
				
				
				if (score <= threshold)
					continue;
				String d1 = temp[0];
				String d2 = temp[1];
				if (!ddis.containsKey(d1))
					ddis.put(d1, new HashSet<String>());
				if (!ddis.containsKey(d2))
					ddis.put(d2, new HashSet<String>());
				ddis.get(d1).add(d2);
				ddis.get(d2).add(d1);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			
		} finally {
			try {
				in.close();
			} catch (Exception e) {
			}
		}
		
		return ddis;
	}
	
	/**
	 * Returns local version of iPfam data, compatibility function
	 */
	public static String getIPfamVersion() {
		return "release 1.0";
	}
	
	/**
	 * Retrieves local iPfam 1.0 interaction data
	 * @return interaction data
	 */
	public static HashMap<String, HashSet<String>> getIPfam() {

		HashMap<String, HashSet<String>> ddis = new HashMap<>();
		
		BufferedReader datastream = null;
		try {
			
			// read
			datastream = new BufferedReader(new InputStreamReader(new GZIPInputStream(DataQuery.class.getResourceAsStream("/data/iPfam_1.0_heterodomain_interaction.csv.gz"))));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty())
					continue;
				if (line.startsWith("family"))
					continue;
				
				String[] temp = line.split("\\s+");
				String d1 = temp[0];
				String d2 = temp[2];
				if (!ddis.containsKey(d1))
					ddis.put(d1, new HashSet<String>());
				if (!ddis.containsKey(d2))
					ddis.put(d2, new HashSet<String>());
				ddis.get(d1).add(d2);
				ddis.get(d2).add(d1);
			}
			
			datastream.close();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(new GZIPInputStream(DataQuery.class.getResourceAsStream("/data/iPfam_1.0_homodomain_interaction.csv.gz"))));
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty())
					continue;
				
				if (line.startsWith("family"))
					continue;
				
				String[] temp = line.split("\\s+");
				String d1 = temp[0];
				if (!ddis.containsKey(d1))
					ddis.put(d1, new HashSet<String>());
				ddis.get(d1).add(d1);
			}
						
		} catch (Exception e) {
			e.printStackTrace();
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}
		
		return ddis;
	}
	
	// iPfam retrieval (if again updated)
//	/**
//	 * Returns version of iPfam data
//	 */
//	public static String getIPfamVersion() {
//		String version_string = "";
//		
//		BufferedReader datastream = null;
//		try {
//			URL server = new URL("ftp://selab.janelia.org/pub/ipfam/Current_Release/relnotes.txt");
//			URLConnection connection = server.openConnection();
//			
//			// read
//			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
//			String line;
//			
//			while ( (line = datastream.readLine()) != null ) {
//				if (line.isEmpty())
//					continue;
//				line = line.trim();
//				
//				if (line.startsWith("RELEASE")) {
//					version_string = line;
//					break;
//				}
//			}
//
//		} catch (Exception e) {
//			if (DataQuery.retries == 10)
//				terminateRetrieval("iPfam");
//			
//			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get release data from iPfam in 10 seconds ..." );
//			try {
//				Thread.sleep(10000);
//			} catch (InterruptedException e1) {
//				e1.printStackTrace();
//			}
//			
//			return getIPfamVersion();
//			
//		} finally {
//			try {
//				datastream.close();
//			} catch (Exception e) {
//			}
//		}
//		
//		return version_string;
//	}
//	
//	/**
//	 * Retrieves current iPfam interaction data
//	 * @return interaction data
//	 */
//	public static HashMap<String, HashSet<String>> getIPfam() {
//
//		HashMap<String, HashSet<String>> ddis = new HashMap<>();
//		
//		BufferedReader datastream = null;
//		try {
//			URL server = new URL("ftp://selab.janelia.org/pub/ipfam/Current_Release/heterodomain_interaction.csv");
//			URLConnection connection = server.openConnection();
//			
//			// read
//			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
//			String line;
//			
//			while ( (line = datastream.readLine()) != null ) {
//				if (line.isEmpty())
//					continue;
//				if (line.startsWith("family"))
//					continue;
//				
//				String[] temp = line.split("\\s+");
//				String d1 = temp[0];
//				String d2 = temp[2];
//				if (!ddis.containsKey(d1))
//					ddis.put(d1, new HashSet<String>());
//				if (!ddis.containsKey(d2))
//					ddis.put(d2, new HashSet<String>());
//				ddis.get(d1).add(d2);
//				ddis.get(d2).add(d1);
//			}
//			
//			datastream.close();
//			
//			server = new URL("ftp://selab.janelia.org/pub/ipfam/Current_Release/homodomain_interaction.csv");
//			connection = server.openConnection();
//			
//			// read
//			datastream = new BufferedReader(new InputStreamReader(connection.getInputStream()));
//			
//			while ( (line = datastream.readLine()) != null ) {
//				if (line.isEmpty())
//					continue;
//				
//				if (line.startsWith("family"))
//					continue;
//				
//				String[] temp = line.split("\\s+");
//				String d1 = temp[0];
//				if (!ddis.containsKey(d1))
//					ddis.put(d1, new HashSet<String>());
//				ddis.get(d1).add(d1);
//			}
//						
//		} catch (Exception e) {
//			if (DataQuery.retries == 10)
//				terminateRetrieval("iPfam");
//			
//			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get interaction data from iPfam in 10 seconds ..." );
//			try {
//				Thread.sleep(10000);
//			} catch (InterruptedException e1) {
//				e1.printStackTrace();
//			}
//			
//			return getIPfam();
//			
//		} finally {
//			try {
//				datastream.close();
//			} catch (Exception e) {
//			}
//		}
//		
//		return ddis;
//	}
	
	/**
	 * Returns version of 3did data
	 */
	public static String get3didVersion() {
		String version_string = "unknown";
		Date date = null;
		
		try {
			URL server = new URL("https://3did.irbbarcelona.org/download/" + DataQuery.specific_3did_release + "/3did_flat.gz");
			URLConnection connection = server.openConnection();
			date = new Date(connection.getLastModified());
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("3did");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get release data from 3did in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
			
			return get3didVersion();
		} 
		
		if (date != null) {
			String[] temp = date.toString().split("\\s+");
			version_string = temp[1] + " " + temp[5];
		}
		
		return version_string;
	}
	
	/**
	 * Retrieval of 3did interaction data
	 * @return interaction data
	 */
	public static HashMap<String, HashSet<String>> get3did() {
		
		HashMap<String, HashSet<String>> ddis = new HashMap<>();
		
		BufferedReader datastream = null;
		try {
			URL server = new URL("https://3did.irbbarcelona.org/download/" + DataQuery.specific_3did_release + "/3did_flat.gz");
			URLConnection connection = server.openConnection();
			
			// read
			datastream = new BufferedReader(new InputStreamReader(new GZIPInputStream(connection.getInputStream())));
			String line;
			
			while ( (line = datastream.readLine()) != null ) {
				if (line.isEmpty())
					continue;
				if (!line.startsWith("#=ID"))
					continue;
				
				String[] temp = line.split("\\s+");
				
				String d1 = temp[3].replace("(", "").split("\\.")[0];
				String d2 = temp[4].split("\\.")[0];
				
				if (!ddis.containsKey(d1))
					ddis.put(d1, new HashSet<String>());
				if (!ddis.containsKey(d2))
					ddis.put(d2, new HashSet<String>());
				ddis.get(d1).add(d2);
				ddis.get(d2).add(d1);
			}
			
		} catch (Exception e) {
			if (DataQuery.retries == 10)
				terminateRetrieval("3did");
			
			err_out.println("Attempting " + (++DataQuery.retries) +". retry to get interaction data from 3did in 10 seconds ..." );
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
			
			return get3did();
			
		} finally {
			try {
				datastream.close();
			} catch (Exception e) {
			}
		}

		return ddis;
	}
	
	/*
	 * concurrent functions
	 */
	
	static class getGenesTranscriptsTask implements Callable<Map<String,List<String>>> {
		final String organism_core_database;
		
		public getGenesTranscriptsTask(String organism_core_database) {
			this.organism_core_database = organism_core_database;
		}
		
		@Override
		public Map<String, List<String>> call() throws Exception {
			Map<String,List<String>> gene_transcripts = new HashMap<>();
			Connection connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/" + this.organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			
			ResultSet rs = st.executeQuery("SELECT gene.stable_id, transcript.stable_id "
					+ "FROM gene, transcript "
					+ "WHERE gene.gene_id=transcript.gene_id");
			
			while (rs.next()) {
				String[] row = {rs.getString(1), rs.getString(2)};
				
				if (!gene_transcripts.containsKey(row[0]))
					gene_transcripts.put(row[0], new LinkedList<String>());
				
				gene_transcripts.get(row[0]).add(row[1]);
			}
			
			connection.close();
			
			return gene_transcripts;
		}
	}
	
	static class getTranscriptProteinsTask implements Callable<Map<String,List<String>>> {
		final String organism_core_database;
		
		public getTranscriptProteinsTask(String organism_core_database) {
			this.organism_core_database = organism_core_database;
		}
		
		@Override
		public Map<String, List<String>> call() throws Exception {
			
			Map<String,List<String>> transcript_proteins = new HashMap<>();
			Connection connection = DriverManager.getConnection("jdbc:mysql://"+ensembl_mysql+"/" + this.organism_core_database, "anonymous", "");
			Statement st = connection.createStatement();
			st.setQueryTimeout(timeout);
			
			ResultSet rs = st.executeQuery("SELECT transcript.stable_id, trans_table.dbprimary_acc FROM transcript, (SELECT translation.translation_id, xref.dbprimary_acc FROM translation, object_xref, xref WHERE translation.translation_id=object_xref.ensembl_id AND object_xref.xref_id=xref.xref_id AND xref.external_db_id='2200') AS trans_table WHERE transcript.canonical_translation_id = trans_table.translation_id");
			
			while (rs.next()) {
				String[] row = {rs.getString(1), rs.getString(2)};
				
				if (!transcript_proteins.containsKey(row[0]))
					transcript_proteins.put(row[0], new LinkedList<String>());
				
				transcript_proteins.get(row[0]).add(row[1]);
			}
			
			connection.close();
			
			return transcript_proteins;
		}
	}
}
