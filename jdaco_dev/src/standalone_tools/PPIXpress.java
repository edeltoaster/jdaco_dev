package standalone_tools;

import java.io.File;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

/**
 * PPIXpress cmd-tool
 * @author Thorsten Will
 */
public class PPIXpress {
	static String version_string = "PPIXpress 1.18";
	
	private static boolean gene_level_only = false;
	private static boolean output_DDINs = false;
	private static boolean output_major_transcripts = false;
	private static String original_network_path;
	private static List<String> input_files = new LinkedList<>();
	private static double threshold = 1.0;
	private static boolean norm_transcripts = false;
	private static double percentile = -1;
	private static String output_folder;
	private static String organism_database;
	private static boolean compress_output = false;
	private static boolean report_reference = false;
	private static boolean remove_decay_transcripts = true;
	private static boolean report_gene_abundance = false;
	
	// stuff that needs to be retrieved
	private static boolean load_UCSC = false;
	private static boolean load_HGNC = false;
	private static boolean up2date_DDIs = true;
	private static boolean update_UniProt = false;
	private static boolean STRING_weights = false;
	private static List<String> matching_files_output = new LinkedList<>();
	
	
	/**
	 * Prints the help message
	 */
	public static void printHelp() {
		System.out.println("usage: java -jar PPIXpress.jar ([OPTIONS]) [INPUT-NETWORK] [OUTPUT-FOLDER] [EXPR-INPUT1] ([EXPR-INPUT2] ...)");
		
		System.out.println();
		
		System.out.println("[OPTIONS] (optional) :");
		System.out.println("	-g : only use gene abundances (default: transcript abundance)");
		System.out.println("	-t=[threshold] : only take transcripts/genes with an expression above [threshold] into account (default: 1.0)");
		System.out.println("	-tp=[percentile] : only take transcripts/genes with an expression above the [percentile]-th percentile into account (default: overrides option above)");
		System.out.println("	-n : normalize counts of transcripts by the length of the transcript (only for major_transcript file, applied after threshold)");
		System.out.println("	-x : do not remove proteins from the network if their coding transcript is subject to degradation (default: remove if tagged as 'nonsense-mediated decay' or 'non-stop decay')");
		System.out.println("	-d : also output underlying domain-domain interaction network(s) (default: no)");
		System.out.println("	-m : also output major transcript per protein (default: no)");
		System.out.println("	-mg : -m option and abundance per protein as the sum of all expressed coding transcripts");
		System.out.println("	-c : use gzip-compression on output (default: no)");
		System.out.println("	-w : add weights using STRING, interactions that are not in STRING are discarded");
		System.out.println("	-u : update outdated UniProt accessions");
		System.out.println("	-l : do not retrieve current 3did data, only use domain-domain interaction data from DOMINE & IDDI & iPfam");
		System.out.println("	-r=[release] : try to use data from a certain Ensembl release, uses newest if specific release not found");
		System.out.println("	-e : extend MySQL timeout to 5 min (default: 3 min)");
		System.out.println("	-s=[US, UK, AS, 'specific URL'] : change initial server, note that US and asian mirrors only store the last two releases (default: UK)");
		
		System.out.println();
		
		System.out.println("[INPUT-NETWORK] :");
		System.out.println("	Any protein-protein interaction network in SIF-format: Protein1 Protein2 (weight).");
		System.out.println("	Proteins are assumed to be given as UniProt or HGNC accessions.");
		System.out.println("	Alternatively: use taxon:[organism taxon] to automatically retrieve current mentha or IntAct data for an organism if the taxon is not found in mentha.");
		System.out.println("	As an example, 'taxon:9606' will retrieve the current mentha PPIN for human.");
		
		System.out.println();
		
		System.out.println("[OUTPUT-FOLDER] :");
		System.out.println("	The outcome is written to this folder. If it does not exist, it is created.");
		
		System.out.println();
		
		System.out.println("[INPUT1] ([INPUT2] ...) :");
		System.out.println("	Arbitrary number of samples to process.");		
		System.out.println("	Usable inputs that are automatically inferred are:");
		System.out.println("	- Cufflinks : isoform.fpkm_tracking files");
		System.out.println("	- Cufflinks : gene.fpkm_tracking files");
		System.out.println("	- GENCODE GTF (transcript) files");
		System.out.println("	- GENCODE GTF (gene) files");
		System.out.println("	- Kallisto (transcript) files");
		System.out.println("	- Kallisto (gene) files");
		System.out.println("	- RSEM (transcript) files");
		System.out.println("	- RSEM (gene) files");
		System.out.println("	- TCGA (RNA-seq V2) : rsem.isoforms.normalized_results files");
		System.out.println("	- TCGA (RNA-seq V2) : rsem.genes.normalized_results files");
		System.out.println("	- simple textfiles with transcript and expression per line");
		System.out.println("	- simple textfile with gene and expression per line");
		System.out.println("	Transcripts/Genes besides TCGA are assumed to be given as Ensembl identifiers.");
		System.out.println("	Typing 'reference' instead will return the network data for the reference PPIN (and full DDIN and transcripts) of the mapping-stage.");
		System.out.println();
		
		System.exit(0);
	}
	
	
	/**
	 * Prints the version of the program
	 */
	public static void printVersion() {
		System.out.println(version_string);
		System.exit(0);
	}
	
	
	/**
	 * Parse arguments
	 * @param args
	 */
	public static void parseInput(String[] args) {
		
		for (String arg:args) {
			
			// help needed?
			if (arg.equals("-h") || arg.equals("-help"))
				printHelp();
			
			// output version
			else if (arg.equals("-version"))
				printVersion();
			
			// gene level only
			else if (arg.equals("-g"))
				gene_level_only = true;
			
			// gene level only
			else if (arg.equals("-x"))
				remove_decay_transcripts = false;
			
			// output domain-domain interaction networks
			else if (arg.equals("-d"))
				output_DDINs = true;
			
			// output domain-domain interaction networks
			else if (arg.equals("-w"))
				STRING_weights = true;
			
			// output major transcript per protein
			else if (arg.equals("-m"))
				output_major_transcripts = true;
			
			// normalize transcripts
			else if (arg.equals("-n"))
				norm_transcripts = true;
			
			// output major transcript per protein with each ones abundance as the sum of all its expressed coding transcripts
			else if (arg.equals("-mg")) {
				output_major_transcripts = true;
				report_gene_abundance = true;
			}
			
			// compress output
			else if (arg.equals("-c"))
				compress_output = true;
			
			// update outdated UniProt accessions
			else if (arg.equals("-u"))
				update_UniProt = true;
			
			// extends MySQL timeout
			else if (arg.equals("-e"))
				DataQuery.extendTimeout(5 * 60);
			
			// turn off retrieval of iPfam domain-domain interaction data
			else if (arg.equals("-l")) {
				up2date_DDIs = false;
				DataQuery.localDDIsOnly();
			}
			
			// set manual threshold
			else if (arg.startsWith("-t="))
				threshold = Double.parseDouble(arg.split("=")[1]);
			
			// set percentile threshold
			else if (arg.startsWith("-tp=")) 
				percentile = Double.parseDouble(arg.split("=")[1]);
				
			// try to enforce a release manually
			else if (arg.startsWith("-r="))
				DataQuery.enforceSpecificEnsemblRelease(arg.split("=")[1]);
			
			// server switching
			else if (arg.startsWith("-s=")){
				String option = arg.split("=")[1];
				if (option.equals("UK"))
					DataQuery.switchServer("ensembldb.ensembl.org:3306");
				else if (option.equals("US"))
					DataQuery.switchServer("useastdb.ensembl.org:3306");
				else if (option.equals("AS"))
					DataQuery.switchServer("asiadb.ensembl.org:3306");
				else
					DataQuery.switchServer(option);
			}
				
			else {
				
				// must be input network then
				if (original_network_path == null)
					original_network_path = arg;
				
				else if (output_folder == null) {
					// must be output-folder
					output_folder = arg;
				}
				// must be input files
				else {
					input_files.add(arg);
				}
			}
		}
		
		// final checks
		
		if (output_folder == null) {
			System.err.println("Please add an output folder or at least one expression file.");
			System.exit(1);
		}
		
		if (input_files.isEmpty()) {
			System.err.println("No transcript/gene abundance files given.");
			System.exit(1);
		}
			
		// check for occuring classes
		for (String path:input_files) {
			
			if (path.equals("reference")) {
				report_reference = true;
				continue;
			}
			
			String type = TranscriptAbundanceReader.inferTranscriptAbundanceFileType(path);
			if (type.equals("TCGA G"))
				load_HGNC = true;
			else if (type.equals("TCGA T"))
				load_UCSC = true;
			else if (type.equals("unknown")) {
				System.err.println("The format of " + path + " is not supported.");
				System.exit(1);
			}
			
			if (type.endsWith("P")) {
				System.err.println("The format of " + path + " is not yet supported.");
				System.exit(1);
			}
		}
		
		// ensure no later processing of this special case
		if (report_reference)
			input_files.remove("reference");
		
		
	}
	
	
	public static void main(String[] args) {
		
		if (args.length == 1 && args[0].equals("-version"))
			printVersion();
		
		if (args.length < 3) {
			printHelp();
		}
		
		// parse cmd-line and set all parameters
		try {
			parseInput(args);
		} catch (Exception e) {
			System.out.println("Something went wrong while reading the command-line, please check your input!");
			System.out.println();
			printHelp();
		}
		
		// check if output-folder exists, otherwise create it
		if (!output_folder.endsWith("/"))
			output_folder += "/";
		File f = new File(output_folder);
		if (!f.exists())
			f.mkdir();
		
		// gathering even more data if necessary
		if (load_UCSC) {
			DataQuery.switchServerGRCh37();
			System.out.println("Retrieving UCSC mapping-data ...");
			DataQuery.getUCSChg19toTranscriptMap();
		}
		
		PPIN original_network = null;
		if (original_network_path.startsWith("taxon:")) {
			String taxon_id = original_network_path.split(":")[1];
			System.out.print("Retrieving mentha interaction data for taxon " + taxon_id + " ... ");
			original_network = DataQuery.getMenthaNetwork(taxon_id, System.out);
			
			if (original_network.getSizes()[0] == 0) {
				System.out.println("no interaction data for taxon " + taxon_id + " available in mentha.");
				System.out.print("Retrieving IntAct interaction data instead ... ");
				original_network = DataQuery.getIntActNetwork(taxon_id, System.out);
			}
			System.out.println("done.");
			
		} else {
			System.out.println("Reading " + original_network_path + " (may take some time if ID conversion is necessary) ... ");
			original_network = new PPIN(original_network_path);
		}
		
		
		System.out.println("Complete network: " + original_network.getSizesStr());
		
		if (update_UniProt) {
			System.out.println("Retrieving data from UniProt to update protein accessions in network ...");
			original_network = original_network.updateUniprotAccessions();
			System.out.println("Size after update of UniProt accessions: " + original_network.getSizesStr());
		}
		
		if (STRING_weights) {
			if (!original_network.isWeighted())
				System.out.println("Retrieving data from STRING to add weights to network ...");
			else
				System.out.println("Retrieving data from STRING to re-weight network ...");
			original_network = original_network.getAsSTRINGWeighted(false);
			System.out.println("Size after adding weights: " + original_network.getSizesStr());
		}
		
		// gathering data that will always be needed
		organism_database = DataQuery.getEnsemblOrganismDatabaseFromProteins(original_network.getProteins());
		String ensembl_version = organism_database.split("_")[organism_database.split("_").length-2];
		System.out.flush();
		
		System.out.print("Retrieving ENSEMBL " + ensembl_version + " data from database " + organism_database + " (may take some minutes) ... ");
		DataQuery.getGenesTranscriptsProteins(organism_database);
		System.out.print("50% ... ");
		System.out.flush();
		DataQuery.getIsoformProteinDomainMap(organism_database);
		System.out.println("100%.");
		System.out.flush();
		
		if (load_HGNC) {
			System.out.println("Retrieving HGNC mapping-data ...");
			DataQuery.getHGNCProteinsGenes();
		}
		
		if (up2date_DDIs) {
			System.out.println("Retrieving current interaction data fom 3did (" + DataQuery.get3didVersion() + ") ...");
			DataQuery.getKnownDDIs();
		}
		
		// start preprocessing
		System.out.print("Initializing PPIXpress with original network ... ");
		System.out.flush();
		NetworkBuilder builder = new NetworkBuilder(original_network);
		System.out.println(Math.round(builder.getMappingDomainPercentage() * 10000)/100.0 +"% of proteins could be annotated with at least one non-artificial domain," );
		System.out.println(Math.round(builder.getMappingPercentage() * 10000)/100.0 +"% of protein interactions could be associated with at least one non-artificial domain interaction." );
		System.out.flush();
		
		/*
		 * optionally output the reference network
		 */
		
		if (report_reference) {
			ConstructedNetworks constr = NetworkBuilder.constructAssociatedIsoformNetworks(original_network);
			System.out.print("Building output data for reference network ");
			
			String file_suffix = "_ppin.txt";
			if (compress_output)
				file_suffix += ".gz";
			constr.getPPIN().writePPIN(output_folder + "reference" + file_suffix);
			
			System.out.println("-> " + constr.getPPIN().getSizesStr());
			System.out.flush();
			
			if (output_DDINs) {
				file_suffix = "_ddin.txt";
				if (compress_output)
					file_suffix += ".gz";
				constr.getDDIN().writeDDIN(output_folder + "reference" + file_suffix);
			}
			
			if (output_major_transcripts) {
				file_suffix = "_major-transcripts.txt";
				if (compress_output)
					file_suffix += ".gz";
				constr.writeProteinToAssumedTranscriptMap(output_folder + "reference" + file_suffix);
			}
		}
		
		/*
		 * process samples
		 */
		
		if (!remove_decay_transcripts)
			System.out.println("Proteins coded by transcripts that are subject to degradation are not removed.");
		
		int sample_no = 1;
		for (String path:input_files) {
			String match_files = path;
			String type = TranscriptAbundanceReader.inferTranscriptAbundanceFileType(path);
			
			if (type.equals("other")) {
				System.out.println("No valid expression format, " + path + " is skipped.");
				continue;
			}
			
			System.out.print("Processing "+ sample_no + ": " + path + " ("+type+") ");
			
			// check threshold / percentile
			if (percentile > 0)
				threshold = TranscriptAbundanceReader.getExprThresholdFromPercentile(path, percentile, gene_level_only, organism_database);
			
			Map<String, Float> abundance = TranscriptAbundanceReader.readSample(path, threshold, gene_level_only, organism_database);
			
			if (norm_transcripts) {
				if (gene_level_only || !type.endsWith("T"))
					System.out.println("Data normalization skipped, currently only applicable to transcript expression data.");
				else
					abundance = TranscriptAbundanceReader.normalizeByTranscriptLength(abundance, organism_database);
			}
			
			// build
			ConstructedNetworks constr;
			if (gene_level_only || !type.endsWith("T")) {
				constr = builder.constructAssociatedNetworksFromGeneAbundance(abundance.keySet(), remove_decay_transcripts);
			} else {
				constr = builder.constructAssociatedNetworksFromTranscriptAbundance(abundance, remove_decay_transcripts, report_gene_abundance);
			}
			
			/*
			 * write output 
			 */
			
			String file_suffix = "_ppin.txt";
			if (compress_output)
				file_suffix += ".gz";
			constr.getPPIN().writePPIN(output_folder + sample_no + file_suffix);
			match_files += " " + sample_no + file_suffix;
			
			String out = "-> " + constr.getPPIN().getSizesStr() + " (threshold: " + threshold;
			
			if (percentile > 0)
				out += ", " + percentile + "-th percentile";
			
			out += ")";
			System.out.println(out);
			System.out.flush();
			
			if (output_DDINs) {
				file_suffix = "_ddin.txt";
				if (compress_output)
					file_suffix += ".gz";
				constr.getDDIN().writeDDIN(output_folder + sample_no + file_suffix);
				match_files += " " + sample_no + file_suffix;
			}
			
			if (output_major_transcripts) {
				file_suffix = "_major-transcripts.txt";
				if (compress_output)
					file_suffix += ".gz";
				constr.writeProteinToAssumedTranscriptMap(output_folder + sample_no + file_suffix);
				match_files += " " + sample_no + file_suffix;
			}
			
			matching_files_output.add(match_files);
			sample_no++;
		}
		
		// only write if necessary
		if (matching_files_output.size() > 0)
			Utilities.writeEntries(matching_files_output, output_folder + "matching_files.txt"); 
	}
}
