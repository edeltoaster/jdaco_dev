package standalone_tools;

import java.awt.EventQueue;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JButton;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JTextArea;
import javax.swing.JSeparator;
import javax.swing.JCheckBox;
import javax.swing.JTextField;
import javax.swing.JLabel;
import javax.swing.SwingConstants;

import framework.ConstructedNetworks;
import framework.DataQuery;
import framework.NetworkBuilder;
import framework.PPIN;
import framework.TranscriptAbundanceReader;
import framework.Utilities;

import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.JProgressBar;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;

import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;

public class PPIXpress_GUI {
	private static boolean gene_level_only = false;
	private static boolean output_DDINs = false;
	private static boolean STRING_weights = false;
	private static String original_network_path;
	private static PPIN original_ppin;
	private static List<String> input_files = new LinkedList<String>();
	private static double threshold = 1.0;
	private static double percentile = -1;
	private static boolean up2date_DDIs = true;
	private static String output_folder;
	private static String organism_database;
	
	// stuff that needs to be retrieved
	private static boolean load_UCSC = false;
	private static boolean load_HGNC = false;
	private static List<String> matching_files_output = new LinkedList<String>();
	
	// GUI stuff
	private JFrame frmPpixpress;
	private JTextArea text_output;
	private JButton btnBuildNetworks;
	private static PrintStream stream_output;
	private static JTextArea text_ppin;
	private JProgressBar progressBar;
	private static Thread compute_thread;
	private static boolean computing = false;
	private List<JComponent> activiy_changing_components = new LinkedList<JComponent>();
	
	// added for compatibility
	private JButton btnLoadNetwork;
	private JTextArea text_expr;
	private JButton btnLoadExpressionData;
	private JCheckBox chckbxSTRING;
	private JCheckBox chckbxOutputDdins;
	private JCheckBox chckbxGenelevelOnly;
	private JTextField text_threshold;
	private JCheckBox chckbxPercentile;
	private JButton btnReset;
	private JComboBox<String> comboBox_server;
	private JButton btnIntAct;
	
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					PPIXpress_GUI window = new PPIXpress_GUI();
					window.frmPpixpress.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public PPIXpress_GUI() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frmPpixpress = new JFrame();
		frmPpixpress.setTitle("PPIXpress");
		frmPpixpress.setResizable(false);
		frmPpixpress.setBounds(100, 100, 900, 750);
		frmPpixpress.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmPpixpress.getContentPane().setLayout(null);
		
		JScrollPane scrollPane_ppin = new JScrollPane();
		scrollPane_ppin.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
		scrollPane_ppin.setBounds(330, 40, 550, 50);
		frmPpixpress.getContentPane().add(scrollPane_ppin);
		
		text_ppin = new JTextArea();
		scrollPane_ppin.setViewportView(text_ppin);
		text_ppin.setEditable(false);
		
		btnLoadNetwork = new JButton("from file");
		activiy_changing_components.add(btnLoadNetwork);
		btnLoadNetwork.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
				fileChooser.setMultiSelectionEnabled(false);
				fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
				fileChooser.setDialogTitle("Open reference protein-protein interaction network");
				int result = fileChooser.showOpenDialog(frmPpixpress);
				if (result == JFileChooser.APPROVE_OPTION) {
				    File selectedFile = fileChooser.getSelectedFile();
				    original_network_path = selectedFile.getAbsolutePath();
				    activiy_changing_components.add(btnBuildNetworks);
				    compute_thread = new Thread() {
				    	public void run() {
				    		try {
				    			load_network();
				    			computing = false;
				    			progressBar.setValue(0);
				    			setActivity(true);
				    			progressBar.setIndeterminate(false);
				    			activiy_changing_components.remove(btnBuildNetworks);
				    			
				    		} catch (Exception ex) {
				    			compute_thread.interrupt();
				    			computing = false;
				    			try {
									Thread.sleep(1000);
								} catch (InterruptedException exc) {
								}
				    			stream_output.println("");
				    			stream_output.println("Loading network aborted.");
				    			progressBar.setValue(0);
				    			setActivity(true);
				    			progressBar.setIndeterminate(false);
				    			activiy_changing_components.remove(btnBuildNetworks);
				    		}
				    	}
				    };
				    try {
				    	computing = true;
				    	setActivity(false);
				    	progressBar.setIndeterminate(true);
				    	compute_thread.start();
		    		} catch (Exception ex) {
		    			computing = false;
		    			compute_thread.interrupt();
		    			compute_thread = null;
		    			try {
							Thread.sleep(1000);
						} catch (InterruptedException exc) {
						}
		    			stream_output.println("");
		    			stream_output.println("Loading network aborted.");
		    			progressBar.setValue(0);
		    			setActivity(true);
		    			progressBar.setIndeterminate(false);
		    			activiy_changing_components.remove(btnBuildNetworks);
		    		}
				     
				    // "reset"
				    if (original_ppin != null && original_ppin.getProteins().size() > 0)
				    	text_output.setText("");
				    else
				    	original_ppin = null;
				}
			}
		});
		btnLoadNetwork.setBounds(19, 40, 124, 40);
		frmPpixpress.getContentPane().add(btnLoadNetwork);
		
		JSeparator separator = new JSeparator();
		separator.setBounds(5, 120, 890, 20);
		frmPpixpress.getContentPane().add(separator);
		
		JScrollPane scrollPane_expr = new JScrollPane();
		scrollPane_expr.setBounds(330, 160, 550, 120);
		frmPpixpress.getContentPane().add(scrollPane_expr);
		
		text_expr = new JTextArea();
		scrollPane_expr.setViewportView(text_expr);
		text_expr.setEditable(false);
		
		btnLoadExpressionData = new JButton("Load expression data");
		activiy_changing_components.add(btnLoadExpressionData);
		btnLoadExpressionData.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
				fileChooser.setMultiSelectionEnabled(true);
				fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
				fileChooser.setDialogTitle("Load one or many compatible expression files");
				int result = fileChooser.showOpenDialog(frmPpixpress);
				if (result == JFileChooser.APPROVE_OPTION) {
				    File[] selectedFiles = fileChooser.getSelectedFiles();
				    
				    int sample_no = input_files.size() + 1;
				    for (File file:selectedFiles) {
				    	
				    	String path = file.getAbsolutePath();
				    	stream_output.println("Reading " + path + " ... ");
						String type = TranscriptAbundanceReader.inferTranscriptAbundanceFileType(path);
						if (type.equals("TCGA G"))
							load_HGNC = true;
						else if (type.equals("TCGA T"))
							load_UCSC = true;
						else if (type.equals("unknown")) {
							System.err.println("The format of " + path + " is not supported.");
							continue;
						}
						
						if (type.endsWith("P")) {
							System.err.println("The format of " + path + " is not yet supported.");
							continue;
						}
						
						input_files.add(path);
						if (text_expr.getText().equals(""))
							text_expr.setText("Sample " + sample_no + ": " + path + " ("+type+")");
						else
							text_expr.setText(text_expr.getText() + System.getProperty("line.separator") + "Sample " + sample_no + ": " + path + " ("+type+")");
						sample_no++;
					}
				}
			}
		});
		btnLoadExpressionData.setBounds(55, 145, 200, 40);
		frmPpixpress.getContentPane().add(btnLoadExpressionData);
		
		JSeparator separator_1 = new JSeparator();
		separator_1.setBounds(5, 304, 890, 29);
		frmPpixpress.getContentPane().add(separator_1);
		
		btnBuildNetworks = new JButton("Set outputfolder and start processing");
		btnBuildNetworks.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (compute_thread == null || !compute_thread.isAlive() || !computing || btnLoadNetwork.getText().equals("Set outputfolder and start processing")) {
					JFileChooser fileChooser = new JFileChooser();
					fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
					fileChooser.setMultiSelectionEnabled(false);
					fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
					fileChooser.setDialogType(JFileChooser.SAVE_DIALOG);
					fileChooser.setDialogTitle("Set output folder");
					int result = fileChooser.showSaveDialog(frmPpixpress);
					if (result == JFileChooser.APPROVE_OPTION) {
					    File selectedFile = fileChooser.getSelectedFile();
					    output_folder = selectedFile.getAbsolutePath();
					    if (!output_folder.endsWith("/"))
							output_folder += "/";
					    compute_thread = new Thread() {
					    	public void run() {
					    		try {
					    			process_calc();
					    			computing = false;
					    			btnBuildNetworks.setText("Set outputfolder and start processing");
					    			progressBar.setValue(0);
					    			setActivity(true);
					    		} catch (Exception ex) {
					    			compute_thread.interrupt();
					    			computing = false;
					    			try {
										Thread.sleep(1000);
									} catch (InterruptedException exc) {
									}
					    			stream_output.println("");
					    			stream_output.println("Processing aborted.");
					    			progressBar.setValue(0);
					    			setActivity(true);
					    		}
					    	}
					    };
					    try {
					    	computing = true;
					    	compute_thread.start();
					    	btnBuildNetworks.setText("Abort processing");
			    		} catch (Exception ex) {
			    			computing = false;
			    			compute_thread.interrupt();
			    			compute_thread = null;
			    			try {
								Thread.sleep(1000);
							} catch (InterruptedException exc) {
							}
			    			stream_output.println("");
			    			stream_output.println("Processing aborted.");
			    			progressBar.setValue(0);
			    			setActivity(true);
			    		}
					}
				} else { // abort function
					computing = false;
					compute_thread.interrupt();
					compute_thread = null;
					try {
						Thread.sleep(1000);
					} catch (InterruptedException exc) {
					}
					stream_output.println("");
					stream_output.println("Processing aborted.");
					progressBar.setValue(0);
					btnBuildNetworks.setText("Set outputfolder and start processing");
					setActivity(true);
				}
			}
		});
		btnBuildNetworks.setBounds(137, 342, 324, 50);
		frmPpixpress.getContentPane().add(btnBuildNetworks);
		
		chckbxSTRING = new JCheckBox("add STRING weights");
		activiy_changing_components.add(chckbxSTRING);
		chckbxSTRING.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxSTRING.isSelected()) 
					STRING_weights = true;
				else
					STRING_weights = false;
			}
		});
		chckbxSTRING.setBounds(74, 90, 193, 30);
		frmPpixpress.getContentPane().add(chckbxSTRING);
		
		chckbxOutputDdins = new JCheckBox("output DDINs");
		activiy_changing_components.add(chckbxOutputDdins);
		chckbxOutputDdins.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxOutputDdins.isSelected())
					output_DDINs = true;
				else
					output_DDINs = false;
			}
		});
		chckbxOutputDdins.setBounds(680, 370, 132, 30);
		frmPpixpress.getContentPane().add(chckbxOutputDdins);
		
		chckbxGenelevelOnly = new JCheckBox("gene-level only");
		activiy_changing_components.add(chckbxGenelevelOnly);
		chckbxGenelevelOnly.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxGenelevelOnly.isSelected())
					gene_level_only = true;
				else
					gene_level_only = false;
			}
		});
		chckbxGenelevelOnly.setBounds(75, 195, 153, 30);
		frmPpixpress.getContentPane().add(chckbxGenelevelOnly);
		
		text_threshold = new JTextField();
		text_threshold.setBounds(160, 231, 50, 30);
		frmPpixpress.getContentPane().add(text_threshold);
		activiy_changing_components.add(text_threshold);
		text_threshold.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					threshold = Double.parseDouble(text_threshold.getText());
				} catch (Exception ex) {
					stream_output.println("Expression threshold no valid floating point value.");
				}
			}
		});
		text_threshold.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				try {
					threshold = Double.parseDouble(text_threshold.getText());
				} catch (Exception ex) {
					stream_output.println("Expression threshold no valid floating point value.");
				}
			}
		});
		text_threshold.setHorizontalAlignment(SwingConstants.RIGHT);
		text_threshold.setText("1.0");
		text_threshold.setColumns(10);
		
		chckbxPercentile = new JCheckBox("percentile-based");
		
		chckbxPercentile.setBounds(75, 259, 153, 30);
		frmPpixpress.getContentPane().add(chckbxPercentile);
		activiy_changing_components.add(chckbxPercentile);
		chckbxPercentile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxPercentile.isSelected()) {
					text_threshold.setText("50.0");
					percentile = 50.0;
				} else {
					percentile = -1;
					text_threshold.setText("1.0");
					threshold = 1.0;
				}
			}
		});
		
		JScrollPane scrollPane_output = new JScrollPane();
		scrollPane_output.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		scrollPane_output.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		scrollPane_output.setBounds(20, 450, 860, 260);
		frmPpixpress.getContentPane().add(scrollPane_output);
		
		text_output = new JTextArea();
		scrollPane_output.setViewportView(text_output);
		text_output.setTabSize(4);
		text_output.setEditable(false);
		stream_output = new PrintStream(new CustomOutputStream(text_output));
		
		btnReset = new JButton("Reset");
		activiy_changing_components.add(btnReset);
		btnReset.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				text_output.setText("");
				
				gene_level_only = false;
				chckbxGenelevelOnly.setSelected(false);
				
				output_DDINs = false;
				chckbxOutputDdins.setSelected(false);
				
				STRING_weights = false;
				chckbxSTRING.setSelected(false);
				
				original_network_path = null;
				text_ppin.setText("");
				original_ppin = null;
				
				input_files.clear();
				text_expr.setText("");
				
				threshold = 1.0;
				text_threshold.setText("1.0");
				percentile = -1;
				chckbxPercentile.setSelected(false);
				
				up2date_DDIs = true;
				output_folder = null;
				organism_database = null;
				load_UCSC = false;
				load_HGNC = false;
				setActivity(true);
				computing = false;
				progressBar.setValue(0);
			}
		});
		btnReset.setBounds(483, 342, 117, 50);
		frmPpixpress.getContentPane().add(btnReset);
		
		progressBar = new JProgressBar();
		progressBar.setBounds(100, 418, 700, 20);
		frmPpixpress.getContentPane().add(progressBar);
		
		comboBox_server = new JComboBox<String>();
		activiy_changing_components.add(comboBox_server);
		comboBox_server.setModel(new DefaultComboBoxModel<String>(new String[] {"US", "UK", "ASIA"}));
		comboBox_server.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				if (comboBox_server.getSelectedItem().toString().equals("UK"))
					DataQuery.switchServer("ensembldb.ensembl.org:3306");
				else if (comboBox_server.getSelectedItem().toString().equals("US"))
					DataQuery.switchServer("useastdb.ensembl.org:3306");
				else if (comboBox_server.getSelectedItem().toString().equals("ASIA"))
					DataQuery.switchServer("asiadb.ensembl.org:3306");
			}
		});
		comboBox_server.setSelectedIndex(0);
		comboBox_server.setBounds(720, 338, 82, 30);
		frmPpixpress.getContentPane().add(comboBox_server);
		
		JLabel lblSever = new JLabel("Server:");
		activiy_changing_components.add(lblSever);
		lblSever.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSever.setBounds(650, 337, 61, 30);
		frmPpixpress.getContentPane().add(lblSever);
		
		JLabel lblProgress = new JLabel("Progress:");
		lblProgress.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProgress.setBounds(17, 412, 75, 30);
		frmPpixpress.getContentPane().add(lblProgress);
		
		JLabel lblNewLabel = new JLabel("threshold:");
		lblNewLabel.setBounds(70, 231, 85, 30);
		frmPpixpress.getContentPane().add(lblNewLabel);
		activiy_changing_components.add(lblNewLabel);
		lblNewLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		
		JLabel lblReferenceProteinproteinInteaction = new JLabel("reference protein-protein interaction network:");
		lblReferenceProteinproteinInteaction.setHorizontalAlignment(SwingConstants.LEFT);
		lblReferenceProteinproteinInteaction.setBounds(330, 10, 414, 30);
		frmPpixpress.getContentPane().add(lblReferenceProteinproteinInteaction);
		
		JLabel lblSpecificExpressionData = new JLabel("specific expression data:");
		lblSpecificExpressionData.setHorizontalAlignment(SwingConstants.LEFT);
		lblSpecificExpressionData.setBounds(330, 130, 256, 30);
		frmPpixpress.getContentPane().add(lblSpecificExpressionData);
		
		btnIntAct = new JButton("from IntAct");
		btnIntAct.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				final String taxon_id = (String)JOptionPane.showInputDialog(frmPpixpress,"organism taxon:", "IntAct retrieval", JOptionPane.PLAIN_MESSAGE, null, null, "9606");
				
				if (taxon_id == null)
					return;
				
				activiy_changing_components.add(btnBuildNetworks);
			    compute_thread = new Thread() {
			    	public void run() {
			    		try {
			    			load_intact(taxon_id);;
			    			computing = false;
			    			progressBar.setValue(0);
			    			setActivity(true);
			    			progressBar.setIndeterminate(false);
			    			activiy_changing_components.remove(btnBuildNetworks);
			    			
			    		} catch (Exception ex) {
			    			compute_thread.interrupt();
			    			computing = false;
			    			try {
								Thread.sleep(1000);
							} catch (InterruptedException exc) {
							}
			    			stream_output.println("");
			    			stream_output.println("IntAct retrieval aborted.");
			    			progressBar.setValue(0);
			    			setActivity(true);
			    			progressBar.setIndeterminate(false);
			    			activiy_changing_components.remove(btnBuildNetworks);
			    		}
			    	}
			    };
			    try {
			    	computing = true;
			    	setActivity(false);
			    	progressBar.setIndeterminate(true);
			    	compute_thread.start();
	    		} catch (Exception ex) {
	    			computing = false;
	    			compute_thread.interrupt();
	    			compute_thread = null;
	    			try {
						Thread.sleep(1000);
					} catch (InterruptedException exc) {
					}
	    			stream_output.println("");
	    			stream_output.println("IntAct retrieval aborted.");
	    			progressBar.setValue(0);
	    			setActivity(true);
	    			progressBar.setIndeterminate(false);
	    			activiy_changing_components.remove(btnBuildNetworks);
	    		}
			    
			    // "reset"
			    if (original_ppin != null && original_ppin.getProteins().size() > 0)
			    	text_output.setText("");
			    else
			    	original_ppin = null;
			}
		});
		btnIntAct.setBounds(166, 40, 124, 40);
		frmPpixpress.getContentPane().add(btnIntAct);
		activiy_changing_components.add(btnIntAct);
		
		JLabel lblLoadProteinInteraction = new JLabel("Load protein interaction data");
		lblLoadProteinInteraction.setBounds(55, 5, 200, 40);
		frmPpixpress.getContentPane().add(lblLoadProteinInteraction);

		// most general
		System.setErr(stream_output);
		System.setOut(stream_output);
	}
	
	// start calculations
	public void process_calc() {
		// check for input data
		if (original_ppin == null) {
			stream_output.println("No protein interaction network given.");
			return;
		}
		
		if (input_files.size() == 0) {
			stream_output.println("No expression data given.");
			return;
		}
		
		if ( !(new File(output_folder).exists())) {
			stream_output.println("Folder " + output_folder + " does not exist.");
			stream_output.println("On OSX this is most likely a localization artefact. In that case please select another folder or subfolder.");
			return;
		}
		
		// processing
		setActivity(false);
		if (STRING_weights) {
			if (!original_ppin.isWeighted())
				stream_output.println("Retrieving data from STRING to add weights to network ...");
			else
				stream_output.println("Retrieving data from STRING to re-weight network ...");
			
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			
			original_ppin = original_ppin.getAsSTRINGWeighted(false);
			
			stream_output.println("New size: " + original_ppin.getSizesStr());
			text_ppin.setText("Complete network (+STRING): " + original_network_path + System.getProperty("line.separator") + "Size: " + original_ppin.getSizesStr());
		}
		
		
		progressBar.setValue(5);
		
		// gathering data that will always be needed
		organism_database = DataQuery.getEnsemblOrganismDatabaseFromProteins(original_ppin.getProteins());
		String ensembl_version = organism_database.split("_")[organism_database.split("_").length-2];
		
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		stream_output.print("Retrieving ENSEMBL " + ensembl_version + " data from database " + organism_database + " (may take some minutes) ... ");
		
		DataQuery.getGenesTranscriptsProteins(organism_database);
		
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		stream_output.print("33% ... ");
		
		progressBar.setValue(30);
		
		DataQuery.getIsoformProteinDomainMap(organism_database);
		
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		stream_output.print("66% ... ");
		
		progressBar.setValue(40);
		
		DataQuery.getTranscriptsDomains(organism_database);
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		stream_output.println("100%");
		
		progressBar.setValue(45);
		
		// gathering even more data if necessary
		if (load_UCSC) {
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			stream_output.println("Retrieving UCSC mapping-data ...");
			DataQuery.getUSCStoTranscriptMap(DataQuery.getEnsemblOrganismDatabaseFromName("homo sapiens"));
		}
		
		if (load_HGNC) {
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			stream_output.println("Retrieving HGNC mapping-data ...");
			DataQuery.getHGNCProteinsGenes();
		}
		
		if (up2date_DDIs) {
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			stream_output.println("Retrieving current interaction data fom 3did (" + DataQuery.get3didVersion() + ") and iPfam (" + DataQuery.getIPfamVersion() + ") ...");
			DataQuery.getKnownDDIs();
		}
		
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		progressBar.setValue(50);
		
		// start preprocessing
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		stream_output.println("Initializing PPIXpress with original network ... ");
		NetworkBuilder builder = new NetworkBuilder(original_ppin);
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		stream_output.println(Math.round(builder.getMappingDomainPercentage() * 10000)/100.0 +"% of proteins could be annotated with at least one non-artificial domain," );
		stream_output.println(Math.round(builder.getMappingPercentage() * 10000)/100.0 +"% of protein interactions could be associated with at least one non-artificial domain interaction." );
		
		int step_size = 50 / input_files.size();
		
		// process samples
		int sample_no = 1;
		matching_files_output.clear();
		for (String path:input_files) {
			String match_files = path;
			String type = TranscriptAbundanceReader.inferTranscriptAbundanceFileType(path);
			
			if (type.equals("other")) {
				if (!computing) {
					progressBar.setValue(0);
					return;
				}
				stream_output.println("No valid expression format, " + path + " is skipped.");
				continue;
			}
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			stream_output.println("");
			stream_output.println("Processing "+ sample_no + ": " + path + " ("+type+") ");
			
			// check threshold / percentile
			if (percentile > 0)
				threshold = TranscriptAbundanceReader.getExprThresholdFromPercentile(path, percentile, gene_level_only, organism_database);
			
			Map<String, Float> abundance = TranscriptAbundanceReader.readSample(path, threshold, gene_level_only, organism_database);
			
			// build
			ConstructedNetworks constr;
			if (gene_level_only || !type.endsWith("T")) {
				constr = builder.constructAssociatedNetworksFromGeneAbundance(abundance.keySet());
			} else {
				constr = builder.constructAssociatedNetworksFromTranscriptAbundance(abundance);
			}
			
			// output
			constr.getPPIN().writePPIN(output_folder + sample_no + "_ppin.tsv");
			match_files += " " + sample_no + "_ppin.tsv";
			
			String out = "-> " + constr.getPPIN().getSizesStr() + " (threshold: " + threshold;
			
			if (percentile > 0)
				out += ", " + percentile + "-th percentile";
			
			out += ")";
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			stream_output.println(out);
			
			if (output_DDINs) {
				constr.getDDIN().writeDDIN(output_folder + sample_no + "_ddin.tsv");
				match_files += " " + sample_no + "_ddin.tsv";
			}
			
			matching_files_output.add(match_files);
			sample_no++;
			progressBar.setValue(Math.min(progressBar.getValue() + step_size, 100));
		}
		
		Utilities.writeEntries(matching_files_output, output_folder + "matching_files.txt");
		
		stream_output.println("Finished.");
	}
	
	private void setActivity(boolean active) {
		for (JComponent c:activiy_changing_components)
			c.setEnabled(active);
	}
	
	public void load_intact(String taxon_id) {
		setActivity(false);
		stream_output.print("Retrieving IntAct interaction network for taxon " + taxon_id + " ... ");
		original_ppin = DataQuery.getIntActNetwork(taxon_id, stream_output);
		original_network_path = "IntAct, taxon:" + taxon_id;
		stream_output.println("done.");
	    text_ppin.setText("Complete network: " + original_network_path + System.getProperty("line.separator") + "Size: " + original_ppin.getSizesStr());
	}
	
	public void load_network() {
		setActivity(false);
		stream_output.println("Reading " + original_network_path + " (may take some time if ID conversion is necessary) ... ");
	    original_ppin = new PPIN(original_network_path);
	    text_ppin.setText("Complete network: " + original_network_path + System.getProperty("line.separator") + "Size: " + original_ppin.getSizesStr());
	}
	
	/**
	 * GUI helper stuff
	 */
	
	// wrapper of stream to JTextArea,  www.codejava.net
	public class CustomOutputStream extends OutputStream {
	    private JTextArea textArea;
	     
	    public CustomOutputStream(JTextArea textArea) {
	        this.textArea = textArea;
	    }

	    public void write(final int b) throws IOException {
	    	EventQueue.invokeLater(new Runnable() {
				public void run() {
					synchronized (textArea) {
			    		textArea.append(String.valueOf((char)b));
				        textArea.setCaretPosition(textArea.getDocument().getLength());
					}
				}
			});
	    	
	    }
	}
}
