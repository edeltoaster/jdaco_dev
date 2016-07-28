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
	private static boolean output_major_transcripts = false;
	private static String original_network_path;
	private static PPIN original_ppin;
	private static List<String> input_files = new LinkedList<>();
	private static double threshold = 1.0;
	private static double percentile = -1;
	private static String output_folder;
	private static String organism_database;
	private static boolean compress_output = false;
	private static boolean report_reference = false;
	private static boolean remove_decay_transcripts = true;
	
	// stuff that needs to be retrieved
	private static boolean load_UCSC = false;
	private static boolean load_HGNC = false;
	private static boolean STRING_weights = false;
	private static boolean up2date_DDIs = true;
	private static boolean update_UniProt = false;
	private static List<String> matching_files_output = new LinkedList<>();
	
	// GUI stuff
	private JFrame frmPpixpress;
	private JTextArea text_output;
	private JButton btnBuildNetworks;
	private static PrintStream stream_output;
	private static JTextArea text_ppin;
	private JProgressBar progressBar;
	private static Thread compute_thread;
	private static boolean computing = false;
	private List<JComponent> activiy_changing_components = new LinkedList<>();
	
	// added for compatibility
	private JButton btnLoadNetwork;
	private JTextArea text_expr;
	private JButton btnLoadExpressionData;
	private JCheckBox chckbxSTRING;
	private JCheckBox chckbxOutputDdins;
	private JCheckBox chckbxOnlyLocalDdi;
	private JCheckBox chckbxGenelevelOnly;
	private JTextField text_threshold;
	private JCheckBox chckbxPercentile;
	private JButton btnReset;
	private JComboBox<String> comboBox_server;
	private JButton btnRetrieve;
	private JCheckBox chckbxReference;
	private JCheckBox chckbxOutputMajorTranscripts;
	private JCheckBox chckbxUpdateUniprotAccs;
	private JCheckBox chckboxCompressOutput;
	
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
		frmPpixpress.setBounds(100, 100, 900, 850);
		frmPpixpress.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmPpixpress.getContentPane().setLayout(null);
		
		JScrollPane scrollPane_ppin = new JScrollPane();
		scrollPane_ppin.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
		scrollPane_ppin.setBounds(330, 40, 550, 52);
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
		separator.setBounds(5, 163, 890, 20);
		frmPpixpress.getContentPane().add(separator);
		
		JScrollPane scrollPane_expr = new JScrollPane();
		scrollPane_expr.setBounds(330, 200, 550, 125);
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
		btnLoadExpressionData.setBounds(55, 180, 200, 40);
		frmPpixpress.getContentPane().add(btnLoadExpressionData);
		
		JSeparator separator_1 = new JSeparator();
		separator_1.setBounds(5, 335, 890, 29);
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
		btnBuildNetworks.setBounds(100, 400, 324, 50);
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
		chckbxOutputDdins.setBounds(630, 409, 132, 30);
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
		chckbxGenelevelOnly.setBounds(75, 222, 153, 30);
		frmPpixpress.getContentPane().add(chckbxGenelevelOnly);
		
		text_threshold = new JTextField();
		text_threshold.setBounds(160, 270, 50, 30);
		frmPpixpress.getContentPane().add(text_threshold);
		activiy_changing_components.add(text_threshold);
		text_threshold.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					if (chckbxPercentile.isSelected())
						percentile = Double.parseDouble(text_threshold.getText());
					else
						threshold = Double.parseDouble(text_threshold.getText());
				} catch (Exception ex) {
					if (chckbxPercentile.isSelected())
						stream_output.println("Percentile no valid floating point value.");
					else
						stream_output.println("Expression threshold no valid floating point value.");
				}
			}
		});
		text_threshold.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				try {
					if (chckbxPercentile.isSelected())
						percentile = Double.parseDouble(text_threshold.getText());
					else
						threshold = Double.parseDouble(text_threshold.getText());
				} catch (Exception ex) {
					if (chckbxPercentile.isSelected())
						stream_output.println("Percentile no valid floating point value.");
					else
						stream_output.println("Expression threshold no valid floating point value.");
				}
			}
		});
		text_threshold.setHorizontalAlignment(SwingConstants.RIGHT);
		text_threshold.setText("1.0");
		text_threshold.setColumns(10);
		
		chckbxPercentile = new JCheckBox("percentile-based");
		
		chckbxPercentile.setBounds(75, 300, 153, 30);
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
		scrollPane_output.setBounds(19, 520, 860, 300);
		frmPpixpress.getContentPane().add(scrollPane_output);
		
		text_output = new JTextArea();
		scrollPane_output.setViewportView(text_output);
		text_output.setTabSize(4);
		text_output.setEditable(false);
		stream_output = new PrintStream(new CustomOutputStream(text_output));
		
		chckbxOutputMajorTranscripts = new JCheckBox("output major transcripts");
		chckbxOutputMajorTranscripts.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxOutputMajorTranscripts.isSelected())
					output_major_transcripts = true;
				else
					output_major_transcripts = false;
			}
		});
		chckbxOutputMajorTranscripts.setBounds(630, 430, 205, 30);
		frmPpixpress.getContentPane().add(chckbxOutputMajorTranscripts);
		activiy_changing_components.add(chckbxOutputMajorTranscripts);
		
		chckbxUpdateUniprotAccs = new JCheckBox("update UniProt accessions");
		chckbxUpdateUniprotAccs.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxUpdateUniprotAccs.isSelected())
					update_UniProt = true;
				else
					update_UniProt = false;
			}
		});
		chckbxUpdateUniprotAccs.setBounds(74, 112, 220, 30);
		frmPpixpress.getContentPane().add(chckbxUpdateUniprotAccs);
		activiy_changing_components.add(chckbxUpdateUniprotAccs);
		
		chckboxCompressOutput = new JCheckBox("compress output");
		chckboxCompressOutput.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckboxCompressOutput.isSelected())
					compress_output = true;
				else
					compress_output = false;
			}
		});
		chckboxCompressOutput.setBounds(630, 451, 205, 30);
		frmPpixpress.getContentPane().add(chckboxCompressOutput);
		activiy_changing_components.add(chckboxCompressOutput);
		
		chckbxOnlyLocalDdi = new JCheckBox("only local DDI data");
		chckbxOnlyLocalDdi.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxOnlyLocalDdi.isSelected()) {
					up2date_DDIs = false;
					DataQuery.localDDIsOnly();
				} else {
					up2date_DDIs = true;
					DataQuery.defaultDDIs();
				}
			}
		});
		chckbxOnlyLocalDdi.setBounds(74, 134, 220, 30);
		frmPpixpress.getContentPane().add(chckbxOnlyLocalDdi);
		activiy_changing_components.add(chckbxOnlyLocalDdi);
		
		chckbxReference = new JCheckBox("output reference network data");
		chckbxReference.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (chckbxReference.isSelected()) {
					report_reference = true;
				} else {
					report_reference = false;
				}
			}
		});
		chckbxReference.setBounds(630, 386, 232, 30);
		frmPpixpress.getContentPane().add(chckbxReference);
		activiy_changing_components.add(chckbxReference);
		
		btnReset = new JButton("Reset");
		activiy_changing_components.add(btnReset);
		btnReset.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				text_output.setText("");
				
				gene_level_only = false;
				chckbxGenelevelOnly.setSelected(false);
				
				output_DDINs = false;
				chckbxOutputDdins.setSelected(false);
				
				output_major_transcripts = false;
				chckbxOutputMajorTranscripts.setSelected(false);
				
				STRING_weights = false;
				chckbxSTRING.setSelected(false);
				
				report_reference = false;
				chckbxReference.setSelected(false);
				
				update_UniProt = false;
				chckbxUpdateUniprotAccs.setSelected(false);
				
				compress_output = false;
				chckboxCompressOutput.setSelected(false);
				
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
				chckbxOnlyLocalDdi.setSelected(false);
				DataQuery.defaultDDIs();
				
				output_folder = null;
				organism_database = null;
				load_UCSC = false;
				load_HGNC = false;
				setActivity(true);
				computing = false;
				progressBar.setValue(0);
				
				System.gc();
			}
		});
		btnReset.setBounds(455, 400, 117, 50);
		frmPpixpress.getContentPane().add(btnReset);
		
		progressBar = new JProgressBar();
		progressBar.setBounds(100, 485, 700, 20);
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
		comboBox_server.setBounds(720, 355, 82, 30);
		frmPpixpress.getContentPane().add(comboBox_server);
		
		JLabel lblSever = new JLabel("Server:");
		activiy_changing_components.add(lblSever);
		lblSever.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSever.setBounds(650, 354, 61, 30);
		frmPpixpress.getContentPane().add(lblSever);
		activiy_changing_components.add(lblSever);
		
		JLabel lblProgress = new JLabel("Progress:");
		lblProgress.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProgress.setBounds(17, 479, 75, 30);
		frmPpixpress.getContentPane().add(lblProgress);
		
		JLabel lblthreshold = new JLabel("threshold:");
		lblthreshold.setBounds(70, 270, 85, 30);
		frmPpixpress.getContentPane().add(lblthreshold);
		activiy_changing_components.add(lblthreshold);
		lblthreshold.setHorizontalAlignment(SwingConstants.RIGHT);
		activiy_changing_components.add(lblthreshold);
		
		JLabel lblReferenceProteinproteinInteaction = new JLabel("reference protein-protein interaction network:");
		lblReferenceProteinproteinInteaction.setHorizontalAlignment(SwingConstants.LEFT);
		lblReferenceProteinproteinInteaction.setBounds(330, 10, 414, 30);
		frmPpixpress.getContentPane().add(lblReferenceProteinproteinInteaction);
		activiy_changing_components.add(lblReferenceProteinproteinInteaction);
		
		JLabel lblSpecificExpressionData = new JLabel("specific expression data:");
		lblSpecificExpressionData.setHorizontalAlignment(SwingConstants.LEFT);
		lblSpecificExpressionData.setBounds(330, 170, 256, 30);
		frmPpixpress.getContentPane().add(lblSpecificExpressionData);
		activiy_changing_components.add(lblSpecificExpressionData);
		
		btnRetrieve = new JButton("from web");
		btnRetrieve.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				final String taxon_id = (String)JOptionPane.showInputDialog(frmPpixpress,"organism taxon:", "iRefIndex/IntAct retrieval", JOptionPane.PLAIN_MESSAGE, null, null, "9606");
				
				if (taxon_id == null)
					return;
				
				activiy_changing_components.add(btnBuildNetworks);
			    compute_thread = new Thread() {
			    	public void run() {
			    		try {
			    			retrieve_network(taxon_id);;
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
			    			stream_output.println("Network retrieval aborted.");
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
	    			stream_output.println("Network retrieval aborted.");
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
		btnRetrieve.setBounds(166, 40, 124, 40);
		frmPpixpress.getContentPane().add(btnRetrieve);
		activiy_changing_components.add(btnRetrieve);
		
		JLabel lblLoadProteinInteraction = new JLabel("Load protein interaction data");
		lblLoadProteinInteraction.setBounds(55, 5, 220, 40);
		frmPpixpress.getContentPane().add(lblLoadProteinInteraction);
		activiy_changing_components.add(lblLoadProteinInteraction);
		
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
		
		if (input_files.size() == 0 && !report_reference) {
			stream_output.println("No networks to construct.");
			return;
		}
		
		if ( !(new File(output_folder).exists())) {
			stream_output.println("Folder " + output_folder + " does not exist.");
			stream_output.println("On OSX this is most likely a localization artefact. In that case please select another folder or subfolder.");
			return;
		}
		
		// processing
		setActivity(false);
		
		if (update_UniProt) {

			stream_output.println("Retrieving data from UniProt to update protein accessions in network ...");
			
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			
			original_ppin = original_ppin.updateUniprotAccessions();
			
			stream_output.println("New size: " + original_ppin.getSizesStr());
			text_ppin.setText("Complete network: " + original_network_path + System.getProperty("line.separator") + "Size: " + original_ppin.getSizesStr());
		}
		
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
		
		stream_output.print("50% ... ");
		progressBar.setValue(30);
		
		DataQuery.getIsoformProteinDomainMap(organism_database);
		
		if (!computing) {
			progressBar.setValue(0);
			return;
		}
		
		stream_output.println("100%.");
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
			stream_output.println("Retrieving current interaction data fom 3did (" + DataQuery.get3didVersion() + ") ...");
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
		
		int add_count = 0;
		if (report_reference)
			add_count = 1;
		int step_size = 50 / (input_files.size() + add_count);
		
		/*
		 * optionally output the reference network
		 */
		
		if (report_reference) {
			ConstructedNetworks constr = NetworkBuilder.constructAssociatedIsoformNetworks(original_ppin);
			stream_output.print("Building output data for reference network ");
			
			String file_suffix = "_ppin.txt";
			if (compress_output)
				file_suffix += ".gz";
			constr.getPPIN().writePPIN(output_folder + "reference" + file_suffix);
			
			stream_output.println("-> " + constr.getPPIN().getSizesStr());
			
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
				constr = builder.constructAssociatedNetworksFromGeneAbundance(abundance.keySet(), remove_decay_transcripts);
			} else {
				constr = builder.constructAssociatedNetworksFromTranscriptAbundance(abundance, remove_decay_transcripts);
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
			if (!computing) {
				progressBar.setValue(0);
				return;
			}
			stream_output.println(out);
			
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
			progressBar.setValue(Math.min(progressBar.getValue() + step_size, 100));
		}
		
		// only write if there is something to write about
		if (matching_files_output.size() > 0)
			Utilities.writeEntries(matching_files_output, output_folder + "matching_files.txt");
		
		stream_output.println("Finished.");
	}
	
	private void setActivity(boolean active) {
		for (JComponent c:activiy_changing_components)
			c.setEnabled(active);
	}
	
	public void retrieve_network(String taxon_id) {
		text_ppin.setText("");
		setActivity(false);
		stream_output.print("Retrieving mentha interaction data for taxon " + taxon_id + " ... ");
		original_ppin = DataQuery.getMenthaNetwork(taxon_id, stream_output);
		original_network_path = "retrieved from mentha, taxon:" + taxon_id;
		
		if (original_ppin.getSizes()[0] == 0) {
			stream_output.println("no interaction data for taxon " + taxon_id + " available in mentha.");
			stream_output.print("Retrieving IntAct interaction data instead ... ");
			original_ppin = DataQuery.getIntActNetwork(taxon_id, stream_output);
			original_network_path = "retrieved from IntAct, taxon:" + taxon_id;
		}
		
		stream_output.println("done.");
		
	    text_ppin.setText("Complete network: " + original_network_path + System.getProperty("line.separator") + "Size: " + original_ppin.getSizesStr());
	}
	
	public void load_network() {
		text_ppin.setText("");
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
