import java.awt.Desktop;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.HPos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Pane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.scene.control.TextField;
import javafx.scene.control.ComboBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import java.awt.image.BufferedImage;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;
import javafx.scene.control.ScrollPane;
import javax.imageio.ImageIO;

/**
 * Mutual Information of Biological Sequences
 * This module is used to take two multiple sequence alignments and determine, position by position, the mutual information between the two. Alternatively, intramolecular mutual information can be found.
 * @author Devin Camenares, PhD
 *
 * @version 20-05-14
 * Latest update: introduced conservation adjustment score (CAS) to gain a sense for where the most conserved or divergent pairs are, and adjust the MI score accordingly. Also introduced scroll bar!
 * 
 * @since 1-19-16
 */

/*
    Mutual Information of Biological Sequences
    Copyright (C) 2016 Devin Camenares

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
 */
 
public final class MutualInfoBio extends Application {
 
    private Desktop desktop = Desktop.getDesktop();
 
	/**
         * This method was written by CC and modified by DJC, based upon code from StackOverflow user Kip.
         * It takes a string and writes it to a specified file, appending said file.
         * The purpose of this method is to save information and results as they are generated, reducing memory load.
         * 
         * @param content The string to be written
         * @param file The destination file for the string output
         */
	private void writeFile(String content, File file){
                try (FileWriter fw = new FileWriter(file, true)) {
                    fw.write(content);//

                } catch (IOException ex) 
                    {
		     Logger.getLogger(MutualInfoBio.class.getName()).log(Level.SEVERE, null, ex);
		    }
		         
            }
        
	/**
         * This method was written by CC and modified by DJC, based upon code from StackOverflow user Kip.
         * It takes a string and writes it to a specified file, appending said file.
         * The purpose of this method is to save information and results as they are generated, reducing memory load.
         * 
         * @param content The string to be written
         * @param file The destination file for the string output
         */
	private void writeFileN(String content, File file){
                try (FileWriter fw = new FileWriter(file, true)) {
                    content += "" + System.lineSeparator();
                    fw.write(content);//

                } catch (IOException ex) 
                    {
		     Logger.getLogger(MutualInfoBio.class.getName()).log(Level.SEVERE, null, ex);
		    }
		         
            }        
	
        	/**
         * This method was written by CC
         * 
         * @param file The desired file to be opened.
         */
        private void openFile(File file) {
            try {
		 desktop.open(file);
		} catch (IOException ex) {
		Logger.getLogger(
                 MutualInfoBio.class.getName()).log(
                 Level.SEVERE, null, ex
                 );
		 }
            }
        
        // Configures the File Chooser, from Oracle Documentation Example 26-5
	private static void configureFileChooser(
            final FileChooser fileChooser) {      
		fileChooser.setTitle("View Files");
		fileChooser.setInitialDirectory(
			new File(System.getProperty("user.dir"))
			);                 
		fileChooser.getExtensionFilters().addAll(
		new FileChooser.ExtensionFilter("Plain Text", "*.txt"),
		new FileChooser.ExtensionFilter("FASTA", "*.fasta"),
		new FileChooser.ExtensionFilter("Rich Text Format", "*.rtf"),
		new FileChooser.ExtensionFilter("All Files", "*.*")
		);
            }

			
	// Initialize All Global Variables, containers
        
        /**
         * Container for first file information / path
         */
        File file1 = new File("");
    
        /**
         * Container for second file information / path
         */    
        File file2 = new File("");        
	
        /**
         * String to hold name of the first file loaded
         */
	public static String flName1 = "";
        
        /**
         * String to hold name of the second file loaded
         */        
	public static String flName2 = "";
        
        /**
         * String to hold path of the first file loaded
         */        
        public static String flPath1 = "";
        /**
         * String to hold path of the second file loaded
         */        
        public static String flPath2 = "";
	
        /**
         * Display text of the grouped protein category, representing selection of String protgr
         */
	public static String protgrCat = "Protein, Grouped";
        
        /**
         * Display text of the individual amino acid category, representing selection of string aa
         */
	public static String aaCat = "Protein";
        
        /**
         * Display text of the individual DNA residue category, representing selection of string dna
         */
	public static String dnaCat = "DNA";
        
        /**
         * Display text of the individual RNA residue category, representing selection of string rna
         */
	public static String rnaCat = "RNA";
        
        /**
         * Display text of the choice for a custom category
         */
	public static String customCat = "Custom Category";
	
        /**
         * Amino acid residue grouping based on chemical properties, separated by comma
         */
	public static String protgr = "MC, VIAL, DE, KR, FWY, TSNQ, H, P, G, -X";
        
        /**
         * Amino acid residue grouping
         */
	public static String aa = "M, C, V, I, A, L, D, E, K, R, F, W, Y, T, S, N, Q, H, P, G, -X";
        
        /**
         * DNA nucleotide residue grouping
         */
	public static String dna = "A, G, T, C, N-";
        
        /**
         * RNA nucleotide residue grouping
         */
	public static String rna = "A, G, U, C, N-"; 
        
        /**
         * Generic counter for determining number of residue pairs analyzed
         */
	public static int zcount = 0;
	
        /**
         * Initial system time
         */
	final long timeID = System.currentTimeMillis();
        
        /**
         * Initial unique ID, based on system time upon loading
         */
	public String uniqueID = Long.toString(timeID);
        
        /**
         * Counter for number of sequences
         */
	public static int seqCount = 0;  
        
        /**
         * Selection to report the instance of each possible combination at any given residue pair
         */
        public static boolean pairData = false;
        
        /**
         * Selection to signify comparison is intramolecular, and that MI should equal zero when x = y
         */
        public static boolean intraMol = false;
        
        /**
         * Verification that the first file has been loaded
         */
        public boolean chosenFile1 = false;
        
        /**
         * Verification that the second file has been loaded
         */
        public boolean chosenFile2 = false;
	
    @Override
    public void start(final Stage stage) {
        stage.setTitle("Biological Mutual Information");
 
		// Initialize Filechooser, Buttons
        final FileChooser fileChooser = new FileChooser();
        final Button openButton1 = new Button("Open Sequence File");
	final Button openButton2 = new Button("Open Sequence File");
        final Button processingButton = new Button("Process File");  
        final Button idButton = new Button("New Unique ID");  
	
        final ScrollPane sp = new ScrollPane();
		// Initialize Gridpane
        final GridPane inputGridPane = new GridPane();

		// Set Master Style Strings
		String bigText = "-fx-font: 35px Tahoma";
		
		// Initialize Texts. Contributed by Christopher Camenares, with Modification
		Text lbl1 = new Text(" Sequence Collection 1 ");
		lbl1.setStyle(bigText);
		
		Text lbl2 = new Text(System.lineSeparator() + " Sequence Collection 2 ");
		lbl2.setStyle(bigText);
		
		Text lbl3 = new Text("No File Loaded");
		
		Text lbl4 = new Text("No File Loaded");
		
		Text lbl5 = new Text("Presets");
		
		Text lbl6 = new Text("Presets");
		
		Text lbl7 = new Text(System.lineSeparator() + "Output Parameters");
		lbl7.setStyle(bigText);
		
		Text lbl8 = new Text("Job ID#: ");

		Text lbl9 = new Text(" Data Submission ");
		lbl9.setStyle(bigText);                
                
		Text statusLbl = new Text(System.lineSeparator() + "Program Awaiting Input");
		
		// Initialize ComboBoxes. Contributed by Christopher Camenares, with Modification
		final ComboBox<String> comboBox1;
		comboBox1 = new ComboBox<String>();
		comboBox1.getItems().addAll(protgrCat, aaCat, dnaCat, rnaCat, customCat);
		comboBox1.setMinWidth(50);
		comboBox1.setMinHeight(25);
		comboBox1.setValue(protgrCat);
		
		final ComboBox<String> comboBox2;
		comboBox2 = new ComboBox<String>();
		comboBox2.getItems().addAll(protgrCat, aaCat, dnaCat, rnaCat, customCat);
		comboBox2.setMinWidth(50);
		comboBox2.setMinHeight(25);
		comboBox2.setValue(protgrCat);
		
		HBox comboBoxSection1;
			comboBoxSection1 = new HBox();
			comboBoxSection1.setSpacing(10);
			comboBoxSection1.getChildren().addAll(comboBox1);

		HBox comboBoxSection2;
			comboBoxSection2 = new HBox();
			comboBoxSection2.setSpacing(10);
			comboBoxSection2.getChildren().addAll(comboBox2);
		
		// Initialize Text Fields, for adding Categories or JobID. Contributed by Christopher Camenares
		TextField txtField1;
		txtField1 = new TextField();
		txtField1.setText("Enter Custom");
		txtField1.setMinWidth(200);
		
		TextField txtField2;
		txtField2 = new TextField();
		txtField2.setText("Enter Custom");
		txtField2.setMinWidth(200);
		
		TextField txtField3;
		txtField3 = new TextField();
		txtField3.setText(uniqueID);
		txtField3.setMinWidth(100);

                TextField txtField4;
		txtField4 = new TextField();
		txtField4.setText("5");
		txtField4.setMaxWidth(100);  
                
		TextField txtField5;
		txtField5 = new TextField();
		txtField5.setText("1");
		txtField5.setMaxWidth(100);                
                
		// Correction Options: APC, RCW, etc
		Text lbl11 = new Text("Which Data to Output:");
		
		CheckBox chck1;
		chck1 = new CheckBox("Sort Scores - Default");
		chck1.setSelected(true);
		
		CheckBox chck2;
		chck2 = new CheckBox("APC");

		CheckBox chck3;
		chck3 = new CheckBox("RCW");

		CheckBox chck5;
		chck5 = new CheckBox("All Pair Data");                
                
		CheckBox chck4;
		chck4 = new CheckBox("Triplet Calculation");                

		CheckBox chck6;
		chck6 = new CheckBox("Substract Shuffled Sequence");                 

		CheckBox chck7;
		chck7 = new CheckBox("Substract Shuffled Sequence"); 
                
		CheckBox chck8;
		chck8 = new CheckBox("Substract Shuffled Order, Iterative");                

		CheckBox chck9;
		chck9 = new CheckBox("Save Shuffled Sequences");                

                CheckBox chck10;
		chck10 = new CheckBox("Save Shuffled Sequences");                

                CheckBox chck11;
		chck11 = new CheckBox("CAS");
                
		Text lbl12 = new Text("Heatmap Cell Size (Pixels):");                        
                                             
                        
		// Setup the Grid, Populate with items
		
		//Header 1
		inputGridPane.add(lbl1, 0, 0);
		inputGridPane.setColumnSpan(lbl1, 2);              
		inputGridPane.setHalignment(lbl1, HPos.CENTER);
        
		//Open Button 1
		inputGridPane.add(openButton1, 0, 1);
		
		//Open Button 2
		inputGridPane.add(openButton2, 0, 5);
		
		//Did Button 1 work?
		inputGridPane.add(lbl3, 1, 1);
		inputGridPane.setColumnSpan(lbl3, 2);

		//Did Button 2 work?		
		inputGridPane.add(lbl4, 1, 5);
		inputGridPane.setColumnSpan(lbl4, 2);
		
		//Header 2	
		inputGridPane.add(lbl2, 0, 4);
		inputGridPane.setColumnSpan(lbl2, 2);               
		inputGridPane.setHalignment(lbl2, HPos.CENTER);
		
		//Presets
		// inputGridPane.add(lbl5, 0, 2);

		//Categories 1
		inputGridPane.add(comboBoxSection1, 0, 2);

		//Presets
		// inputGridPane.add(lbl6, 0, 6);
		
		//Categories 2		
		inputGridPane.add(comboBoxSection2, 0, 6);		
		
		//Custom Categories
		inputGridPane.add(txtField1, 1, 2);
		// inputGridPane.setColumnSpan(txtField1, 2);		

		// Shuffle DNA 1
		inputGridPane.add(chck6, 0, 3);                
		inputGridPane.add(chck9, 1, 3);
                
                // Shuffle DNA 2
		inputGridPane.add(chck7, 0, 7);
		inputGridPane.add(chck10, 1, 7);
                
		//Custom Categories		
		inputGridPane.add(txtField2, 1, 6);
		// inputGridPane.setColumnSpan(txtField2, 2);	              
                
		//Output		
		inputGridPane.add(lbl7, 0, 8);
		inputGridPane.setColumnSpan(lbl7, 2);
		inputGridPane.setHalignment(lbl7, HPos.CENTER);		
            
		//Correction		
		inputGridPane.add(lbl11, 0, 9);
		inputGridPane.add(chck1, 0, 10);				
		inputGridPane.add(chck2, 0, 11);	
		inputGridPane.add(chck3, 0, 12);			
		inputGridPane.add(chck11, 0, 13);
                inputGridPane.add(chck5, 0, 14);

		//Graphics		
		inputGridPane.add(lbl12, 1, 9);
		inputGridPane.add(txtField4, 1, 10);			
                
		// Triplet Scoring
		inputGridPane.add(chck4, 1, 11);

		// Shuffle Organism 1 (Sufficient vs other collection)
		inputGridPane.add(chck8, 1, 12);
		// Shuffle Iteration		
		inputGridPane.add(txtField5, 1, 14);                
                
		//Submission	
		inputGridPane.add(lbl9, 0, 15);
		inputGridPane.setColumnSpan(lbl9, 2);
		inputGridPane.setHalignment(lbl9, HPos.CENTER);	
                
		//Job ID Label
		inputGridPane.add(lbl8, 0, 16);
		inputGridPane.setHalignment(lbl8, HPos.RIGHT);		

		//Job ID Entry		
		inputGridPane.add(txtField3, 1, 16);	

		//ID Button		
		inputGridPane.add(idButton, 0, 17);
		inputGridPane.setColumnSpan(idButton, 2);                
		inputGridPane.setHalignment(idButton, HPos.CENTER);
                
		//Status Update
		inputGridPane.add(statusLbl, 0, 18);
		inputGridPane.setColumnSpan(statusLbl, 2);		
		inputGridPane.setHalignment(statusLbl, HPos.CENTER);		
		
		//Process Button		
		inputGridPane.add(processingButton, 0, 19);
		inputGridPane.setColumnSpan(processingButton, 2);		
		inputGridPane.setHalignment(processingButton, HPos.CENTER);
		
	inputGridPane.setHgap(25);
        inputGridPane.setVgap(5);
        
        final Pane rootGroup = new VBox(15);
        rootGroup.getChildren().addAll(inputGridPane);
        rootGroup.setPadding(new Insets(25, 25, 25, 25));
 
		// Button to Generate New ID
        idButton.setOnAction(
			new EventHandler<ActionEvent>() {
                        @Override
                        public void handle(final ActionEvent e) {
                            final long timeButton = System.currentTimeMillis();                   
                            txtField3.setText(Long.toString(timeButton));
                    }
                }
            );
 
	    // Hanlde Events for Opening File#1, read into program as String
        openButton1.setOnAction(
            new EventHandler<ActionEvent>() {
                @Override
                public void handle(final ActionEvent e) {
                    configureFileChooser(fileChooser);
                    file1 = fileChooser.showOpenDialog(stage);
                    if (file1 != null) {
				flName1 = file1.getName();
                                flPath1 = file1.getPath();
				lbl3.setText(flName1);
				statusLbl.setText("File 1 Loaded, Awaiting File 2");	
				
                                chosenFile1 = true;
                                
					if (chosenFile2){
						statusLbl.setText(System.lineSeparator() + "Two Files loaded, Ready to Process!");	
                                                }
						
                    }
                }
            });

	    // Hanlde Events for Opening File#2, read into program as String			
		openButton2.setOnAction(
            new EventHandler<ActionEvent>() {
                @Override
                public void handle(final ActionEvent e) {
                    configureFileChooser(fileChooser);
                    file2 = fileChooser.showOpenDialog(stage);
                    if (file2 != null) {
				flName2 = file2.getName();
                                flPath2 = file2.getPath();                                
				lbl4.setText(flName2);   
				statusLbl.setText("File 2 Loaded, Awaiting File 1");
                                
                                chosenFile2 = true;
                                
					if (chosenFile1){
						statusLbl.setText(System.lineSeparator() + "Two Files loaded, Ready to Process!");	
                                                }						
						
                    }
                }
            });
 
 	    // Hanlde Events for Calculations and Saving the File
        processingButton.setOnAction(
            new EventHandler<ActionEvent>() {
                @Override
                public void handle(final ActionEvent e) {
                        boolean reportMI = chck1.isSelected();
                        boolean reportAPC = chck2.isSelected();
                        boolean reportRCW = chck3.isSelected();
                        boolean reportCAS = chck11.isSelected();
                        boolean tripletZ = chck4.isSelected();
                        boolean shuffleSeq1 = chck6.isSelected();
                        boolean shuffleSeq2 = chck7.isSelected();
                        boolean shuffleOrg1 = chck8.isSelected();
                        boolean saveShuf1 = chck9.isSelected();
                        boolean saveShuf2 = chck10.isSelected();
                        pairData = chck5.isSelected();
 
                        int processCount = 1;
                        
                        /**
                         * Beginning timestamp to track processing time
                         */
			final long timeStart = System.currentTimeMillis();                        

                        // Interpret category choice and identity, with user input specified
			String catA1 = catSwitch(comboBox1.getValue(), txtField1.getText());
			String catB1 = catSwitch(comboBox2.getValue(), txtField2.getText());

                        // Interpret user choice regarding shuffling
                        int shuffleIt1 =  Integer.parseInt(txtField5.getText());
                        
                        // Interpret user choice regarding graphical interface
                        int pixelSize =  Integer.parseInt(txtField4.getText());
                        
                        // Initialize temp. string variable
                        String reporter = "";
                                  
                        // Get new name for directory, initialize directory
			String dirName = txtField3.getText();
			File dir = new File(dirName);
			dir.mkdir();
                        
                        File dir_sub = null;
                        String dir_subName = null;
                                                         
                        //Initalize final files
                        /**
                         * Recipient file for runtime and processing information
                         */
			File fileRun = new File(dirName + "\\Runtime_Info_" + dirName + ".txt");                                        
                                  
                        /**
                         * Recipient file for raw MI scores from comparing collection A against collection B
                         */
			File fRawMIxy = new File(dirName + "\\rawMIxy" + ".txt");                        
                        
                        /**
                         * Recipient file for raw MI scores from comparing collection A against itself
                         */
			File fRawMIx = new File(dirName + "\\rawMIx" + ".txt");
         
                        /**
                         * Recipient file for raw MI scores from comparing collection B against itself
                         */
			File fRawMIy = new File(dirName + "\\rawMIy" + ".txt");                                  

                        /**
                         * Recipient file for sorted MI scores from comparing collection A against collection B
                         */          
                        File fSortMIxy = new File(dirName + "\\sortedMIxy" + ".txt");
       
                        /**
                         * Recipient file for sorted MI scores from comparing collection A against itself
                         */
                        File fSortMIx = new File(dirName + "\\sortedMIx" + ".txt");

                        /**
                         * Recipient file for sorted MI scores from comparing collection B against itself
                         */                        
                        File fSortMIy = new File(dirName + "\\sortedMIy" + ".txt");                              
                                  
                        // If pair data is selected, initialize the file

                        if (pairData)
                        {
                            String[] tempCatA = catA1.split(", ");
                            String[] tempCatB = catB1.split(", ");
                                      
                            ArrayList<String> tempCatCombo = new ArrayList<>();
                                      
                            for (String tempCatA1 : tempCatA) {
                                for (String tempCatB1 : tempCatB) {
                                    tempCatCombo.add(tempCatA1 + "%" + tempCatB1);
                                }
                            }
                                      
                        reporter = "\t" + "\t" + "\t";
                                      
                        for (String tempCatCombo1 : tempCatCombo) {
                            reporter += tempCatCombo1 + "\t";
                        }
                                    
                        writeFileN(reporter, fRawMIxy);
                        }
                        
			// Convert category string into array
			char[][] catA2 = formatCat(catA1);
			char[][] catB2 = formatCat(catB1);	
                 
                        statusLbl.setText(System.lineSeparator() + "Processing Incomplete: Please Standby or Restart");
				                          
                        // End of common code for processing run
                        
                 // If only file 2 is selected
                 if   (file1 == null & file2 != null)
                 {
                     file1 = file2;
                     file2 = null;
                     chosenFile1 = true;
                     chosenFile2 = false;
                     
                     if (shuffleSeq2 && !shuffleSeq1)
                     {
                         shuffleSeq1 = true;
                         shuffleSeq2 = false;
                     }
                     if (saveShuf2 && !saveShuf1)
                     {
                         saveShuf1 = true;
                         saveShuf2 = false;
                     }
                 }
                 
                 // If only 1 file is selected
                 if (chosenFile1 && !chosenFile2)
                  {
                        intraMol = true;
                        
                        char[][] seqA1 = null;
                        int[][][] seqA2 = null;                                    
                        
                        seqA1 = formatFasta(file1);
                        
                        int countSeqA = seqCount;
                       
                        seqA2 = ccSeq(seqA1, catA2);
			
                        double[][][] arrayMIx =  mutualInfo(seqA2, seqA2, true, true, fRawMIx);
                        
                        drawHeatMap(arrayMIx[0], dirName + "\\heatmapMIx.PNG", pixelSize, "gray");                           
                            
                        if(reportMI)
                         {					  
                            sortReport(arrayMIx[0], fSortMIx);
                            
                         }
                          
                        String res1 = "";
                        
                            
                        if(reportAPC)
                         {
                            dir_subName = dirName + "\\APCx";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            apcProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);
			 }

                        if(reportRCW)
			 {
                            dir_subName = dirName + "\\RCWx";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            rcwProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize); 
                         }
                    
                        
                    // Initiate shuffle / background routine    
                    if (shuffleSeq1)
                    {       

                        //Initialize Files
                        /**
                         * Recipient file for raw MI scores from comparing collection A against itself, shuffle subtracted
                         */
                        dir_subName = dirName + "\\shuffleX";
                        dir_sub = new File(dir_subName);
                        dir_sub.mkdir();
                        
                        File fSortMIxShSub = new File(dirName + "\\sortedMIxShuffle" + ".txt");
                        File fRawMIxShSub = new File(dirName + "\\rawMIxShuffle" + ".txt");
                        
                        countSeqA = seqCount;
                        
                        char[][] seqA1s = null;  
                        
                        double shuffleIt1D = (double)shuffleIt1;
                        
                        double[][][] arrayMIxS1 = arrayMIx;                        
                      
                      for (int sh = 0; sh < shuffleIt1; sh++)
                      {
                          
                        File fRawMIxShN = new File(dir_sub + "\\rawMIxShuffle_" + sh + ".txt");

                        File fSortMIxShN = new File(dir_sub + "\\sortedMIxShuffle_" + sh + ".txt");                        
                       
                        // Generate Shuffle, substract shuffle
                        seqA1s = new char[seqA1.length][seqA1[0].length];
                        
                        for (int i = 0; i < seqA1.length; i++)
                        {
                            int shuffleArr[] = shuffleRandom(seqA1[i].length);
                            
                            for (int j = 0; j < shuffleArr.length; j++)
                            {
                                int newPosition = shuffleArr[j];
                                
                                seqA1s[i][j] = seqA1[i][newPosition];
                            }
                            
                        }
                        
                        seqA2 = ccSeq(seqA1s, catA2);                         
                        
                        double[][][] arrayMIxS2 = mutualInfo(seqA2, seqA2, true, saveShuf1, fRawMIxShN);

                        if(reportMI && saveShuf1)
                         {					  
                            sortReport(arrayMIxS2[0], fSortMIxShN);
                            
                         }
                        
                        arrayMIxS1 = average3Darr(arrayMIxS1, arrayMIxS2, sh);
                        
                      }
                        arrayMIx = subtract3Darr(arrayMIx, arrayMIxS1);                      
                      
                        write2Darr(arrayMIx[0], fRawMIxShSub);
                        
                        drawHeatMap(arrayMIx[0], dirName + "\\heatmapMIx_shuffle_subtracted.PNG", pixelSize, "gray");
                            
                        if(reportMI)
                         {					  
                            sortReport(arrayMIx[0], fSortMIxShSub);
                            
                         }
                          
                        res1 = "";
                            
                            
                        if(reportAPC)
                         {
                            dir_subName = dirName + "\\APCx_" + "_subtract_shuffle";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            apcProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);
			 }

                        if(reportRCW)
			 {
                            dir_subName = dirName + "\\RCWx_" + "_subtract_shuffle";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            rcwProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);                           
                         }
                    }
                    
                    // Get final runtime statistics
			final long timeEnd = System.currentTimeMillis();
			long runMinutesL = (timeEnd - timeStart) / 60000;
			int runMinutesI = (int) runMinutesL;
			String runRep = "Job ID#: " + dirName;
			runRep += System.lineSeparator();
			runRep += System.lineSeparator();
			runRep += Long.toString(timeEnd - timeStart) + " milliseconds of processing time";
			runRep += System.lineSeparator();
			runRep += "- or about " + Integer.toString(runMinutesI) + " minutes";
			runRep += System.lineSeparator();
			runRep += System.lineSeparator();
			runRep += "Source file 1: " + flPath1 + ", contained " + countSeqA + " sequences, aligned to a length of " + Integer.toString(seqA2[0].length);		
			runRep += System.lineSeparator();
			runRep += "Analyzed with category identities: " + catA1;						  
			runRep += System.lineSeparator();				  				  
			runRep += System.lineSeparator();				  
			runRep += zcount + " total residue pairs analyzed";
                        if(shuffleSeq1)
                        {
                        runRep += System.lineSeparator();				  
			runRep += shuffleIt1 + " iterations of shuffling subtracted from sequence";						  
                        }	  
			// Save runtime information

			writeFile(runRep, fileRun);
				  
                        openFile(dir);                                  
                                  
			statusLbl.setText(System.lineSeparator() + "Processing Complete! Ready for new input.");				  
                                      
                        
                    
                  }                 
                 // If both files are selected, proceed with full routine.
                 else if (chosenFile1 && chosenFile2)
		   {
				  
                        // Initialize all variables
                        char[][] seqA1 = null;
                        char[][] seqB1 = null;
                        int[][][] seqA2 = null;
                        int[][][] seqB2 = null;
                        
                        // Read and process file 1
                        seqA1 = formatFasta(file1);
                        int countSeqA = seqCount;                        
                        seqA2 = ccSeq(seqA1, catA2);

                        // Read and process file 2                        
                        seqB1 = formatFasta(file2);
                        int countSeqB = seqCount;                        
                        seqB2 = ccSeq(seqB1, catB2);
                                                
                          
                        double[][][] arrayMIx =  mutualInfo(seqA2, seqA2, true, true, fRawMIx);                                                       
                        
                        double[][][] arrayMIy =  mutualInfo(seqB2, seqB2, true, true, fRawMIy);                        
                        
                        double[][][] arrayMIxy =  mutualInfo(seqA2, seqB2, false, true, fRawMIxy);                        
                        
                        
                        if(tripletZ)
                        {
                        /**
                         * Recipient File for Triplet Calculations
                         */
			File fRawMIxyz = new File(dirName + "\\rawMIxyz" + ".txt");
                        
                        /**
                         * Recipient File for Triplet Calculations
                         */
			File fSortMIxyz = new File(dirName + "\\sortedMIxyz" + ".txt");    
                            
                        double[][] arrayMIxyZ = new double[seqA2[0].length][seqB2[0].length];
                        arrayMIxyZ = dumdoub2Darr(arrayMIxyZ);

                            for (int i = 0; i < arrayMIxyZ.length; i++)
                            {
                                for (int j = 0; j < arrayMIxyZ[i].length; j++)
                                {
                                    arrayMIxyZ[i][j] = arrayMIxy[0][i][j] / (arrayMIx[1][0][i] * arrayMIy[2][0][j]);
                                    if(arrayMIx[1][0][i] == 0 || arrayMIy[2][0][j] == 0)
                                    {
                                        arrayMIxyZ[i][j] = 0.0;
                                    }
                                    
                                    String zRep = Integer.toString(i + 1) + "\t" + Integer.toString(j + 1) + "\t" +  String.format("%.10f", arrayMIxyZ[i][j]);
                                    writeFileN(zRep, fRawMIxyz);
                                    
                                }
                            }
                            
                            if(reportMI)
                            {
                            sortReport(arrayMIxyZ, fSortMIxyz);
                            }
                            
                            drawHeatMap(arrayMIxyZ, dirName + "\\heatmapMIxyZ.PNG", pixelSize, "red");                            
                        }
                            
                            // sortReport, Draw Heatmaps
                            if(reportMI)
                            {
                                sortReport(arrayMIx[0], fSortMIx);
                                sortReport(arrayMIy[0], fSortMIy);
                                sortReport(arrayMIxy[0], fSortMIxy);                                  
                            }
                              
                                drawHeatMap(arrayMIx[0], dirName + "\\heatmapMIx.PNG", pixelSize, "gray");
                                drawHeatMap(arrayMIy[0], dirName + "\\heatmapMIy.PNG", pixelSize, "gray");
                                drawHeatMap(arrayMIxy[0], dirName + "\\heatmapMIxy.PNG", pixelSize, "green");
                                
				  if(reportAPC)
                                  {
                                    dir_subName = dirName + "\\APCx";
                                    dir_sub = new File(dir_subName);
                                    dir_sub.mkdir();
                                    apcProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);				  

                                    dir_subName = dirName + "\\APCy";
                                    dir_sub = new File(dir_subName);
                                    dir_sub.mkdir();
                                    apcProcess(arrayMIy, dir_subName, seqB2[0].length, seqB2[0].length, pixelSize);
                                    
                                    dir_subName = dirName + "\\APCxy";
                                    dir_sub = new File(dir_subName);
                                    dir_sub.mkdir();
                                    apcProcess(arrayMIxy, dir_subName, seqA2[0].length, seqB2[0].length, pixelSize);
                                  
                                  }
				  
				  if(reportRCW)
				  {
                                    dir_subName = dirName + "\\RCWx";
                                    dir_sub = new File(dir_subName);
                                    dir_sub.mkdir();
                                    rcwProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);				  

                                    dir_subName = dirName + "\\RCWy";
                                    dir_sub = new File(dir_subName);
                                    dir_sub.mkdir();
                                    rcwProcess(arrayMIy, dir_subName, seqB2[0].length, seqB2[0].length, pixelSize);
                                    
                                    dir_subName = dirName + "\\RCWxy";
                                    dir_sub = new File(dir_subName);
                                    dir_sub.mkdir();
                                    rcwProcess(arrayMIxy, dir_subName, seqA2[0].length, seqB2[0].length, pixelSize);
                                  
                                  }

                                  if(reportCAS)
                                {
                                dir_subName = dirName + "\\CASx";
                                dir_sub = new File(dir_subName);
                                dir_sub.mkdir();
                                casProcess(arrayMIx, dir_subName, pixelSize, seqA2, seqA2);                           

                                dir_subName = dirName + "\\CASy";
                                dir_sub = new File(dir_subName);
                                dir_sub.mkdir();                                
                                casProcess(arrayMIy, dir_subName, pixelSize, seqB2, seqB2);
                                
                                dir_subName = dirName + "\\CASxy";
                                dir_sub = new File(dir_subName);
                                dir_sub.mkdir();
                                casProcess(arrayMIxy, dir_subName, pixelSize, seqA2, seqB2); 
                                }                                  

                    // Initiate shuffle / background routine    
                    if (shuffleSeq1 || shuffleSeq2 || shuffleOrg1)
                    {       

                        //Initialize Files
                        /**
                         * Recipient file for raw MI scores from comparing collection A against itself, shuffle subtracted
                         */
                        dir_subName = dirName + "\\shuffleX_XY";
                        dir_sub = new File(dir_subName);
                        dir_sub.mkdir();
                        
                        File fRawMIxShSub = new File(dirName + "\\RawMIxShuffleSub" + ".txt");
                        File fRawMIyShSub = new File(dirName + "\\RawMIyShuffleSub" + ".txt");                        
                        File fRawMIxyShSub = new File(dirName + "\\RawMIxyShuffleSub" + ".txt");
                        
                        File fSortMIxShSub = new File(dirName + "\\sortedMIxShuffleSub" + ".txt");
                        File fSortMIyShSub = new File(dirName + "\\sortedMIyShuffleSub" + ".txt");                        
                        File fSortMIxyShSub = new File(dirName + "\\sortedMIxyShuffleSub" + ".txt");
                        
                        countSeqA = seqCount;
                        
                        char[][] seqA1s = null;  
                        char[][] seqB1s = null; 
                        
                        double shuffleIt1D = (double)shuffleIt1;

                        double[][][] arrayMIxS1 = new double[arrayMIx.length][][];
                        arrayMIxS1 = clone3Darr(arrayMIx);
                        double[][][] arrayMIyS1 = new double[arrayMIy.length][][];
                        arrayMIyS1 = clone3Darr(arrayMIy);
                        double[][][] arrayMIxyS1 = new double[arrayMIxy.length][][];
                        arrayMIxyS1 = clone3Darr(arrayMIxy);                        

                      for (int sh = 0; sh < shuffleIt1; sh++)
                      {
                        

                        File fRawMIxShN = new File(dir_sub + "\\rawMIxShuffle_" + sh + ".txt");
                        File fSortMIxShN = new File(dir_sub + "\\sortedMIxShuffle_" + sh + ".txt");                          

                        File fRawMIyShN = new File(dir_sub + "\\rawMIyShuffle_" + sh + ".txt");
                        File fSortMIyShN = new File(dir_sub + "\\sortedMIyShuffle_" + sh + ".txt");                        

                        File fRawMIxyShN = new File(dir_sub + "\\rawMIxyShuffle_" + sh + ".txt");
                        File fSortMIxyShN = new File(dir_sub + "\\sortedMIxyShuffle_" + sh + ".txt");                             

                        // Generate Shuffle, substract shuffle
                        seqA1s = new char[seqA1.length][seqA1[0].length];
                        seqB1s = new char[seqB1.length][seqB1[0].length];
                        
                        if(shuffleSeq1)
                        {
                         for (int i = 0; i < seqA1.length; i++)
                         {
                            int shuffleArr[] = shuffleRandom(seqA1[i].length);
                            
                            for (int j = 0; j < shuffleArr.length; j++)
                            {
                                int newPosition = shuffleArr[j];
                                
                                seqA1s[i][j] = seqA1[i][newPosition];
                            }
                            
                         }                            
                        }
                        else
                        {
                            seqA1s = seqA1;
                        }

                        if(shuffleSeq2)
                        {
                         for (int i = 0; i < seqB1.length; i++)
                         {
                            int shuffleArr[] = shuffleRandom(seqB1[i].length);
                            
                            for (int j = 0; j < shuffleArr.length; j++)
                            {
                                int newPosition = shuffleArr[j];
                                
                                seqB1s[i][j] = seqB1[i][newPosition];
                            }
                            
                         }                            
                        }
                        else
                        {
                            seqB1s = seqB1;
                        }

                        int shuffSize = Math.min(seqA1.length, seqB1.length);
                        
                        if(shuffleOrg1)
                        {
                            int shuffleArr[] = shuffleRandom(shuffSize);
                            
                            for (int i = 0; i < shuffleArr.length; i++)
                            {
                                int newPosition = shuffleArr[i];
                                
                                seqA1s[i] = seqA1s[newPosition];
                                seqB1s[i] = seqB1s[newPosition];
                            }                            
                        }
                        
                        seqA2 = ccSeq(seqA1s, catA2);                         
                        seqB2 = ccSeq(seqB1s, catB2);                         
                        
                        double[][][] arrayMIxS2 = mutualInfo(seqA2, seqA2, true, saveShuf1, fRawMIxShN);
                        double[][][] arrayMIyS2 = mutualInfo(seqB2, seqB2, true, saveShuf2, fRawMIyShN);                        
                        double[][][] arrayMIxyS2 = mutualInfo(seqA2, seqB2, false, (saveShuf1 || saveShuf2), fRawMIxyShN);
                        
                        if(reportMI && saveShuf1)
                         {					  
                            sortReport(arrayMIxS2[0], fSortMIxShN);                         
                         }
                        
                        if(reportMI && saveShuf2)
                         {					  
                            sortReport(arrayMIyS2[0], fSortMIyShN);                            
                         }
                        
                        if(reportMI && (saveShuf1 || saveShuf2))
                         {					  
                            sortReport(arrayMIxyS2[0], fSortMIxyShN);                            
                         }                        
                        
                        arrayMIxS1 = average3Darr(arrayMIxS1, arrayMIxS2, (sh + 1));
                        arrayMIyS1 = average3Darr(arrayMIyS1, arrayMIyS2, (sh + 1));                        
                        arrayMIxyS1 = average3Darr(arrayMIxyS1, arrayMIxyS2, (sh + 1));
                      }
                      
                        arrayMIx = subtract3Darr(arrayMIx, arrayMIxS1);
                        arrayMIy = subtract3Darr(arrayMIy, arrayMIyS1);                        
                        arrayMIxy = subtract3Darr(arrayMIxy, arrayMIxyS1);                        

                        write2Darr(arrayMIx[0], fRawMIxShSub);
                        write2Darr(arrayMIy[0], fRawMIyShSub);
                        write2Darr(arrayMIxy[0], fRawMIxyShSub);                        
                        
                        drawHeatMap(arrayMIx[0], dirName + "\\heatmapMIx_shuffle_subtracted.PNG", pixelSize, "gray");
                        drawHeatMap(arrayMIy[0], dirName + "\\heatmapMIy_shuffle_subtracted.PNG", pixelSize, "gray");                        
                        drawHeatMap(arrayMIxy[0], dirName + "\\heatmapMIxy_shuffle_subtracted.PNG", pixelSize, "blue");

                        if(reportMI)
                         {					  
                            sortReport(arrayMIx[0], fSortMIxShSub);
                            sortReport(arrayMIy[0], fSortMIyShSub);                            
                            sortReport(arrayMIxy[0], fSortMIxyShSub); 
                         }                       
  
                        if(reportAPC)
                         {
                            dir_subName = dirName + "\\APC" + "_subtract_shuffle";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            apcProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);
                            apcProcess(arrayMIy, dir_subName, seqB2[0].length, seqB2[0].length, pixelSize);                            
                            apcProcess(arrayMIxy, dir_subName, seqA2[0].length, seqB2[0].length, pixelSize);
                         }

                        if(reportRCW)
			 {
                            dir_subName = dirName + "\\RCW" + "_subtract_shuffle";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            rcwProcess(arrayMIx, dir_subName, seqA2[0].length, seqA2[0].length, pixelSize);                           
                            rcwProcess(arrayMIy, dir_subName, seqB2[0].length, seqB2[0].length, pixelSize);
                            rcwProcess(arrayMIxy, dir_subName, seqA2[0].length, seqB2[0].length, pixelSize); 
                         }
                   
                        if(reportCAS)
			 {
                            dir_subName = dirName + "\\CAS" + "_subtract_shuffle";
                            dir_sub = new File(dir_subName);
                            dir_sub.mkdir();
                            casProcess(arrayMIx, dir_subName, pixelSize, seqA2, seqA2);                           
                            casProcess(arrayMIy, dir_subName, pixelSize, seqB2, seqB2);
                            casProcess(arrayMIxy, dir_subName, pixelSize, seqA2, seqB2); 
                         }
                        
                    }                                      
                                  
				    // Get final runtime statistics
				  final long timeEnd = System.currentTimeMillis();
				  long runMinutesL = (timeEnd - timeStart) / 60000;
				  int runMinutesI = (int) runMinutesL;
				  String runRep = "Job ID#: " + dirName;
				  runRep += System.lineSeparator();
				  runRep += System.lineSeparator();
				  runRep += Long.toString(timeEnd - timeStart) + " milliseconds of processing time";
				  runRep += System.lineSeparator();
				  runRep += "- or about " + Integer.toString(runMinutesI) + " minutes";
				  runRep += System.lineSeparator();
				  runRep += System.lineSeparator();
				  runRep += "Source file 1: " + flPath1 + ", contained " + countSeqA + " sequences, aligned to a length of " + Integer.toString(seqA2[0].length);		
				  runRep += System.lineSeparator();
				  runRep += "Analyzed with category identities: " + catA1;						  
				  runRep += System.lineSeparator();				  
				  runRep += System.lineSeparator();				  
				  runRep += "Source file 2: " + flPath2 + ", contained " + countSeqB + " sequences, aligned to a length of " + Integer.toString(seqB2[0].length);				  
				  runRep += System.lineSeparator();
				  runRep += "Analyzed with category identities: " + catB1;	
				  runRep += System.lineSeparator();				  
				  runRep += System.lineSeparator();				  
				  runRep += zcount + " total residue pairs analyzed";
                                  if(shuffleSeq1 || shuffleSeq2 || shuffleOrg1)
                                  {
                                  runRep += System.lineSeparator();				  
				  runRep += shuffleIt1 + " iterations of shuffling subtracted from sequence";					  
                                  }
					// Save runtime information

				  writeFile(runRep, fileRun);
				  
                                  openFile(dir);                                  
                                  
				  statusLbl.setText(System.lineSeparator() + "Processing Complete! Ready for new input.");				  
                                  
				 }
				 else
				 {
				  txtField1.setText("Input parameters undefined");
				 }
				}
				});
 
        sp.setContent(rootGroup);
        stage.setScene(new Scene(sp));
        stage.show();
    }
 
	// Main Method, Launch Program
    public static void main(String[] args) {
        Application.launch(args);
    }
	
	// Group String Array into String, line Separator (Used for Debugging)
	private static String arrStr(String[] strArr){
		StringBuilder strBuilder = new StringBuilder();
		for (int i = 0; i < strArr.length; i++) 
		{
		strBuilder.append(strArr[i]);
		strBuilder.append(System.lineSeparator());
		}
		String newString = strBuilder.toString();
		return newString;
	}
	
	// Group 2D char Array into Strings, line Separator for 1st Dimension (Used for Debugging)
	private static String arr2Dstr(char[][] strArr){
		String res = new String();
		for (int i = 0; i < strArr.length; i++) 
		{
		 String inter = new String(strArr[i]);
		 res += inter;
		 res += System.lineSeparator();
		}
		return res;
	}

	// Format Sequences and Create Array; splitting by sequence; then, flip position & organism dimensions
	private static char[][] formatFasta(File fileA){
            /**
             * Boolean to determine if new sequence body should be added
             */
            boolean addBody = false;
            
            /**
             * Temporary container to hold or build the sequence body from lines
             */
            String seqBody = "";
            
            seqCount = 0;
            
            ArrayList<String> seqArr = new ArrayList<>();
            
            try (BufferedReader br = new BufferedReader(new FileReader(fileA))) {
                String line;
                
                while ((line = br.readLine()) != null) {
                // process the line.

                if (line.contains(">") && addBody)    
                {
                    //Process preceeding sequence
                    seqArr.add(seqBody);
                    
                    // Reset Parameters
                    addBody = false;
                    seqBody = "";
                    
                }
                
                if (line.contains(">"))
                {
                    addBody = true;
                    seqCount++;
                }
                else if (addBody)
                {
                    seqBody += line.toUpperCase();
                }
                
              }
                // This code is for the last sequence block in collection 2
                if (addBody)    
                {
                    //Process preceeding sequence
                    seqArr.add(seqBody);
                }                
                
            } catch (FileNotFoundException ex) {
              Logger.getLogger(MutualInfoBio.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
              Logger.getLogger(MutualInfoBio.class.getName()).log(Level.SEVERE, null, ex);
            } 
		  
	char[][] seqC = new char[seqArr.size()][];
		  
	int[] lenB = new int[seqArr.size()];
		  
	for (int i = 0; i < seqArr.size(); i++)
            {
		seqC[i] = seqArr.get(i).toCharArray();
		lenB[i] = seqC[i].length;
            }
		  
        int maxLen = 0;
		  
	for (int i = 0; i < lenB.length; i++)
            {
                if (lenB[i] > maxLen)
		{
                    maxLen = lenB[i];
		}
            }
		  
        // Flip array
	char[][] seqD = new char[maxLen][seqC.length];
		  
        for (int i = 0; i < seqC.length; i++)
            {
                for (int j = 0; j < seqC[i].length; j++)
                    {
                        seqD[j][i] = seqC[i][j];
                    }
            }
		  
        return seqD;
        
	}
	
	// Process Category String into Groups and Identity Double Array.
	private static char[][] formatCat(String catA){
		// Find the number of categories or groups to be analyzed, Collection A
		int catlenA = catA.length() - catA.replace(",", "").length() + 1;

		// Split the category A string, by category, into separate arrays
		String[] catB = catA.split(", ");	
	
		// Within category A array, sub-split any identities that are to be considered equal
		char[][] catC = new char[catB.length][];
		for (int i = 0; i < catB.length; i++)
		{
			catC[i] = catB[i].toCharArray();
		}
		
		return catC;
	}
	
	// Determine Category String to use, Based on User Choise
	private static String catSwitch(String choice, String userI){
		String res = new String();
		switch (choice)
		{
			case "Protein, Grouped": res = protgr; break;
			case "Protein": res = aa; break;
			case "DNA": res = dna; break;
			case "RNA": res = rna; break;
			case "Custom Category": res = userI; break;
		}
		return res;
	}
	
	// Subroutine for converting sequence string based on categories
	// Label sequence identity according to category groups, Count instances
	private static int[][][] ccSeq(char[][] seq, char[][] cat){
		int[][] seq1 = new int[seq.length][seq[0].length];
		int[][] obs1 = new int[seq.length][cat.length];
		
		 for (int i = 0; i < seq.length; i++)
		{	
			for (int j = 0; j < seq[i].length; j++)
			{
				for (int k = 0; k < cat.length; k++)
				{    
					for (int m = 0; m < cat[k].length; m++)
					{
						if (seq[i][j] == cat[k][m])
						{
							seq1[i][j] = k;
							obs1[i][k] += 1;
						}
					}
				}
			}
		}
     
     
     int[][][] res = new int[2][seq1.length][obs1.length];
     res[0] = seq1;
     res[1] = obs1;

	 return res;
		
	}
	
	/**
         * Takes a sequence set and an observed count and determines the frequency
         * @param combo A 3D dimensional array: the first dimension splits between the sequence identity and the observed count 2D arrays.
         * @return 
         */
	private static double[][] oFreq(int[][][] combo){
		double[][] res = new double[combo[1].length][combo[1][0].length];
		   for (int i = 0; i < combo[1].length; i++)
			{
				for (int j = 0; j < combo[1][i].length; j++)
				{   
					int calc = (combo[1][i][j] * 1000) / (combo[0][i].length);
					 double calc2 = (double) calc;
					 calc2 = calc2 / 1000.0;
					res[i][j] = calc2;
                                        
				}
			};
		return res;
	}
	
        /**
         * MutualInfo Method: Calculate MI for all pairs across two sets of biological data
         * @param seq1 The 2nd sequence collection
         * @param seq2 The first sequence collection
         * @param intraMol Trigger switch to determine if diagonals (i = j) are ignored.
         * @param saveIt Trigger switch to determine if the results are saved
         * @param resFile File to which the information is saved
         * @return A 3D array; In first dimension, 0 = MI scores(2D), 1,0 = Row sums(1D), 2,0 = Col sums(1D), 3,0,0 = total sum(0D)
         */
	private double[][][] mutualInfo(int[][][] seq1, int[][][] seq2, boolean intraMol, boolean saveIt, File resFile)
        {

                 // BEGIN BLOCK FOR MI METHOD
                        double[][] pX = null;
                        double[][] pY = null;
                 
                        pX = oFreq(seq1);
                        pY = oFreq(seq2);
                        
                        double[] sumMIi = new double[seq1[0].length];
                        double[] sumMIj = new double[seq2[0].length];                                  
                        
                        sumMIi = dumdoub1Darr(sumMIi);
                        sumMIj = dumdoub1Darr(sumMIj);                        
                        
                        double sumMIij = 0.0;
                        
                        double[][] arrayMImet = new double[seq1[0].length][seq2[0].length];
                                    
                            for (int i = 0; i < seq1[0].length; i++)
                                {
                                    for (int j = 0; j < seq2[0].length; j++)
                                        {
                                            double ijMI = 0.0;
                                            // Find and Report raw MI

                                            ijMI = miIJ(seq1[0][i], seq2[0][j], pX[i], pY[j], i, j, saveIt, resFile);                                               

                                            if(intraMol && i == j)
                                            {
                                                ijMI = 0.0;
                                            }
                                            
                                            // Add MI to Arrays for Later Processing					
                                            arrayMImet[i][j] += ijMI;
						
                                            // Add to sumMIi array
                                            sumMIi[i] += ijMI;
						
                                            // Add to sumMIj array
                                            sumMIj[j] += ijMI;
						
                                            // Add to overall total
                                            sumMIij += ijMI;

                                            // Add to counter
                                            zcount++;
                                        }
                                }
                            
                        double[][][] resArray = new double[4][][];
                        
                        resArray[1] = new double[1][];
                        resArray[2] = new double[1][];
                        
                        resArray[0] = arrayMImet;
                        
                        resArray[1][0] = sumMIi;
                        
                        resArray[2][0] = sumMIj;
                        
                        resArray[3] = new double[1][1];
                        
                        resArray[3][0][0] = sumMIij;                              
                        
                        return resArray;
        }
        
        /**
         * miIJ: Calculates the Mutual Information between two residue pairs, I & J
         * @param seq1 The collection of identities from the MSA for the 1st residue
         * @param seq2 The collection of identities from the MSA for the 2nd residue
         * @param pX The probabilities for the identities for the 1st residue
         * @param pY The probabilities for the identities for the 2nd residue
         * @param iA The 1st residue position, used for reporting purposes
         * @param jB The 1st residue position, used for reporting purposes
         * @param resFile The file location to which pair data information is written
         * @return The mutual information score, give as a whole number or fraction
         */
        private double miIJ(int[] seq1, int[] seq2, double[] pX, double[] pY, int iA, int jB, boolean saveIt, File resFile)
        {
                int[][] obsXY = new int[pX.length][pY.length];
		
                //Initialize the obsXY array
		for (int i = 0; i < obsXY.length; i++)
		{
			for (int j = 0; j < obsXY[i].length; j++)
			{
				obsXY[i][j] = 0;
			}
		}
		
		double ijMI = 0.0;
		int minLen = Math.min(seq1.length, seq2.length);
		
		double[][] pXY = new double[pX.length][pY.length];
		
		for (int i = 0; i < minLen; i++)
		 {
			int x = seq1[i];
			int y = seq2[i];
			obsXY[x][y] += 1;
		 };
  
		for (int k = 0; k < obsXY.length; k++)
		 {
			for (int m = 0; m  < obsXY[k].length; m++)
			 { 
				double calc1 = (double) obsXY[k][m];
				double calc2 = (double) minLen;
				pXY[k][m] = calc1 / calc2;
			 }
		 }
  
		for (int x = 0; x < pX.length; x++)
		 {
			for (int y = 0; y < pY.length; y++)
			 {
                             
			 double xyMI = 0.0;
			 double dpX = (double) pX[x];
			 double dpY = (double) pY[y];
			 double pXpY = dpX * dpY;
			 
			 if (pXpY > 0.0 && pXY[x][y] > 0.0)
			  {
				xyMI = pXY[x][y] * Math.log(pXY[x][y] / pXpY);
                                
			  }
			 else 
			  {
				xyMI = 0.0;
			  }    
                                
                         if (intraMol && x == y)
                         {
                             xyMI = 0.0;
                         }

                         ijMI += xyMI;
			
                        }
		}

                
		String reporter = Double.toString(ijMI);
                String res1 = Integer.toString(iA + 1) + "\t" + Integer.toString(jB + 1) + "\t" +  String.format("%.10f", ijMI);
                if (pairData)
                {
                    res1 += "\t"; 
                    for (int q = 0; q < obsXY.length; q++)
                    {
                        for (int s = 0; s < obsXY[q].length; s++)
                        {
                            res1 += obsXY[q][s] + "\t";
                        }
                    }
                }
                
                if(saveIt)
                {
                writeFileN(res1, resFile);
                }
                
		return ijMI;
		
	}
	
	private static double[] dumdoub1Darr (double[] arr){
				for (int i = 0; i < arr.length; i++)
				   {
						arr[i] = 0.0;
				   }
			return arr;
	}

	
	private static double[][] dumdoub2Darr (double[][] arr){
				  for (int i = 0; i < arr.length; i++)
				   {
						for (int j = 0; j < arr[i].length; j++)
						 {
						  arr[i][j] = 0.0;
						 }
				   }
			return arr;
	}
	
        /**
         * A method for taking an array and sorting, from highest value to lowest, and returning a sorted matrix
         * @param arr The array to be sorted
         * @return A multidimensional array, with the first dimension containing the scores, i-th positions, and j-th positions, in that order as doubles
         */
	private double[][] sortArr (double[][] arr){
	
                int totalArrSize = 0;
		                
                //Find Min
                double arrMin = 10000;
                
                for (int i = 0; i < arr.length; i++)
                {
                    for (int j = 0; j < arr[i].length; j++)
                    {
                        totalArrSize++;
                        arrMin = Math.min(arrMin, arr[i][j]);
                    }
                }
 
                String[] arrS = new String[totalArrSize];
                
                // Boost by minimum, to ensure all values are positive
                
                for (int i = 0; i < arr.length; i++)
                {
                    for (int j = 0; j < arr[i].length; j++)
                    {
                           arr[i][j] -= arrMin;
                    }
                }                
                
		int z = 0;
		
		for (int i = 0; i < arr.length; i++)
		{
			for (int j = 0; j < arr[i].length; j++)
			{
				
				arrS[z] = String.format("%.10f", arr[i][j]) + "\t" + (i + 1) + "\t" + (j + 1);
				// rep = Double.toString(arr[i][j]);
                                // arrS[z] = rep + "\t" + i + "\t" + j;
                                z++;
			}
		}
		
		Arrays.sort(arrS);
                
                double[][] res = new double[3][arrS.length];           
                
                for (int i = 0; i < arrS.length; i++)
                {
                    String[] arrStemp = arrS[i].split("\t");
                    
                    res[0][i] = Double.parseDouble(arrStemp[0]) + arrMin;
                    res[1][i] = Double.parseDouble(arrStemp[1]);
                    res[2][i] = Double.parseDouble(arrStemp[2]);
                    
                }
		
                return res;
	}

        /**
         * A simple method for sorting and reporting an Array, used to eliminate the need to return large array values just to record them
         * @param arr The array to be sorted
         * @param resFile The file to which the results will be written
         */
	private void sortReport (double[][] arr, File resFile){
	
                double[][] sorted = sortArr(arr);
                
                int totalArrSize = sorted[0].length;
                   
                
                for (int i = totalArrSize - 1; i > 0; i--)
                {
                    String temp = String.format("%.10f", sorted[0][i]) + "\t" + String.format("%.0f", sorted[1][i]) + "\t" + String.format("%.0f", sorted[2][i]);
                    writeFileN(temp, resFile);
                }

		
	}        
        
        /**
         * A method for drawing a filled square
         * This method will take X and Y starting positions and draw a square of a specified size and color
         * @param img The image object to which the square will be drawn
         * @param size The size of the square
         * @param xStart The starting X coordinate
         * @param yStart The starting Y coordinate
         * @param color The square color, as a RBG int
         */
        private static void fillSq (BufferedImage img, int size, int xStart, int yStart, int color)
        {
            size += 1;
            int xEnd = xStart + size;
            int yEnd = yStart + size;
            
               for (int i = xStart; i < xEnd; i++)
                   {
                     for (int j = yStart; j < yEnd; j++)
                          {
                              img.setRGB(i, j, color);
                              
                          }
                   }  
        }

        /**
         * Method for drawing a 2D heatmap from array values
         * @param arr The 2D double array for which the values and size will be used for the heatmap
         * @param fileName The name for the destination file for the heatmap
         * @param scale  The scale for the heatmap, in terms of square pixels per data point.
         */
        private static void drawHeatMap (double[][] arr, String fileName, int scale, String colorChoice)
        {
            
            int totXsize = arr.length;
            
            int totYsize = 0;
            
            double arrMax = 0;
            double arrMin = 10000;
            
            // Find max, min score value, max Y size
            for (int i = 0; i < arr.length; i++)
            {
                totYsize = Math.max(arr[i].length, totYsize);
                for (int j = 0; j < arr[i].length; j++)
                {
                    arrMax = Math.max(arrMax, arr[i][j]);
                    arrMin = Math.min(arrMin, arr[i][j]);
                }
            }
            
            arrMax -= arrMin;
            
            double normFactor = 255 / arrMax;
            
            int[][] normScore = new int[totXsize][totYsize];
            
            // Normalize score
            for (int i = 0; i < arr.length; i++)
            {
                for (int j = 0; j < arr[i].length; j++)
                {
                   double tempScore = arr[i][j];
                   tempScore -= arrMin;
                   tempScore *= normFactor;
                   
                   normScore[i][j] = 255;
                   normScore[i][j] -= (int) tempScore;
                }
            }            
            
            totXsize++;
            totXsize *= scale;
            
            totYsize++;
            totYsize *= scale;
            
            
            // Normalize / categorize sccores
            
            BufferedImage img = new BufferedImage(totXsize, totYsize, BufferedImage.TYPE_INT_RGB);  
            
            int posX = 0;
            
            int posY = 0;

            int red = 65536;
            int green = 256;
            int blue = 1;
            
            for (int i = 0; i < arr.length; i++)
            {
                for (int j = 0; j < arr[i].length; j++)
                {
                    
                // default is gray
                int color = (red * normScore[i][j]) + (green * normScore[i][j]) + (blue * normScore[i][j]);
            
                switch (colorChoice)
                 {
                    case "red":
                    color = (red * 255) + (green * normScore[i][j]) + (blue * normScore[i][j]);
                    break;
                
                    case "green":
                    color = (red * normScore[i][j]) + (green * 255) + (blue * normScore[i][j]);
                    break;
                
                    case "blue":
                    color = (red * normScore[i][j]) + (green * normScore[i][j]) + (blue * 255);
                    break;
                
                    case "gray":
                    color = (red * normScore[i][j]) + (green * normScore[i][j]) + (blue * normScore[i][j]);
                    break;
                    
                    }
                    
                    fillSq(img, scale, posX, posY, color);
                    posY += scale;
                }
                posX += scale;
                posY = 0;
            }
                                  
            try {
                 ImageIO.write(img, "PNG", new File(fileName));
                 } catch (IOException ex) {
                 Logger.getLogger(MutualInfoBio.class.getName()).log(Level.SEVERE, null, ex);
                 }            
        }
        
        /**
         * Method for generating an array that represents a random, shuffled order of another array (that can be used to shuffle several arrays in the same way)
         * @param size The size of the array that needs a shuffle order
         * @return The array with the order listed
         */
        private int[] shuffleRandom (int size)
        {
            int[] resArr = new int[size];
            double[] randomArr1 = new double[size];
            double[] randomArr2 = new double[size];
            
            for (int i = 0; i < randomArr1.length; i++)
            {
                double randomNum = Math.random();
                
                // Check to make sure random number isn't a duplicate of value already present
                for (int j = 0; j < randomArr1.length; j++)
                {
                    if (randomNum == randomArr1[j])
                    {
                        randomNum = Math.random();
                        j = 0;
                    }
                }
                
                randomArr1[i] = randomNum;
                randomArr2[i] = randomNum;
                
            };
            
            Arrays.sort(randomArr2);
            
            for (int i = 0; i < resArr.length; i++)
            {
                for (int j = 0; j < randomArr1.length; j++)
                {
                    if (randomArr1[j] == randomArr2[i])
                    {
                        resArr[i] = j;
                    }
                }
            }
            
            return resArr;
        }
        
        /**
         * Method is designed to apply the APC correction factor on a set of mutual information data
         * @param arrayMI The mutual information array to be corrected
         * @param sumMIi An array with the sum of each row
         * @param sumMIj An array with the sum of each column
         * @param sumMIij The overall sum of the mutual information double array
         * @param dirName The name of the directory to place the files generated
         * @param size The size of the arrayMI
         * @param pixelSize The size of the pixels in the heatmap
         * @param processCount The iteration of this process, for use in naming files
         */
        private void apcProcess (double[][][] arrayMI, String dirName, int size1, int size2, int pixelSize)
        {
            
                        double[] sumMIi = arrayMI[1][0];
                        double[] sumMIj = arrayMI[2][0];
                        double sumMIij = arrayMI[3][0][0];     
                        
                        /**
                         * Recipient file for Raw APC scores, comparing collection B and A
                         */
                        File fRawAPC = new File(dirName + "\\APC_raw" + ".txt");

                        /**
                         * Recipient file for sorted APC scores, comparing collection A against itself
                         */
                        File fSortAPC = new File(dirName + "\\APC_sorted" + ".txt");
                        
                        
                            // Calculate raw APC values and Report as String
                            String res1 = "";
                                  
                            double[][] arrayAPCx = new double[size1][size2];
                                  
                            // Fill with dummy values
                            arrayAPCx = dumdoub2Darr(arrayAPCx);
                                  
                            for (int i = 0; i < arrayMI[0].length; i++)
                                    {
					  for (int j = 0; j < arrayMI[0][i].length; j++)
					  {
						
						// Make APC Adjustments
						double adjAPC = (sumMIi[i] * sumMIj[j]) / sumMIij;
                                                if(sumMIij == 0.0)
                                                {
                                                    adjAPC = 0.0;
                                                }
						arrayAPCx[i][j] = arrayMI[0][i][j] - adjAPC;
                                                
						res1 = Integer.toString(i) + "\t" + Integer.toString(j) + "\t" +  String.format("%.10f", arrayAPCx[i][j]);
						writeFileN(res1, fRawAPC);
					  }
                                    }				   
								  
				  sortReport(arrayAPCx, fSortAPC);

                                  drawHeatMap(arrayAPCx, dirName + "\\heatmap" + ".PNG", pixelSize, "blue");  
        }

        private void rcwProcess(double[][][] arrayMI, String dirName, int size1, int size2, int pixelSize)
        {
                        
                        double[] sumMIi = arrayMI[1][0];
                        double[] sumMIj = arrayMI[2][0];
                        double sumMIij = arrayMI[3][0][0];    
                        
                        /**
                         * Recipient file for Raw APC scores, comparing collection B and A
                         */
                        File fRawRCW = new File(dirName + "\\RCW_raw" + ".txt");

                        /**
                         * Recipient file for sorted APC scores, comparing collection A against itself
                         */
                        File fSortRCW = new File(dirName + "\\RCW_sorted" + ".txt");
            
                        double[][] arrayRCW = new double[size1][size2];				                             
                            
                            // Fill with dummy values
                            arrayRCW = dumdoub2Darr(arrayRCW);
                                  
                            // Calculate raw RCW values and Report as String
                            String res1 = "";
                            for (int i = 0; i < arrayMI[0].length; i++)
                                {
                                    for (int j = 0; j < arrayMI[0][i].length; j++)
					  {
						
						// Make RCW Adjustments
						double adjRCW = (sumMIi[i] + sumMIj[j]) / 2.0;
                                                if(adjRCW != 0.0)
                                                {
						arrayRCW[i][j] = arrayMI[0][i][j] / adjRCW;
                                                }
                                                else
                                                {
                                                    arrayRCW[i][j] = 0.0;
                                                }
						res1 = Integer.toString(i) + "\t" + Integer.toString(j) + "\t" +  String.format("%.10f", arrayRCW[i][j]);
						writeFileN(res1, fRawRCW);
					  }
                                }		

			  
                            sortReport(arrayRCW, fSortRCW);
                            
                            drawHeatMap(arrayRCW, dirName + "\\heatmap" + ".PNG", pixelSize, "red");
        }
        
        /*
        A method which takes each MI score and multiples it by the conservation frequency, to give a score that is reflective of the level of conservation versus the strenght of co-evolution
         */
        private void casProcess(double[][][] arrayMI, String dirName, int pixelSize, int[][][] seq1, int[][][] seq2)
        {
                        
                        // BEGIN BLOCK FOR MI METHOD
                        double[][] pX = null;
                        double[][] pY = null;
                 
                        int size1 = seq1[0].length;
                        int size2 = seq2[0].length;
                        
                        pX = oFreq(seq1);
                        pY = oFreq(seq2); 
                        
                        /**
                         * Recipient file for Raw APC scores, comparing collection B and A
                         */
                        File fRawCAS = new File(dirName + "\\CAS_noMI_raw" + ".txt");
                        
                        /**
                         * Recipient file for sorted APC scores, comparing collection A against itself
                         */
                        File fSortCAS = new File(dirName + "\\CAS_noMI_sorted" + ".txt");
                        

                        /**
                         * Recipient file for Raw APC scores, comparing collection B and A
                         */
                        File fRawCAS2 = new File(dirName + "\\CAS_raw" + ".txt");

                        /**
                         * Recipient file for sorted APC scores, comparing collection A against itself
                         */
                        File fSortCAS2 = new File(dirName + "\\CAS_sorted" + ".txt");
            
                        double[][] arrayCAS = new double[size1][size2];				                             
                            
                            // Fill with dummy values
                            arrayCAS = dumdoub2Darr(arrayCAS);

                                double[] maxPX = new double[pX.length];
                                double[] maxPY = new double [pY.length];
                                
                                dumdoub1Darr(maxPX);
                                dumdoub1Darr(maxPY);
                                
                                for (int i = 0; i < pX.length; i++)
                                {
                                    for (int k = 0; k < pX[i].length; k++)
                                    {
                                        maxPX[i] = Math.max(maxPX[i], pX[i][k]);
                                    }                                                        
                                }

                                for (int i = 0; i < pY.length; i++)
                                {
                                    for (int k = 0; k < pY[i].length; k++)
                                    {
                                        maxPY[i] = Math.max(maxPY[i], pY[i][k]);
                                    }                                                        
                                }
                                
                                double[][] adjCAS = new double[pX.length][pY.length];

                                dumdoub2Darr(adjCAS);

                            // Calculate raw CAS values and report as a string
                            String res1 = "";
                            for (int i = 0; i < arrayMI[0].length; i++)
                                {
                                    for (int j = 0; j < arrayMI[0][i].length; j++)
					  {     
						// Make RCW Adjustments
						adjCAS[i][j] = maxPX[i] * maxPY[j];
                                                
                                                res1 = Integer.toString(i) + "\t" + Integer.toString(j) + "\t" +  String.format("%.10f", adjCAS[i][j]);
						writeFileN(res1, fRawCAS);
					  }
                                }		

			  
                            sortReport(adjCAS, fSortCAS);
                            
                            drawHeatMap(adjCAS, dirName + "\\heatmap_noMI" + ".PNG", pixelSize, "red");

                            // Calculate raw CAS values and report as a string
                            res1 = "";
                            for (int i = 0; i < arrayMI[0].length; i++)
                                {
                                    for (int j = 0; j < arrayMI[0][i].length; j++)
					  {     
                                                arrayCAS[i][j] = arrayMI[0][i][j] * adjCAS[i][j];
                                                
                                                res1 = Integer.toString(i) + "\t" + Integer.toString(j) + "\t" +  String.format("%.10f", arrayCAS[i][j]);
						writeFileN(res1, fRawCAS2);
					  }
                                }		

			  
                            sortReport(arrayCAS, fSortCAS2);
                            
                            drawHeatMap(arrayCAS, dirName + "\\heatmap" + ".PNG", pixelSize, "red");
                            
        }
        
        /**
         * Subtract one 3D double array by another.
         * @param arr1 The 1st array
         * @param arr2 The 2nd array
         * @return The transformed array
         */
        private static double[][][] subtract3Darr(double[][][] arr1, double[][][] arr2)
        {
            for (int i = 0; i < arr1.length; i++)
            {
                for (int j = 0; j < arr1[i].length; j++)
                {
                    for (int k = 0; k < arr1[i][j].length; k++)
                    {
                        arr1[i][j][k] -= arr2[i][j][k];                       
                    }
                }
            }
            
            return arr1;
        }
        
        /**
         * Method for taking two 3D double arrays and averaging them together, using a weight factor to allow iteration in code.
         * @param arr1 The first array
         * @param arr2 The second array
         * @param weight The weighting to apply to the first array
         * @return The averaged values
         */
        private static double[][][] average3Darr(double[][][] arr1, double[][][] arr2, int weight)
        {
            for (int i = 0; i < arr1.length; i++)
            {                
                for (int j = 0; j < arr1[i].length; j++)
                {
                    for (int k = 0; k < arr1[i][j].length; k++)
                    {   
                        double weightD = (double)weight;
                        arr1[i][j][k] = ((arr1[i][j][k] * weightD) + arr2[i][j][k]) / (weightD + 1);
                    }
                }
            }
            
            return arr1;
        }        

        private static double[][][] clone3Darr(double[][][] arr)
        {
            double[][][] arr2 = new double[arr.length][][];
            for (int i = 0; i < arr.length; i++)
            {
                double[][] arr2i = new double[arr[i].length][];
                for (int j = 0; j < arr[i].length; j++)
                {
                    double[] arr2j = new double[arr[i][j].length];
                    for (int k = 0; k < arr[i][j].length; k++)
                    {
                        double arr2k = arr[i][j][k];
                        arr2j[k] = arr2k;
                    }
                    arr2i[j] = arr2j;
                }
                arr2[i] = arr2i;
            }
            return arr2;
        }
        
        private void write2Darr(double[][] arr, File file1)
        {
                        for (int i = 0; i < arr.length; i++)
                        {
                            for (int j = 0; j < arr[i].length; j++)
                            {                                 
                                    String zRep = Integer.toString(i + 1) + "\t" + Integer.toString(j + 1) + "\t" +  String.format("%.10f", arr[i][j]);
                                    writeFileN(zRep, file1);
                            }
                        }
        }
}