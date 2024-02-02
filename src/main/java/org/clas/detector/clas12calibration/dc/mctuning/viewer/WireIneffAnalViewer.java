/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.viewer;

import java.awt.BorderLayout;
import java.awt.Dialog;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.TreeMap;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import org.clas.detector.clas12calibration.dc.mctuning.analysis.configButtonPanel;
import org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.WireIneffAnal;
import org.jlab.clas.swimtools.MagFieldsEngine;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.decode.CodaEventDecoder;
import org.jlab.detector.decode.DetectorEventDecoder;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.detector.view.DetectorListener;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;
import org.jlab.rec.dc.Constants;
/**
 *
 * @author ziegler
 */

    

public class WireIneffAnalViewer implements IDataEventListener, DetectorListener, ActionListener, ChangeListener {
    
    List<DetectorPane2D> AnalysisPanels 	= new ArrayList<DetectorPane2D>();
    JTabbedPane tabbedpane           		= null;
    JPanel mainPanel 				= null;
    JSplitPane splitpane                        = null;
    JMenuBar menuBar                            = null;
    DataSourceProcessorPane processorPane 	= null;
    EmbeddedCanvasTabbed CLAS12Canvas           = null;
    JFrame  innerConfigFrame = new JFrame("Configure calibration settings");
    JDialog configFrame = new JDialog(innerConfigFrame, "Configure calibration settings");
    JTabbedPane configPane = new JTabbedPane();
    
    CodaEventDecoder               decoder = new CodaEventDecoder();
    DetectorEventDecoder   detectorDecoder = new DetectorEventDecoder();
    public static DCGeant4Factory dcDetector;   
    public static ConstantsManager ccdb = new ConstantsManager();
    TreeMap<String, List<H2F>>  histos = new TreeMap<String,List<H2F>>();
    
    private int canvasUpdateTime = 10000;
    private int analysisUpdateTime = 10000;
    private int runNumber  = 0;
    private String Dir = System.getProperty("user.dir");
    
    
    private JLabel[] superlayer = {new JLabel("", JLabel.CENTER),new JLabel("", JLabel.CENTER),new JLabel("", JLabel.CENTER),new JLabel("", JLabel.CENTER),new JLabel("", JLabel.CENTER),new JLabel("", JLabel.CENTER)};
    
     // detector monitors
    AnalysisMonitor[] monitors ; 
    private WireIneffAnal _wireIneffAnal ;     
    public WireIneffAnalViewer() throws FileNotFoundException {  
        this._wireIneffAnal = new WireIneffAnal("Wire Inefficiency Analysis",ccdb);
        this.monitors = new AnalysisMonitor[]{this._wireIneffAnal};		
	// create menu bar
        menuBar = new JMenuBar();
        JMenuItem menuItem;
        JMenu file = new JMenu("File");
        file.getAccessibleContext().setAccessibleDescription("File options");
        menuItem = new JMenuItem("Open histograms file...");
        menuItem.getAccessibleContext().setAccessibleDescription("Open histograms file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuItem = new JMenuItem("Print histograms to file...");
        menuItem.getAccessibleContext().setAccessibleDescription("Print histograms to file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuItem = new JMenuItem("Save histograms to file...");
        menuItem.getAccessibleContext().setAccessibleDescription("Save histograms to file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuBar.add(file);
        JMenu settings = new JMenu("Settings");
        settings.getAccessibleContext().setAccessibleDescription("Choose monitoring parameters");
        menuItem = new JMenuItem("Set GUI update interval...");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GUI update interval");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        
        
        menuItem = new JMenuItem("Set Style to default");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GROOT style to default");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuItem = new JMenuItem("Set Style for performance plots");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GROOT style for performance plots");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuBar.add(settings);
        
        JMenu fits = new JMenu("Fits");
        fits.getAccessibleContext().setAccessibleDescription("Choose parameters");
        menuItem = new JMenuItem("Refit");
        menuItem.getAccessibleContext().setAccessibleDescription("...");
        menuItem.addActionListener(this);
        fits.add(menuItem);
        menuBar.add(fits);
           
        // create main panel
        mainPanel = new JPanel();	
	mainPanel.setLayout(new BorderLayout());
        splitpane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
      	tabbedpane 	= new JTabbedPane();

        processorPane = new DataSourceProcessorPane();
        processorPane.setUpdateRate(analysisUpdateTime);

        mainPanel.add(splitpane);
        splitpane.setTopComponent(tabbedpane);
        
        splitpane.setBottomComponent(monitors[0].getCcview());
        splitpane.setDividerLocation(0.75);
        splitpane.setResizeWeight(0.75);
        mainPanel.add(processorPane,BorderLayout.PAGE_END);
        
    
        GStyle.getAxisAttributesX().setTitleFontSize(24);
        GStyle.getAxisAttributesX().setLabelFontSize(18);
        GStyle.getAxisAttributesY().setTitleFontSize(24);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        
 
        tabbedpane.addChangeListener(this);
       
        for(int k =0; k<this.monitors.length; k++) {
                this.tabbedpane.add(this.monitors[k].getAnalysisPanel(), this.monitors[k].getAnalysisName());
        	this.monitors[k].getAnalysisView().getView().addDetectorListener(this);
        }
        this.processorPane.addEventListener(this);
        
        this.setCanvasUpdate(canvasUpdateTime);
        // init constants manager
        ccdb.init(Arrays.asList(new String[]{
            "/geometry/dc/superlayer/wpdist",
            "/calibration/dc/time_to_distance/t2d_pressure", 
            "/hall/weather/pressure",
            "/calibration/dc/time_to_distance/ref_pressure",
            "/calibration/dc/time_jitter"}));
        ccdb.setVariation("default");
        ConstantProvider provider = GeometryFactory.getConstants(DetectorType.DC, 11, "default");
        for(int l=0; l<6; l++) {
            Constants.getInstance().wpdist[l] = provider.getDouble("/geometry/dc/superlayer/wpdist", l);
        }
        Constants.getInstance().setT2D(1);
        dcDetector = new DCGeant4Factory(provider, DCGeant4Factory.MINISTAGGERON, true);

        // set directory to local
        this.Dir = System.getProperty("user.dir");
        System.out.println("Work directory set to " + this.Dir);
        
        enf.init();
        //tm.init(); 
    }
    public static MagFieldsEngine enf = new MagFieldsEngine();
        
    //public static LayerEfficiencyAnalyzer tm = new LayerEfficiencyAnalyzer();
        
    public void actionPerformed(ActionEvent e) {
        System.out.println(e.getActionCommand());
        if(e.getActionCommand()=="Set GUI update interval...") {
            this.chooseUpdateInterval();
        }
        if(e.getActionCommand()=="Set Style to default") {
            this.setStyle(0);
        }
        if(e.getActionCommand()=="Set Style for performance plots") {
            this.setStyle(1);
        }
        if(e.getActionCommand()=="Open histograms file...") {
            String fileName = null;
            JFileChooser fc = new JFileChooser();
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            File workingDirectory = new File(this.Dir + "/cal-histos");
            fc.setCurrentDirectory(workingDirectory);
            int option = fc.showOpenDialog(null);
            if (option == JFileChooser.APPROVE_OPTION) {
                fileName = fc.getSelectedFile().getAbsolutePath();            
            }
            if(fileName != null) this.loadHistosFromFile(fileName);
        }        
        if(e.getActionCommand()=="Print histograms to file...") {
            this.printHistosToFile();
        }
        if(e.getActionCommand()=="Refit") {
            this._wireIneffAnal.Refit();
        }
        if(e.getActionCommand()=="Save histograms to file...") {
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "clas12rec_run_" + this.runNumber + "_" + df.format(new Date()) + ".hipo";
            JFileChooser fc = new JFileChooser();
            File workingDirectory = new File(this.Dir + "/kpp-histos");
            fc.setCurrentDirectory(workingDirectory);
            File file = new File(fileName);
            fc.setSelectedFile(file);
            int returnValue = fc.showSaveDialog(null);
            if (returnValue == JFileChooser.APPROVE_OPTION) {
               fileName = fc.getSelectedFile().getAbsolutePath();            
            }
            this.saveHistosToFile(fileName);
        }
                    if (e.getActionCommand().compareTo("Next")==0) {
                    int currentTab = configPane.getSelectedIndex();
                    for (int i=currentTab+1; i<configPane.getTabCount(); i++) {
                        if (configPane.isEnabledAt(i)) {
                                configPane.setSelectedIndex(i);
                                break;
                        }
                    }
		}
		if (e.getActionCommand().compareTo("Back")==0) {
			int currentTab = configPane.getSelectedIndex();
			for (int i=currentTab-1; i>=0; i--) {
				if (configPane.isEnabledAt(i)) {
					configPane.setSelectedIndex(i);
					break;
				}
			}        
		}
		if (e.getActionCommand().compareTo("Cancel")==0) {
			System.exit(0);
		}

		if (e.getActionCommand().compareTo("Finish")==0) {
			configFrame.setVisible(false);
                        //read field in config panel

                }
    }

    public void chooseUpdateInterval() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "GUI update interval (ms)",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "1000");
        if(s!=null){
            int time = 1000;
            try { 
                time= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            if(time>0) {
                this.setCanvasUpdate(time);
            }
            else {
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
        }
    }
        
    
    public JPanel  getPanel(){
        return mainPanel;
    }

    private int getRunNumber(DataEvent event) {
        int rNum = this.runNumber;
        DataBank bank = event.getBank("RUN::config");
        if(bank!=null) {
            rNum = bank.getInt("run", 0);
        }
        return rNum;
    }
    
    @Override
    public void dataEventAction(DataEvent event) {
    	
	if(event!=null ){
//            event.show();
            if (event.getType() == DataEventType.EVENT_START) {
                this.runNumber = this.getRunNumber(event);
            }
            if(this.runNumber != this.getRunNumber(event)) {
//                this.saveToFile("mon12_histo_run_" + runNumber + ".hipo");
                this.runNumber = this.getRunNumber(event);
//                resetEventListener();
            }
            for(int k=0; k<this.monitors.length; k++) {
                this.monitors[k].dataEventAction(event);
            }      
            if (event.getType() == DataEventType.EVENT_STOP) {
                for(int k=0; k<this.monitors.length; k++) {
                    this.monitors[k].analysis();
                }                
            }
	}
   }

    public void loadHistosFromFile(String fileName) {
        // TXT table summary FILE //
        System.out.println("Opening file: " + fileName);
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        System.out.println(dir.getDirectoryList());
        dir.cd();
        dir.pwd();
        
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].readDataGroup(dir);
        }
    }

    
    public void printHistosToFile() {
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String data = this.Dir + "/clas12rec_run_" + this.runNumber + "_" + df.format(new Date());        
        File theDir = new File(data);
        // if the directory does not exist, create it
        if (!theDir.exists()) {
            boolean result = false;
            try{
                theDir.mkdir();
                result = true;
            } 
            catch(SecurityException se){
                //handle it
            }        
            if(result) {    
            System.out.println("Created directory: " + data);
            }
        }
        String fileName = data + "/clas12_canvas.png";
        System.out.println(fileName);
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].printCanvas(data);
        }
    }

     @Override
    public void processShape(DetectorShape2D shape) {
        System.out.println("SHAPE SELECTED = " + shape.getDescriptor());
    }
    
    @Override
    public void resetEventListener() {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].resetEventListener();
            this.monitors[k].timerUpdate();
        } 
    }
    
   public void saveHistosToFile(String fileName) {
        // TXT table summary FILE //
        TDirectory dir = new TDirectory();
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].writeDataGroup(dir);
        }
        System.out.println("Saving histograms to file " + fileName);
        dir.writeFile(fileName);
    }
            
    public void setCanvasUpdate(int time) {
        System.out.println("Setting " + time + " ms update interval");
        this.canvasUpdateTime = time;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setCanvasUpdate(time);
        }
    }

    public void setStyle(int mode) {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setStyle(mode);
            this.monitors[k].plotHistos();
        }    
    }
    
    public void stateChanged(ChangeEvent e) {
        this.timerUpdate();
    }
    
    @Override
    public void timerUpdate() {
//        System.out.println("Time to update ...");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].timerUpdate();
        }
    }


    
    public void configure() {

        configFrame.setSize(900, 830);
        //configFrame.setSize(1000, 600); // vnc size
        configFrame.setLocationRelativeTo(mainPanel);
        configFrame.setDefaultCloseOperation(configFrame.DO_NOTHING_ON_CLOSE);

        // Which steps    
        JPanel stepOuterPanel = new JPanel(new BorderLayout());
        JPanel stepPanel = new JPanel(new GridBagLayout());
        stepOuterPanel.add(stepPanel, BorderLayout.NORTH);
        GridBagConstraints c = new GridBagConstraints();
        JPanel butPage1 = new configButtonPanel(this, false, "Next");
        stepOuterPanel.add(butPage1, BorderLayout.SOUTH);

        //configPane.add("Select steps", stepOuterPanel);    

        // Previous calibration values
        JPanel confOuterPanel = new JPanel(new BorderLayout());
        Box confPanel = new Box(BoxLayout.Y_AXIS);
        
        

        JPanel butPage2 = new configButtonPanel(this, false, "Next");
        confOuterPanel.add(confPanel, BorderLayout.NORTH);
        confOuterPanel.add(butPage2, BorderLayout.SOUTH);

        //configPane.add("Previous calibration values", confOuterPanel);


	}
    
    
    public static void main(String[] args) throws FileNotFoundException {
        
        JFrame frame = new JFrame("DC Calibration");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        WireIneffAnalViewer viewer = new WireIneffAnalViewer();
        frame.add(viewer.mainPanel);
        frame.setJMenuBar(viewer.menuBar);
        frame.setSize(1400, 800);
        frame.setVisible(true);
        viewer.configFrame.setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
	viewer.configure();

    }
   
}