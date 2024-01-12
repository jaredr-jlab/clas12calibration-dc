/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt0;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.calt2d.SegmentProperty;
import org.clas.detector.clas12calibration.dc.calt2d.Utilities;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.clas.detector.clas12calibration.viewer.T0Viewer;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
/**
 *
 * @author KPAdhikari, ziegler
 */
public class T00Calib extends AnalysisMonitor{
    //public HipoDataSync writer = null;
    //private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    PrintWriter pw = null;
    File outfile = null;
    private int runNumber;
    private String analTabs = "Fully Corrected Time";
    //private Utilities util = new Utilities();
   
    public T00Calib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames(analTabs);
        this.init(false, "T00");
        T00Array = new double[nsec][nsl];
        outfile = new File("Files/ccdbConstantstT00.txt");
        pw = new PrintWriter(outfile);
        pw.printf("#& Sector Superlayer T0Correction T0Error\n");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
       
        if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
            System.out.println(" BANK FOUND........");
        } else {
            System.out.println(" BANK NOT FOUND........");
        }
        //writer = new HipoDataSync(schemaFactory);
        //writer.setCompressionType(2);
        //hipoEvent = (HipoDataEvent) writer.createEvent();
        //writer.open("TestOutPut.hipo");
        //writer.writeEvent(hipoEvent);

        
    }
    int nsl  = 6;
    int nsec = 6;
    
    boolean[][] Fitted = new boolean[nsec][nsl];
    int[] nTdcBins =
    { 50, 50, 50, 50, 50, 50 };
    int[] nTimeBins =
    { 50, 50, 50, 50, 50, 50 };
    double[] tLow =
    { 80.0, 80.0, 80.0, 80.0, 80.0, 80.0 };
    
    public static final double[] tLow4T0Fits  = {-40.0, -40.0, -40.0, -40.0, -40.0, -40.0};
    public static final double[] tHigh4T0Fits  = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0}; 

    public static  double[][] fitMax ;


    //H1F[][][][] h = new H1F[6][6][nSlots7][nCables6];
    private Map<Coordinate, H1F> TDCHis        = new HashMap<Coordinate, H1F>();    
    public  Map<Coordinate, FitLine> TDCFits   = new HashMap<Coordinate, FitLine>();
    public  Map<Coordinate, Double> T0s        = new HashMap<Coordinate, Double>();
    
    @Override
    public void createHistos() {
        //histo max range for the fit
        fitMax = new double[nsec][nsl]; 
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup hgrps = new DataGroup();
        String hNm;
        String hTtl;
        for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                hNm = String.format("timeS%dS%d", i + 1, j + 1);

                TDCHis.put(new Coordinate(i,j), new H1F(hNm, 80, tLow4T0Fits[j], tHigh4T0Fits[j])); 
                                                                                                                                                                                // HBHits
                hTtl = String.format("time (Sec%d SL%d)", i + 1, j + 1);
                TDCHis.get(new Coordinate(i,j)).setTitleX(hTtl);
                TDCHis.get(new Coordinate(i,j)).setLineColor(1);
                TDCFits.put(new Coordinate(i,j), new FitLine());
                hgrps.addDataSet(TDCHis.get(new Coordinate(i, j)), 0);

                T0s.put(new Coordinate(i,j), 0.0);

                Fitted[i][j] = false;
                
            }
                
            this.getDataGroup().add(hgrps, i+1,0,0);
        }

        this.getDataGroup().add(hgrps,0,0,0);
        
        for (int i = 0; i < nsec; i++) {
                this.getCalib().addEntry(i+1,0,0);
        }
        this.getCalib().setName("T00 Table");
        this.getCalib().fireTableDataChanged();
    }
     
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas(analTabs).setGridX(false);
        this.getAnalysisCanvas().getCanvas(analTabs).setGridY(false);
        this.getAnalysisCanvas().getCanvas(analTabs).divide(3, 2);
        this.getAnalysisCanvas().getCanvas(analTabs).update();
        
        
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        this.plotFits();
    }
    public void plotFits() {
        
        //pw.close();
        File file2 = new File("");
        file2 = outfile;
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String fileName = "Files/ccdb_T00Corr_run" + this.runNumber + "time_" 
                + df.format(new Date())  + ".txt";
        file2.renameTo(new File(fileName));
        
        for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                if(this.fitThisHisto(this.TDCHis.get(new Coordinate(i,j)))==true) {
                    this.runFit(i, j);
                    int binmax = this.TDCHis.get(new Coordinate(i,j)).getMaximumBin();
                    fitMax[i][j] = this.TDCHis.get(new Coordinate(i,j)).getDataX(binmax);

                }
            }
        }
        pw.close();
        this.getCalib().fireTableDataChanged();  
        
    }
    
    public int NbRunFit = 0;
    public static double[][] T00Array; 
    public void runFit(int i, int j) {
        
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** "); 
	
        double[] Tminmax = this.getT0(i, j);
        T00Array[i][j] = Tminmax[0];
        //Sector Superlayer Slot Cable T0Correction T0Error
        pw.printf("%d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), 
            Tminmax[0], 
            Tminmax[1]);
        System.out.printf("%d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1),
            Tminmax[0], 
            Tminmax[1]);
        this.updateTable(i, Tminmax[0]);
        Fitted[i][j] = true;
        System.out.println(" FITTED ? "+Fitted[i][j]);
    }
     private void updateTable(int i, double t0) {
       this.getCalib().setDoubleValue(t0, "T00", i+1, 0, 0);
     }
    
    int counter = 0;
    public  HipoDataSource reader = new HipoDataSource();
    

    int count = 0;
    public static int polarity =-1;
    public List<FittedHit> hits = new ArrayList<>();
    Map<Integer, ArrayList<Integer>> segMapTBHits = new HashMap<Integer, ArrayList<Integer>>();
    Map<Integer, SegmentProperty> segPropMap = new HashMap<Integer, SegmentProperty>();
    @Override
    public void processEvent(DataEvent event) {
        
        if (!event.hasBank("RUN::config")) {
            return ;
        }

        if (!event.hasBank("REC::Particle")) {
            return ;
        }
        
        DataBank bank = event.getBank("RUN::config");
        DataBank bank2 = event.getBank("REC::Particle");
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) {
           return ;
        } else {
           count++;
        }
        
        if(count==1) {
            //Constants.getInstance().initialize("DCCAL");
            Driver.init();
            TableLoader.FillT0Tables(newRun, "default");

            TableLoader.Fill(T0Viewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/time2dist"));  
            ReadTT.Load(newRun, "default"); 
            runNumber = newRun; 
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
        // get segment property
        
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        
        for (int j = 0; j < bnkHits.rows(); j++) {
            
            int sec = bnkHits.getInt("sector", j);
            int sl = bnkHits.getInt("superlayer", j);
            
            double time = (double) bnkHits.getFloat("time", j);
            double calibtime = time + bnkHits.getFloat("tBeta",j); //time without tbeta correction
            //tbeta correction function is in Utilities.java named calcDeltaTimeBetaNew(calibtime,distbeta,beta);
            //double beta = bnkHits.getFloat("beta",j);
            //distbeta is from ccdb for superlayer and sector 
            //double distbeta = TableLoader.distbeta[sec-1][sl-1];
           //scale factor on how much of the new beta correction should be applied. values should be from 0 to 1
            //double tbetascalefactor = 0; 
            //calculate tbeta correction
            //double newtbetacorr = tbetascalefactor*util.calcDeltaTimeBetaNewFCN(calibtime, distbeta, beta);            
           //apply tbeta correction
            //double timecorrected = calibtime - newtbetacorr;
            //debug output
           // System.out.println("Sector " + sec + " sl " + sl + " distbeta is " + distbeta);
           // System.out.println("time " + time + " calibtime " + calibtime + " calib-newtbeta " + timecorrected + " tbetafromfile " + bnkHits.getFloat("tBeta",j) + " tbetarecal " + newtbetacorr);
            int trkID = bnkHits.getByte("trkID", j);
            if(trkID!=-1 && 11 == bank2.getInt("pid", 0)) {
                 //int pid = this.readPID(event, trkID);            
            	 //pid 11 is electron, 2212 proton and 211 pion
                //pass if the track is identified as an electron or as a hadron
                 //if(pid==11) {
                	 this.TDCHis.get(new Coordinate(sec-1, sl-1))
                     //.fill(timecorrected);
                     .fill(calibtime);
                 //}
                }
            }
        } 
    
    public void Plot(int i ) {
        for(int j = 0; j<nsl; j++) {
            this.getAnalysisCanvas().getCanvas(analTabs).cd(j);
            this.getAnalysisCanvas().getCanvas(analTabs)
                    .draw(this.TDCHis.get(new Coordinate(i, j)));

            if(Fitted[i][j]==true) {
                this.getAnalysisCanvas().getCanvas(analTabs).cd(j);
                            this.getAnalysisCanvas().getCanvas(analTabs)
                                .draw(this.TDCFits.get(new Coordinate(i, j)), "same");
            }
        }
    }
    
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        
        System.out.println("sector" +str_sector );
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);

       if(group.hasItem(sector,0,0)==true){
           this.Plot(sector-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
   
    }

    private boolean fitThisHisto(H1F h) {
        boolean pass = false;
        int nevent = 0;
        int maxbin = h.getMaximumBin();
        for (int ix =0; ix< maxbin; ix++) {
            double y = h.getBinContent(ix);
            double err = h.getBinError(ix);
            
            if(err>0 && y>0) {
                nevent+=y;
            }
        }
        
        if(nevent > 99)
            pass = true;
        
        return pass;
    }

    private double[] getT0(int i, int j) {
        System.out.println("Getting t0 for i,j = "+i+" "+j);
        H1F h = this.TDCHis.get(new Coordinate(i,j));
        
        double [] T0val = new double[2];

        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+[p0]", 
                h.getDataX(0), h.getDataX(h.getMaximumBin())); 

        gausFunc.setParameter(0, h.getMax());
        gausFunc.setParameter(1, h.getDataX(h.getMaximumBin()));
        gausFunc.setParameter(2, 25);
        gausFunc.setParameter(3, 0);
                
        DataFitter.fit(gausFunc, h, "Q");
        
        F1D erfcFunc = new F1D("erfcFunc", "[amp]*erf(-x, -[mean], [sigma])+[p0]", 
                        h.getDataX(0), h.getDataX(h.getMaximumBin()));
                
        erfcFunc.setParameter(0, gausFunc.getParameter(0));
        erfcFunc.setParameter(1, gausFunc.getParameter(1));
        erfcFunc.setParameter(2, gausFunc.getParameter(2));
        erfcFunc.setParameter(3, gausFunc.getParameter(3));
        
        DataFitter.fit(erfcFunc, h, "Q");

        double T0 = erfcFunc.getParameter(1)-1.5*erfcFunc.getParameter(2);
        double T0Err = Math.sqrt(erfcFunc.parameter(1).error()*erfcFunc.parameter(1).error()+erfcFunc.parameter(2).error()*erfcFunc.parameter(2).error());
        if(Double.isNaN(T0)|| Double.isNaN(T0Err)){
            T0 = 0;
            T0Err = 1.42;
        }
        T0val[1] =T0Err;
        T0val[0] = T0;

        h.setOptStat(0);
        String t = "T00 = "+(float)T0;
        h.setTitle(t);
        T0s.put(new Coordinate(i,j), T0);
        
        return T0val;
    }
    
    /*private int readPID(DataEvent event, int trkId) {
        int pid = 0;
        //fetch the track associated pid from the REC tracking bank
        if (!event.hasBank("REC::Particle") || !event.hasBank("REC::Track"))
            return pid;
        DataBank bank = event.getBank("REC::Track");
        //match the index and return the pid
        int rows = bank.rows();
        for (int i = 0; i < rows; i++) {
            if (bank.getByte("detector", i) == 6 &&
                    bank.getShort("index", i) == trkId - 1) {
                DataBank bank2 = event.getBank("REC::Particle");
                if(bank2.getByte("charge", bank.getShort("pindex", i))!=0) {
                    pid = bank2.getInt("pid", bank.getShort("pindex", i));
                }
            }
        }
        
        return pid;
    } */
}