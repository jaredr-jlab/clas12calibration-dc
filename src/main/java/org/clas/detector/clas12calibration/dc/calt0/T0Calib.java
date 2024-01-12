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

//import javax.sound.sampled.DataLine;

import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import static org.clas.detector.clas12calibration.dc.calt0.ReadTT.t0min;
import org.clas.detector.clas12calibration.dc.calt2d.SegmentProperty;
import org.clas.detector.clas12calibration.dc.calt2d.Utilities;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.DataVector;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.jlab.groot.data.DataLine;
/**
 *
 * @author KPAdhikari, ziegler
 */
public class T0Calib extends AnalysisMonitor{
    
    //public HipoDataSync writer = null;
    //private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    PrintWriter pw = null;
    PrintWriter pw2 = null;
    PrintWriter pw3 = null;
    File outfile = null;
    File outfile2 = null;
    File outfile3 = null;
    private int runNumber;
    private String analTabs = "Corrected Time + T0";
    private Utilities util = new Utilities();
    public T0Calib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames(analTabs);
        this.init(false, "T0");
        outfile = new File("Files/ccdbConstantst0.txt");
        outfile2 = new File("Files/ccdbConstantst00.txt");
        outfile3 = new File("Files/ccdbConstantst00todb.txt");
        pw = new PrintWriter(outfile);
        pw.printf("#& Sector Superlayer Slot Cable T0Correction T0Error\n");
        pw2 = new PrintWriter(outfile2);
        pw2.printf("#& Sector Superlayer Slot Cable T0Correction T0Error\n");
        pw3 = new PrintWriter(outfile3);
        pw3.printf("#& Sector Superlayer Slot Cable T0Correction T0Error\n");
        
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
    int nCrates = 18;// Goes from 41 to 58 (one per chamber)
    int nSlots = 20; // Total slots in each crate (only 14 used)
    int nChannels = 96;// Total channels per Slot (one channel per wire)
    int nLayers0to35 = 36;// Layers in each sector (0th is closest to CLAS
    int nCables = 84;
    int nCables6 = 6; // # of Cables per DCRB or STB.
    int nSlots7 = 7; // # of STBs or occupied DCRB slots per SL.
    boolean[][][][] Fitted = new boolean[nsec][nsl][nSlots7][nCables6];
    int[] nTdcBins =
    { 50, 50, 50, 50, 50, 50 };
    int[] nTimeBins =
    { 50, 50, 50, 50, 50, 50 };
    double[] tLow =
    { 80.0, 80.0, 80.0, 80.0, 80.0, 80.0 };
    
    public static final double[] tLow4T0Fits  = {-140.0, -140.0, -140.0, -140.0, -140.0, -140.0};
    public static final double[] tHigh4T0Fits  = {380.0, 380.0, 680.0, 780.0, 1080.0, 1080.0}; 
    public static final int[] tBinnum4T0Fits = {(380+140)/2, (380+140)/2, (680+140)/2, (780+140)/2, (1080+140)/3, (1080+140)/3};
    
    public static  double[][][][] fitMax ;


    //H1F[][][][] h = new H1F[6][6][nSlots7][nCables6];
    private Map<Coordinate, H1F> TDCHis        = new HashMap<Coordinate, H1F>();    
    public  Map<Coordinate, FitLine> TDCFits   = new HashMap<Coordinate, FitLine>();
    public  Map<Coordinate, F1D> ERFFits   = new HashMap<Coordinate, F1D>();
    public  Map<Coordinate, F1D> DERFits   = new HashMap<Coordinate, F1D>();
    public  Map<Coordinate, H1F> DerHis        = new HashMap<Coordinate, H1F>();
    public  Map<Coordinate, Double> T0s        = new HashMap<Coordinate, Double>();
    public  Map<Coordinate, Double> T0CCDBs        = new HashMap<Coordinate, Double>();
    public  Map<Coordinate, Double> T0ERFs        = new HashMap<Coordinate, Double>();
    public  Map<Coordinate, Double> T0Errs        = new HashMap<Coordinate, Double>();
    public  Map<Coordinate, Double> T0CCDBErrs        = new HashMap<Coordinate, Double>();
    public  Map<Coordinate, Double> T0ERFErrs        = new HashMap<Coordinate, Double>();
    
    @Override
    public void createHistos() {
        //histo max range for the fit
        fitMax = new double[nsec][nsl][nSlots][nCables]; 
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup hgrps = new DataGroup(6,7);
        String hNm;
        String hTtl;
        int ijk=-1;
        for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                for (int k = 0; k < nSlots7; k++)
                {
                    for (int l = 0; l < nCables6; l++)
                    {
                        hNm = String.format("timeS%dS%dS%dCbl%d", i + 1, j + 1, k + 1, l + 1);                        
                        TDCHis.put(new Coordinate(i,j,k, l), new H1F(hNm, tBinnum4T0Fits[j], tLow4T0Fits[j], tHigh4T0Fits[j])); 

                        ERFFits.put(new Coordinate(i,j,k, l), new F1D("erfcFunc","0.5*[amp]*erf(-x, -[mean], [sigma])+[p0]", TDCHis.get(new Coordinate(i,j,k, l)).getDataX(0), TDCHis.get(new Coordinate(i,j,k, l)).getDataX(TDCHis.get(new Coordinate(i,j,k, l)).getMaximumBin())));

                        DerHis.put(new Coordinate(i,j,k, l), new H1F(hNm, tBinnum4T0Fits[j], tLow4T0Fits[j], tHigh4T0Fits[j]));

                        DERFits.put(new Coordinate(i,j,k, l), new F1D("derFunc", "[amp]*gaus(x,[mean],[sigma])", TDCHis.get(new Coordinate(i,j,k, l)).getDataX(0), TDCHis.get(new Coordinate(i,j,k, l)).getDataX(TDCHis.get(new Coordinate(i,j,k, l)).getMaximumBin())));
                        // HBHits
                        hTtl = String.format("time (Sec%d SL%d Slot%d Cable%d)", i + 1, j + 1, k + 1, l + 1);
                        TDCHis.get(new Coordinate(i,j,k, l)).setTitleX(hTtl);
                        TDCHis.get(new Coordinate(i,j,k, l)).setLineColor(1);
                        DerHis.get(new Coordinate(i,j,k, l)).setLineColor(8);
                        DERFits.get(new Coordinate(i,j,k, l)).setLineColor(8);
                        TDCFits.put(new Coordinate(i,j,k, l), new FitLine());
                        hgrps.addDataSet(TDCHis.get(new Coordinate(i, j, k, l)), 0);
                        
                        T0s.put(new Coordinate(i,j,k, l), ReadTT.T0[i][j][k][l]);
                        Fitted[i][j][k][l] = false;
                    }
                    this.getDataGroup().add(hgrps, i+1, j+1, k+1);
                }
                
            }
        }

        this.getDataGroup().add(hgrps,0,0,0);
        for (int i = 0; i < nsec; i++) {
            for (int j = 0; j < nsl; j++) {
                for (int k = 0; k < nSlots7; k++) {
                    this.getCalib().addEntry(i+1,j+1,k+1);
                }
            }
        }
        this.getCalib().setName("T0 Table (slot 7)");
        this.getCalib().fireTableDataChanged();
    }
     
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas(analTabs).setGridX(false);
        this.getAnalysisCanvas().getCanvas(analTabs).setGridY(false);
        this.getAnalysisCanvas().getCanvas(analTabs).divide(nCables6/2, 2);
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

        File directory = new File("./T0");
        if (! directory.exists()){
            directory.mkdirs();

        }
        
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("./T0/T0s_output.txt", false));
            writer.write("T0 values generated at " + df.format(new Date()) + "\n");
            writer.append("DerT0\tDerT0Err\tCCDBT0\tCCDBT0Err\tErfT0\tErfT0Err\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        //pw.close();
        File file2 = new File("");
        file2 = outfile;
        File file20 = new File("");
        file20 = outfile2;
        File file20db = new File("");
        file20db = outfile3;
        String fileName = "Files/ccdb_T0Corr_run" + this.runNumber + "time_" 
                + df.format(new Date())  + ".txt";
        file2.renameTo(new File(fileName));
        String fileName2 = "Files/ccdb_T0CorrT00Sub_run" + this.runNumber + "time_" 
                + df.format(new Date())  + ".txt";
        file20.renameTo(new File(fileName2));
        String fileName3 = "Files/ccdb_T0CorrT00SubT0DB_run" + this.runNumber + "time_" 
                + df.format(new Date())  + ".txt";
        file20db.renameTo(new File(fileName3));
        for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                for (int k = 0; k < nSlots7; k++)
                {
                    for (int l = 0; l < nCables6; l++)
                    { 
                        if(this.fitThisHisto(this.TDCHis.get(new Coordinate(i,j,k,l)))==true) {
                            this.runFit(i, j, k, l);
                            int binmax = this.TDCHis.get(new Coordinate(i,j,k,l)).getMaximumBin();
                            fitMax[i][j][k][l] = this.TDCHis.get(new Coordinate(i,j,k,l)).getDataX(binmax);
                            
                        } else {
                            this.mkTableT0(i, j, k, l);
                            this.mkTableT0Sub1(i, j, k, l);
                        }
                        this.mkTableT0Sub(i, j, k, l);

                        printHistoToFile(i, j, k, l);
                    }
                }
            }
        }
        pw.close();
        pw2.close();
        pw3.close();
        this.getCalib().fireTableDataChanged();  
        
    }
    
    public int NbRunFit = 0;
    int countFits = 0;
    public void runFit(int i, int j, int k, int l) {
            
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** "); 
	
        double[] Tminmax = this.getT0(i, j, k, l);
        
        //Sector Superlayer Slot Cable T0Correction T0Error
        pw.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            Tminmax[0], 
            Tminmax[1]);
        System.out.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            Tminmax[0], 
            Tminmax[1]);
        
        pw2.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            (Tminmax[0]+T00Calib.T00Array[i][j]), 
            Tminmax[1]);
        
        Fitted[i][j][k][l] = true;
        System.out.println((countFits++) +") FITTED ? "+Fitted[i][j][k][l]);
    }
    
    public void mkTableT0Sub1(int i, int j, int k, int l) {
       
        pw2.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            (ReadTT.T0[i][j][k][l]+T00Calib.T00Array[i][j]), 
            ReadTT.T0ERR[i][j][k][l]);
    }
    
    public void mkTableT0Sub(int i, int j, int k, int l) {
       
        pw3.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            (ReadTT.T0[i][j][k][l]+T00Calib.T00Array[i][j]), 
            ReadTT.T0ERR[i][j][k][l]);
    }
    public void mkTableT0(int i, int j, int k, int l) {
       
        pw3.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            ReadTT.T0[i][j][k][l], 
            ReadTT.T0ERR[i][j][k][l]);
    }
     
    private void updateTable(int i, int j,  int k, double t0) {
       this.getCalib().setDoubleValue(t0, "T0", i+1, j+1, k+1);
    }
    
    int counter = 0;
    public  HipoDataSource reader = new HipoDataSource();
    

    int count = 0;
    public static int polarity =-1;
    public List<FittedHit> hits = new ArrayList<>();
    Map<Integer, ArrayList<Integer>> segMapTBHits = new HashMap<Integer, ArrayList<Integer>>();
    Map<Integer, SegmentProperty> segPropMap = new HashMap<Integer, SegmentProperty>();
    //List<FittedHit> hitlist = new ArrayList<>();
    private ReadTT cableMap = new ReadTT();
    @Override
    public void processEvent(DataEvent event) {
        
        if (!event.hasBank("RUN::config")) {
            return ;
        }
        
        DataBank bank = event.getBank("RUN::config");
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
            ReadTT.Load(newRun, "default"); 
            runNumber = newRun; 
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
        // get segment property
        
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");

        //DataBank bank2 = event.getBank("REC::Particle");
        
        for (int j = 0; j < bnkHits.rows(); j++) {
            
            int sec = bnkHits.getInt("sector", j);
            int sl = bnkHits.getInt("superlayer", j);
            int lay = bnkHits.getInt("layer", j);// layer goes from 1 to 6 in data
            int wire = bnkHits.getInt("wire", j);// wire goes from 1 to 112 in data
            //int lay0to35 = (sl - 1) * 6 + lay - 1;
            //int region0to2 = (int) ((sl - 1) / 2);
            int slot1to7  = (int) ((wire - 1) / 16) + 1;
            int wire1to16 = (int) ((wire - 1) % 16 + 1);
            int cable1to6 = ReadTT.CableID[lay - 1][wire1to16 - 1];

            double time = (double) bnkHits.getFloat("time", j)
                    + (double) bnkHits.getFloat("T0", j)
                    + (double)  bnkHits.getFloat("tBeta",j);

            double beta = bnkHits.getFloat("beta",j);

            double TProp = bnkHits.getFloat("TProp",j);

            if(bnkHits.getByte("trkID", j)!=-1 /*&& 11 == bank2.getInt("pid", 0)*/) {

                final int trkID = bnkHits.getByte("trkID", j);

                final int pid = this.readPID(event, trkID);

                if(beta > 0.95 && (pid == 11 || pid == 211 || pid == -211))

                {

                this.TDCHis.get(new Coordinate(sec-1, sl-1, slot1to7-1, cable1to6-1))
                    .fill(time);

                }
            }
            }
        }
    
    public void Plot(int i , int j, int k) {
        
        for (int l = 0; l < nCables6; l++){
            this.getAnalysisCanvas().getCanvas(analTabs).cd(l);
            this.getAnalysisCanvas().getCanvas(analTabs)
                    .draw(this.TDCHis.get(new Coordinate(i, j, k, l)));
            
            if(Fitted[i][j][k][l]==true) {
                this.getAnalysisCanvas().getCanvas(analTabs).cd(l);
                            this.getAnalysisCanvas().getCanvas(analTabs)
                                .draw(this.TDCFits.get(new Coordinate(i, j, k, l)), "same");
            }
        }
    }
    
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_slot     = (String) cc.getValueAt(row, 2);
        System.out.println(str_sector + " " + str_layer + " " );
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);
       int slot      = Integer.parseInt(str_slot);

       if(group.hasItem(sector,layer,slot)==true){
           this.Plot(sector-1, layer-1, slot-1);
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
        
        if(nevent >= 0)
            pass = true;
        
        return pass;
    }

    private double getThreshold(H1F h) {
        // find the bin at which the integral corresponds to 1% of the full integral to max 
        // this is the max range to obtain a threshold from a flat line fit
        double integral = 0;
        double partintegral = 0;
        
        GraphErrors gr = new GraphErrors(); 
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            integral+= h.getBinContent(ix);
        }
        double x = 0;
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            x = h.getDataX(ix);
            double y = h.getBinContent(ix);
            double err = h.getBinError(ix);
            if(err<1) {
                err = 1.4142;
            }
            // fill graph
            gr.addPoint(x, y, 0, err);
            partintegral += h.getBinContent(ix);
            
            if(partintegral>0.5*integral)
                break;
        }
        // fit the graph 
        F1D f0 = new F1D("f0","[p0]", h.getDataX(0), x);
        f0.setParameter(0, 0);
        DataFitter.fit(f0, gr, "QN"); 
        return f0.getParameter(0);
    }
    
    private double[] getT0(int i, int j, int k, int l) {
        System.out.println("Getting t0 for i,j,k,l = "+i+" "+j+" "+k+" "+l );
        H1F h = this.TDCHis.get(new Coordinate(i,j,k,l));

        double thres = 0;//this.getThreshold(h);
        double [] T0val = new double[2];
        double [] T0valCCDB = new double[2];
        double [] T0valERF = new double[2];
        F1D f1 = new F1D("f1","[a]*x+[b]", h.getDataX(0), h.getDataX(20));
        
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+[p0]", 
            h.getDataX(0), h.getDataX(h.getMaximumBin())); 
        
        gausFunc.setParameter(0, h.getMax());
        gausFunc.setParameter(1, h.getDataX(h.getMaximumBin()));
        gausFunc.setParameter(2, 0.05);
        gausFunc.setParameter(3, 0);
        
        DataFitter.fit(gausFunc, h, "QNR");

        System.out.println("gaus");

        //erf fit

        F1D erfcFunc = new F1D("erfcFunc", "0.5*[amp]*erf(-x, -[mean], [sigma])+[p0]", 
                h.getDataX(0), h.getDataX(h.getMaximumBin()));
        
        erfcFunc.setParameter(0, gausFunc.getParameter(0));
        erfcFunc.setParameter(1, gausFunc.getParameter(1));
        erfcFunc.setParameter(2, gausFunc.getParameter(2));
        erfcFunc.setParameter(3, gausFunc.getParameter(3));
        erfcFunc.setParLimits(3, 0., Double.POSITIVE_INFINITY);

        DataFitter.fit(erfcFunc, h, "QNR");

        double T0 = erfcFunc.getParameter(1)-1.5*erfcFunc.getParameter(2);
        double T0Err = Math.sqrt(erfcFunc.parameter(1).error()*erfcFunc.parameter(1).error()+erfcFunc.parameter(2).error()*erfcFunc.parameter(2).error());
        if(Double.isNaN(T0)|| Double.isNaN(T0Err)){
            T0 = 0;
            T0Err = 1.42;
        }
        T0valERF[1] =T0Err;
        T0valERF[0] = T0;

        //second derivative

        double tmidY = gausFunc.getParameter(0)/2;
        double tminY = gausFunc.getParameter(3);
        double del_min_halfmaxY = tmidY-tminY;
        
        double minRangeY = tmidY-del_min_halfmaxY/2;
        double maxRangeY = tmidY;
        if(h.getMax()>tmidY && tmidY+(h.getMax()-tmidY)/3 < h.getMax() ) {
            maxRangeY+=(h.getMax()-tmidY)/3;
        }

        int t0idx  = -1;
        int t0midx = -1;
        double t0 = Double.NEGATIVE_INFINITY; 
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            if(h.getBinContent(ix)>=maxRangeY) {
                t0midx= ix;
                break;
            }
        }
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            if(h.getBinContent(ix)>=minRangeY) {
                t0idx= ix;
                break;
            }
        }

        double[] sec_der = new double[h.getMaximumBin()];
        double[] sec_der_avg = new double[h.getMaximumBin()];
        double[] times = new double[h.getMaximumBin()];

        int index = 0;

        for (int ix =1; ix< h.getMaximumBin()-1; ix++) {
            times[index]=h.getDataX(ix);
            sec_der[index++]=h.getBinContent(ix+1)-2*h.getBinContent(ix)+h.getBinContent(ix-1); //no division by bin size - only a constant scaling factor
        }

        for(int t = 0; t < h.getMaximumBin(); ++t)
        {
            sec_der_avg[t]=(((t-1)>=0)?(sec_der[t-1]):(0.)+sec_der[t]+(((t+1)<=h.getMaximumBin())?(sec_der[t+1]):(0.)))/3.;
        }

        sec_der = sec_der_avg; //use averaging for 2nd derivative

        DataVector vec1 = new DataVector(times);

        DataVector vec2 = new DataVector(sec_der);

        H1F derhis = H1F.create("der", index, vec1, vec2);

        DerHis.put(new Coordinate(i,j,k, l), derhis);

        F1D gausderFunc = new F1D("gausderFunc", "[amp]*gaus(x,[mean],[sigma])", 
        derhis.getDataX(derhis.getMaximumBin()-2), derhis.getMaximumBin()+2);

        gausderFunc.setParameter(0, derhis.getDataY(derhis.getMaximumBin()));
        gausderFunc.setParLimits(0, 0., Double.POSITIVE_INFINITY);
        gausderFunc.setParameter(1, erfcFunc.getParameter(1)-1.5*erfcFunc.getParameter(2));
        gausderFunc.setParLimits(0, 0., 2.*(derhis.getDataX(derhis.getMaximumBin())-derhis.getDataX(derhis.getMaximumBin()-1)));
        gausderFunc.setParameter(2, erfcFunc.getParameter(2)/2);

        DataFitter.fit(gausderFunc, derhis, "QR");

        this.DERFits.get(new Coordinate(i,j,k,l)).setParameter(0, gausderFunc.getParameter(0));
        this.DERFits.get(new Coordinate(i,j,k,l)).setParameter(1, gausderFunc.getParameter(1));
        this.DERFits.get(new Coordinate(i,j,k,l)).setParameter(2, gausderFunc.getParameter(2));

        T0val[0] = ReadTT.T0[i][j][k][l];
        T0val[1] = ReadTT.T0ERR[i][j][k][l];

        T0valCCDB[0] = T0val[0];
        T0valCCDB[1] = T0val[1];

        T0CCDBs.put(new Coordinate(i,j,k, l), T0valCCDB[0]);
        T0CCDBErrs.put(new Coordinate(i,j,k, l), T0valCCDB[1]);

        T0ERFs.put(new Coordinate(i,j,k, l), T0valERF[0]);
        T0ERFErrs.put(new Coordinate(i,j,k, l), T0valERF[1]);

        T0val[0] = gausderFunc.getParameter(1);
        T0val[1] =gausderFunc.parameter(1).error();

        T0s.put(new Coordinate(i,j,k, l), T0val[0]);
        T0Errs.put(new Coordinate(i,j,k, l), T0val[1]);
        this.updateTable(i, j, k, T0val[0]);

        File directory = new File("./T0");
        if (! directory.exists()){
            directory.mkdirs();

        }

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("./T0/T0s_output.txt", true));
            writer.append("" + T0val[0] + "\t");
            writer.append("" + T0val[1] + "\t");
            writer.append("" + T0valCCDB[0] + "\t");
            writer.append("" + T0valCCDB[1] + "\t");
            writer.append("" + T0valERF[0] + "\t");
            writer.append("" + T0valERF[1] + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return T0val;
    }

    private double calcError(double n, double en, double d, double ed) {
        return Math.sqrt((en/d)*(en/d)+(ed*n/(d*d))*(ed*n/(d*d)));
    }

    private int readPID(DataEvent event, int trkId) {
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
    } 

    private void printHistoToFile(int i, int j, int k, int l) {
        File directory = new File("./T0");
        if (! directory.exists()){
            directory.mkdirs();

        }

        EmbeddedCanvas canvas = new EmbeddedCanvas(1200, 800);
        canvas.draw(this.TDCHis.get(new Coordinate(i, j, k, l)));

        DataLine T0Line = new DataLine(this.T0s.get(new Coordinate(i, j, k, l)), 0, this.T0s.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin()));

        T0Line.setLineColor(3);

        DataLine T0LinepE = new DataLine(this.T0s.get(new Coordinate(i, j, k, l))-this.T0Errs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin())/3., this.T0s.get(new Coordinate(i, j, k, l))+this.T0Errs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin())/3.);

        T0LinepE.setLineColor(3);

        canvas.draw(T0Line);

        canvas.draw(T0LinepE);


        DataLine T0ERFLine = new DataLine(this.T0ERFs.get(new Coordinate(i, j, k, l)), 0, this.T0ERFs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin()));

        T0ERFLine.setLineColor(2);

        DataLine T0ERFLinepE = new DataLine(this.T0ERFs.get(new Coordinate(i, j, k, l))-this.T0ERFErrs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin())/3.*2., this.T0ERFs.get(new Coordinate(i, j, k, l))+this.T0ERFErrs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin())/3.*2.);

        T0ERFLinepE.setLineColor(2);

        canvas.draw(T0ERFLine);

        canvas.draw(T0ERFLinepE);

        DataLine T0CCDBLine = new DataLine(this.T0CCDBs.get(new Coordinate(i, j, k, l)), 0, this.T0CCDBs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin()));

        T0CCDBLine.setLineColor(1);

        DataLine T0LineCCDBpE = new DataLine(this.T0CCDBs.get(new Coordinate(i, j, k, l))-this.T0CCDBErrs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin()), this.T0CCDBs.get(new Coordinate(i, j, k, l))+this.T0CCDBErrs.get(new Coordinate(i, j, k, l)), this.TDCHis.get(new Coordinate(i, j, k, l)).getDataY(this.TDCHis.get(new Coordinate(i, j, k, l)).getMaximumBin()));

        T0LineCCDBpE.setLineColor(1);

        canvas.draw(T0CCDBLine);

        canvas.draw(T0LineCCDBpE);
        
        canvas.save("./T0/T0s_" + i + "_"+ j + "_"+ k + "_"+ l + ".png", org.jlab.groot.data.SaveType.PNG);

        EmbeddedCanvas canvas2 = new EmbeddedCanvas(1200, 800);
        canvas2.draw(this.TDCHis.get(new Coordinate(i, j, k, l)));

        canvas2.draw(this.TDCFits.get(new Coordinate(i, j, k, l)), "same");

        canvas2.draw(this.DerHis.get(new Coordinate(i, j, k, l)), "same");

        canvas2.draw(this.DERFits.get(new Coordinate(i, j, k, l)), "same");

        canvas2.save("./T0/T0fit_" + i + "_"+ j + "_"+ k + "_"+ l + ".png", org.jlab.groot.data.SaveType.PNG);
    }
}