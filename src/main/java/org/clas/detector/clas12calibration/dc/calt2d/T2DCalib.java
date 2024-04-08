/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.analysis.FitPanel;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import static org.clas.detector.clas12calibration.dc.t2d.TableLoader.calc_Time;
import org.clas.detector.clas12calibration.dc.t2d.TimeToDistanceEstimator;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import static org.clas.detector.clas12calibration.viewer.T2DViewer.ccdb;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
/**
 *
 * @author ziegler
 */
public class T2DCalib extends AnalysisMonitor{

    private static double DeltaTimeCut = 50;
    public HipoDataSync calwriter = null;
    public HipoDataSync writer = null;
    private HipoDataEvent calhipoEvent = null;
    private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    public FitPanel fp;
    PrintWriter pw = null;
    PrintWriter pw2 = null;
    public PrintWriter pw3 = null;
    File outfile = null;
    private int runNumber;
    private Utilities util = new Utilities();
    private int numberprocessedevents;
    private static double betaAve = 1;
    
    public T2DCalib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames("TrackDoca vs T","TrackDoca vs T Graphs","TrackDoca vs T Fit Resi", 
                "CalcDoca vs T","Time Residuals","Parameters", "Fit Function");
        this.init(false, "v0:vmid:R:tmax:distbeta:delBf:b1:b2:b3:b4");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
       
        if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
            System.out.println(" BANK FOUND........");
        } else {
            System.out.println(" BANK NOT FOUND........");
        }
        calwriter = new HipoDataSync(schemaFactory);
        calwriter.setCompressionType(2);
        writer = new HipoDataSync(schemaFactory);
        writer.setCompressionType(2);
        calhipoEvent = (HipoDataEvent) calwriter.createEvent();
        hipoEvent = (HipoDataEvent) writer.createEvent();
        this.checkFile("TestCalOutPut.hipo");
        this.checkFile("TestOutPut.hipo");
        calwriter.open("TestCalOutPut.hipo");
        calwriter.writeEvent(calhipoEvent);
        writer.open("TestOutPut.hipo");
        writer.writeEvent(hipoEvent);
        
        //init BBin Centers
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k < this.BBins; k++) {
                    BfieldValuesUpd[i][j][k] = BfieldValues[k];
                }
            }
        }
        //init ABin Centers
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k <= this.BBins; k++) {
                    AlphaValuesUpd[i][j][k] = AlphaValues[j];
                }
            }
        }
        
    }
    private Map<Coordinate, H2F> Tvstrkdocas                    = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tvscalcdocas                   = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tresvstrkdocas                 = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, GraphErrors> TvstrkdocasProf        = new HashMap<Coordinate, GraphErrors>();
    private final Map<Coordinate, GraphErrors> TvstrkdocasInit  = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> TvstrkdocasFit         = new HashMap<Coordinate, FitFunction>();
    public Map<Coordinate, MnUserParameters> TvstrkdocasFitPars = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> TvstrkdocasFits            = new HashMap<Coordinate, FitLine>();
    private Map<Coordinate, H1F> timeResi                       = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiFromFile               = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiNew                    = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> fitResi                        = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> B                              = new HashMap<Coordinate, H1F>(); //histogram to get B values centroids
    private Map<Coordinate, H1F> A                              = new HashMap<Coordinate, H1F>(); //histogram to get Alpha values centroids
    private Map<Coordinate, H1F> BAlphaBins                     = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> ParsVsIter                     = new HashMap<Coordinate, H1F>();
    public int nsl = 6;

    public static double[] BfieldValues = new double[]{0.707106781,1.224744871,1.58113883,1.87082869,2.121320344,2.34520788,2.549509757,2.738612788};   
    public static double[] AlphaValues = new double[]{-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40};
    public static double AlphaBinHalfWidth = 2;
    public static int alphaBins = AlphaValues.length; 
    public static int BBins = BfieldValues.length;
    //update middle of the B bins
    public static double[][][] BfieldValuesUpd = new double[2][alphaBins][BBins];
    public static double[][][] AlphaValuesUpd = new double[6][alphaBins][BBins+1];
    int nbinx[] = new int[6];
    int nbiny[] = new int[6];
    double docaBinWidth = 0.025;
    double timeBinWidth = 1.0;
    double maxx[] = new double[]{0.85,0.97,1.34,1.38,1.95,2.0};
    double maxy[] = new double[]{400,410,1200,1500,1000,1200};
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup td = new DataGroup(7,2);
        DataGroup tdp = new DataGroup(14,8);
        DataGroup tdp2 = new DataGroup(14,8);
        DataGroup cd = new DataGroup(7,2);
        DataGroup tr = new DataGroup(6,1);
        DataGroup fr = new DataGroup(6,1);
        
        for(int l=0; l<6; l++) {
            nbinx[l] = (int) Math.ceil(maxx[l]/docaBinWidth);
            nbiny[l] = (int) Math.ceil(maxy[l]/timeBinWidth);
        }
        int ijk = 0;
        int ij = 0;
        for (int i = 0; i < nsl; i++) {
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            timeResi.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            timeResiFromFile.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            timeResiNew.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            fitResi.put(new Coordinate(i), new H1F("fit residual for sly " + (i+1), 100, -0.5, 0.5));
            
            tr.addDataSet(timeResi.get(new Coordinate(i)), i);
            tr.addDataSet(timeResiFromFile.get(new Coordinate(i)), i);
            tr.addDataSet(timeResiNew.get(new Coordinate(i)), i);
            fr.addDataSet(fitResi.get(new Coordinate(i)), i);
            
            for (int j = 0; j < alphaBins; j++) {
                DataGroup trkdvst = new DataGroup(1,1);
                DataGroup dvst = new DataGroup(1,1);
                
                for (int k = 0; k < BBins+1; k++) {
                    DataGroup prfdvst = new DataGroup(1,1);
                    Tvstrkdocas.put(new Coordinate(i,j,k), new H2F("trkDocavsT" + (i + 1)*1000+(j+1)+26, "superlayer" + (i + 1)
                            + ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")"
                            +", B "+k, nbinx[i], 0, maxx[i], nbiny[i], 0, maxy[i]));
                    Tresvstrkdocas.put(new Coordinate(i,j,k), new H2F("trkDocavsTres" + (i + 1)*1000+(j+1)+26, "superlayer" + (i + 1)
                            + ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")"
                            +", B "+k, nbinx[i], 0, maxx[i], nbiny[i], -50, 50));
                    
                    TvstrkdocasProf.put(new Coordinate(i,j,k), new GraphErrors());
                    TvstrkdocasInit.put(new Coordinate(i,j,k), new GraphErrors());
                    TvstrkdocasProf.get(new Coordinate(i,j,k)).setMarkerColor(k+1);
                    TvstrkdocasInit.get(new Coordinate(i,j,k)).setMarkerColor(k+1);
                    TvstrkdocasInit.get(new Coordinate(i,j,k)).setMarkerStyle(2);
                    TvstrkdocasInit.get(new Coordinate(i,j,k)).setTitle( "superlayer" + (i + 1)
                            + ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")");
                    Tvscalcdocas.put(new Coordinate(i,j,k), new H2F("calcDocavsT" + (i + 1)*1000+(j+1)+26, "superlayer" + (i + 1)
                            + ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")"
                            +", B "+k, nbinx[i], 0, maxx[i], nbiny[i], 0, maxy[i]));
                    tdp.addDataSet(TvstrkdocasProf.get(new Coordinate(i,j,k)), ijk);
                    tdp.addDataSet(Tresvstrkdocas.get(new Coordinate(i,j,k)), 0);
                    //trkdvst.addDataSet(Tvstrkdocas.get(new Coordinate(i,j,k)), 0);
                    prfdvst.addDataSet(TvstrkdocasProf.get(new Coordinate(i,j,k)), 0);
                    prfdvst.addDataSet(TvstrkdocasInit.get(new Coordinate(i,j,k)), 0);
                    
                    TvstrkdocasFits.put(new Coordinate(i,j,k), new FitLine());
                    prfdvst.addDataSet(TvstrkdocasFits.get(new Coordinate(i,j,k)), 0);
                    this.getDataGroup().add(prfdvst, 0, i+1, j+1);
                    
                    ijk++;
                }
            }
        }
         //Alpha centroids
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k <= this.BBins+1; k++) {
                    A.put(new Coordinate(i,j,k), new H1F("A centroid " +(i + 1)*1000+(j+1)+26, 100, -36, 36));
                }
            }
        }
        //B centroids
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                BAlphaBins.put(new Coordinate(i,j), new H1F("B  " , 100, 0.0, 3.0));
                for (int k = 0; k < this.BBins; k++) {
                    B.put(new Coordinate(i,j,k), new H1F("B centroid " +(i + 1)*1000+(j+1)+26, 100, 0.0, 3.0));
                }
            }
        }
        this.getDataGroup().add(td, 0,0,0);
        this.getDataGroup().add(tdp,1,0,0);
        this.getDataGroup().add(cd, 2,0,0);
        this.getDataGroup().add(tr, 3,0,0);
        this.getDataGroup().add(fr, 4,0,0);
        
        for (int i = 0; i < nsl; i++) {
            for (int j = 0; j < alphaBins; j++) {
                this.getCalib().addEntry(0,i+1,j+1);
                //blank out
                this.getCalib().setDoubleValue((double)999, "v0", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "vmid", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "R", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "tmax", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "distbeta", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "delBf", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "b1", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "b2", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "b3", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "b4", 0, i+1, j+1);
            }
        }
        
        this.getCalib().fireTableDataChanged();
    }
    private void updateTable(int i, int j) {
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(0), "v0", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(1), "vmid", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(2), "R", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(3), "tmax", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(4), "distbeta", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(5), "delBf", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(6), "b1", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(7), "b2", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(8), "b3", 0, i+1, j+1);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(9), "b4", 0, i+1, j+1);
    }    
    private void resetTable(int i, int j) {
        this.getCalib().setDoubleValue(999., "v0", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "vmid", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "R", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "tmax", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "distbeta", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "delBf", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "b1", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "b2", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "b3", 0, i+1, j+1);
        this.getCalib().setDoubleValue(999., "b4", 0, i+1, j+1);
    }    
    @Override
    public void plotHistos() {
        String[] Names = {"TrackDoca vs T","TrackDoca vs T Graphs","TrackDoca vs T Fit Resi","CalcDoca vs T","Time Residuals","Parameters",
                            "Fit Function"};
        for(int s = 0; s<3; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
            int NumPads = 
            this.getAnalysisCanvas().getCanvas(Names[s]).getCanvasPads().size();
            for (int n = 0; n < NumPads; n++) {
                this.getAnalysisCanvas().getCanvas(Names[s]).getPad(n).getAxisZ().setLog(true);
            }
        }
        this.getAnalysisCanvas().getCanvas(Names[2]).getPad(0).getAxisZ().setLog(true);
        for(int s = 3; s<6; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
        }
        this.getAnalysisCanvas().getCanvas(Names[4]).divide(this.nsl, 2);
        this.getAnalysisCanvas().getCanvas(Names[5]).divide(6, 6);
        
        this.getAnalysisCanvas().getCanvas(Names[6]).divide(4, 3);
        
        this.getAnalysisCanvas().getCanvas( "TrackDoca vs T Fit Resi").getPad().getAxisY().setRange(-150, 150);
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").update();
        this.getAnalysisCanvas().getCanvas( "TrackDoca vs T Fit Resi").update();
        this.getAnalysisCanvas().getCanvas("CalcDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("Time Residuals").update();
        this.getAnalysisCanvas().getCanvas("Parameters").update();
        this.getAnalysisCanvas().getCanvas("Fit Function").update();
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        
        try {
            calwriter.close();
            writer.close();
            this.UpdateBBinCenters();
            this.UpdateAlphaBinCenters();
            
            for (int i = 0; i < this.nsl; i++) {
                for (int j = 0; j < this.alphaBins; j++) {
                    this.filltrkDocavsTGraphs(i,j);
                }
                //runFit(i);
            }
            
            reLoadFitPars();
            
            //fp.refit();
            //pw.close();
            this.getAnalysisCanvas().getCanvas("Time Residuals").divide(nsl, 3);
            //
            for(int i = 0; i<this.nsl; i++) {
                this.runInitFit(i);
                this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i);
                this.fitTimeResPlot(timeResiFromFile.get(new Coordinate(i)),
                        this.getAnalysisCanvas().getCanvas("Time Residuals"));
            }
           
            //fp.setGreenFitButton();
            this.plotFits(true);
            
            System.out.println("PLOTS WITH CCDB CONSTANTS DONE");
            this.plotHistos();
            
            for (int i = 0; i < this.nsl; i++) {
                for (int j = 0; j < this.alphaBins; j++) {
                    this.Plot(i,j); 
                }
            }
            eventProcessingDone = true;
            System.out.println("ANALYSIS Done ....");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(T2DCalib.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void plotFits(boolean fitted) throws FileNotFoundException {
        if(fitted==true) {
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "Files/ccdb_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
           
            String fileName2 = "Files/parameteranderror_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
           
            
            
            pw = new PrintWriter(fileName);
            pw.printf("#& sector superlayer component v0 deltanm tmax distbeta delta_bfield_coefficient b1 b2 b3 b4 delta_T0 c1 c2 c3\n");
            pw2 = new PrintWriter(fileName2);
            pw2.printf("#& sector superlayer component v0 +/-v0 tmax +/-tmax vmid +/-vmid delta_bf +/-delta_bf distbeta +/-distbeta \n");
            
            int ij =0;
            int ip =0;
            NbRunFit++;
            for (int i = 0; i < this.nsl; i++) {
               
                for(int p = 0; p<6; p++) {
                    ParsVsIter.get(new Coordinate(i,p)).setBinContent(NbRunFit, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
                    ParsVsIter.get(new Coordinate(i,p)).setBinError(NbRunFit, TvstrkdocasFitPars.get(new Coordinate(i)).error(p));
                    ParsVsIter.get(new Coordinate(i,p)).setOptStat(0);
                    this.getAnalysisCanvas().getCanvas("Parameters").cd(ip);
                    GraphErrors gr = new GraphErrors();
                    double min = ParsVsIter.get(new Coordinate(i,p)).getMin()-2*ParsVsIter.get(new Coordinate(i,p)).getBinError(1);
                    double max = ParsVsIter.get(new Coordinate(i,p)).getMax()+2*ParsVsIter.get(new Coordinate(i,p)).getBinError(1);
                    if(Math.abs(min)<1.e-06 && Math.abs(max)<1.e-06) {
                        min = -0.1;
                        max = 0.1;
                    }
                    
                    this.getAnalysisCanvas().getCanvas("Parameters").
                            draw(ParsVsIter.get(new Coordinate(i,p)));

                    //this.getAnalysisCanvas().getCanvas("Parameters").getPad(ip).getAxisX().setRange(0.5, 11.5);
                    //this.getAnalysisCanvas().getCanvas("Parameters").getPad(ip).getAxisY().setRange(min, max);
                    ip++;
                }
            }

            
            for (int i = 0; i < this.nsl; i++) {

                for (int j = 0; j < this.alphaBins; j++) {
                    if(i<2 || i>3) {
                        if(TvstrkdocasProf.get(new Coordinate(i, j, BBins)).getVectorX().size()>0) {
                            this.updateTable(i,j);
                            TvstrkdocasFits.put(new Coordinate(i,j,BBins), new FitLine("f"+""+i+""+j+"0", i, j, BBins, 
                            TvstrkdocasFitPars.get(new Coordinate(i))));
                            TvstrkdocasFits.get(new Coordinate(i, j, BBins)).setLineStyle(4);
                            TvstrkdocasFits.get(new Coordinate(i, j, BBins)).setLineWidth(5);
                            TvstrkdocasFits.get(new Coordinate(i, j, BBins)).setLineColor(8);
                        } else {
                            //this.resetTable(i,j);
                        }

                    } else {
                        for(int k = 0; k < this.BBins; k++) { 
                            if(TvstrkdocasProf.get(new Coordinate(i, j, k)).getVectorX().size()>0){
                                this.updateTable(i,j);
                                TvstrkdocasFits.put(new Coordinate(i,j,k), new FitLine("f"+""+i+""+j+""+k, i, j, k, 
                                TvstrkdocasFitPars.get(new Coordinate(i))));
                                TvstrkdocasFits.get(new Coordinate(i, j, k)).setLineStyle(4);
                                TvstrkdocasFits.get(new Coordinate(i, j, k)).setLineWidth(5);
                                TvstrkdocasFits.get(new Coordinate(i, j, k)).setLineColor(k+1);
                            } else {
                                //this.resetTable(i,j);
                            }
                        }
                    }
                    ij++;
                }
            }
            this.getCalib().fireTableDataChanged();    
            for(int isec = 0; isec < 6; isec++) {
                for(int i = 0; i<6; i++) {
                    pw.printf("%d\t %d\t %d\t %.6f\t %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %d\t %.6f\t %.6f\t %d\n",
                    (isec+1), (i+1), 0,
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(0),
                    0,
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(3),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(4),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(5),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(6),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(7),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(8),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(9),
                    0,
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(2),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(1),
                    0);
                    
                    pw2.printf("%d\t %d\t %d\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\n",
                            (isec+1), (i+1), 0,
                            //v0
                            TvstrkdocasFitPars.get(new Coordinate(i)).value(0),
                            TvstrkdocasFitPars.get(new Coordinate(i)).error(0),
                            //tmax
                            TvstrkdocasFitPars.get(new Coordinate(i)).value(3),
                            TvstrkdocasFitPars.get(new Coordinate(i)).error(3),
                            //vmid
                            TvstrkdocasFitPars.get(new Coordinate(i)).value(1),
                            TvstrkdocasFitPars.get(new Coordinate(i)).error(1),
                            //deltaBf
                            TvstrkdocasFitPars.get(new Coordinate(i)).value(5),
                            TvstrkdocasFitPars.get(new Coordinate(i)).error(5),
                            //distbeta
                            TvstrkdocasFitPars.get(new Coordinate(i)).value(4),
                            TvstrkdocasFitPars.get(new Coordinate(i)).error(4));
                            
                    
                }
            }
            pw.close(); 
            pw2.close(); 
            //this.rePlotResi();
        }
    }
    int maxIter = 10;
    private int maxNfits = 10;
    double[][] fixPars = new double[6][10];
    public boolean eventProcessingDone=false;
    
    public int NbRunFit = 0;
    
    public void runInitFit(int i) {
        TvstrkdocasFit.put(new Coordinate(i), 
                new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
        
    }
    public void initFitParsToFile() throws FileNotFoundException {
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String fileName3 = "Files/minuit_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
        pw3 = new PrintWriter(fileName3);
    }
    
    double v0limits[][]         = new double[][]{{0.002, 0.009},{0.002, 0.009},{0.002, 0.009},{0.002, 0.009},{0.002, 0.009},{0.002, 0.009}};     //  limits for each superlayer
    double vmidlimits[][]       = new double[][]{{0.002, 0.009},{0.002, 0.009},{0.001, 0.009},{0.001, 0.009},{0.001, 0.009},{0.001, 0.009}};     
    double Rlimits[][]          = new double[][]{{0.55,0.75},{0.55,0.75},{0.55,0.75},{0.55,0.75},{0.55,0.75},{0.55,0.75}};                            
    double distbetalimits[][]   = new double[][]{{0.02, 0.1},{0.02, 0.1},{0.02, 0.1},{0.02, 0.1},{0.02, 0.1},{0.02, 0.1}};                       
    double delBflimits[][]      = new double[][]{{0.10, 0.4},{0.10, 0.4},{0.10, 0.4},{0.10, 0.4},{0.10, 0.4},{0.10, 0.4}};                       //  limits for each superlayer
    
    double limits[][][] = new double[][][]{v0limits, vmidlimits, Rlimits, distbetalimits, delBflimits};
    boolean useFixedBoundsMethod = true;
    double[][] chi2FitNum = new double[100][6];
    int[] fitNum =new int[6];
    
    //get R, distbeta, delBf, for best combine ch2 for all superlayers
    
    public double[][] estimateFixedParsPerRegion(boolean fixFit[][], MnMigrad scanner[], MnMigrad fitter[]) { 
        double errs2 = errs[2];
        double errs4 = errs[4];
        int nR = (int) Math.ceil((Rlimits[0][1]-Rlimits[0][0])/(double)errs2);
        int ndbeta = (int) Math.ceil((distbetalimits[0][1]-distbetalimits[0][0])/(double)errs4);
       
        double[] bestR = new double[3];
        double[] bestDistbeta = new double[3];

        double bestchi2[] = new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};

        double R ;
        double distbeta ;
        System.out.println("Estimating PARS "+nR+" "+ndbeta);
        int cnt =0;
        
        for(int ri =1; ri<nR-1; ri++) {
            for(int di =1; di<ndbeta-1; di++) { 
                R=Rlimits[0][0]+(double)ri*errs2;
                distbeta=distbetalimits[0][0]+(double)di*errs4;
                double[] c2=new double[3];
                cnt++;
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta);
                for(int i =0; i<2; i++) {
                    String s="";
                    s+=(" ******************************************");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1));
                    s+=(" ******************************************");
                    fMin fm = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2[0]+= fm.getChi2();
                }
                for(int i =2; i<4; i++) {
                    String s="";
                    s+=(" ******************************************");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1));
                    s+=(" ******************************************");
                    fMin fm = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2[1]+= fm.getChi2();
                }
                for(int i =4; i<6; i++) {
                    String s="";
                    s+=(" ******************************************");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1));
                    s+=(" ******************************************");
                    fMin fm = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2[2]+= fm.getChi2();
                }
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta +" c2 "+c2[0]+" "+c2[1]+" "+c2[2]);
                for(int j = 0; j < 3; j++) {
                    if(c2[j]<bestchi2[j]) {
                        bestR[j] = R;
                        bestDistbeta[j] = distbeta;
                        bestchi2[j] = c2[j];
                    }
                }
            }
        }
        double[][] result = new double[2][3];
        for(int k =0; k<3; k++) {
            result[0][k] = bestR[k];
            result[1][k] = bestDistbeta[k];
        }
        return result;
    }
    
    public double[] estimateFixedPars(boolean fixFit[][], MnMigrad scanner[], MnMigrad fitter[]) { 
        double errs2 = errs[2];
        double errs4 = errs[4];
        int nR = (int) Math.ceil((Rlimits[0][1]-Rlimits[0][0])/(double)errs2);
        int ndbeta = (int) Math.ceil((distbetalimits[0][1]-distbetalimits[0][0])/(double)errs4);
       
        double bestR = Rlimits[0][0];
        double bestDistbeta = distbetalimits[0][0];

        double bestchi2 = Double.POSITIVE_INFINITY;

        double R = Rlimits[0][0];
        double distbeta = distbetalimits[0][0];
        System.out.println("Estimating PARS "+nR+" "+ndbeta);
        int cnt =0;
        
        for(int ri =1; ri<nR-1; ri++) {
            for(int di =1; di<ndbeta-1; di++) { 
                R=Rlimits[0][0]+(double)ri*errs2;
                distbeta=distbetalimits[0][0]+(double)di*errs4;
                double c2=0;
                cnt++;
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta);
                for(int i =0; i<6; i++) {
                    String s="";
                    s+=(" ******************************************");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1));
                    s+=(" ******************************************");
                    fMin fm = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2+= fm.getChi2();
                }
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta +" c2 "+c2);
                if(c2<bestchi2) {
                    bestR = R;
                    bestDistbeta = distbeta;
                    bestchi2 = c2;
                    System.out.println(cnt+"] best R "+R+" disbeta "+distbeta +" c2 "+c2);
                }
            }
        }
        
        return new double[] {bestR, bestDistbeta};
    }
    
    public fMin getfMinFixedRDPars(int i, boolean fixFit[][], MnMigrad scanner, MnMigrad fitter, 
            double R, double distbeta, boolean reset, String s) {
        
        double edm = Double.POSITIVE_INFINITY;
        double edm2 = Double.POSITIVE_INFINITY;
        double bestchi2 = Double.POSITIVE_INFINITY;
        double bestMchi2 = Double.POSITIVE_INFINITY;
        
        
        System.out.println(s); 
        FunctionMinimum min = null ;
        FunctionMinimum bestmin = null ;
        FunctionMinimum min2 = null ;
        FunctionMinimum bestmin2 = null ;
        
        
        for(int pi = 0; pi<3; pi++) {
            scanner.setLimits(pi, limits[pi][i][0], limits[pi][i][1]);
            fitter.setLimits(pi, limits[pi][i][0], limits[pi][i][1]);
        }

        for(int pi = 4; pi<6; pi++) {
            scanner.setLimits(pi, limits[pi-1][i][0], limits[pi-1][i][1]);
            fitter.setLimits(pi, limits[pi-1][i][0], limits[pi-1][i][1]);
        }
        scanner.setLimits(3, TvstrkdocasFitPars.get(new Coordinate(i)).value(3)-50, TvstrkdocasFitPars.get(new Coordinate(i)).value(3)+50);
        fitter.setLimits(3, TvstrkdocasFitPars.get(new Coordinate(i)).value(3)-50, TvstrkdocasFitPars.get(new Coordinate(i)).value(3)+50);
        
        
        try {
            min = scanner.minimize();
        } catch (Exception e) {
            // Handle the exception appropriately
            System.err.println("An error occurred during minimization: " + e.getMessage());
            // You may want to log the exception or take other actions depending on your application
        }
        if(fixFit[2][i]==false) {
            min.userParameters().setValue(2,R);
            scanner.fix(2);
            fitter.fix(2);
        }
        
        if(fixFit[4][i]==false) {
            min.userParameters().setValue(4,distbeta);
            scanner.fix(4);
            fitter.fix(4);
        }
        int itercnt=0;
        for(int it = 0; it<maxIter; it++) {
                try {
                        min = scanner.minimize();
                    } catch (Exception e) {
                        // Handle the exception appropriately
                        System.err.println("An error occurred during minimization: " + e.getMessage());
                        // You may want to log the exception or take other actions depending on your application
                    }
                itercnt++;
                if(FitFunction.chi2<bestchi2) {
                    bestchi2 = FitFunction.chi2;
                    bestmin = min;
                    
                }
                if(edm-FitFunction.chi2<0.1 || FitFunction.chi2+10>edm) break;
                edm = FitFunction.chi2;
        } 
        for (int p = 0; p < 10; p++) {
            fitter.setValue(p, bestmin.userParameters().value(p));
        }
        
        int itercnt2=0;
        for(int it = 0; it<maxIter; it++) {
            
            try {
                    min2 = fitter.minimize();
                } catch (Exception e) {
                    // Handle the exception appropriately
                    System.err.println("An error occurred during minimization: " + e.getMessage());
                    // You may want to log the exception or take other actions depending on your application
                }
            itercnt2++;
            
            if(FitFunction.chi2<bestMchi2) {
                bestMchi2 = FitFunction.chi2;
                bestmin2 = min2;
                if(edm2-FitFunction.chi2<0.01) break;
                edm2 = FitFunction.chi2;
            }
        }
        if(bestmin2==null || bestMchi2>bestchi2) {
            bestMchi2 = bestchi2;
            bestmin2 = bestmin;
        }
        
        if(fixFit[2][i]==false) {
            scanner.release(2);
            fitter.release(2);
        }
        if(fixFit[4][i]==false) {
            scanner.release(4);
            fitter.release(4);
        }
       
        if(reset) {
            for (int p = 0; p < 10; p++) {
                scanner.setValue(p, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
                fitter.setValue(p, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
            }
        }
        System.out.println(+itercnt+"] SCAN CHI2 "+bestchi2);
        System.out.println(itercnt2+"] MIGRAD CHI2 "+bestMchi2);
        
        System.gc();
        
        return new fMin(min2, bestMchi2);
    }
    
    
    
    public void fitWithFixedParsPerRegion(int ridx,boolean fixFit[][], double pars[], MnMigrad scanner[], MnMigrad fitter[]) { 
        fMin results [] = new fMin[6];
       
        
        for(int i =2*ridx; i<2*ridx+2; i++) {
            
            String s2="";
            s2+=(" ******************************************");
            s2+=("   RUNNING THE PARAMETER FIT FOR SUPERLAYER "+(i+1));
            s2+=(" ******************************************");
            fMin fm2 = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], pars[0], pars[1], false, s2);
            FunctionMinimum fmin=null;
            if(fm2.getFcnMin().isValid()) {
                results[i] = fm2;
                fmin = fm2.getFcnMin();
                System.out.println("UPDATED "+fmin.toString());
                TvstrkdocasFitPars.put(new Coordinate(i),fmin.userParameters()); 
            } 
        }
    }
    public void fitWithFixedPars(boolean fixFit[][], double pars[], MnMigrad scanner[], MnMigrad fitter[]) { 
        fMin results [] = new fMin[6];
        
        for(int i =0; i<6; i++) {
            
            String s2="";
            s2+=(" ******************************************");
            s2+=("   RUNNING THE PARAMETER FIT FOR SUPERLAYER "+(i+1));
            s2+=(" ******************************************");
            fMin fm2 = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], pars[0], pars[1], false, s2);
            FunctionMinimum fmin=null;
            if(fm2.getFcnMin().isValid()) {
                results[i] = fm2;
                fmin = fm2.getFcnMin();
                System.out.println("UPDATED "+fmin.toString());
                TvstrkdocasFitPars.put(new Coordinate(i),fmin.userParameters()); 
            } 
        }
    }
    public boolean useBProf = false;
    public void runFit(boolean fixFit[][]) {
        
        MnMigrad scanner[] = new MnMigrad[6];
        MnMigrad fitter[] = new MnMigrad[6];
        double[] pars = new double[2];
        pars[0] = TvstrkdocasFitPars.get(new Coordinate(0)).value(2);
        pars[1] = TvstrkdocasFitPars.get(new Coordinate(0)).value(4);
        
        for(int i =0; i<6; i++) {
            TvstrkdocasFitPars.get(new Coordinate(i)).fix(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).fix(p);
                }
            }
            TvstrkdocasFit.put(new Coordinate(i), 
                                         new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
            scanner[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                TvstrkdocasFitPars.get(new Coordinate(i)),0);
            fitter[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                TvstrkdocasFitPars.get(new Coordinate(i)),1);
            
        }
        
        this.fitWithFixedPars(fixFit, pars, scanner, fitter);
        for(int i =0; i<6; i++) {
            TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
                }
            }
        }
    }
    
    public void runParamScan(boolean fixFit[][]) {
        MnMigrad scanner[] = new MnMigrad[6];
        MnMigrad fitter[] = new MnMigrad[6];
        for(int i =0; i<6; i++) {
            TvstrkdocasFitPars.get(new Coordinate(i)).fix(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).fix(p);
                }
            }
            TvstrkdocasFit.put(new Coordinate(i), 
                                         new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
            scanner[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                TvstrkdocasFitPars.get(new Coordinate(i)),0);
            fitter[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                TvstrkdocasFitPars.get(new Coordinate(i)),1);
            
        }
        
        //double[] pars = this.estimateFixedPars(fixFit, scanner, fitter);
        //System.out.println("PARAMETERS R "+pars[0]+" distbeta "+pars[1]);
        double[][] pars2 = this.estimateFixedParsPerRegion(fixFit, scanner, fitter);
        for(int i = 0; i<3; i++) {
            double[] pars = new double[2];
            pars[0] = pars2[0][i];
            pars[1] = pars2[1][i];
            System.out.println(i+"] PARAMETERS R "+pars2[0][i]+" distbeta "+pars2[1][i]);
            this.fitWithFixedParsPerRegion(i,fixFit, pars, scanner, fitter);
        }
        for(int i =0; i<6; i++) {
            TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
                }
            }
        }
    }
    
    int counter = 0;
    private int iterationNum = 0;
    public  HipoDataSource calreader = new HipoDataSource();
    public  HipoDataSource reader = new HipoDataSource();
    private TimeToDistanceEstimator t2d = new TimeToDistanceEstimator();
    public void reCook() {
        iterationNum++;
        fp.setRedFitButton();
        //reset histos to refill
        for (int i = 0; i < this.nsl; i++) {
            timeResi.get(new Coordinate(i)).reset();
            timeResiNew.get(new Coordinate(i)).reset();
            for (int j = 0; j < alphaBins; j++) {
                for (int k = 0; k < BBins+1; k++) {
                    Tvstrkdocas.get(new Coordinate(i,j,k)).reset();
                    Tvscalcdocas.get(new Coordinate(i,j,k)).reset();
                    Tresvstrkdocas.get(new Coordinate(i,j,k)).reset();
                    if(A.containsKey(new Coordinate(i,j,k))) A.get(new Coordinate(i,j,k)).reset();
                    if(k<BBins && B.containsKey(new Coordinate(i,j,k))) B.get(new Coordinate(i,j,k)).reset();
                }
            }
        }
        System.out.println("***********************************************");       
        System.out.println("****** Reprocessing TestCalOutPut.hipo ********");
        System.out.println("***********************************************");
        calreader = new HipoDataSource();
        calreader.open("TestCalOutPut.hipo");
        System.out.println("Events in hipofile " +  calreader.getSize() );  
        int numberofeventsinfile = calreader.getSize();
        int eventcounter = 0;
        while (calreader.hasEvent()) { 
            hits.clear();
            DataEvent event = calreader.getNextEvent();
            eventcounter++;
            if ((eventcounter%10000 == 0) && (eventcounter < numberofeventsinfile) ) {
             	 System.out.println("Processed " + eventcounter + " events from " + numberofeventsinfile);  
            }
            if(event.hasBank("TimeBasedTrkg::TBHits")) { 
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        hits.add(this.getCalHit(bnkHits, i));
                }
                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit, true); //update the hit with previous cal results
                }
                //refit with new constants
                Refit rf = new Refit();
                rf.reFit(hits, false);    //refit to get the parameters
                
                for(FittedHit hit : hits) {
                    //filling the timeResi for the previously calibrated hits
                    timeResi.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                }
            } 
        }
        
        // the newly calibrated hits
        System.out.println("*************************************************");       
        System.out.println("*** Done Reprocessing with initial parameters ***");
        System.out.println("*************************************************");
        
        
        reLoadFitPars();
        //Parameters are now fit values
        System.out.println("************  Fit Parameters Reloaded! ************");
        calreader.gotoEvent(0);
        eventcounter = 0;
        while (calreader.hasEvent()) { 
            calhits.clear();
            DataEvent event = calreader.getNextEvent();
            eventcounter++;
            if ((eventcounter%10000 == 0) && (eventcounter < numberofeventsinfile) ) {
            	 System.out.println("Processed " + eventcounter + " events from " + numberofeventsinfile);  
            }
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        calhits.add(this.getCalHit(bnkHits, i));
                }
                for(FittedHit hit : calhits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit, true);
                }
                //refit with new constants
                Refit rf = new Refit();
                rf.reFit(calhits, true);    //uaw only not out of time hits to get the calib cst
                for(FittedHit hit : calhits) {
                    if(hit.get_OutOfTimeFlag()) continue;
                    double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
                    double alphaUncor = hit.getAlpha()+(double)T2DCalib.polarity*theta0;
                    int alphaBin = this.getAlphaBin(alphaUncor);
                    double bFieldVal = (double) hit.getB();
                    
                    if(alphaBin!=-1 && !hit.get_OutOfTimeFlag()) {
                        double calibTime = (double) (hit.get_TDC() - hit.getTProp()
                                            - hit.getTFlight() - hit.getTStart()
                                            - hit.getT0()); 
                        double yf = TvstrkdocasFits.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, this.BBins)).evaluate(hit.get_ClusFitDoca());
                        
                            Tvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, this.BBins))
                                        .fill(hit.get_ClusFitDoca(), calibTime);
                            Tvscalcdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, this.BBins))
                                        .fill(hit.get_Doca(), calibTime);
                        if(!Double.isNaN(yf)) {
                            Tresvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, this.BBins)).fill(hit.get_ClusFitDoca(), calibTime-yf);
                        }
                        //Fill region 2 for different b-field values
                        if(hit.get_Superlayer() >2 && hit.get_Superlayer() <5) { 
                            int bBin = this.getBBin(bFieldVal);
                            double r2yf = TvstrkdocasFits.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, bBin)).evaluate(hit.get_ClusFitDoca());
                            
                                Tvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, bBin))
                                            .fill(hit.get_ClusFitDoca(), calibTime);
                                Tvscalcdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, bBin))
                                            .fill(hit.get_Doca(),calibTime);
                            if(!Double.isNaN(r2yf)) {
                                Tresvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, bBin))
                                        .fill(hit.get_ClusFitDoca(), calibTime-r2yf);
                            }
                        }
                        
                        if(hit.get_Superlayer()<3 || hit.get_Superlayer()>4) { 
                            A.get(new Coordinate(hit.get_Superlayer()-1, alphaBin, BBins))
                                        .fill(alphaUncor);
                            //System.out.println("Filling alpha for spl "+superlayer+" alpha bin "+alphaBin+" alpha "+alphaUncor);
                        }

                        // fill B values histograms
                        if(hit.get_Superlayer() ==3 || hit.get_Superlayer() ==4) {
                            int bBin = this.getBBin(bFieldVal);
                            B.get(new Coordinate(hit.get_Superlayer()-3, alphaBin, bBin))
                                    .fill(bFieldVal);
                            A.get(new Coordinate(hit.get_Superlayer()-1, alphaBin, this.getBBin(bFieldVal)))
                                    .fill(alphaUncor);
                            //System.out.println("Filling alpha for spl "+superlayer+" alpha bin "+alphaBin+" Bin "+this.getBBin(bFieldVal)+" alpha "+alphaUncor);
                        }
                        
                    }
                }
                
                rf.reFit(calhits, false);    //use all hits
                for(FittedHit hit : calhits) {
                    //filling the timeResi for the newly calibrated hits
                        timeResiNew.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                    
                }
            }
        }
        this.UpdateAlphaBinCenters();
        this.UpdateBBinCenters();
        System.out.println("REMAKING PROFILES");
        for (int i = 0; i < this.nsl; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                this.filltrkDocavsTGraphs(i,j);
                System.out.println("PROFILE "+i+" "+j+" is OK!");
            }
        }
        System.out.println("RECOOKING DONE WITH THE NEW CONSTANTS!");
        System.out.println("CHECK THE RESIDUALS!");
        calreader.close();
        fp.setGreenFitButton();
        
    }
    
    public void rePlotResi() {
        
        this.getAnalysisCanvas().getCanvas("Time Residuals").clear();
        this.getAnalysisCanvas().getCanvas("Time Residuals").divide(nsl, 3);
        //
        for(int i = 0; i<this.nsl; i++) {
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i);
            this.fitTimeResPlot(timeResiFromFile.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i+6);
            this.fitTimeResPlot(timeResi.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i+12);
            //this.getAnalysisCanvas().getCanvas("Time Residuals").draw(timeResiNew.get(new Coordinate(i)));
            this.fitTimeResPlot(timeResiNew.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
        }
        //for(int i = 0; i<this.nsl; i++) {
        //    this.getAnalysisCanvas().getCanvas("Fit Residuals").cd(i);
        //    this.getAnalysisCanvas().getCanvas("Fit Residuals").draw(fitResi.get(new Coordinate(i)));
        //}
        System.out.println("REPROCESSING DONE!");
        reader.close();
        this.getCalib().fireTableDataChanged();  
    }
         
    private void fitTimeResPlot(H1F h1, EmbeddedCanvas canvasRes) {
        if (h1==null) return;
        F1D gaus1Func = new F1D("gaus1Func", "[amp]*gaus(x,[mean],[sigma])", -0.3, 0.3); 
        gaus1Func.setParameter(0, h1.getMax());
        gaus1Func.setParameter(1, 0.0);
        gaus1Func.setParameter(2, 0.05);
        DataFitter.fit(gaus1Func, h1, "Q");
        
        int binMax = h1.getMaximumBin(); 
        double hm = h1.getDataY(binMax)/3; 
        double xlo = 0; 
        double xhi =0; 
        for(int b = 0; b<binMax-1; b++) {
            if(h1.getDataY(b)<=hm && h1.getDataY(b+1) >hm) {
                xlo = h1.getDataX(b);
                break;
            }
                
        }
        
        for(int b = binMax; b<99; b++) {
            if(h1.getDataY(b+1)<=hm && h1.getDataY(b) >hm) {
                xhi = h1.getDataX(b);
                break;
            }
                
        }
        gaus1Func.setRange(xlo, xhi);
        DataFitter.fit(gaus1Func, h1, "Q");
        
        //refit using a double gaussian 
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+0.25*[amp]*gaus(x,[mean],[sigma2])", -0.3, 0.3); 
        gausFunc.setLineColor(4);
        gausFunc.setLineStyle(1);
        gausFunc.setLineWidth(2);
        gausFunc.setParameter(0, gaus1Func.getParameter(0));
        gausFunc.setParameter(1, gaus1Func.getParameter(1));
        gausFunc.setParameter(2, gaus1Func.getParameter(2)*0.75);
        gausFunc.setParameter(3, gaus1Func.getParameter(2));
        gausFunc.setOptStat(11100);
        h1.setOptStat(11); //only number of entries
        //canvasRes.clear();
        
        DataFitter.fit(gausFunc, h1, "Q");
        //gausFunc.setOptStat(101100);
        //gausFunc.setOptStat(101100); //mean and both sigmas
        
        
        int effSig = (int) Math.round(10000*Math.sqrt(gausFunc.getParameter(2)*gausFunc.getParameter(2)+
                                            0.25*0.25*gausFunc.getParameter(3)*gausFunc.getParameter(3))/
                                            Math.sqrt(1+0.25*0.25));
        String t = "Eff. Sig. "+effSig + "microns";
        h1.setTitle(t);
        canvasRes.draw(h1, "E1");
        //canvasRes.draw(gausFunc, "same");
        
        
    }
    public int getAlphaBin(double alpha) {
        int v = -1;
        for(int i = 0; i<T2DCalib.AlphaValues.length; i++) {
            
            if(Math.abs(alpha-T2DCalib.AlphaValues[i])<this.AlphaBinHalfWidth)
                v = i;
        } 
        
        return v;
    }

    private int getBBin(double bFieldVal) {
        
        int v = BfieldValues.length-1;
        //BfieldValues = new double[]{0.0000, 1.0000, 1.4142, 1.7321, 2.0000, 2.2361, 2.4495, 2.6458};
        //BfieldValues^2 = new double[]{0.0000, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000};
        double BSqrBinHalfWidth = 0.5;
        for(int i = 0; i<BfieldValues.length; i++) {
            if(Math.abs(bFieldVal*bFieldVal-T2DCalib.BfieldValues[i]*T2DCalib.BfieldValues[i])<BSqrBinHalfWidth)
                v = i;
        }      
        
        //return bbinIdx ;
        return v ;
    }
    private int MINENTRIES = 100;
    private double CHI2CUT = 5;
    F1D f3 = new F1D("f3","[amp1]*gaus(x,[mean],[sigma1])+[amp2]*landau(x,[mean],[sigma2])", 0, 1.8);
    F1D f2 = new F1D("f2","[amp1]*gaus(x,[mean1],[sigma1])+[amp2]*gaus(x,[mean2],[sigma2])+[p02]", 0, 1.8);
    F1D f1 = new F1D("f1","[amp1]*gaus(x,[mean1],[sigma1])", 0, 2);
    
    private void filltrkDocavsTGraphs(int i, int j, int k) {
        Logger.getLogger("org.freehep.math.minuit").setLevel(Level.WARNING);
        
        if(TvstrkdocasProf.get(new Coordinate(i, j, k))!=null) {
            TvstrkdocasProf.get(new Coordinate(i, j, k)).reset();
        }
        if(Tvstrkdocas.get(new Coordinate(i, j, k))!=null) {
            H2F h2 = Tvstrkdocas.get(new Coordinate(i, j, k));
//            ParallelSliceFitter psf = new ParallelSliceFitter(h2);
//            psf.setBackgroundOrder(0);
//            psf.fitSlicesX();
//            GraphErrors meangr  = psf.getMeanSlices();
            double integ = h2.getEntries();
            ArrayList<H1F> hslice = h2.getSlicesX();
            int n=0;
            boolean saveG =false;
            for(int si=0; si<hslice.size(); si++) {
                double x = h2.getXAxis().getBinCenter(si);
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                double dmax = 2.*Constants.getInstance().wpdist[i];
                if(hslice.get(si).getMean()==0 || integ<200 || x>dmax) 
                    continue;
                H1F hn = hslice.get(si).histClone("h");
                int nbrebin=0;
                while( nbrebin<4 && amp<MINENTRIES) {
                    hn= this.rebin(hn.histClone("hnc"));
                    amp = hn.getBinContent(hn.getMaximumBin());
                    nbrebin++;
                }
                if(amp>MINENTRIES)
                    n=this.fitLandau(x, dmax, hn, TvstrkdocasProf.get(new Coordinate(i, j, k)));
                
                if(n>0 && x/dmax >0.8) 
                    saveG = true;  //ensure the large doca entries are filled to avoid biases at small docas
            }
            if(!useBProf && (i==2 || i==3) && k>0)
                saveG = false;
            if(!saveG) {
                TvstrkdocasProf.get(new Coordinate(i, j, k)).reset();
            }
        }
    }

    
    private int fitLandau(double x, double dmax, H1F hi, GraphErrors ge) {
        
        H1F h = hi.histClone("hc");
        double meanh = h.getDataX(h.getMaximumBin());
        double amp   = h.getBinContent(h.getMaximumBin());
        double sigma = h.getRMS(); 
        double mean = h.getMean();
        double binSize = h.getDataX(1)-h.getDataX(0);
        double min = meanh-2*sigma;
        double max1 = meanh+2*sigma;
        double max2 = meanh+3*sigma;
        
        if(min<0) min=0;
        f1.reset(); 
        f1.setRange(min, max1);
        f1.setParameter(0, amp);
        f1.setParameter(1, mean);
        f1.setParameter(2, sigma);
        f1.setParLimits(1, mean-2*binSize, mean-2*binSize);
        DataFitter.fit(f1, h,"LQ"); //No options uses error for sigma 
        double f1mu = f1.getParameter(1);
        double ch2f1 = f1.getChiSquare();
        f1.setParameter(1, meanh);
        f1.setParLimits(1, meanh-2*binSize, meanh+2*binSize);
        DataFitter.fit(f1, h,"LQ"); //No options uses error for sigma 
        if(f1.getChiSquare()<ch2f1)
           f1mu = f1.getParameter(1);
        
        f3.reset(); 
        f3.setRange(min, max2);
        f3.setParameter(0, amp/3.0);
        f3.setParameter(1, f1mu);
        f3.setParameter(2, sigma);
        f3.setParameter(3, amp);
        f3.setParameter(4, sigma);
        f3.setParLimits(0, 0, amp+10000);
        f3.setParLimits(3, 0, amp+10000);
        //f3.setParLimits(1, f1mu-sigma, f1mu+sigma);
        
        DataFitter.fit(f3, h,"LQ"); //No options uses error for sigma 
       
        int cnt =0;
        if(f3.getParameter(0)>0 && (f3.getChiSquare()/(double)f3.getNDF())<CHI2CUT && Math.abs(f3.getParameter(1)-f1mu)<4) {
            double mu = f3.getParameter(1);
            double emu = f3.parameter(1).error();
            if(x<0.05) emu = Math.max(mu, emu);
            if(x/dmax>0.9) emu = sigma;
            ge.addPoint(x,  mu, 0, emu);
            cnt++;
        } else {
            f1.setParameter(1, f1mu);
            f1.setRange(min, max1);
            DataFitter.fit(f1, h,"LQ"); //No options uses error for sigma 
            f1mu = f1.getParameter(1);
            f3.setRange(min, max1);
            f3.setParameter(1, f1mu);
            DataFitter.fit(f3, h,"LQ"); //No options uses error for sigma 
            if(f3.getParameter(0)>0 && (f3.getChiSquare()/(double)f3.getNDF())<CHI2CUT && Math.abs(f3.getParameter(1)-f1mu)<4) {
                double mu = f3.getParameter(1);
                double emu = f3.parameter(1).error();
                if(x<0.05) emu = Math.max(mu, emu);
                if(x/dmax>0.9) emu = sigma;
                ge.addPoint(x,  mu, 0, emu);
                cnt++;
            } else {
                ge.addPoint(x,  mean, 0, sigma);
                cnt++;
            }
            
        }
       
        return cnt;    
    }
    private H1F rebin(H1F h) {
        int nbins = h.getData().length;
        double binWidth = h.getDataX(1)-h.getDataX(0);
        double lowEdge = h.getDataX(0)-binWidth/2.0;
        double highEdge = h.getDataX(nbins-1)+binWidth/2.0;
        int j =0;
        for(int i = 0; i<nbins-1; i++) {
            if(i%2==0) j++;
        } 
       
        H1F hn = new H1F("hn", "","", (int) j, lowEdge, highEdge);
        j =0;
        for(int i = 0; i<nbins-1; i++) {
            if(i%2==0) {
                double y1 = h.getBinContent(i);
                double y2 = h.getBinContent(i+1);
                double e1 = h.getBinError(i);
                double e2 = h.getBinError(i+1);
                double y = y1+y2;
                double e = Math.sqrt(e1*e1+e2*e2);
                hn.setBinContent(j, y);
                hn.setBinError(j, e);
                j++;
            }
            
        }
        
        return hn;
    }
    private void filltrkDocavsTGraphs(int i, int j) {
        
        if(i<2 || i>3) { //region 1 and 3
            filltrkDocavsTGraphs(i, j, BBins);
        } else {
            //for(int k = 0; k < this.BBins; k++) {
            for(int k = 0; k < this.BBins; k++) {    
                filltrkDocavsTGraphs(i, j, k);
            }
        }       
    }


    int count = 0;
    public static int polarity =-1;
    public static double field =1;
    public List<FittedHit> calhits = new ArrayList<>();
    public List<FittedHit> hits = new ArrayList<>();
    Map<Integer, ArrayList<Integer>> segMapTBHits = new HashMap<Integer, ArrayList<Integer>>();
    Map<Integer, SegmentProperty> segPropMap = new HashMap<Integer, SegmentProperty>();
    Map<Integer, FittedHit> calhitmap = new HashMap<>();
    List<FittedHit> calhitlist = new ArrayList<>();
    Map<Integer, FittedHit> hitmap = new HashMap<>();
    List<FittedHit> hitlist = new ArrayList<>();
    Refit rf = new Refit();
    boolean refitSegs = false;
    @Override
    public void processEvent(DataEvent event) {
        hitmap.clear();
        calhitmap.clear();
        hitlist.clear();
        calhitlist.clear();
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
         
        //if(count>20000) return;
        if(count==1) {
            //Constants.getInstance().initialize("DCCAL");
            Driver.init();
            String newVar = String.valueOf(T2DViewer.calVariation.getSelectedItem());
            System.out.println("* VARIATION *"+newVar);
            ccdb.setVariation(newVar);
            TableLoader.t2dc=this;
            TableLoader.FillT0Tables(newRun, newVar);
            
            TableLoader.Fill(T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/t2d_pressure"),
                    T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/ref_pressure"),
                    T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/ref_pressure"));  
            
            this.loadFitPars(); 
            polarity = (int)Math.signum(event.getBank("RUN::config").getFloat("torus",0));
            field = event.getBank("RUN::config").getFloat("torus",0);
            runNumber = newRun;
             
            numberprocessedevents = Integer.parseInt(T2DViewer.enternofevents.getText());
            if (numberprocessedevents==-1)
            {
            	System.out.println("All events will be processed!!!!");
            }
            else {
            	System.out.println(numberprocessedevents + " events will be processed!!!!");
            }
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
     
        if (count > numberprocessedevents && numberprocessedevents!=-1) {
        	return;
        }
        
        
        // get segment property     
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        this.getSegProperty(bnkHits);
        
        for (int i = 0; i < bnkHits.rows(); i++) {
            double bFieldVal = (double) bnkHits.getFloat("B", i);
            int superlayer = bnkHits.getInt("superlayer", i);
            // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
            double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
            double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
            
            int alphaBin = this.getAlphaBin(alphaUncor); 
            if(alphaBin==-1) 
                continue;
            
            FittedHit theHit = this.getHit(bnkHits, i);
            if(this.passResiCuts(event, bnkHits, i)){//no previous entries
                if(hitmap.get(theHit.get_Id())==null) {
                    hitmap.put(theHit.get_Id(), theHit);
                }
            }
            
            if(this.passCalibCuts(event,bnkHits, i) && calhitmap.get(theHit.get_Id())==null){ 
                calhitmap.put(theHit.get_Id(), theHit);
                if(!this.refitSegs) {
                    
                    double calibTime = (double) (bnkHits.getInt("TDC", i) - bnkHits.getFloat("TProp", i)
                                            - bnkHits.getFloat("TFlight", i) - bnkHits.getFloat("TStart", i) 
                                            - bnkHits.getFloat("T0", i));
                    
                    Tvstrkdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, this.BBins))
                                    .fill(bnkHits.getFloat("trkDoca", i), calibTime);
                    
                    Tvscalcdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, this.BBins))
                                    .fill(bnkHits.getFloat("doca", i), calibTime);
                    
                    double yf = TvstrkdocasFits.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, this.BBins)).evaluate(bnkHits.getFloat("trkDoca", i));
                    if(!Double.isNaN(yf))
                        Tresvstrkdocas.get(new Coordinate(bnkHits.getInt("superlayer", i)  - 1, alphaBin, this.BBins))
                                        .fill(bnkHits.getFloat("trkDoca", i), calibTime-yf);
                    
                    //Fill region 2 for different b-field values
                    if(superlayer>2 && superlayer<5) { 
                        int bBin = this.getBBin(bFieldVal);
                        Tvstrkdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, bBin))
                                    .fill(bnkHits.getFloat("trkDoca", i), calibTime);
                        Tvscalcdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, bBin))
                                    .fill(bnkHits.getFloat("doca", i), calibTime);
                        double r2yf = TvstrkdocasFits.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, bBin)).evaluate(bnkHits.getFloat("trkDoca", i));
                        if(!Double.isNaN(r2yf))
                        Tresvstrkdocas.get(new Coordinate(bnkHits.getInt("superlayer", i)  - 1, alphaBin, bBin))
                                        .fill(bnkHits.getFloat("trkDoca", i), calibTime-r2yf);
                    }
                }
                // fill uncalibrated plot
                timeResiFromFile.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1))
                                .fill(bnkHits.getFloat("timeResidual", i));
                //Fill region 2 for different b-field values
                if(superlayer<3 || superlayer>4) { 
                    A.get(new Coordinate(superlayer-1, alphaBin, BBins))
                                .fill(alphaUncor);
                    //System.out.println("Filling alpha for spl "+superlayer+" alpha bin "+alphaBin+" alpha "+alphaUncor);
                }
                    
                    // fill B values histograms
                if(superlayer ==3 || superlayer ==4) {
                    int bBin = this.getBBin(bFieldVal);
                    B.get(new Coordinate(superlayer-3, alphaBin, bBin))
                            .fill(bFieldVal);
                    A.get(new Coordinate(superlayer-1, alphaBin, this.getBBin(bFieldVal)))
                            .fill(alphaUncor);
                    //System.out.println("Filling alpha for spl "+superlayer+" alpha bin "+alphaBin+" Bin "+this.getBBin(bFieldVal)+" alpha "+alphaUncor);
                }
            }
        }
        hitmap.forEach((k,v) -> hitlist.add(v));
        calhitmap.forEach((k,v) -> calhitlist.add(v));
        calhipoEvent = (HipoDataEvent) calwriter.createEvent();
        hipoEvent = (HipoDataEvent) writer.createEvent();
        //selected.show();
        calhipoEvent.appendBank(this.fillTBHitsBank(event, calhitlist));
        calhipoEvent.appendBank(event.getBank("RUN::config"));
        hipoEvent.appendBank(this.fillTBHitsBank(event, hitlist));
        calwriter.writeEvent(calhipoEvent);
        writer.writeEvent(hipoEvent);
        if(betacnt!=0) {
            setBetaAve(avebeta/(double)betacnt);
            avebeta=0;
            betacnt=0;
        }
        
        calhitlist.clear();
        hitlist.clear();
        
        
    }
    
    private double[][] resetPars = new double[6][11];
    private String[] parNames = {"v0", "vmid", "R", "tmax", "distbeta", "delBf", 
        "b1", "b2", "b3", "b4", "dmax"};
    
    double[] errs = {0.00001,0.00001,0.01,0.1,0.01,0.001,0.00001,0.00001,0.00001,0.00001,0.0000001};
    //private double[] errs = {0.00001,0.00001,0.001,0.1,0.001,0.001,0.00001,0.00001,0.00001,0.00001,0.0000001};
  //private double[] errs = {0.000001,0.000001,0.00001,0.01,0.0001,0.00001,0.00001,0.00001,0.00001,0.00001,0.0000001};
    
    public void resetPars() {
        TvstrkdocasFitPars.clear();
        for (int i = 0; i < this.nsl; i++) {
            double[] pars = resetPars[i];
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            for(int p = 0; p < 10; p++) {
                //TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p]);
                //create graphs of parameters for various iterations
               ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p+": " +parNames[p]+", superlayer "+(i+1),this.maxNfits+1, 0.5,this.maxNfits+1.5));
            }
            //TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10], errs[10]);
            TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10]);
        }
        fp.openFitPanel("fit panel", TvstrkdocasFitPars);
        
        reLoadFitPars();
        
        
        
    }
    public void loadFitPars() {
        for (int i = 0; i < this.nsl; i++) {
            double[] pars = new double[11];
            //T2DFunctions.polyFcnMac(x, alpha, bfield, v0[s][r], vmid[s][r], FracDmaxAtMinVel[s][r], 
            //tmax, dmax, delBf, Bb1, Bb2, Bb3, Bb4, superlayer) ;
            
            pars[0] = TableLoader.v0[0][i];
            pars[1] = TableLoader.vmid[0][i];
            pars[2] = TableLoader.FracDmaxAtMinVel[0][i];
            pars[3] = TableLoader.Tmax[0][i];
            pars[4] = TableLoader.distbeta[0][i];
            pars[5] = TableLoader.delta_bfield_coefficient[0][i];
            pars[6] = TableLoader.b1[0][i];
            pars[7] = TableLoader.b2[0][i];
            pars[8] = TableLoader.b3[0][i];
            pars[9] = TableLoader.b4[0][i];
            pars[10] = 2.*Constants.getInstance().wpdist[i];//fix dmax
            
            resetPars[i] = pars;
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            for(int p = 0; p < 10; p++) {
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                //create graphs of parameters for various iterations
                ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p, "superlayer "+(i+1)+" par "+p,this.maxNfits+1, 0.5,this.maxNfits+1.5));
            }
            TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10], errs[10]);
        }   
        // Fit panel
        fp = new FitPanel(this);
        fp.openFitPanel("fit panel", TvstrkdocasFitPars);
        
    }
    private void reLoadFitPars() {
        for (int s =0; s < 6; s++) {
            for (int i = 0; i < this.nsl; i++) {
                TableLoader.v0[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(0);
                TableLoader.vmid[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(1);
                TableLoader.FracDmaxAtMinVel[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(2);
                System.out.println("TMAX "+TableLoader.Tmax[s][i]+" ==> ");
                TableLoader.Tmax[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(3);
                System.out.println(".................... "+TableLoader.Tmax[s][i]+"");
                TableLoader.distbeta[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(4);
                TableLoader.delta_bfield_coefficient[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(5);
                TableLoader.b1[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(6);
                TableLoader.b2[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(7);
                TableLoader.b3[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(8);
                TableLoader.b4[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(9);
                
            }
        }
        for (int i = 0; i < this.nsl; i++) {
            for(int p = 0; p<6; p++) {
                ParsVsIter.get(new Coordinate(i,p)).setBinContent(0, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
                ParsVsIter.get(new Coordinate(i,p)).setBinError(0, TvstrkdocasFitPars.get(new Coordinate(i)).error(p));
            }
        }
        TableLoader.ReFill();
        
    }

    public void Plot(int i , int j) {
        DataLine l = new DataLine(0, 0, 2.5, 0);
        l.setLineStyle(2);
        l.setLineColor(2);
        System.out.println("Iteration num "+iterationNum);
        if(i<2 || i>3) { // regions 1 and 3 --> no b-field
            if(TvstrkdocasProf.get(new Coordinate(i, j, BBins)).getVectorX().size()>0) {
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T").cd(0);
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T").draw(Tvstrkdocas.get(new Coordinate(i, j, BBins)));
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").cd(0);
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").draw(Tvstrkdocas.get(new Coordinate(i, j, BBins)));
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").cd(0);
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").draw(Tvscalcdocas.get(new Coordinate(i, j, BBins)));
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").
                        draw(TvstrkdocasProf.get(new Coordinate(i, j, BBins)), "same");
                //TvstrkdocasFits.get(new Coordinate(i, j, BBins)).setRange(0, maxx[i]);
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").
                        draw(TvstrkdocasFits.get(new Coordinate(i, j, BBins)), "same");
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").
                        draw(TvstrkdocasFits.get(new Coordinate(i, j, BBins)), "same");
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").cd(0);
                GraphErrors g1 = new GraphErrors();
                GraphErrors g2 = new GraphErrors();
                g2.setTitle(TvstrkdocasInit.get(new Coordinate(i, j, BBins)).getTitle());
                g1.copy(TvstrkdocasProf.get(new Coordinate(i, j, BBins)));
                for(int ip =0; ip<g1.getVectorX().getSize(); ip++) {
                    if(g1.getDataEY(ip)!=0) {
                        double yf = TvstrkdocasFits.get(new Coordinate(i, j, BBins)).evaluate(g1.getDataX(ip));
                        double y = g1.getDataY(ip);
                        g2.addPoint(g1.getDataX(ip), y-yf, 0, g1.getDataEY(ip));
                        if(iterationNum==0) {
                            TvstrkdocasInit.get(new Coordinate(i, j, BBins)).addPoint(g1.getDataX(ip), y-yf, 0, 0);
                        }
                    }       
                }
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").clear();
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(Tresvstrkdocas.get(new Coordinate(i, j, BBins)));
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(g2, "Esame");
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(TvstrkdocasInit.get(new Coordinate(i, j, BBins)), "Esame");                   
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(l);
            }
        } else {   
            //plot the profiles for the various B-field components
            
            this.getAnalysisCanvas().getCanvas("TrackDoca vs T").cd(0);
            this.getAnalysisCanvas().getCanvas("TrackDoca vs T").draw(Tvstrkdocas.get(new Coordinate(i, j, BBins)));
            this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").cd(0);
            this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").draw(Tvstrkdocas.get(new Coordinate(i, j, BBins)));
            this.getAnalysisCanvas().getCanvas("CalcDoca vs T").cd(0);
            this.getAnalysisCanvas().getCanvas("CalcDoca vs T").draw(Tvscalcdocas.get(new Coordinate(i, j, BBins)));    
            
            for(int k = 0; k < this.BBins; k++) {
                if(TvstrkdocasProf.get(new Coordinate(i, j, k)).getVectorX().size()>0){
                    this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").
                            draw(TvstrkdocasProf.get(new Coordinate(i, j, k)), "same");
                    //TvstrkdocasFits.get(new Coordinate(i, j, k)).setRange(0, maxx[i]);
                    this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").
                            draw(TvstrkdocasFits.get(new Coordinate(i, j, k)), "same");
                    this.getAnalysisCanvas().getCanvas("CalcDoca vs T").
                    draw(TvstrkdocasFits.get(new Coordinate(i, j, k)), "same");
                }
            }
            GraphErrors g0 = new GraphErrors();
            g0.addPoint(0, 0, 0, 400);
            g0.setMarkerColor(0);
            g0.setLineColor(0);
            g0.setTitle(TvstrkdocasInit.get(new Coordinate(i, j, 1)).getTitle());
            this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").cd(0);
            this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(g0, "E");
            for(int k = 0; k < this.BBins; k++) {
                if(TvstrkdocasProf.get(new Coordinate(i, j, k)).getVectorX().size()>0){
                    GraphErrors g1 = new GraphErrors();
                    GraphErrors g2 = new GraphErrors();
                    g2.setMarkerColor(k+1);
                    g1.copy(TvstrkdocasProf.get(new Coordinate(i, j, k)));
                    for(int ip =0; ip<g1.getVectorX().getSize(); ip++) {
                        if(g1.getDataEY(ip)!=0) {
                            double yf = TvstrkdocasFits.get(new Coordinate(i, j, k)).evaluate(g1.getDataX(ip));
                            double y = g1.getDataY(ip);
                            g2.addPoint(g1.getDataX(ip), y-yf, 0, g1.getDataEY(ip));
                            if(iterationNum==0) {
                                TvstrkdocasInit.get(new Coordinate(i, j, k)).addPoint(g1.getDataX(ip), y-yf, 0, 0);
                            }
                        } 
                    }
                    
                    this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(g2, "Esame");
                    this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(TvstrkdocasInit.get(new Coordinate(i, j, k)), "Esame");
                    this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(l);
                }
            }
        }
    }
    
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_component = (String) cc.getValueAt(row, 2);
        System.out.println(str_sector + " " + str_layer + " " + str_component);
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);
       int component = Integer.parseInt(str_component);

       if(group.hasItem(sector,layer,component)==true){
           this.Plot(layer-1, component-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       
   
    }

    private FittedHit getCalHit(DataBank bnkHits, int i) {
        FittedHit hit = null;
        int id = bnkHits.getShort("id", i);;
        int sector = bnkHits.getByte("sector", i);
        int superlayer = bnkHits.getByte("superlayer", i);
        int layer = bnkHits.getByte("layer", i);
        int wire = bnkHits.getShort("wire", i);
        int TDC = bnkHits.getInt("TDC", i);
        double doca = bnkHits.getFloat("doca", i);
        double docaError = bnkHits.getFloat("docaError", i);
        double trkDoca = bnkHits.getFloat("trkDoca", i);
        int LR = bnkHits.getByte("LR", i);
        double X = bnkHits.getFloat("X", i);
        double Z = bnkHits.getFloat("Z", i);
        double B = bnkHits.getFloat("B", i);
        double TProp = bnkHits.getFloat("TProp", i);
        double TFlight = bnkHits.getFloat("TFlight", i);
        double T0 = bnkHits.getFloat("T0", i);
        double TStart = bnkHits.getFloat("TStart", i);
        int clusterID = bnkHits.getShort("clusterID", i);
        int trkID = bnkHits.getByte("trkID", i);
        double time = bnkHits.getFloat("time", i);
        double beta = bnkHits.getFloat("beta", i);
        double tBeta = bnkHits.getFloat("tBeta", i);
        double resiTime = bnkHits.getFloat("timeResidual", i);
        double resiFit = bnkHits.getFloat("fitResidual", i);
        int jitter = (int) bnkHits.getByte("jitter", i);
        double dDoca = bnkHits.getFloat("dDoca", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        //double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        //double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double alpha = bnkHits.getFloat("Alpha", i);

        hit = new FittedHit(sector, superlayer, layer, wire, TDC, jitter, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
        hit.setB(B);
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(T2DViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(T2DViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
        hit.setAlpha(alpha);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        hit.set_DeltaDocaBeta(dDoca);
        //this.getSegProperty(bnkHits);
        //if(this.passCalibCuts(bnkHits, i)) {            
            return hit;
        //} else {
        //    return null;
        //}
    }
    private FittedHit getHit(DataBank bnkHits, int i) {
        FittedHit hit = null;
        int id = bnkHits.getShort("id", i);;
        int sector = bnkHits.getByte("sector", i);
        int superlayer = bnkHits.getByte("superlayer", i);
        int layer = bnkHits.getByte("layer", i);
        int wire = bnkHits.getShort("wire", i);
        int TDC = bnkHits.getInt("TDC", i);
        double doca = bnkHits.getFloat("doca", i);
        double docaError = bnkHits.getFloat("docaError", i);
        double trkDoca = bnkHits.getFloat("trkDoca", i);
        int LR = bnkHits.getByte("LR", i);
        double X = bnkHits.getFloat("X", i);
        double Z = bnkHits.getFloat("Z", i);
        double B = bnkHits.getFloat("B", i);
        double TProp = bnkHits.getFloat("TProp", i);
        double TFlight = bnkHits.getFloat("TFlight", i);
        double T0 = bnkHits.getFloat("T0", i);
        double TStart = bnkHits.getFloat("TStart", i);
        int clusterID = bnkHits.getShort("clusterID", i);
        int trkID = bnkHits.getByte("trkID", i);
        double time = bnkHits.getFloat("time", i);
        double beta = bnkHits.getFloat("beta", i);
        double tBeta = bnkHits.getFloat("tBeta", i);
        double resiTime = bnkHits.getFloat("timeResidual", i);
        double resiFit = bnkHits.getFloat("fitResidual", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        //double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        //double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double alpha = bnkHits.getFloat("Alpha", i);
        int jitter = (int) bnkHits.getByte("jitter", i);
        double dDoca = bnkHits.getFloat("dDoca", i);
        hit = new FittedHit(sector, superlayer, layer, wire, TDC, jitter, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
        hit.setB(B);
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(T2DViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(T2DViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
        hit.setAlpha(alpha);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        hit.set_DeltaDocaBeta(dDoca);
        //this.getSegProperty(bnkHits);
        //double dmax = 2.*Constants.getInstance().wpdist[superlayer-1]; 
        //double cos30minusalpha = Math.cos(Math.toRadians(30.-util.getReducedAngle(alpha)));
        return hit;
    }
    private DataBank fillTBHitsBank(DataEvent event, List<FittedHit> hitlist) {
        //if(event.hasBank("TimeBasedTrkg::TBHits")) { // for second pass tracking
        //     event.removeBank("TimeBasedTrkg::TBHits");
        //}
        DataBank bank = event.createBank("TimeBasedTrkg::TBHits", hitlist.size());

        for (int i = 0; i < hitlist.size(); i++) {
            bank.setShort("id", i, (short) hitlist.get(i).get_Id());
            bank.setShort("status", i, (short) hitlist.get(i).get_QualityFac());
            bank.setByte("superlayer", i, (byte) hitlist.get(i).get_Superlayer());
            bank.setByte("layer", i, (byte) hitlist.get(i).get_Layer());
            bank.setByte("sector", i, (byte) hitlist.get(i).get_Sector());
            bank.setShort("wire", i, (short) hitlist.get(i).get_Wire());

            bank.setFloat("X", i, (float) hitlist.get(i).get_X());
            bank.setFloat("Z", i, (float) hitlist.get(i).get_Z());
            bank.setByte("LR", i, (byte) hitlist.get(i).get_LeftRightAmb());
            bank.setFloat("time", i, (float) (hitlist.get(i).get_Time())); //time is the fully corrected time
            bank.setFloat("tBeta", i, (float) hitlist.get(i).get_DeltaTimeBeta());
            bank.setFloat("fitResidual", i, (float) hitlist.get(i).get_Residual());
            bank.setFloat("Alpha", i, (float) hitlist.get(i).getAlpha());
            
            bank.setFloat("doca", i, (float) hitlist.get(i).get_Doca());
            bank.setFloat("docaError", i, (float) hitlist.get(i).get_DocaErr());
            bank.setFloat("trkDoca", i, (float) hitlist.get(i).get_ClusFitDoca());

            bank.setShort("clusterID", i, (short) hitlist.get(i).get_AssociatedClusterID());
            bank.setByte("trkID", i, (byte) hitlist.get(i).get_AssociatedHBTrackID());
            bank.setFloat("timeResidual", i, (float) hitlist.get(i).get_TimeResidual());
            
            bank.setInt("TDC",i,hitlist.get(i).get_TDC());
            bank.setFloat("B", i, (float) hitlist.get(i).getB());
            bank.setFloat("TProp", i, (float) hitlist.get(i).getTProp());
            bank.setFloat("TFlight", i, (float) hitlist.get(i).getTFlight());
            bank.setFloat("T0", i, (float) hitlist.get(i).getT0());
            bank.setFloat("TStart", i, (float) hitlist.get(i).getTStart());
            bank.setFloat("dDoca", i, (float) hitlist.get(i).get_DeltaDocaBeta());
            bank.setFloat("beta", i, (float) hitlist.get(i).get_Beta());
            
        }
        //System.out.println(" Created Bank "); bank.show();
        return bank;

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
                    double chi2pid = bank2.getFloat("chi2pid", bank.getShort("pindex", i));
                    if(Math.abs(chi2pid)>10) pid =0;
                }
            }
        }
        
        return pid;
    } 
    //method by Raffaella
    private boolean isElectronZeroField(DataEvent event, int trackId) {
       
        if(!event.hasBank("REC::Particle") ||
           !event.hasBank("REC::Calorimeter") ||
           !event.hasBank("REC::Cherenkov") ||
           !event.hasBank("REC::Track"))
            return false;
        
        DataBank particleBank    = event.getBank("REC::Particle");
        DataBank calorimeterBank = event.getBank("REC::Calorimeter");
        DataBank cherenkovBank   = event.getBank("REC::Cherenkov");
        DataBank trackBank       = event.getBank("REC::Track");
         
        int pindex = -1;
        for (int i=0; i<trackBank.rows(); i++) {
            if(trackBank.getShort("index", i)==trackId-1) {
                pindex = trackBank.getShort("pindex", i);
                break;
            }
        }    
        
        if(pindex>=0) {
            double beta    = particleBank.getFloat("beta", pindex);
            double vtx     = particleBank.getFloat("vz", pindex);

            double nphe = 0;
            for(int i=0; i<cherenkovBank.rows(); i++) {
                if(cherenkovBank.getShort("pindex", i)==pindex) {
                    nphe = cherenkovBank.getFloat("nphe", i);
                    break;
                }
            }

            double energy = 0;
            for (int i=0; i<calorimeterBank.rows(); i++) {
                if(calorimeterBank.getShort("pindex", i)==pindex) {
                    energy+=calorimeterBank.getFloat("energy", i);
                }
            }

            if(beta>0 && nphe>2 && energy>0.5) { 
                return true;
            }
        }    
        return false;
    }    
    
    private void updateHit(FittedHit hit, boolean flagOT) { 
        double distbeta = TvstrkdocasFitPars.get(new Coordinate(hit.get_Superlayer()-1)).value(4);
        double v0 = TvstrkdocasFitPars.get(new Coordinate(hit.get_Superlayer()-1)).value(0);
        double d = hit.get_ClusFitDoca();
        double beta = hit.get_Beta();
        if(beta>1.0) {
            beta=1.0;
        }
        //beta = T2DCalib.getBetaAve();
        double calibTime = (double) (hit.get_TDC() - hit.getTProp()
                                            - hit.getTFlight() - hit.getTStart()
                                            - hit.getT0());
  
        double deltatime_beta = util.getDeltaTimeBeta(d, beta, distbeta, v0);
        hit.set_DeltaTimeBeta(deltatime_beta);
        double deltadoca_beta = 0;
       
        hit.set_Doca(this.timeToDistance(hit.get_Sector(), hit.get_Superlayer(), 
                    hit.getAlpha(), hit.get_Beta(), hit.getB(), calibTime, 0));
        double ralpha = (double) util.getReducedAngle(hit.getAlpha());
        if(flagOT) {
            double calctime = calc_Time( hit.get_Doca(),  ralpha, hit.getB(), hit.get_Sector(), hit.get_Superlayer()) ;
            double deltatimebeta = util.getDeltaTimeBeta(hit.get_Doca(),beta,distbeta, v0);
            calctime+=deltatimebeta;
            if(calibTime-calctime>T2DCalib.DeltaTimeCut) {
                //System.out.println(hit.printInfo()+" alpha "+hit.getAlpha()+" ralpha "+ralpha+ " tDOCA "+d+ " DOCA "+hit.get_Doca()+ " calibTime "+calibTime+" calctime "+calctime);
                hit.set_OutOfTimeFlag(true);
            }
        }
        double x = hit.get_XWire(); 
        double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
        //fix alpha to get the local angle
        double alphaRadUncor = Math.toRadians(hit.getAlpha()+(double)T2DCalib.polarity*theta0);
        double trkAngle = Math.tan(alphaRadUncor);
        double cosTkAng = 1./Math.sqrt(trkAngle*trkAngle + 1.);
        hit.set_X(x + hit.get_LeftRightAmb() * (hit.get_Doca() / cosTkAng) );
    }
    
    
    private TimeToDistanceEstimator tdee = new TimeToDistanceEstimator();
    private double timeToDistance(int sector, int superlayer, double alpha, double beta, double B, double time,
                                double deltadoca_beta) {
        double distance = 0;
        
        //reduce the corrected angle
        double ralpha = (double) util.getReducedAngle(alpha);
         
        if(beta>1.0)
            beta = 1.0;
        
        double correctedTime = time;
        if(correctedTime<=0)
            return 0; // fixes edge effects ... to be improved
        
        distance = tdee.interpolateOnGrid(B, ralpha, beta,
                    correctedTime, 
                    sector-1, superlayer-1) ; 
        
        return distance;
    }
    private void UpdateAlphaBinCenters() {
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k <= this.BBins; k++) {
                    AlphaValuesUpd[i][j][k] = AlphaValues[j];
                    if(A.get(new Coordinate(i,j,k)).getBinContent(A.get(new Coordinate(i,j,k)).getMaximumBin())>0) {
                        AlphaValuesUpd[i][j][k] = A.get(new Coordinate(i,j,k)).getMean();
                    System.out.println("ijk" +i+" "+j+" "+k+" Alpha Bin UPdated "+AlphaValues[j]+" --> "+AlphaValuesUpd[i][j][k] );
                    }
                }
            }
        }
    }
    boolean filledBspectra = false;
    private void UpdateBBinCenters() {
        if(field==0) return;
         if(filledBspectra) return;
        TCanvas can1 = new TCanvas("superlayer3 B", 800, 800);
        TCanvas can2 = new TCanvas("superlayer4 B", 800, 800);
        
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k < this.BBins; k++) {
                    BfieldValuesUpd[i][j][k] = BfieldValues[k];
                    if(B.get(new Coordinate(i,j,k)).getBinContent(B.get(new Coordinate(i,j,k)).getMaximumBin())>10) {
                        BfieldValuesUpd[i][j][k] = B.get(new Coordinate(i,j,k)).getMean();
                    System.out.println("ijk" +i+" "+j+" "+k+" BBin UPdated "+BfieldValues[k]+" --> "+BfieldValuesUpd[i][j][k] );
                    }
                }
            }
        }
        can1.divide(7, 3);
        can2.divide(7, 3);
        int can1idx=0;
        int can2idx=0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                BAlphaBins.get(new Coordinate(i,j)).setOptStat(0);
                BAlphaBins.get(new Coordinate(i,j)).setTitle("alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")");
                
                for (int k = 0; k < this.BBins; k++) {
                    if(B.get(new Coordinate(i,j,k)).getBinContent(B.get(new Coordinate(i,j,k)).getMaximumBin())>1)
                        BAlphaBins.get(new Coordinate(i,j)).add(B.get(new Coordinate(i,j,k)));
                }    
                if(i==0) {
                        can1.cd(can1idx);
                        can1.draw(BAlphaBins.get(new Coordinate(i,j)));
                        can1idx++;
                }
                if(i==1) {
                    can2.cd(can2idx);
                    can2.draw(BAlphaBins.get(new Coordinate(i,j)));
                    can2idx++;
                }
            }
        }
        filledBspectra=true;
    }
    
    

    private boolean selectOnAlpha(int superlayer, double alphaRadUncor) {
        boolean pass = false;
        int i = superlayer - 1;
        if(alphaRadUncor>Double.parseDouble(T2DViewer.alphaCuts1[i].getText()) &&
                alphaRadUncor<Double.parseDouble(T2DViewer.alphaCuts2[i].getText())) {
            pass = true;
        }
        return pass;        
    }

    private void getSegProperty(DataBank bnkHits) {
        
        segMapTBHits.clear();
        segPropMap.clear();
               
        for (int j = 0; j < bnkHits.rows(); j++){
            Integer cID = bnkHits.getInt("clusterID", j);
            if(segMapTBHits.containsKey(cID)==false) {
                segMapTBHits.put(cID, new ArrayList<Integer>());
            }
            
            segMapTBHits.get(cID).add((int)bnkHits.getShort("wire", j));
            
        }
        
        Iterator<Map.Entry<Integer, ArrayList<Integer>>> itr = segMapTBHits.entrySet().iterator(); 
          
        while(itr.hasNext()) { 
            Map.Entry<Integer, ArrayList<Integer>> entry = itr.next(); 
            segPropMap.put(entry.getKey() , 
                    new SegmentProperty(entry.getKey(),entry.getValue(),Integer.parseInt(T2DViewer.deltaWire.getText())));
        } 

    }

    

    private boolean passPID(DataEvent event, DataBank bnkHits, int rowIdxinTrkBank) {
        boolean pass = false;
        int trkID = bnkHits.getByte("trkID", rowIdxinTrkBank);
        if(polarity==0) return this.isElectronZeroField(event, trkID);
        int pid = this.readPID(event, trkID);
        //pass if the track is identified as an electron or as a hadron
        //if(pid==11 || Math.abs(pid)==211 || Math.abs(pid)==2212 || Math.abs(pid)==321) {
        if(Integer.parseInt(T2DViewer.pid.getText())==-1) {
            pass=true;
        }
        if(pid==Integer.parseInt(T2DViewer.pid.getText())) {    
            pass = true;
        }
        
        return pass;
    }
    
    private boolean passResiCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        
        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("TFlight", i)>0 
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && this.passPID(event, bnkHits, i)==true)
            {
                pass = true;
            }
        return pass;
    }
    boolean betaLoaded=false;
    double avebeta = 0;
    int betacnt =0;
    private boolean passCalibCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        if(betaLoaded==false) {
            double betaLow  = Double.parseDouble(T2DViewer.betaCut.getText());
            double betaHigh = Double.parseDouble(T2DViewer.betaCut2.getText());   
            //setBetaAve(0.99);
           // System.out.println("AVERAGE BETA = "+getBetaAve());
            betaLoaded = true;
        }
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double beta = (double) bnkHits.getFloat("beta", i);
        
        if (bnkHits.getByte("trkID", i) >0 
                    && beta> Double.parseDouble(T2DViewer.betaCut.getText()) 
                    && beta< Double.parseDouble(T2DViewer.betaCut2.getText()) 
                    && this.selectOnAlpha(superlayer, alphaRadUncor)==true
                    && bnkHits.getFloat("TFlight", i)>0 
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getNumWireWithinDW()<=Integer.parseInt(T2DViewer.npassWires.getText())
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()>Integer.parseInt(T2DViewer.nWires.getText())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && this.passPID(event, bnkHits, i)==true
                )
            {
                avebeta+=beta;
                betacnt++;
                pass = true;
            }
        return pass;
    } 

    /**
     * @return the betaAve
     */
    public static double getBetaAve() {
        return betaAve;
    }

    /**
     * @param aBetaAve the betaAve to set
     */
    public static void setBetaAve(double aBetaAve) {
        betaAve = aBetaAve;
    }
    
    private void checkFile(String testCalOutPuthipo) {
        File file = new File(testCalOutPuthipo);
        if (file.exists()) {
            // Delete the file
            if (file.delete()) {
                System.out.println("File deleted successfully.");
            } else {
                System.out.println("Failed to delete the file.");
            }
        } else {
            System.out.println("File does not exist.");
        }
        
    }

    public void resetPars(int i, boolean[][] fixFit) {
        // reset
        TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
        for (int p = 0; p < 10; p++) {
            if(fixFit[p][i]==true) {
                TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
            }
        }
    }
    
}

