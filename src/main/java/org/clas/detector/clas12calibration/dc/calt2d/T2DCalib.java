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
import org.freehep.math.minuit.MnScan;
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

    private static double DeltaTimeCut = 20;
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
    private Map<Coordinate, H2F> Tvstrkdocas                = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tvscalcdocas               = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, GraphErrors> TvstrkdocasProf    = new HashMap<Coordinate, GraphErrors>();
    private final Map<Coordinate, GraphErrors> TvstrkdocasInit    = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> TvstrkdocasFit             = new HashMap<Coordinate, FitFunction>();
    public Map<Coordinate, MnUserParameters> TvstrkdocasFitPars    = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> TvstrkdocasFits                = new HashMap<Coordinate, FitLine>();
    private Map<Coordinate, H1F> timeResi               = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiFromFile       = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiNew            = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> fitResi        = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> B              = new HashMap<Coordinate, H1F>(); //histogram to get B values centroids
    private Map<Coordinate, H1F> A              = new HashMap<Coordinate, H1F>(); //histogram to get Alpha values centroids
    private Map<Coordinate, H1F> BAlphaBins              = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> ParsVsIter    = new HashMap<Coordinate, H1F>();
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
    double docaBinWidth = 0.015;
    double timeBinWidth = 2.0;
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
                BAlphaBins.put(new Coordinate(i,j), new H1F("A  " , 100, -36, 36));
                for (int k = 0; k <= this.BBins; k++) {
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
    private int maxIter = 10;
    double[][] fixPars = new double[6][10];
    
    private MnScan  scanner = null;
    private MnMigrad migrad = null;
    
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
    public void runFit(int i, boolean fixFit[][]) {
        // i = superlayer - 1;
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** ");
        TvstrkdocasFit.put(new Coordinate(i), 
                new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
        
        scanner = new MnScan((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                TvstrkdocasFitPars.get(new Coordinate(i)),2);
	
        scanner.fix(10);
        //v0
        scanner.setLimits(0, 0.002, 0.009);
        //vmid
        scanner.setLimits(1, 0.001, 0.005);
        //R
        scanner.setLimits(2, 0.6, 0.72);
        //distbeta
        scanner.setLimits(4, 0.02, 0.1);
        //delBf
        scanner.setLimits(5, 0.10, 0.4);
        
        for (int p = 0; p < 10; p++) {
            if(fixFit[p][i]==true) {
                scanner.fix(p); 
            }
        }
        FunctionMinimum scanmin = scanner.minimize();
        if(scanmin.isValid())
            TvstrkdocasFitPars.put(new Coordinate(i),scanmin.userParameters());
        
        migrad = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                TvstrkdocasFitPars.get(new Coordinate(i)),1);
        migrad.setCheckAnalyticalDerivatives(true);
        
        FunctionMinimum min = null ;
        
        
        for(int it = 0; it<maxIter; it++) {
            
            min = migrad.minimize();
           // System.err.println("****************************************************");
           // System.err.println("*   FIT RESULTS  FOR SUPERLAYER  "+(i+1)+" at iteration "+(it+1)+"  *");
           // System.err.println("****************************************************");  
           // for(int pi = 0; pi<6; pi++) 
           //     System.out.println("par["+pi+"]="+(TvstrkdocasFitPars.get(new Coordinate(i)).value(pi)-min.userParameters().value(pi)));
            
            if(min.isValid()) {
                TvstrkdocasFitPars.put(new Coordinate(i),min.userParameters());  
            }
            
            
            //System.err.println(min);
            
            
        }
        if(min!=null) {
            String t = "FINAL MINUIT FIT RESULT FOR SUPERLAYER ";
            t+=(i+1);
            t+=" :";
            System.out.println(t);
            String s = String.valueOf(min);  System.out.println(s);
            pw3.println(t);
            pw3.println(s);
        }

        //release all so they can be fixed again
        TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
        for (int p = 0; p < 10; p++) {
            if(fixFit[p][i]==true) {
                TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
            }
        }
        if(i==5) {
            System.err.println("***************************************************");
            System.err.println("*                  FITS ARE DONE!                 *");
            System.err.println("***************************************************");  
            
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
                    //TvstrkdocasProf.get(new Coordinate(i,j,k)).reset();
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
                    updateHit(hit, false); //update the hit with previous cal results
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
                    double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
                    double alphaUncor = hit.getAlpha()+(double)T2DCalib.polarity*theta0;
                    int alphaBin = this.getAlphaBin(alphaUncor);
                    double bFieldVal = (double) hit.getB();
                    
                    if(alphaBin!=-1 && !hit.get_OutOfTimeFlag()) {
                        double calibTime = (double) (hit.get_TDC() - hit.getTProp()
                                            - hit.getTFlight() - hit.getTStart()
                                            - hit.getT0()); 
                        
                        Tvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, this.BBins))
                                        .fill(hit.get_ClusFitDoca(), calibTime);
                        Tvscalcdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, this.BBins))
                                        .fill(hit.get_Doca(), calibTime);
                        
                        //Fill region 2 for different b-field values
                        if(hit.get_Superlayer() >2 && hit.get_Superlayer() <5) { 
                            int bBin = this.getBBin(bFieldVal);
                            Tvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, bBin))
                                        .fill(hit.get_ClusFitDoca(), calibTime);
                            Tvscalcdocas.get(new Coordinate(hit.get_Superlayer() - 1, alphaBin, bBin))
                                        .fill(hit.get_Doca(),calibTime);
                            
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
        
        
        int effSig = (int) (10000*(gaus1Func.getParameter(2)));
        String t = "Eff. Sig. "+effSig + "microns";
        h1.setTitle(t);
        canvasRes.draw(h1, "E1");
        //canvasRes.draw(gausFunc, "same");
        
        
    }
    public int getAlphaBin(double alpha) {
        int v = -1;
        for(int i = 0; i<this.AlphaValues.length; i++) {
            
            if(Math.abs(alpha-this.AlphaValues[i])<this.AlphaBinHalfWidth)
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
            if(Math.abs(bFieldVal*bFieldVal-this.BfieldValues[i]*this.BfieldValues[i])<BSqrBinHalfWidth)
                v = i;
        }      
        
        //return bbinIdx ;
        return v ;
    }
    private int MINENTRIES = 5;
    
    F1D f2 = new F1D("f2","[amp1]*gaus(x,[mean1],[sigma1])+[amp2]*gaus(x,[mean2],[sigma2])+[p02]", 0, 1.8);
    F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])+[p0]", 0, 1.8);
    private void filltrkDocavsTGraphs(int i, int j, int k) {

        Logger.getLogger("org.freehep.math.minuit").setLevel(Level.WARNING);

        if(TvstrkdocasProf.get(new Coordinate(i, j, k))!=null) {
            TvstrkdocasProf.get(new Coordinate(i, j, k)).reset();
            H2F h2 = Tvstrkdocas.get(new Coordinate(i, j, k));
//            ParallelSliceFitter psf = new ParallelSliceFitter(h2);
//            psf.setBackgroundOrder(0);
//            psf.fitSlicesX();
//            GraphErrors meangr  = psf.getMeanSlices();
            ArrayList<H1F> hslice = h2.getSlicesX();
            for(int si=0; si<hslice.size(); si++) {
                double x = h2.getXAxis().getBinCenter(si);
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                double dmax = 2.*Constants.getInstance().wpdist[i];
                if(hslice.get(si).getMean()==0 || amp<this.MINENTRIES) 
                    continue;
               
                if( x/dmax<0.1 ) { // for large docas (80% of dmax take fit with a double gauss )
                    this.fitSingleGaus(si, x, hslice.get(si), TvstrkdocasProf.get(new Coordinate(i, j, k)));
                } else {
                    this.fitDoubleGaus(si, x,hslice.get(si), TvstrkdocasProf.get(new Coordinate(i, j, k)));
                }
            }
        }
    }

    private void fitSingleGaus(int si, double x, H1F h, GraphErrors ge) {
        double mean = h.getMean(); 
        double meanh = h.getDataX(h.getMaximumBin());
        double amp   = h.getBinContent(h.getMaximumBin());
        double sigma = h.getRMS(); 

        double min = mean-3*sigma;
        double max = mean+3*sigma;
        
        f1.setRange(min, max);
        f1.setParameter(0, amp);
        f1.setParameter(1, mean);
        f1.setParameter(2, sigma);
        f1.setParameter(3, 0);
        f1.setParLimits(0, this.MINENTRIES, amp+10);
        f1.setParLimits(1, min, max);
        
        DataFitter.fit(f1, h,"LQ"); //No options uses error for sigma 
        double mu = meanh;
        double emu= sigma;
        if(f1.isFitValid() && f1.getParameter(1)>0
                && f1.parameter(1).error()<sigma
                && Math.abs(meanh-f1.getParameter(1))<sigma && f1.getParameter(1)>0) {
            
            if(Math.abs(mean-mu)<sigma && f1.getChiSquare()<500)
                mu = f1.getParameter(1);
            emu = f1.parameter(1).error();
            if(si==0 || si==1) 
                emu = 0;
        }
        ge.addPoint(x, mu, 0, emu);
    }

    private void fitDoubleGaus(int si, double x, H1F h, GraphErrors ge) {
        double mean = h.getMean();
        double meanh = h.getDataX(h.getMaximumBin());
        double amp   = h.getBinContent(h.getMaximumBin());
        double sigma = h.getRMS(); 

        double min = mean-3*sigma;
        double max1 = mean+3*sigma;
        double max2 = mean+5*sigma;
        if(min<0) min=0;
        f2.setRange(min, max2);
        f2.setParameter(0, amp);
        f2.setParameter(1, mean);
        f2.setParameter(2, sigma);
        f2.setParameter(3, amp/3.0);
        f2.setParameter(4, mean);
        f2.setParameter(5, sigma);
        f2.setParameter(6, 0);
        f2.setParLimits(0, this.MINENTRIES, amp+10);
        f2.setParLimits(1, min, max1);
        f2.setParLimits(4, min, max1);

        DataFitter.fit(f2, h,"LQ"); //No options uses error for sigma 
        double mu = meanh;
        double emu= sigma;
        if(f2.isFitValid() && f2.getParameter(1)>0 && f2.getParameter(4)>0
                && f2.parameter(1).error()<sigma && f2.parameter(4).error()<sigma
                && f1.getParameter(1)>0) {
           
            double mutest = (f2.getParameter(0)*f2.getParameter(1)+f2.getParameter(3)*f2.getParameter(4))/(f2.getParameter(0)+f2.getParameter(3));
            
            if(Math.abs(meanh-mutest)<sigma && f1.getChiSquare()<500)
                mu = mutest;
            emu = Math.sqrt(f2.parameter(1).error()*f2.parameter(1).error()+f2.parameter(4).error()*f2.parameter(4).error());
        }
        ge.addPoint(x, mu, 0, emu);
    }
    
    private void filltrkDocavsTGraphs(int i, int j) {
        if(i<2 || i>3) { //region 1 and 3
            filltrkDocavsTGraphs(i, j, BBins);
        } else {
            for(int k = 0; k < this.BBins; k++) {
                filltrkDocavsTGraphs(i, j, k);
            }
        }       
    }


    int count = 0;
    public static int polarity =-1;
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
                    //Fill region 2 for different b-field values
                    if(superlayer>2 && superlayer<5) { 
                        int bBin = this.getBBin(bFieldVal);
                        Tvstrkdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, bBin))
                                    .fill(bnkHits.getFloat("trkDoca", i), calibTime);
                        Tvscalcdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1, alphaBin, bBin))
                                    .fill(bnkHits.getFloat("doca", i), calibTime);
                        
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
    private double[] errs = {0.00001,0.000001,0.0001,1.0,0.001,0.0001,0.001,0.001,0.001,0.001,0.00001};
    
    public void resetPars() {
        TvstrkdocasFitPars.clear();
        for (int i = 0; i < this.nsl; i++) {
            double[] pars = resetPars[i];
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            for(int p = 0; p < 10; p++) {
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                
                //create graphs of parameters for various iterations
               ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p+": " +parNames[p]+", superlayer "+(i+1),this.maxIter+1, 0.5,this.maxIter+1.5));
            }
            TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10], errs[10]);
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
                ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p, "superlayer "+(i+1)+" par "+p,this.maxIter+1, 0.5,this.maxIter+1.5));
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
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Fit Resi").draw(g2, "E");
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
        if(flagOT) {
            double calctime = calc_Time( hit.get_Doca(),  hit.getAlpha(), hit.getB(), hit.get_Sector(), hit.get_Superlayer()) ;
            double deltatimebeta = util.getDeltaTimeBeta(hit.get_Doca(),beta,distbeta, v0);
            calctime+=deltatimebeta;
            if(calibTime-calctime>T2DCalib.DeltaTimeCut) 
                hit.set_OutOfTimeFlag(true);
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
            correctedTime=0.01; // fixes edge effects ... to be improved
        
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
    private void UpdateBBinCenters() {
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
        System.out.println("Number of alpha bins "+this.alphaBins);
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
        if(polarity==0) return true;
        int trkID = bnkHits.getByte("trkID", rowIdxinTrkBank);
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

    
}

