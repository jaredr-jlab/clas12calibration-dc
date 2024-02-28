/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
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
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.analysis.FitPanel;
import org.clas.detector.clas12calibration.dc.t2d.TimeToDistanceEstimator;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
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
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.graphics.EmbeddedCanvas;
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
    public HipoDataSync calwriter = null;
    public HipoDataSync writer = null;
    private HipoDataEvent calhipoEvent = null;
    private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    FitPanel fp;
    PrintWriter pw = null;
    File outfile = null;
    private int runNumber;
    private Utilities util = new Utilities();
    private int numberprocessedevents;
    
    public T2DCalib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames("TrackDoca vs T","TrackDoca vs T Graphs","CalcDoca vs T","Time Residuals","Parameters");
        this.init(false, "v0:vmid:R:tmax:distbeta:delBf:b1:b2:b3:b4");
        //outfile = new File("Files/ccdbConstants.txt");
        //pw = new PrintWriter(outfile);
        //pw.printf("#& sector superlayer component v0 deltanm tmax distbeta delta_bfield_coefficient b1 b2 b3 b4 delta_T0 c1 c2 c3\n");
        
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
        
    }
    private Map<Coordinate, H2F> Tvstrkdocas                = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tvscalcdocas               = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, GraphErrors> TvstrkdocasProf    = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> TvstrkdocasFit             = new HashMap<Coordinate, FitFunction>();
    private Map<Coordinate, MnUserParameters> TvstrkdocasFitPars    = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> TvstrkdocasFits                = new HashMap<Coordinate, FitLine>();
    private Map<Coordinate, H1F> timeResi               = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiFromFile       = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiNew            = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> fitResi        = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> B              = new HashMap<Coordinate, H1F>(); //histogram to get B values centroids
    private Map<Coordinate, H1F> ParsVsIter    = new HashMap<Coordinate, H1F>();
    int nsl = 6;

    public static double[] BfieldValues = new double[]{0.707106781,1.224744871,1.58113883,1.87082869,2.121320344,2.34520788,2.549509757,2.738612788};   
    public static double[] AlphaValues = new double[]{-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26};
    double AlphaBinHalfWidth = 2;
    public static int alphaBins = AlphaValues.length;
    public static int BBins = BfieldValues.length;
    //update middle of the B bins
    public static double[][][] BfieldValuesUpd = new double[2][alphaBins][BBins];
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup td = new DataGroup(7,2);
        DataGroup tdp = new DataGroup(14,8);
        DataGroup cd = new DataGroup(7,2);
        DataGroup tr = new DataGroup(6,1);
        DataGroup fr = new DataGroup(6,1);
        
        int ijk = 0;
        int ij = 0;
        for (int i = 0; i < nsl; i++) {
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            timeResi.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            timeResiFromFile.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            timeResiNew.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            fitResi.put(new Coordinate(i), new H1F("fit residual for sly " + (i+1), 100, -0.3, 0.3));
            
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
                            +", B "+k, 200, 0, 2.0, 200, 0, 500.0+(int)(i/2)*450.0));
                    
                    TvstrkdocasProf.put(new Coordinate(i,j,k), new GraphErrors());
                    TvstrkdocasProf.get(new Coordinate(i,j,k)).setMarkerColor(k+1);
                    
                    Tvscalcdocas.put(new Coordinate(i,j,k), new H2F("calcDocavsT" + (i + 1)*1000+(j+1)+26, "superlayer" + (i + 1)
                            + ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")"
                            +", B "+k, 200, 0, 2.0, 200, 0, 500.0+(int)(i/2)*450.0));
                    tdp.addDataSet(TvstrkdocasProf.get(new Coordinate(i,j,k)), ijk);
                    //trkdvst.addDataSet(Tvstrkdocas.get(new Coordinate(i,j,k)), 0);
                    prfdvst.addDataSet(TvstrkdocasProf.get(new Coordinate(i,j,k)), 0);
                    TvstrkdocasFits.put(new Coordinate(i,j,k), new FitLine());
                    prfdvst.addDataSet(TvstrkdocasFits.get(new Coordinate(i,j,k)), 0);
                    this.getDataGroup().add(prfdvst, 0, i+1, j+1);
                    
                    ijk++;
                }
            }
        }
        //B centroids
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
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
    @Override
    public void plotHistos() {
        String[] Names = {"TrackDoca vs T","TrackDoca vs T Graphs","CalcDoca vs T","Time Residuals","Parameters"};
        for(int s = 0; s<3; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
            int NumPads = 
            this.getAnalysisCanvas().getCanvas(Names[s]).getCanvasPads().size();
            for (int n = 0; n < NumPads; n++) {
                this.getAnalysisCanvas().getCanvas(Names[s]).getPad(n).getAxisZ().setLog(true);
            }
        }
        
        for(int s = 3; s<5; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
        }
        this.getAnalysisCanvas().getCanvas(Names[3]).divide(this.nsl, 2);
        this.getAnalysisCanvas().getCanvas(Names[4]).divide(6, 6);
        
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").update();
        this.getAnalysisCanvas().getCanvas("CalcDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("Time Residuals").update();
        this.getAnalysisCanvas().getCanvas("Parameters").update();
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        
        calwriter.close();
        writer.close();
        this.UpdateBBinCenters();
        for (int i = 0; i < this.nsl; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                this.filltrkDocavsTGraphs(i,j);
            }
            //runFit(i); 
        }
        //fp.refit();
        //pw.close();
        this.getAnalysisCanvas().getCanvas("Time Residuals").divide(nsl, 3);
        //
        for(int i = 0; i<this.nsl; i++) {
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i);
            this.fitTimeResPlot(timeResiFromFile.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
        }
        fp.setGreenFitButton();
    }
    public void plotFits(boolean fitted) throws FileNotFoundException {
        if(fitted==true) {
            //pw.close();
            //File file2 = new File("");
            //file2 = outfile;
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "Files/ccdb_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
            //file2.renameTo(new File(fileName));
           
            pw = new PrintWriter(fileName);
            pw.printf("#& sector superlayer component v0 deltanm tmax distbeta delta_bfield_coefficient b1 b2 b3 b4 delta_T0 c1 c2 c3\n");
        
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
                    
                }
            }
            pw.close();
            //this.rePlotResi();
        }
    }
    private int maxIter = 10;
    double[][] fixPars = new double[6][10];
    
    private MnScan  scanner = null;
    private MnMigrad migrad = null;
    
    public int NbRunFit = 0;
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
        
        FunctionMinimum min ;
        
        
        for(int it = 0; it<maxIter; it++) {
            
            min = migrad.minimize();
            System.err.println("****************************************************");
            System.err.println("*   FIT RESULTS  FOR SUPERLAYER  "+(i+1)+" at iteration "+(it+1)+"  *");
            System.err.println("****************************************************");  
            for(int pi = 0; pi<6; pi++) 
                System.out.println("par["+pi+"]="+(TvstrkdocasFitPars.get(new Coordinate(i)).value(pi)-min.userParameters().value(pi)));
            
            if(min.isValid()) {
                TvstrkdocasFitPars.put(new Coordinate(i),min.userParameters());  
            }
            
            
            System.err.println(min);
            
        }
        
//        for(int isec = 0; isec < 6; isec++) {
//           
//            pw.printf("%d\t %d\t %d\t %.6f\t %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %d\t %.6f\t %.6f\t %d\n",
//                (isec+1), (i+1), 0,
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(0),
//                0,
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(3),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(4),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(5),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(6),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(7),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(8),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(9),
//                0,
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(2),
//                TvstrkdocasFitPars.get(new Coordinate(i)).value(1),
//                0);
//        }
        
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
    
    private final Map<Integer, Integer> hitid2tid = new HashMap<Integer, Integer>();
    private final Map<Integer, Double> hitid2B = new HashMap<Integer, Double>();
    private final Map<Integer, Double> hitid2tProp = new HashMap<Integer, Double>();
    private final Map<Integer, Double> hitid2tFlight = new HashMap<Integer, Double>();
    
    public void reCook() {
        iterationNum++;
        fp.setRedFitButton();
        
        hitid2tid.clear();
        hitid2B.clear();
        hitid2tProp.clear();
        hitid2tFlight.clear();
        
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
                
          //      DataBank bnkHBHits = event.getBank("HitBasedTrkg::HBHits");
                DataBank bnkHBHitTrkID = event.getBank("HitBasedTrkg::HBHitTrkId"); 
              //  if (eventcounter==2) bnkHBHitTrkID.show();
                
                for (int idrow = 0; idrow < bnkHBHitTrkID.rows(); idrow++ ) {
                	hitid2tid.put((int)bnkHBHitTrkID.getShort("id", idrow), (int)bnkHBHitTrkID.getShort("tid", idrow));
                	hitid2B.put((int)bnkHBHitTrkID.getShort("id", idrow), (double)bnkHBHitTrkID.getFloat("B", idrow));
                	hitid2tFlight.put((int)bnkHBHitTrkID.getShort("id", idrow), (double)bnkHBHitTrkID.getFloat("TFlight", idrow));
                	hitid2tProp.put((int)bnkHBHitTrkID.getShort("id", idrow), (double)bnkHBHitTrkID.getFloat("TProp", idrow));                
                }
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        hits.add(this.getCalHit(bnkHits, i));
                }
                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    int tempid = hit.get_Id();
                   // if (eventcounter==2) {
                   // 	System.out.println(hit.printInfo());
                   // 	System.out.println("Hit before TFlight " + hit.getTFlight());
                   // }
                    hit.setTFlight(hitid2tFlight.get(tempid));
                    hit.setTProp(hitid2tProp.get(tempid));
                    hit.setB(hitid2B.get(tempid));
                    
                    
                  //  if (eventcounter==2) {
                   // 	System.out.println(hit.printInfo());
                  //  	System.out.println("Hit after TFlight " + hit.getTFlight());
                   // }
                    
                    
                    updateHit(hit); //update the hit with previous cal results
                    
                   // if (eventcounter==2) {
                   // 	System.out.println("after updateHit function");
                   // 	System.out.println(hit.printInfo());
                   // }
                }
              //  refit with new constants
                Refit rf = new Refit();
                rf.reFit(hits);    //refit to get the parameters
                
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
                    int tempid = hit.get_Id();
                    hit.setTFlight(hitid2tFlight.get(tempid));
                    hit.setTProp(hitid2tProp.get(tempid));
                    hit.setB(hitid2B.get(tempid));
                    updateHit(hit);
                }
               // refit with new constants
                Refit rf = new Refit();
                rf.reFit(calhits);    
                for(FittedHit hit : calhits) {
                    //filling the timeResi for the newly calibrated hits
                    timeResiNew.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                    double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
                    double alphaUncor = hit.getAlpha()+(double)T2DCalib.polarity*theta0;
                    int alphaBin = this.getAlphaBin(alphaUncor);
                    double bFieldVal = (double) hit.getB();
                         
                    if(alphaBin!=-1) {
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
            }
        }
        for (int i = 0; i < this.nsl; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                this.filltrkDocavsTGraphs(i,j);
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
        //refit using a double gaussian 
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+[amp2]*gaus(x,[mean],[sigma2])", -0.3, 0.3); 
        gausFunc.setLineColor(4);
        gausFunc.setLineStyle(1);
        gausFunc.setLineWidth(2);
        gausFunc.setParameter(0, gaus1Func.getParameter(0));
        gausFunc.setParameter(1, gaus1Func.getParameter(1));
        gausFunc.setParameter(2, gaus1Func.getParameter(2)*0.75);
        gausFunc.setParameter(3, gaus1Func.getParameter(0)*0.15);
        gausFunc.setParameter(4, gaus1Func.getParameter(2));
        gausFunc.setOptStat(1110);
        h1.setOptStat(1110);
        //canvasRes.clear();
        
        DataFitter.fit(gausFunc, h1, "Q");
        //gausFunc.setOptStat(101100);
        gausFunc.setOptStat(0);
        
        int effSig = (int) (10000*this.getError(gausFunc.getParameter(0), gausFunc.getParameter(1), gausFunc.getParameter(2), 
                gausFunc.getParameter(3), gausFunc.getParameter(1), gausFunc.getParameter(4)));
        String t = "eff. Sig. "+effSig + "microns";
        h1.setTitle(t);
        canvasRes.draw(h1, "E1");
        //canvasRes.draw(gausFunc, "same");
        
        
    }
    private int getAlphaBin(double alpha) {
        int v = -1;
                //double[] AlphaValues = new double[]{-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26};
        for(int i = 0; i<AlphaValues.length; i++) {
            if(Math.abs(alpha-AlphaValues[i])<this.AlphaBinHalfWidth)
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
     private int MINENTRIES = 10;
    F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])+[p0]", 0, 1.8);
    F1D f2 = new F1D("f2","[amp1]*gaus(x,[mean1],[sigma1])+[amp2]*gaus(x,[mean2],[sigma2])+[p02]", 0, 1.8);
    
    private void filltrkDocavsTGraphs(int i, int j, int k) {
        
        if(TvstrkdocasProf.get(new Coordinate(i, j, k))!=null) {
            
            TvstrkdocasProf.get(new Coordinate(i, j, k)).reset();
            H2F h2 = Tvstrkdocas.get(new Coordinate(i, j, k));
            ArrayList<H1F> hslice = h2.getSlicesX();
            
            for(int si=0; si<hslice.size(); si++) {
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                double mean = hslice.get(si).getDataX(hslice.get(si).getMaximumBin());
                if(amp<this.MINENTRIES) {
                    
                } else {
                    
                    double x = h2.getXAxis().getBinCenter(si);
                    double y = hslice.get(si).getMean();
                    double sigma = hslice.get(si).getRMS();
                   
                    if(x/(2.*Constants.getInstance().wpdist[i])>0.9) { // for large docas (90% of dmax take fit with a double gauss and take the peak of the one at smaller distance as y
                        f2.setParameter(0, amp);
                        f2.setParameter(1, mean);
                        f2.setParameter(4, y+10);
                        f2.setParameter(2, sigma);
                        f2.setParameter(5, sigma+10);
                        f2.setParameter(6, 0);
                    
                        DataFitter.fit(f2, hslice.get(si), "Q"); //No options uses error for sigma 
                        if(f2.getChiSquare()<200 && f2.getParameter(1)>0 && f2.parameter(1).error()<50) {
                            //TvstrkdocasProf.get(new Coordinate(i, j, k)).
                            //        addPoint(x, f2.getParameter(1), 0, this.getError(f2.getParameter(0),
                            //                f2.getParameter(1),f2.getParameter(2),f2.getParameter(3),
                            //                f2.getParameter(4), f2.getParameter(5) ));
                            TvstrkdocasProf.get(new Coordinate(i, j, k)).
                                    addPoint(x, f2.getParameter(1), 0, sigma);
                        }
                    } else {
                        f1.setParameter(0, amp);
                        f1.setParameter(1, mean);
                        f1.setParameter(2, sigma);
                        f1.setParameter(3, 0);
                    
                        DataFitter.fit(f1, hslice.get(si), "Q"); //No options uses error for sigma 
                        if(f1.getChiSquare()<200 && f1.getParameter(1)>0 && f1.parameter(1).error()<50) {
                            TvstrkdocasProf.get(new Coordinate(i, j, k)).
                                    addPoint(x, f1.getParameter(1), 0, f1.getParameter(2));
                        }
                    }
                }
            }
        }
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
    boolean refitSegs = true;
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
            Constants.getInstance().initialize("DCCAL");
            String newVar = String.valueOf(T2DViewer.calVariation.getSelectedItem());
            System.out.println("* VARIATION *"+newVar);
            ccdb.setVariation(newVar);
            TableLoader.FillT0Tables(newRun, newVar);
            TableLoader.Fill(T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/time2dist"));  
            this.loadFitPars(); 
            polarity = (int)Math.signum(event.getBank("RUN::config").getFloat("torus",0));
            runNumber = newRun;
            String distBetaFCNString = String.valueOf(T2DViewer.distBetaFCN.getSelectedItem());
            if(distBetaFCNString.equalsIgnoreCase("old")) {
                Utilities.NEWDELTATBETAFCN = false;
                System.out.println("USING OLD DISTBETA FCN!!!!");
            } else {
                Utilities.NEWDELTATBETAFCN = true;
                System.out.println("USING NEW DISTBETA FCN!!!!");
            }
            numberprocessedevents = Integer.parseInt(T2DViewer.enternofevents.getText());
            if (numberprocessedevents==-1)
            {
            	System.out.println("All events are processed!!!!");
            }
            else {
            	System.out.println(numberprocessedevents + " events will be processed!!!!");
            }
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits") || !event.hasBank("HitBasedTrkg::HBHits") || !event.hasBank("HitBasedTrkg::HBHitTrkId") ) {
            return;
        } 
     
        if (count > numberprocessedevents && numberprocessedevents!=-1) {
        	return;
        }
        
        
        // get segment property     
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        DataBank hitbasedbank = event.getBank("HitBasedTrkg::HBHits");
        DataBank hitbasedhittrkid = event.getBank("HitBasedTrkg::HBHitTrkId");
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
                if(superlayer>2 && superlayer<5) { 
                    int bBin = this.getBBin(bFieldVal);
                    // fill B values histograms
                    if(superlayer ==3 || superlayer ==4)
                        B.get(new Coordinate(superlayer-3, alphaBin, bBin))
                                .fill(bFieldVal);
                }
            }
        }
        hitmap.forEach((k,v) -> hitlist.add(v));
        calhitmap.forEach((k,v) -> calhitlist.add(v));
        calhipoEvent = (HipoDataEvent) calwriter.createEvent();
        hipoEvent = (HipoDataEvent) writer.createEvent();
        //selected.show();
        calhipoEvent.appendBank(this.fillTBHitsBank(event, calhitlist));
        hipoEvent.appendBank(this.fillTBHitsBank(event, hitlist));
        calhipoEvent.appendBank(hitbasedbank);
        hipoEvent.appendBank(hitbasedbank);
        calhipoEvent.appendBank(hitbasedhittrkid);
        hipoEvent.appendBank(hitbasedhittrkid);
        calwriter.writeEvent(calhipoEvent);
        writer.writeEvent(hipoEvent);
        
        //refit with new constants
        if(this.refitSegs) {
            rf.reFit(calhitlist);    
            for(FittedHit hit : calhitlist) {// get the new track docas after refitting using the TB tFlight and tProp
                double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
                double alphaUncor = hit.getAlpha()+(double)T2DCalib.polarity*theta0;
                int alphaBin = this.getAlphaBin(alphaUncor);
                double bFieldVal = (double) hit.getB();

                if(alphaBin!=-1) {
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
        }
        calhitlist.clear();
        hitlist.clear();
        
    }
    
    private double[][] resetPars = new double[6][11];
    private String[] parNames = {"v0", "vmid", "R", "tmax", "distbeta", "delBf", 
        "b1", "b2", "b3", "b4", "dmax"};
    private double[] errs = {0.001,0.001,0.01,1.0,0.01,0.001,0.001,0.001,0.001,0.001,0.00001};
    
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
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").
                        draw(TvstrkdocasFits.get(new Coordinate(i, j, BBins)), "same");
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").
                        draw(TvstrkdocasFits.get(new Coordinate(i, j, BBins)), "same");
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
                    this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").
                            draw(TvstrkdocasFits.get(new Coordinate(i, j, k)), "same");
                    this.getAnalysisCanvas().getCanvas("CalcDoca vs T").
                    draw(TvstrkdocasFits.get(new Coordinate(i, j, k)), "same");
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
        double Alpha = bnkHits.getFloat("Alpha", i);
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
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        //double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        //double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double alpha = bnkHits.getFloat("Alpha", i);

        hit = new FittedHit(sector, superlayer, layer, wire, TDC, id);
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
        double Alpha = bnkHits.getFloat("Alpha", i);
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
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        //double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        //double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double alpha = bnkHits.getFloat("Alpha", i);
        
        hit = new FittedHit(sector, superlayer, layer, wire, TDC, id);
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
        //this.getSegProperty(bnkHits);
        
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
                }
            }
        }
        
        return pid;
    } 
    
   // private void updatewithHBInfo(FittedHit hit, double TFlight, double TProp) {
   // 	
   // }
    
    
    private void updateHit(FittedHit hit) { 
        //double distbeta = TvstrkdocasFitPars.get(new Coordinate(hit.get_Superlayer()-1)).value(4);
        double distbeta = TableLoader.distbeta[hit.get_Sector()-1][hit.get_Superlayer()-1];
        double d = hit.get_ClusFitDoca();
        double beta = hit.get_Beta();
        double calibTime = (double) (hit.get_TDC() - hit.getTProp()
                                            - hit.getTFlight() - hit.getTStart()
                                            - hit.getT0());
        double deltatime_beta = 0;
        
        if(Utilities.NEWDELTATBETAFCN==false) {
            deltatime_beta = util.calcDeltaTimeBeta(d, distbeta, beta);
        } else {
            deltatime_beta = util.calcDeltaTimeBetaNewFCN(calibTime, distbeta, beta);
        }
        
        hit.set_DeltaTimeBeta(deltatime_beta);
        hit.set_Doca(this.timeToDistance(hit.get_Sector(), hit.get_Superlayer(), 
                hit.getAlpha(), hit.get_Beta(), hit.getB(), calibTime, deltatime_beta));
        
        double x = hit.get_XWire(); 
        double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
        //fix alpha to get the local angle
        double alphaRadUncor = Math.toRadians(hit.getAlpha()+(double)T2DCalib.polarity*theta0);
        double trkAngle = Math.tan(alphaRadUncor);
        double cosTkAng = 1./Math.sqrt(trkAngle*trkAngle + 1.);
        hit.set_X(x + hit.get_LeftRightAmb() * (hit.get_Doca() / cosTkAng) );
    }
    
    
    
    private TimeToDistanceEstimator tde = new TimeToDistanceEstimator();
    private double timeToDistance(int sector, int superlayer, double alpha, double beta, double B, double time,
                                double deltatime_beta) {
        double distance = 0;
        
        //reduce the corrected angle
        double ralpha = (double) util.getReducedAngle(alpha);
         
        if(beta>1.0)
            beta = 1.0;
        
        double correctedTime = time -deltatime_beta;
        if(correctedTime<=0)
            correctedTime=0.01; // fixes edge effects ... to be improved
        distance = tde.interpolateOnGrid(B, ralpha, 
                        correctedTime, 
                        sector-1, superlayer-1) ; 
        return distance;
    }

    private void UpdateBBinCenters() {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k < this.BBins; k++) {
                    BfieldValuesUpd[i][j][k] = BfieldValues[k];
                    if(B.get(new Coordinate(i,j,k)).getBinContent(B.get(new Coordinate(i,j,k)).getMaximumBin())>10)
                        BfieldValuesUpd[i][j][k] = B.get(new Coordinate(i,j,k)).getMean();
                    //System.out.println("ijk" +i+" "+j+" "+k+" BBin UPdated "+BfieldValues[k]+" --> "+BfieldValuesUpd[i][j][k] );
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

    private double getError(double a1, double mu1, double sig_1, double a2, double mu2, double sig_2) {
       // if(a2<0.95*a1)
       //     return sig_1;
       // if(a1<0.95*a2)
       //     return sig_2;
        
        double sig1 = Math.abs(sig_1);
        double sig2 = Math.abs(sig_2);
        
        double A1 = Math.sqrt(2*Math.PI)*sig1*a1;
        double A2 = Math.sqrt(2*Math.PI)*sig2*a2;
        
        double Momt1 = A1/(A1+A2) * mu1 + A2/(A1+A2) * mu2;
        double Momt2 = A1/(A1+A2) * (mu1*mu1 + sig1*sig1) + A2/(A1+A2) * (mu2*mu2 + sig2*sig2);
        
        double var = Momt2 - Momt1*Momt1;
        
        return Math.sqrt(var);
    }

    private boolean passPID(DataEvent event, DataBank bnkHits, int rowIdxinTrkBank) {
        boolean pass = false;
        int trkID = bnkHits.getByte("trkID", rowIdxinTrkBank);
        int pid = this.readPID(event, trkID);
        //pass if the track is identified as an electron or as a hadron
        if(pid==11 || Math.abs(pid)==211 || Math.abs(pid)==2212 || Math.abs(pid)==321) {
            pass = true;
        }
        return pass;
    }
    
    private boolean passResiCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;

        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("TFlight", i)>0 
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && this.passPID(event, bnkHits, i)==true)
            {
                pass = true;
            }
        return pass;
    }

    private boolean passCalibCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;

        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("beta", i)> Double.parseDouble(T2DViewer.betaCut.getText()) 
                    && this.selectOnAlpha(superlayer, alphaRadUncor)==true
                    && bnkHits.getFloat("TFlight", i)>0 
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getNumWireWithinDW()<=Integer.parseInt(T2DViewer.npassWires.getText())
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()>Integer.parseInt(T2DViewer.nWires.getText())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && this.passPID(event, bnkHits, i)==true)
            {
                pass = true;
            }
        return pass;
    } 
    //Function is not used, could be deleted F.H. 2/2/22
  /*  private boolean passCalibCuts(DataBank bnkHits, int i) {
        boolean pass = false;
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;

        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("beta", i)> Double.parseDouble(T2DViewer.betaCut.getText()) 
                    && this.selectOnAlpha(superlayer, alphaRadUncor)==true
                    && bnkHits.getFloat("TFlight", i)>0 
                    //&& segPropMap.get(bnkHits.getInt("clusterID", i)).getNumWireWithinDW()<=Integer.parseInt(T2DViewer.npassWires.getText())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    //&& segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()>Integer.parseInt(T2DViewer.nWires.getText())
                )
            {
                pass = true;
            }
        return pass;
    } */
}

