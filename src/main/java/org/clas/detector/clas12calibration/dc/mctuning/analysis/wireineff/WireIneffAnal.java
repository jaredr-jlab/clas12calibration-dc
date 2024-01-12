/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;
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
import org.clas.detector.clas12calibration.dc.mctuning.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.AnalysisMonitor;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.WireIneffAnalViewer;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.clas.clas.math.FastMath;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.cluster.Cluster;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.segment.Segment;
//import org.jlab.rec.dc.segment.SegmentFinder;
import org.jlab.rec.dc.trajectory.SegmentTrajectory;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
/**
 *
 * @author ziegler
 */
public class WireIneffAnal extends AnalysisMonitor{
    private SchemaFactory schemaFactory = new SchemaFactory();
    PrintWriter pw = null;
    File outfile = null;
    private int runNumber;
    
    public WireIneffAnal(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames("Inefficiency vs TrackDoca");
        this.init(false, "p0:p1:p2:p3:p4");
        outfile = new File("Files/wireineffConstantsIneff.txt");
        pw = new PrintWriter(outfile);
        pw.printf("#& Superlayer parameter1 parameter2 parameter3 parameter4 scale\n");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
       
        if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
            System.out.println(" BANK FOUND........");
        } else {
            System.out.println(" BANK NOT FOUND........");
        }
       
    }
    private Map<Coordinate, H1F> Ineffvstrkdocas                = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, GraphErrors> IneffvstrkdocasProf    = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> IneffvstrkdocasFit             = new HashMap<Coordinate, FitFunction>();
    private Map<Coordinate, MnUserParameters> IneffvstrkdocasFitPars    = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> IneffvstrkdocasFits                = new HashMap<Coordinate, FitLine>();
    
    int nsl = 6;

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup tr = new DataGroup(6,1);
        for(int i =0; i<nsl; i++) {
            Ineffvstrkdocas.put(new Coordinate(i), 
                    new H1F("layer inefficiencies vs trkDoca" + (i + 1), "superlayer" + (i + 1), 40, 0.0, 1.0));
            Ineffvstrkdocas.get(new Coordinate(i)).setTitleX("Track Doca (cm)" );
            Ineffvstrkdocas.get(new Coordinate(i)).setTitleY("Inefficiency for superlayer "+(i+1) );
            tr.addDataSet(Ineffvstrkdocas.get(new Coordinate(i)), 0);
            this.getDataGroup().add(tr, 0, i+1, 0);
            IneffvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            IneffvstrkdocasProf.put(new Coordinate(i), new GraphErrors());
            IneffvstrkdocasProf.get(new Coordinate(i)).setMarkerColor(1);
            IneffvstrkdocasFits.put(new Coordinate(i), new FitLine());
        }
        
        this.getDataGroup().add(tr, 0,0,0);
        this.loadFitPars() ;
        for (int i = 0; i < nsl; i++) {
            this.getCalib().addEntry(0,i+1,0);
            this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(0), "p0", 0, i+1, 0);
            this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(1), "p1", 0, i+1, 0);
            this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(2), "p2", 0, i+1, 0);
            this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(3), "p3", 0, i+1, 0);
            this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(4), "p4", 0, i+1, 0);
        }
        
        this.getCalib().fireTableDataChanged();
    }
    private void updateTable(int i) {
        this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(0), "p0", 0, i+1, 0);
        this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(1), "p1", 0, i+1, 0);
        this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(2), "p2", 0, i+1, 0);
        this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(3), "p3", 0, i+1, 0);
        this.getCalib().setDoubleValue(IneffvstrkdocasFitPars.get(new Coordinate(i)).value(4), "p4", 0, i+1, 0);
    }    
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").update();
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        int count = 0;
        int n = 0;
        for (int i = 0; i < this.nsl; i++) {
            for(int bb =0; bb<40; bb++) {
                n+= totLayA[i][bb];
                count++;
            }
        }
        for (int i = 0; i < this.nsl; i++) {
            for(int bb =0; bb<40; bb++) {
                float ddc = effLayA[i][bb];
                float ndc = totLayA[i][bb];
                float errdc = (float) Math.sqrt(ndc*(ndc+ddc)/40/(ddc*ddc*ddc));
                if(ndc>0) {
                    Ineffvstrkdocas
                    .get(new Coordinate(i))
                    .setBinContent(bb,100*(1. - (float)ddc/ (float)ndc));
                    Ineffvstrkdocas
                    .get(new Coordinate(i))
                    .setBinError(bb,(float)100*errdc);
                }
            }
            System.out.println("Filling graph for superlayer "+(i+1));
            this.fillGraphs(i);
            this.runFit(i);
        }
        this.plotFits();
    }
    
    public  void Refit() {
        for (int i = 0; i < this.nsl; i++) {
                this.runFit(i);
        }
        this.plotFits();
    }
    
    public void plotFits() {
        {
            pw.close();
            File file2 = new File("");
            file2 = outfile;
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "Files/wireineff_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
            file2.renameTo(new File(fileName));
            int ij =0;
            int ip =0;
            NbRunFit++;
            
            for (int i = 0; i < this.nsl; i++) {
                if(IneffvstrkdocasProf.get(new Coordinate(i)).getVectorX().size()>0) {
                    this.updateTable(i);
                    IneffvstrkdocasFits.put(new Coordinate(i), new FitLine("f"+""+i, i, 
                    IneffvstrkdocasFitPars.get(new Coordinate(i))));
                    IneffvstrkdocasFits.get(new Coordinate(i)).setLineStyle(4);
                    IneffvstrkdocasFits.get(new Coordinate(i)).setLineWidth(5);
                    IneffvstrkdocasFits.get(new Coordinate(i)).setLineColor(8);
                    this.Plot(i);
                }
            }
            this.getCalib().fireTableDataChanged();  
        }
    }
    private int maxIter = 5;
    private MnScan  scanner = null;
    private MnMigrad migrad = null;
    
    public int NbRunFit = 0;
    public void runFit(int i) {
        // i = superlayer - 1;
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** ");
        IneffvstrkdocasFit.put(new Coordinate(i), 
                new FitFunction(i, (Map<Coordinate, GraphErrors>) IneffvstrkdocasProf));
        
        scanner = new MnScan((FCNBase) IneffvstrkdocasFit.get(new Coordinate(i)), 
                IneffvstrkdocasFitPars.get(new Coordinate(i)),2);
	
        for (int p = 0; p < 1; p++) {
            scanner.fix(p); 
        }
        FunctionMinimum scanmin = scanner.minimize();
        //if(scanmin.isValid())
            IneffvstrkdocasFitPars.put(new Coordinate(i),scanmin.userParameters());
        
        migrad = new MnMigrad((FCNBase) IneffvstrkdocasFit.get(new Coordinate(i)), 
                IneffvstrkdocasFitPars.get(new Coordinate(i)),1);
        migrad.setCheckAnalyticalDerivatives(true);
        
        FunctionMinimum min ;
        
        
        for(int it = 0; it<4; it++) {
           // scanner.release(it+1);
            min = migrad.minimize();
            System.err.println("****************************************************");
            System.err.println("*   FIT RESULTS  FOR SUPERLAYER  "+(i+1)+" at iteration "+(it+1)+"  *");
            System.err.println("****************************************************");  
            
            //if(min.isValid()) {
                IneffvstrkdocasFitPars.put(new Coordinate(i),min.userParameters());  
            //}
            
            
            System.err.println(min);
            
        }
        
        for(int isec = 0; isec < 6; isec++) {
           
            pw.printf("%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n",
                (i+1), 
                IneffvstrkdocasFitPars.get(new Coordinate(i)).value(1),
                IneffvstrkdocasFitPars.get(new Coordinate(i)).value(2),
                IneffvstrkdocasFitPars.get(new Coordinate(i)).value(3),
                IneffvstrkdocasFitPars.get(new Coordinate(i)).value(4),
                IneffvstrkdocasFitPars.get(new Coordinate(i)).value(0));
        }
        
        for (int p = 0; p < 1; p++) {
            IneffvstrkdocasFitPars.get(new Coordinate(i)).release(p);
        }
        
    }
    
    int counter = 0;
    private int iterationNum = 1;
    private void fillGraphs(int i) {
        if(Ineffvstrkdocas.get(new Coordinate(i))!=null) {
            IneffvstrkdocasProf.get(new Coordinate(i)).reset(); 
            for(int si=0; si<Ineffvstrkdocas.get(new Coordinate(i)).getDataSize(0); si++) {
                double x = Ineffvstrkdocas.get(new Coordinate(i)).getDataX(si);
                double y = Ineffvstrkdocas.get(new Coordinate(i)).getDataY(si);
                double ey = Ineffvstrkdocas.get(new Coordinate(i)).getDataEY(si);
                
                if(ey>0 ) {
                    IneffvstrkdocasProf.get(new Coordinate(i)).
                                addPoint(x, y, 0, ey);
                }
            }
        }
    }
                
    int count = 0; 
    
    int[][] totLay = new int[6][6];
    int[][] effLay = new int[6][6];
    int[][] totLayA = new int[6][40];
    int[][] effLayA = new int[6][40];
    float trkDBinning = (float) ((float) 1.0/40.0);
    static ConstantsManager ccdb = WireIneffAnal.ccdb;	
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
        int jitter = (int) bnkHits.getByte("jitter", i);
        
        hit = new FittedHit(sector, superlayer, layer, wire, TDC, jitter, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
        hit.setB(B);
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        //double T0Sub = (TDC - TProp - TFlight - T0);
        //hit.set_Time(T0Sub-TStart);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(WireIneffAnalViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(WireIneffAnalViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
        hit.setAlpha(Alpha);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        
        if (bnkHits.getShort("clusterID", i) >0 ) {            
            return hit;
        } else {
            return null;
        }
    }
    
    List<FittedHit> hits = new ArrayList<FittedHit>();
    List<FittedCluster> clusters = new ArrayList<FittedCluster>();
    List<Segment> segments = new ArrayList<Segment>();
    //SegmentFinder segFinder = new SegmentFinder(); 
    
    public DataBank getLayIneffBank(DataEvent event) {
        hits.clear();
        clusters.clear();
        segments.clear();
        
        if (!event.hasBank("RUN::config")) {
            return null;
        }
        
        DataBank bank = event.getBank("RUN::config");
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) {
           return null;
        } else {
           count++;
        }
       
        if(count==1) {
            //Constants.getInstance().initialize("DCCAL");
            Driver.init();
            TableLoader.FillT0Tables(newRun, "default");
            TableLoader.Fill(ccdb.getConstants(newRun, Constants.TIME2DIST));  
            runNumber = newRun;
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return null;
        } 
        // get segment property
        
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        FittedHit hit = null;
        for (int i = 0; i < bnkHits.rows(); i++) {
            hit = this.getHit(bnkHits, i);
            if(hit!=null)
                hits.add(hit);
        }
        clusters = this.recomposeClusters(hits);
        segments =  this.get_Segments(clusters, event, WireIneffAnalViewer.dcDetector);
        DataBank bankE = null;
        if(segments!=null && segments.size()>0) {
            bankE = event.createBank("TimeBasedTrkg::TBSegmentTrajectory", segments.size() * 6);  
            int index = 0;
            for (Segment aSeglist : segments) {
                if (aSeglist.get_Id() == -1) {
                    continue;
                }
                SegmentTrajectory trj = aSeglist.get_Trajectory();
                for (int l = 0; l < 6; l++) {
                    bankE.setShort("segmentID", index, (short) trj.get_SegmentId());
                    bankE.setByte("sector", index, (byte) trj.get_Sector());
                    bankE.setByte("superlayer", index, (byte) trj.get_Superlayer());
                    bankE.setByte("layer", index, (byte) (l + 1));
                    bankE.setShort("matchedHitID", index, (short) trj.getMatchedHitId()[l]);
                    bankE.setFloat("trkDoca", index, (float) trj.getTrkDoca()[l]);
                    index++;
                }
            }
        }
        //bankE.show();
        return bankE;
    }
    
    public List<FittedCluster> recomposeClusters(List<FittedHit> fhits) {
        ClusterFitter cf = new ClusterFitter();
        Map<Integer, ArrayList<FittedHit>> grpHits = new HashMap<Integer, ArrayList<FittedHit>>();
        List<FittedCluster> clusters = new ArrayList<FittedCluster>();
        for (FittedHit hit : fhits) { 
            if (hit.get_AssociatedClusterID() == -1 || hit.get_AssociatedHBTrackID() == -1) {
                continue;
            }
            if (hit.get_AssociatedClusterID() != -1 &&
                    hit.get_AssociatedHBTrackID() != -1) {
                int index = hit.get_AssociatedHBTrackID()*10000+hit.get_AssociatedClusterID();
               
                if(grpHits.get(index)==null) { // if the list not yet created make it
                    grpHits.put(index, new ArrayList<FittedHit>()); 
                    grpHits.get(index).add(hit); // append hit
                    //System.out.println("appended first hit "+hit.get_Sector()+","+hit.get_Superlayer()+", "+hit.get_Layer()+","+hit.get_Wire()+
                    //        " to cluster "+index);
                } else {
                    grpHits.get(index).add(hit); // append hit
                   // System.out.println("appended subs hit "+hit.get_Sector()+","+hit.get_Superlayer()+", "+hit.get_Layer()+","+hit.get_Wire()+
                    //        " to cluster "+index);
                }
            }
        }
        Iterator<Map.Entry<Integer, ArrayList<FittedHit>>> itr = grpHits.entrySet().iterator(); 
          
        while(itr.hasNext()) {
            Map.Entry<Integer, ArrayList<FittedHit>> entry = itr.next(); 
             
            if(entry.getValue().size()>3) {
                Cluster cluster = new Cluster(entry.getValue().get(0).get_Sector(), 
                        entry.getValue().get(0).get_Superlayer(), entry.getValue().get(0).get_AssociatedClusterID());
                FittedCluster fcluster = new FittedCluster(cluster);
                fcluster.addAll(entry.getValue());
                for(FittedHit h : fcluster) {
                    //System.out.println(h.printInfo());
                    //System.out.println("          .........."+h.get_X()+" d "+h.get_Doca()+ " +/- "+h.get_DocaErr());
                    h.set_TrkgStatus(1);
                }
                cf.SetFitArray(fcluster, "TSC");
                cf.Fit(fcluster, true);
                cf.SetResidualDerivedParams(fcluster, true, false, WireIneffAnalViewer.dcDetector); //calcTimeResidual=false, resetLRAmbig=false 
                cf.Fit(fcluster, false);
                cf.SetSegmentLineParameters(fcluster.get(0).get_Z(), fcluster);
              
                clusters.add(fcluster);
            }
        }
        return clusters;
    }
    
    @Override
    public void processEvent(DataEvent event) {
        
        //WireIneffAnalViewer.tm.processDataEvent(event);
        //if(event.hasBank("TimeBasedTrkg::TBSegmentTrajectory")==false)
        //    return;
        //DataBank Bank = event.getBank("TimeBasedTrkg::TBSegmentTrajectory") ;
        
        DataBank Bank = this.getLayIneffBank(event);
        
        int nrows =  Bank.rows();
        //Bank.show(); System.out.println(" NUMBER OF ENTRIES IN BANK = "+nrows);
        for (int i = 0; i < nrows; i++) {
            totLay[Bank.getByte("superlayer", i)-1][Bank.getByte("layer", i)-1]++;

            //System.out.println(" Layer eff denom for ["+Bank.getByte("sector", i)+"]["+ Bank.getByte("superlayer", i)+"]["+Bank.getByte("layer", i)+"] = "+totLay[Bank.getByte("sector", i)-1][Bank.getByte("superlayer", i)-1][Bank.getByte("layer", i)-1]);
            if(Bank.getShort("matchedHitID", i)!=-1) {
                effLay[Bank.getByte("superlayer", i)-1][Bank.getByte("layer", i)-1]++;
            }
            //System.out.println(" Layer eff num for ["+Bank.getByte("sector", i)+"]["+ Bank.getByte("superlayer", i)+"]["+Bank.getByte("layer", i)+"] = "+effLay[Bank.getByte("sector", i)-1][Bank.getByte("superlayer", i)-1][Bank.getByte("layer", i)-1]);

            int bb = this.getTrkDocBin(Math.abs(Bank.getFloat("trkDoca", i))/(2.0*Constants.getInstance().wpdist[Bank.getByte("superlayer", i)-1]));
            
            if(bb<40) {
                //(int)((Math.floor(Math.abs(Bank.getFloat("trkDoca", i))/(2.*Constants.wpdist[Bank.getByte("superlayer", i)-1])/trkDBinning)))
                totLayA[Bank.getByte("superlayer", i)-1][bb]++;
                if(Bank.getShort("matchedHitID", i)!=-1) {
                        effLayA[Bank.getByte("superlayer", i)-1][bb]++;
                }
            }
        }
    }
    private String[] parNames = {"p0", "p1", "p2", "p3", "p4"};
    public double[][] resetPars = new double[6][parNames.length];
    
    private double[] errs = {0.01,0.001,0.001,0.001,0.001};
    
    
    public void loadFitPars() {
        DatabaseConstantProvider dbprovider = new DatabaseConstantProvider(11, "default");
        dbprovider.loadTable("/calibration/dc/signal_generation/intrinsic_inefficiency");
        //disconnect from database. Important to do this after loading tables.
        dbprovider.disconnect();
//        for (int i = 0; i < dbprovider.length("/calibration/dc/signal_generation/intrinsic_inefficiency/Superlayer"); i++) {
//            int iSly = dbprovider.getInteger("/calibration/dc/signal_generation/intrinsic_inefficiency/Superlayer", i);
//            double parameter1 = dbprovider.getDouble("/calibration/dc/signal_generation/intrinsic_inefficiency/parameter1", i);
//            double parameter2 = dbprovider.getDouble("/calibration/dc/signal_generation/intrinsic_inefficiency/parameter2", i);
//            double parameter3 = dbprovider.getDouble("/calibration/dc/signal_generation/intrinsic_inefficiency/parameter3", i);
//            double parameter4 = dbprovider.getDouble("/calibration/dc/signal_generation/intrinsic_inefficiency/parameter4", i);
//            double scale = dbprovider.getDouble("/calibration/dc/signal_generation/intrinsic_inefficiency/scale", i);
//            scale = 100;
//            IneffvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[0], scale, errs[0]);
//            IneffvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[1], parameter1, errs[1]);
//            IneffvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[2], parameter2, errs[2]);
//            IneffvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[3], parameter3, errs[3]);
//            IneffvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[4], parameter4, errs[4]);
//           
//            System.out.println("LOADED PARS "+scale+" : "+parameter1+" : "+parameter2+" : "+parameter3+" : "+parameter4);
//        }
        //renit
            IneffvstrkdocasFitPars.get(new Coordinate(0)).add(parNames[0], 100, errs[0]);
            IneffvstrkdocasFitPars.get(new Coordinate(0)).add(parNames[1], 0.00001, errs[1]);
            IneffvstrkdocasFitPars.get(new Coordinate(0)).add(parNames[2], 0.0282, errs[2]);
            IneffvstrkdocasFitPars.get(new Coordinate(0)).add(parNames[3], 0.0413, errs[3]);
            IneffvstrkdocasFitPars.get(new Coordinate(0)).add(parNames[4], 0.0578, errs[4]);
            IneffvstrkdocasFitPars.get(new Coordinate(1)).add(parNames[0], 100, errs[0]);
            IneffvstrkdocasFitPars.get(new Coordinate(1)).add(parNames[1], 0.00001, errs[1]);
            IneffvstrkdocasFitPars.get(new Coordinate(1)).add(parNames[2], 0.0282, errs[2]);
            IneffvstrkdocasFitPars.get(new Coordinate(1)).add(parNames[3], 0.0413, errs[3]);
            IneffvstrkdocasFitPars.get(new Coordinate(1)).add(parNames[4], 0.0578, errs[4]);
            IneffvstrkdocasFitPars.get(new Coordinate(2)).add(parNames[0], 100, errs[0]);
            IneffvstrkdocasFitPars.get(new Coordinate(2)).add(parNames[1], 0.000015, errs[1]);
            IneffvstrkdocasFitPars.get(new Coordinate(2)).add(parNames[2], 0.0421, errs[2]);
            IneffvstrkdocasFitPars.get(new Coordinate(2)).add(parNames[3], 0.011, errs[3]);
            IneffvstrkdocasFitPars.get(new Coordinate(2)).add(parNames[4], 0.01413, errs[4]);
            IneffvstrkdocasFitPars.get(new Coordinate(3)).add(parNames[0], 100, errs[0]);
            IneffvstrkdocasFitPars.get(new Coordinate(3)).add(parNames[1], 0.000015, errs[1]);
            IneffvstrkdocasFitPars.get(new Coordinate(3)).add(parNames[2], 0.0421, errs[2]);
            IneffvstrkdocasFitPars.get(new Coordinate(3)).add(parNames[3], 0.011, errs[3]);
            IneffvstrkdocasFitPars.get(new Coordinate(3)).add(parNames[4], 0.01413, errs[4]);
            IneffvstrkdocasFitPars.get(new Coordinate(4)).add(parNames[0], 100, errs[0]);
            IneffvstrkdocasFitPars.get(new Coordinate(4)).add(parNames[1], 0.000017, errs[1]);
            IneffvstrkdocasFitPars.get(new Coordinate(4)).add(parNames[2], 0.03, errs[2]);
            IneffvstrkdocasFitPars.get(new Coordinate(4)).add(parNames[3], 0.004, errs[3]);
            IneffvstrkdocasFitPars.get(new Coordinate(4)).add(parNames[4], 0.0, errs[4]);
            IneffvstrkdocasFitPars.get(new Coordinate(5)).add(parNames[0], 100, errs[0]);
            IneffvstrkdocasFitPars.get(new Coordinate(5)).add(parNames[1], 0.000017, errs[1]);
            IneffvstrkdocasFitPars.get(new Coordinate(5)).add(parNames[2], 0.03, errs[2]);
            IneffvstrkdocasFitPars.get(new Coordinate(5)).add(parNames[3], 0.004, errs[3]);
            IneffvstrkdocasFitPars.get(new Coordinate(5)).add(parNames[4], 0.0, errs[4]);
        
    }
    
    
    public void Plot(int i ) {
    
        if(IneffvstrkdocasProf.get(new Coordinate(i)).getVectorX().size()>0) {
            this.getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(0);
            this.getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(IneffvstrkdocasProf.get(new Coordinate(i)));
                    
            this.getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").
                    draw(IneffvstrkdocasFits.get(new Coordinate(i)), "same");
            this.getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").getPad(0).getAxisY().setRange(0, 20);
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
           this.Plot(layer-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       
   
    }

    private int getTrkDocBin(double d) {
        int bin = 39;
        double binWidth = 1./40.;
        double lo = 0.;
        double hi = 40.;
        int n = 40;
        
        for (int i = 0; i < n; i++) {
            double blo = i*binWidth;
            double bhi = (i+1)*binWidth;
            if(d>blo && d<=bhi)
                bin = i;
        }
        return bin;
    }
    
     /**
     * 
     * @param seg Segment
     * @param event HipoDataEvent
     * @param DcDetector DC detector utility
     */
    public void get_LayerEfficiencies(Segment seg, DataEvent event, DCGeant4Factory DcDetector) {
        if (seg!=null) { 
            // get all the hits to obtain layer efficiency
            if (event.hasBank("DC::tdc") != false || event.hasBank("TimeBasedTrkg::TBHits") != false) {
              
                DataBank bankDGTZ ;
                String wireStg;
                Boolean layer1_36;
                if(event.hasBank("DC::tdc") != false) {
                    bankDGTZ = event.getBank("DC::tdc");
                    wireStg = "component";
                    layer1_36 = true;
                } else {
                    bankDGTZ = event.getBank("TimeBasedTrkg::TBHits");
                    wireStg = "wire";
                    layer1_36 = false;
                }

                int rows = bankDGTZ.rows();
                int[] sector = new int[rows];
                int[] layer = new int[rows];
                int[] wire = new int[rows];
                int[] tdc = new int[rows];

                for (int i = 0; i < rows; i++) {
                    sector[i] = bankDGTZ.getByte("sector", i);
                    layer[i] = bankDGTZ.getByte("layer", i);
                    wire[i] = bankDGTZ.getShort(wireStg, i);
                    tdc[i] = bankDGTZ.getInt("TDC", i);
                }

                int size = layer.length;
                int[] layerNum = new int[size];
                int[] superlayerNum = new int[size];
                //double[] smearedTime = new double[size];

                for (int i = 0; i < size; i++) { 
                    if(layer1_36==true) {
                        superlayerNum[i] = (layer[i] - 1) / 6 + 1;
                        layerNum[i] = layer[i] - (superlayerNum[i] - 1) * 6;
                    } else {
                        superlayerNum[i] = bankDGTZ.getByte("superlayer", i);
                        layerNum[i] = layer[i];
                    }
                }

                // Get the Segment Trajectory
                SegmentTrajectory trj = new SegmentTrajectory();
                trj.set_SegmentId(seg.get_Id());
                trj.set_Superlayer(seg.get_Superlayer());
                trj.set_Sector(seg.get_Sector());
                double[] trkDocas = new double[6];
                int[] matchHits = new int[6];

                int[][] matchedHits = new int[1][6]; // first arrays = how many wires off
                for (int i1 = 0; i1 < 1; i1++) {
                    for (int i2 = 0; i2 < 6; i2++) {
                        matchedHits[i1][i2] = -1;
                    }
                }

                for (int l = 0; l < 6; l++) {
                    double z = DcDetector.getWireMidpoint(seg.get_Sector() - 1, seg.get_Superlayer() - 1, l, 0).z;
                    double trkXMP = seg.get_fittedCluster().get_clusterLineFitSlopeMP() * z + seg.get_fittedCluster().get_clusterLineFitInterceptMP();
                    double trkX = seg.get_fittedCluster().get_clusterLineFitSlope() * z + seg.get_fittedCluster().get_clusterLineFitIntercept();

                    if (trkX == 0) {
                        continue; // should always get a cluster fit
                    }
                    int trjWire = this.getWireOnTrajectory(seg, l + 1, trkXMP, DcDetector);
                    if(trjWire == -1) {
                        seg.set_Status(-1);
                    }
                    double x = DcDetector.getWireMidpoint(seg.get_Sector() - 1, seg.get_Superlayer() - 1, l, trjWire - 1).x;
                    double cosTkAng = FastMath.cos(Math.toRadians(6.)) * Math.sqrt(1. + seg.get_fittedCluster().get_clusterLineFitSlope() * seg.get_fittedCluster().get_clusterLineFitSlope());
                    double calc_doca = (x - trkX) * cosTkAng;
                    trkDocas[l] = calc_doca;
                    
                    double xWire = 
                            DcDetector.getWireMidpoint(seg.get_Sector() - 1, seg.get_Superlayer() - 1, l, trjWire - 1).x;
                    double zWire = 
                            DcDetector.getWireMidpoint(seg.get_Sector() - 1, seg.get_Superlayer() - 1, l, trjWire - 1).z;
                   
                    Line3D FitLine = new Line3D();
                    double z0 = DcDetector.getWireMidpoint(seg.get_Sector() - 1, seg.get_Superlayer() - 1, 0, 0).z;
                    double x0 = seg.get_fittedCluster().get_clusterLineFitSlope() * z0 + seg.get_fittedCluster().get_clusterLineFitIntercept();

                    Point3D pointOnTrk = new Point3D(z0, x0, 0);
                    Vector3D trkDir = new Vector3D(1, seg.get_fittedCluster().get_clusterLineFitSlope(), 0);
                    trkDir.unit();
                    FitLine.set(pointOnTrk, trkDir);
                    Point3D Wire = new Point3D(zWire, xWire, 0);

                    //double trkDocaMP = -xWire + (FitPars.slope()*FitArray[0][i]+FitPars.intercept());
                    double trkDocaMP = FitLine.distance(Wire).length();
                    double trkDoca = trkDocaMP * stereo;
                    int matchWireTightness = 1;
                    //if(trkDoca/cellSize[l]>0.8)
                    //    matchWireTightness =2;
                    
                    //
                    trkDocas[l] = trkDoca;
                    for(FittedHit h : seg) {
                        if(h.get_Layer()==(l+1)) {
                            matchedHits[0][l] = h.get_Id();
                            trkDoca = h.get_Doca();
                        }
                    }     
                    for (int j = 0; j < sector.length; j++) {
                        if (sector[j] == seg.get_Sector() && superlayerNum[j] == seg.get_Superlayer()) {
                            if (layerNum[j] == l + 1) {
                                
                                for (int wo = 0; wo < matchWireTightness; wo++) {
                                    if (Math.abs(trjWire - wire[j]) == wo) {
                                        matchedHits[wo][l] = (j + 1);
                                    } else {
                                    }
                                }
                            }
                        }
                    }
                    matchHits[l] = -1;
                    for (int wo = 0; wo < matchWireTightness; wo++) {
                        if (matchedHits[wo][l] != -1) {
                            matchHits[l] = matchedHits[wo][l];
                            wo = matchWireTightness; 
                        }
                    }
                }
                trj.setTrkDoca(trkDocas);
                trj.setMatchedHitId(matchHits);
                if(seg.get_Status()!=-1)
                    seg.set_Trajectory(trj);
            }
        }
    }
    
    private double[] cellSize = new double[]{0.38616*2, 0.40422*2, 0.621906*2,0.658597*2,
        0.9351482*2, 0.977982*2};
    private double stereo = FastMath.cos(Math.toRadians(6.));
    public List<Segment> get_Segments(List<FittedCluster> allClusters, DataEvent event, DCGeant4Factory DcDetector) {
        List<Segment> segList = new ArrayList<Segment>();
        int nCountPerLayer=0;
        for (FittedCluster fClus : allClusters) {
            
            if (fClus.size() > 12) {
                continue;
            }
            if (fClus.get_TrkgStatus() == -1) {
                return segList;
            }

            Segment seg = new Segment(fClus);
            seg.set_fitPlane(DcDetector);
            
            this.get_LayerEfficiencies(seg, event, DcDetector);
            
            double sumRes=0;
            double sumTime=0;
            
            for(FittedHit h : seg) {
                sumRes+=h.get_TimeResidual();
                sumTime+=h.get_Time();
            }
            seg.set_ResiSum(sumRes);
            seg.set_TimeSum(sumTime);
            
            segList.add(seg);
        }

        //this.setAssociatedID(segList);
        return segList;

    }

    public int getWireOnTrajectory(Segment seg, int layer, double trkX, DCGeant4Factory DcDetector) {
        int sector = seg.get_Sector();
        int superlayer = seg.get_Superlayer();
        //double x1 = GeometryLoader.dcDetector.getSector(0).getSuperlayer(superlayer-1).getLayer(layer-1).getComponent(1).getMidpoint().x();
        //double x0 = GeometryLoader.dcDetector.getSector(0).getSuperlayer(superlayer-1).getLayer(layer-1).getComponent(0).getMidpoint().x();
        double x1 = DcDetector.getWireMidpoint(sector-1, superlayer-1, layer-1, 1).x;
        double x0 = DcDetector.getWireMidpoint(sector-1, superlayer-1, layer-1, 0).x;

        double deltax = Math.abs(x1-x0);

        //double xFirstCell = GeometryLoader.dcDetector.getSector(0).getSuperlayer(superlayer-1).getLayer(layer-1).getComponent(0).getMidpoint().x();
        double xFirstCell = DcDetector.getWireMidpoint(sector-1, superlayer-1, layer-1, 0).x;
        double res = 0.035;
        int nearestWire = (int)Math.ceil((trkX-xFirstCell+deltax/2.)/deltax) ;
        int nearestWireM = (int)Math.ceil((trkX-res-xFirstCell+deltax/2.)/deltax) ;
        int nearestWireP = (int)Math.ceil((trkX+res-xFirstCell+deltax/2.)/deltax) ;
        
        if(nearestWire>112)
                nearestWire = 112;
        if(nearestWire<1 && nearestWire!=-1)
                nearestWire = 1;
        return nearestWire;
    }

}

