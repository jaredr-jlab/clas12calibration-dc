package org.clas.detector.clas12calibration.dc.t2d;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map; 
import org.clas.detector.clas12calibration.dc.calt2d.FitFunction;
import org.clas.detector.clas12calibration.dc.calt2d.T2DCalib;
import org.clas.detector.clas12calibration.dc.calt2d.Utilities;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.groot.data.GraphErrors;
import org.jlab.rec.dc.Constants;
import org.jlab.utils.groups.IndexedTable;


public class TableLoader {

    public TableLoader() {
            // TODO Auto-generated constructor stub
    }
    static final protected int nBinsT=2000;
    //public static double[][][][][] DISTFROMTIME = new double[6][6][6][6][850]; // sector slyr alpha Bfield time bins
    
    static boolean T2DLOADED = false;
    static boolean T0LOADED = false;
    
    public static double[] BfieldValues = new double[]{0.0000, 1.0000, 1.4142, 1.7321, 2.0000, 2.2361, 2.4495, 2.6458};
    public static double[] betaValues = new double[] {0.6, 0.7, 0.8, 0.9, 1.0};
    static int minBinIdxB = 0;
    static int maxBinIdxB = BfieldValues.length-1;
    static int minBinIdxAlpha = 0;
    static int maxBinIdxAlpha = 5;
    public static double[] AlphaMid = new double[6];
    public static double[][] AlphaBounds = new double[6][2];
    static int minBinIdxT  = 0;
    static int[][][][] maxBinIdxT  = new int[6][6][8][6];
    public static double[][][][][][] DISTFROMTIME = new double[6][6][maxBinIdxB+1][maxBinIdxAlpha+1][betaValues.length][nBinsT]; // sector slyr alpha Bfield time bins [s][r][ibfield][icosalpha][tbin]
    
    
    private static double[][][][] T0 ;
    private static double[][][][] T0ERR ;
    public static synchronized void FillT0Tables(int run, String variation) {
        if (T0LOADED) return;
        System.out.println(" T0 TABLE FILLED..... for Run "+run+" with VARIATION "+variation);
        DatabaseConstantProvider dbprovider = new DatabaseConstantProvider(run, variation);
        dbprovider.loadTable("/calibration/dc/time_corrections/T0Corrections");
        //disconnect from database. Important to do this after loading tables.
        dbprovider.disconnect();
        // T0-subtraction
        double[][][][] T0 ;
        double[][][][] T0ERR ;
        //T0s
        T0 = new double[6][6][7][6]; //nSec*nSL*nSlots*nCables
        T0ERR = new double[6][6][7][6]; //nSec*nSL*nSlots*nCables
        for (int i = 0; i < dbprovider.length("/calibration/dc/time_corrections/T0Corrections/Sector"); i++) {
            int iSec = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Sector", i);
            int iSly = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Superlayer", i);
            int iSlot = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Slot", i);
            int iCab = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Cable", i);
            double t0 = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Correction", i);
            double t0Error = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Error", i);
            T0[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0; 
            T0ERR[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0Error;
            TableLoader.setT0(T0);
            TableLoader.setT0Err(T0ERR);
            //System.out.println("T0 = "+t0);
        }
        T0LOADED = true;
    }

    /**
     * @return the T0
     */
    public static double[][][][] getT0() {
        return T0;
    }

    /**
     * @param aT0 the T0 to set
     */
    public static void setT0(double[][][][] aT0) {
        T0 = aT0;
    }

    /**
     * @return the T0ERR
     */
    public static double[][][][] getT0Err() {
        return T0ERR;
    }

    /**
     * @param aT0ERR the T0ERR to set
     */
    public static void setT0Err(double[][][][] aT0ERR) {
        T0ERR = aT0ERR;
    }
    
    public static int getAlphaBin(double Alpha) {
        int bin = 0;
        for(int b =0; b<6; b++) {
            if(Alpha>=AlphaBounds[b][0] && Alpha<=AlphaBounds[b][1] )
                bin = b;
        }
        return bin;
    }
    public static int maxTBin = -1;
    public static synchronized void FillAlpha() {
        for(int icosalpha =0; icosalpha<maxBinIdxAlpha+1; icosalpha++) {

            double cos30minusalphaM = Math.cos(Math.toRadians(30.)) + (double) 
                    (icosalpha)*(1. - Math.cos(Math.toRadians(30.)))/5.;
            double alphaM = -(Math.toDegrees(Math.acos(cos30minusalphaM)) - 30);
            AlphaMid[icosalpha]= alphaM;
            double cos30minusalphaU = Math.cos(Math.toRadians(30.)) + (double) 
                    (icosalpha+0.5)*(1. - Math.cos(Math.toRadians(30.)))/5.;
            double alphaU = -(Math.toDegrees(Math.acos(cos30minusalphaU)) - 30);
            AlphaBounds[icosalpha][1] = alphaU;
            double cos30minusalphaL = Math.cos(Math.toRadians(30.)) + (double) 
                    (icosalpha-0.5)*(1. - Math.cos(Math.toRadians(30.)))/5.;
            double alphaL = -(Math.toDegrees(Math.acos(cos30minusalphaL)) - 30);
            AlphaBounds[icosalpha][0] = alphaL;
        }
        AlphaMid[0] = 0;
        AlphaMid[5] = 30;
        AlphaBounds[0][0] = 0;
        AlphaBounds[5][1] = 30;
    }
    static Map<Integer, GraphErrors> sup1 = new HashMap<Integer, GraphErrors>();
    static Map<Integer, GraphErrors> sup2 = new HashMap<Integer, GraphErrors>();
    static Map<Integer, GraphErrors> sup3 = new HashMap<Integer, GraphErrors>();
    static Map<Integer, GraphErrors> sup4 = new HashMap<Integer, GraphErrors>();
    static Map<Integer, GraphErrors> sup5 = new HashMap<Integer, GraphErrors>();
    static Map<Integer, GraphErrors> sup6 = new HashMap<Integer, GraphErrors>();
    public static synchronized void Fill(IndexedTable t2dPressure, IndexedTable t2dPressRef, IndexedTable pressure) {
        for(int a = 0; a<6; a++ ){ // loop over alpha
            sup1.put(a, new GraphErrors()); 
            sup1.get(a).setMarkerColor(a+1);
            sup2.put(a, new GraphErrors()); 
            sup2.get(a).setMarkerColor(a+1);
            sup3.put(a, new GraphErrors()); 
            sup3.get(a).setMarkerColor(a+1);
            sup4.put(a, new GraphErrors()); 
            sup4.get(a).setMarkerColor(a+1);
            sup5.put(a, new GraphErrors()); 
            sup5.get(a).setMarkerColor(a+1);
            sup6.put(a, new GraphErrors()); 
            sup6.get(a).setMarkerColor(a+1);
        }
        sup1.get(0).setTitleX("time (ns)");
        sup1.get(0).setTitleY("calc doca (cm)");
        sup2.get(0).setTitleX("time (ns)");
        sup2.get(0).setTitleY("calc doca (cm)");
        sup3.get(0).setTitleX("time (ns)");
        sup3.get(0).setTitleY("calc doca (cm)");
        sup4.get(0).setTitleX("time (ns)");
        sup4.get(0).setTitleY("calc doca (cm)");
        sup5.get(0).setTitleX("time (ns)");
        sup5.get(0).setTitleY("calc doca (cm)");
        sup6.get(0).setTitleX("time (ns)");
        sup6.get(0).setTitleY("calc doca (cm)");
        //CCDBTables 0 =  "/calibration/dc/signal_generation/doca_resolution";
        //CCDBTables 1 =  "/calibration/dc/time_to_distance/t2d";
        //CCDBTables 2 =  "/calibration/dc/time_corrections/T0_correction";	
        if (T2DLOADED) return;
        
        double stepSize = 0.0010;
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);
        
        FillAlpha();
        double p_ref = t2dPressRef.getDoubleValue("pressure", 0,0,0);
        double p = pressure.getDoubleValue("value", 0,0,3);
        double dp = p - p_ref;
        double dp2scale = 0;
        double dpscale = 1;
        boolean useP = Boolean.parseBoolean(T2DViewer.usePressureTerm.getText());
        if(!useP) 
            dpscale = 0;
        for(int s = 0; s<6; s++ ){ // loop over sectors
            for(int r = 0; r<6; r++ ){ //loop over slys
                // Fill constants
                FracDmaxAtMinVel[s][r] = t2dPressure.getDoubleValue("c1_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("c1_a1", s+1,r+1,0)*dp*dpscale;
                v0[s][r] = t2dPressure.getDoubleValue("v0_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("v0_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("v0_a2", s+1,r+1,0)*dp*dp*dp2scale;
                vmid[s][r] = t2dPressure.getDoubleValue("vmid_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("vmid_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("vmid_a2", s+1,r+1,0)*dp*dp*dp2scale;
                distbeta[s][r] = t2dPressure.getDoubleValue("distbeta_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("distbeta_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("distbeta_a2", s+1,r+1,0)*dp*dp*dp2scale;
                delta_bfield_coefficient[s][r] = t2dPressure.getDoubleValue("delta_bfield_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("delta_bfield_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("delta_bfield_a2", s+1,r+1,0)*dp*dp*dp2scale;
                b1[s][r] = t2dPressure.getDoubleValue("b1_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("b1_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("b1_a2", s+1,r+1,0)*dp*dp*dp2scale;
                b2[s][r] = t2dPressure.getDoubleValue("b2_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("b2_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("b2_a2", s+1,r+1,0)*dp*dp*dp2scale;
                b3[s][r] = t2dPressure.getDoubleValue("b3_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("b3_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("b3_a2", s+1,r+1,0)*dp*dp*dp2scale;
                b4[s][r] = t2dPressure.getDoubleValue("b4_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("b4_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("b4_a2", s+1,r+1,0)*dp*dp*dp2scale;
                Tmax[s][r] = t2dPressure.getDoubleValue("tmax_a0", s+1,r+1,0)
                        +t2dPressure.getDoubleValue("tmax_a1", s+1,r+1,0)*dp*dpscale
                        +t2dPressure.getDoubleValue("tmax_a2", s+1,r+1,0)*dp*dp*dp2scale;
             
                // end fill constants
                //System.out.println("sector "+(s+1)+" sly "+(r+1)+" v0 "+v0[s][r]+" vmid "+vmid[s][r]+" R "+FracDmaxAtMinVel[s][r]);
                double dmax = 2.*Constants.getInstance().wpdist[r]; 
                //double tmax = CCDBConstants.getTMAXSUPERLAYER()[s][r];
                for(int ibfield =0; ibfield<maxBinIdxB+1; ibfield++) {
                    double bfield = BfieldValues[ibfield];
                    for(int ibeta=0; ibeta<betaValues.length; ibeta++) {

                        for(int icosalpha =0; icosalpha<maxBinIdxAlpha+1; icosalpha++) {
                            maxBinIdxT[s][r][ibfield][icosalpha] = nBinsT; 
                            double cos30minusalpha = Math.cos(Math.toRadians(30.)) + (double) (icosalpha)*(1. - Math.cos(Math.toRadians(30.)))/5.;
                            double alpha = -(Math.toDegrees(Math.acos(cos30minusalpha)) - 30);
                            int nxmax = (int) (dmax*cos30minusalpha/stepSize)+1; 

                            for(int idist =0; idist<nxmax; idist++) {

                                double x = (double)(idist+1)*stepSize;
                                double timebfield = calc_Time( x,  alpha, bfield, s+1, r+1) ;
                                double deltatime_beta = util.getDeltaTimeBeta(x,betaValues[ibeta],distbeta[s][r],v0[s][r]);
                                timebfield+=deltatime_beta;

                                int tbin = Integer.parseInt(df.format(timebfield/2.) ) -1;

                                if(tbin<0 || tbin>nBinsT-1) {
                                    //System.err.println("Problem with tbin");
                                    continue;
                                }
                                if(tbin>maxTBin)
                                    maxTBin = tbin;
                                if(DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]==0) {
                                    DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]=x;
                                } else {
                                    DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]+=stepSize;
                                }
                                if(DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]>dmax) {
                                    DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]=dmax*cos30minusalpha;
                                    idist=nxmax;
                                }
                            }
                        }
                    }
                }
            }
        }	
        TableLoader.fillMissingTableBins();
        fillT2DGraphs();
        System.out.println(" T2D TABLE FILLED.....");
        //testBeq1();
        //test();
        T2DLOADED = true;
     }
    
    
    
    public static synchronized void ReFill() {
        //reset
        DISTFROMTIME = new double[6][6][maxBinIdxB+1][maxBinIdxAlpha+1][betaValues.length][nBinsT]; // sector slyr alpha Bfield time bins [s][r][ibfield][icosalpha][tbin]
        minBinIdxT  = 0;
        maxBinIdxT  = new int[6][6][8][6];
        maxTBin = -1;
        double stepSize = 0.0010;
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);
        
        for(int s = 0; s<6; s++ ){ // loop over sectors
                for(int r = 0; r<6; r++ ){ //loop over slys
                    // end fill constants
                    //System.out.println(v0[s][r]+" "+vmid[s][r]+" "+FracDmaxAtMinVel[s][r]);
                    double dmax = 2.*Constants.getInstance().wpdist[r]; 
                    //double tmax = CCDBConstants.getTMAXSUPERLAYER()[s][r];
                    for(int ibfield =0; ibfield<maxBinIdxB+1; ibfield++) {
                        double bfield = BfieldValues[ibfield];

                        for(int icosalpha =0; icosalpha<maxBinIdxAlpha+1; icosalpha++) {
                                maxBinIdxT[s][r][ibfield][icosalpha] = nBinsT; 
                                double cos30minusalpha = Math.cos(Math.toRadians(30.)) + (double) (icosalpha)*(1. - Math.cos(Math.toRadians(30.)))/5.;
                                double alpha = -(Math.toDegrees(Math.acos(cos30minusalpha)) - 30);
                                int nxmax = (int) (dmax*cos30minusalpha/stepSize); 
                                for(int ibeta=0; ibeta<betaValues.length; ibeta++) {
                                    for(int idist =0; idist<nxmax; idist++) {

                                        double x = (double)(idist+1)*stepSize;
                                        double timebfield = calc_Time( x,  alpha, bfield, s+1, r+1) ;
                                        double deltatime_beta = util.getDeltaTimeBeta(x,betaValues[ibeta],distbeta[s][r],v0[s][r]);
                                        timebfield+=deltatime_beta;
                                        
                                        int tbin = Integer.parseInt(df.format(timebfield/2.) ) -1;

                                        if(tbin<0 || tbin>nBinsT-1) {
                                            //System.err.println("Problem with tbin");
                                            continue;
                                        }
                                        if(tbin>maxTBin)
                                            maxTBin = tbin;
                                        if(DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]==0) {
                                            DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]=x;
                                        } else {
                                            DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]+=stepSize;
                                        }

                                    }
                                }
                            }
                        }
                }
        }	
        TableLoader.fillMissingTableBins();
        //TableLoader.test();
        fillT2DGraphs();
        System.out.println(" T2D TABLE RE-FILLED.....");
     }
    

    private static void fillMissingTableBins() {
        
        for(int s = 0; s<6; s++ ){ // loop over sectors

            for(int r = 0; r<6; r++ ){ //loop over slys
                
                for(int ibfield =0; ibfield<maxBinIdxB+1; ibfield++) {
                    
                    for(int icosalpha =0; icosalpha<maxBinIdxAlpha+1; icosalpha++) {
                        
                        for(int ibeta=0; ibeta<betaValues.length; ibeta++) {
                            
                            for(int tbin = 0; tbin<maxTBin; tbin++) {
                                
                                if(DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin]!=0 && DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin+1]==0) {
                                    DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin+1] = DISTFROMTIME[s][r][ibfield][icosalpha][ibeta][tbin];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    public static void fillT2DGraphs() {
        TimeToDistanceEstimator t2de = new TimeToDistanceEstimator();
        for(int r = 0; r<6; r++ ){ //loop over slys
            for(int ibfield =0; ibfield<1; ibfield++) {
                for(int icosalpha =0; icosalpha<maxBinIdxAlpha+1; icosalpha++) {
                    double cos30minusalpha = Math.cos(Math.toRadians(30.)) + (double) (icosalpha)*(1. - Math.cos(Math.toRadians(30.)))/5.;
                    double alpha = -(Math.toDegrees(Math.acos(cos30minusalpha)) - 30);
                    for(int ibeta=4; ibeta<5; ibeta++) {
                        for(int tbin = 0; tbin<maxTBin; tbin++) {
                            double time = 2*(tbin+1);
                            double doca = DISTFROMTIME[0][r][0][icosalpha][ibeta][tbin];
                            double calctime = calc_Time( doca,  alpha, 0, 1, r+1) ;
                            double deltatime_beta = util.getDeltaTimeBeta(doca,betaValues[ibeta],distbeta[0][r],v0[0][r]);
                            calctime+=deltatime_beta;
                            double calcdoca = t2de.interpolateOnGrid(0, alpha, 1, time,  0, r);
                            if(r==0 && time-calctime <4)System.out.println("alpha "+alpha+" time "+time+" calctime "+calctime+" doca "+doca+" calcdoca "+calcdoca);
                        }
                    }
                }
            }
        }
    }
    /**
     * 
     * @param x distance to wire in cm
     * @param alpha local angle in deg
     * @param bfield B field value a x in T
     * @param sector sector  
     * @param superlayer superlayer 
     * @return returns time (ns) when given inputs of distance x (cm), local angle alpha (degrees) and magnitude of bfield (Tesla).  
     */
    private static Utilities util = new Utilities();
    
    private static int getAlphaBinT2DC(double alpha) {
        int v = -1;
        for(int i = 0; i<T2DCalib.AlphaValues.length; i++) {
            
            if(Math.abs(alpha-T2DCalib.AlphaValues[i])<T2DCalib.AlphaBinHalfWidth)
                v = i;
        } 
        
        return v;
    }
    
    public static synchronized double calc_Time(double x, double alpha, double bfield, int sector, int superlayer) {
        int s = sector - 1;
        int r = superlayer - 1;
        double dmax = 2.*Constants.getInstance().wpdist[r]; 
        double tmax = Tmax[s][r];
        double delBf = delta_bfield_coefficient[s][r]; 
        double Bb1 = b1[s][r];
        double Bb2 = b2[s][r];
        double Bb3 = b3[s][r];
        double Bb4 = b4[s][r];
        if(x>dmax)
            x=dmax;
       return FitFunction.polyFcnMac(x, alpha, bfield, v0[s][r], vmid[s][r], FracDmaxAtMinVel[s][r], 
                tmax, dmax, delBf, Bb1, Bb2, Bb3, Bb4, superlayer) ;
        
    }
    
    public static double[][] delta_T0 = new double[6][6];
    public static double[][] delta_bfield_coefficient = new double[6][6];
    public static double[][] distbeta = new double[6][6];
    public static double[][] vmid = new double[6][6];
    public static double[][] v0 = new double[6][6];
    public static double[][] b1 = new double[6][6];
    public static double[][] b2 = new double[6][6];
    public static double[][] b3 = new double[6][6];
    public static double[][] b4 = new double[6][6];
    public static double[][] Tmax = new double[6][6];
    public static double[][] FracDmaxAtMinVel = new double[6][6];		// fraction of dmax corresponding to the point in the cell where the velocity is minimal

}
