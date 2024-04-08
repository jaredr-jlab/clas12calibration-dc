/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.io.FileNotFoundException;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.jlab.detector.calib.utils.ConstantsManager;

/**
 *
 * @author ziegler
 */
public class T2DFitter extends T2DCalib {
    
    public T2DFitter(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
    }
    
    private void fitWithFixedRDPars(boolean fixFit[][], double pars[], MnMigrad scanner[], MnMigrad fitter[]) { 
        
        for(int i =0; i<6; i++) {
            
            String s2="";
            s2+=(" ******************************************");
            s2+=("   RUNNING THE PARAMETER FIT FOR SUPERLAYER "+(i+1));
            s2+=(" ******************************************");
            //fMin fm = this.getfMinFixedRDPars(i, fixFit, scanner[i], fitter[i], pars[0], pars[1], true, false, s1);
            fMin fm2 = this.fMinFixedRDPars(i, fixFit, scanner[i], fitter[i], pars[0], pars[1], false, s2);
            FunctionMinimum fmin=null;
            if(fm2.getFcnMin().isValid()) {
                    fmin = fm2.getFcnMin();
                    System.out.println("UPDATED "+fmin.toString());
                TvstrkdocasFitPars.put(new Coordinate(i),fmin.userParameters()); 
            } 
        }
    }
    
    private double[] estimateRDPars(boolean fixFit[][], MnMigrad scanner[], MnMigrad fitter[]) { 
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
                    fMin fm = this.fMinFixedRDPars(i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
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
        for(int i =0; i<6; i++) {
            TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
                }
            }
        }
        return new double[] {bestR, bestDistbeta};
    }
    
    private fMin fMinFixedRDPars(int i, boolean fixFit[][], MnMigrad scanner, MnMigrad fitter, 
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
        
        
        scanner.fix(10);
        for (int p = 0; p < 10; p++) {
            if(fixFit[p][i]==true) {
                scanner.fix(p); 
            }
        }
       
        min = scanner.minimize();
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
                min = scanner.minimize();
                itercnt++;
                if(FitFunction.chi2<bestchi2) {
                    bestchi2 = FitFunction.chi2;
                    bestmin = min;
                    
                }
                if(edm-FitFunction.chi2<0.1 || FitFunction.chi2+10>edm) break;
                edm = FitFunction.chi2;
        } 
        System.out.println("MIN "+min.isValid());
        for (int p = 0; p < 10; p++) {
            fitter.setValue(p, bestmin.userParameters().value(p));
        }
        
        int itercnt2=0;
        for(int it = 0; it<maxIter; it++) {
            min2 = fitter.minimize();
            itercnt2++;
            System.out.println("Iteration "+itercnt2+" valid "+min2.isValid());
            if(min2.isValid()) {
                bestMchi2 = FitFunction.chi2;
                bestmin2 = min2;
                System.out.println(" iteration "+(it+1));
                System.out.println("CHI2   "+FitFunction.chi2);
                System.out.println("******SCAN RESULT******");
                System.out.println(String.valueOf(bestmin));
                System.out.println("******MIGRAD RESULT******");
                System.out.println(String.valueOf(min2));
            } else {
                min = scanner.minimize();
                for (int p = 0; p < 10; p++) {
                    fitter.setValue(p, min.userParameters().value(p));
                }
                min2 = fitter.minimize();
                bestMchi2 = FitFunction.chi2;
                bestmin2 = min2;
                System.out.println("******MIGRAD refit RESULT******");
                System.out.println("******M0 refit RESULT******");
                System.out.println(String.valueOf(min));
                System.out.println("******M1 refit RESULT******");
                System.out.println(String.valueOf(min2));
            }
            if(min2.isValid() && edm2-FitFunction.chi2<0.1) break;
            edm2 = FitFunction.chi2;
        }
        
        
        if(fixFit[2][i]==false) {
            scanner.release(2);
            fitter.release(2);
        }
        if(fixFit[4][i]==false) {
            scanner.release(4);
            fitter.release(4);
        }
        //release all so they can be fixed again
        scanner.release(10);
        fitter.release(10);
        System.out.println("INDEX "+fitter.index("dmax"));
        for (int p = 0; p < 10; p++) {
            if(fixFit[p][i]==true) {
                scanner.release(p); 
                fitter.release(p); 
            }
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
    
}
