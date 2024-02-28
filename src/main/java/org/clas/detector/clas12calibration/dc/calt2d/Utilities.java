/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

//import org.jlab.rec.dc.Constants;

/**
 *
 * @author ziegler
 */
public class Utilities {
    
    public double calcDeltaDocaBeta(double doca, double distbeta, double scale, double bbeta){
        double beta=bbeta;
        if(beta>1)
            beta=1;
        double delta_doca = scale*0.5 * (distbeta *beta*beta) *(distbeta *beta*beta) *(distbeta *beta*beta) * doca / 
                ( (distbeta *beta*beta) *(distbeta *beta*beta)*(distbeta *beta*beta) + doca *doca *doca );
        
        return delta_doca;
    }
    
    
    
    public double getReducedAngle(double alpha) {
        double ralpha = 0;

        ralpha = Math.abs(Math.toRadians(alpha));

        while (ralpha > Math.PI / 3.) {
            ralpha -= Math.PI / 3.;
        }
        if (ralpha > Math.PI / 6.) {
            ralpha = Math.PI / 3. - ralpha;
        }

        return Math.toDegrees(ralpha);
    }  

    public double getDeltaTimeBeta(double x, double beta, double distbeta, double v_0) {
        /*
        double dtbscl = 0.5*(x/Math.sqrt(x*x+T2DCalib.DBF)+1);
       // dtbscl = 0.1;
        double dtb0 =(Math.sqrt(x * x + (distbeta * beta * beta) * (distbeta* beta * beta)) - x) / v_0;
        double dtb = dtbscl*dtb0;
        return dtb;
        */
        
        ////double a = beta*beta*distbeta*Math.pow(0.5, 1/3); 
       
        double value = (0.5*Math.pow(beta*beta*distbeta,3)*x/(Math.pow(beta*beta*distbeta,3)+x*x*x))/v_0;
        
        return value;
    }
    
    
}
