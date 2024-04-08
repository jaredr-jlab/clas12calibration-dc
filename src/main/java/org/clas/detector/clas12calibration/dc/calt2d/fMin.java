/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import org.freehep.math.minuit.FunctionMinimum;

/**
 *
 * @author ziegler
 */
public class fMin {

    public fMin(FunctionMinimum fM, double chi2) {
        _fcnMin = fM;
        _chi2 = chi2;
    }
    /**
     * @return the _fcnMin
     */
    public FunctionMinimum getFcnMin() {
        return _fcnMin;
    }

    /**
     * @param _fcnMin the _fcnMin to set
     */
    public void setFcnMin(FunctionMinimum _fcnMin) {
        this._fcnMin = _fcnMin;
    }

    /**
     * @return the _chi2
     */
    public double getChi2() {
        return _chi2;
    }

    /**
     * @param _chi2 the _chi2 to set
     */
    public void setChi2(double _chi2) {
        this._chi2 = _chi2;
    }
    private FunctionMinimum _fcnMin;
    private double _chi2;
}
