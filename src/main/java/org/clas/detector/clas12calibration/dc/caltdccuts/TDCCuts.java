/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.caltdccuts;

import java.util.HashMap;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

/**
 *
 * @author ziegler
 */
public class TDCCuts extends AnalysisMonitor {
    private int runNumber;
    private String analTabs = "Raw TDC";
    public TDCCuts(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames(analTabs);
        this.init(false, "TDC");
    }
    private Map<Coordinate, H1F> TDCHis        = new HashMap<Coordinate, H1F>(); 
    int nsl  = 6;
    int nsec = 6;
    public static final double[] tLow4T0Fits  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public static final double[] tHigh4T0Fits  = {380.0, 380.0, 680.0, 780.0, 1080.0, 1080.0}; 

     @Override
    public void createHistos() { 
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
                
                hNm = String.format("tdcS%dS%d", i + 1, j + 1);

                TDCHis.put(new Coordinate(i,j), new H1F(hNm, 150, tLow4T0Fits[j], tHigh4T0Fits[j])); 
                hgrps.addDataSet(TDCHis.get(new Coordinate(i, j)), 0);
                this.getDataGroup().add(hgrps, i+1, j+1,0);
            }
        }

        this.getDataGroup().add(hgrps,0,0,0);
        for (int i = 0; i < nsec; i++) {
            for (int j = 0; j < nsl; j++) {
                this.getCalib().addEntry(i+1,j+1,0);
                
            }
        }
        this.getCalib().fireTableDataChanged();
    }
    int count = 0;
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
            
            runNumber = newRun; 
        }
        
        // get hits property
        if(event.hasBank("DC::tdc") && !(event.hasBank("TimeBasedTrkg::TBHits") )
                && !(event.hasBank("HitBasedTrkg::HBHits") )) {
            DataBank bnkHits = event.getBank("DC::tdc");

            for (int j = 0; j < bnkHits.rows(); j++) {

                int sec = bnkHits.getInt("sector", j);
                int sl = (bnkHits.getInt("layer", j) - 1)/6 + 1;// layer goes from 1 to 36 in data
                int tdc = bnkHits.getInt("TDC", j);// wire goes from 1 to 112 in data

                    this.TDCHis.get(new Coordinate(sec-1, sl-1))
                        .fill((float)tdc);
            } 
        }
        if(event.hasBank("TimeBasedTrkg::TBHits")) {
            DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");

            for (int j = 0; j < bnkHits.rows(); j++) {

                int sec = bnkHits.getInt("sector", j);
                int sl = bnkHits.getInt("superlayer", j);// layer goes from 1 to 36 in data
                int tdc = bnkHits.getInt("TDC", j);// wire goes from 1 to 112 in data

                    this.TDCHis.get(new Coordinate(sec-1, sl-1))
                        .fill((float)tdc);
            } 
        }
    }
    
    public void Plot(int i , int j) {
        
        this.getAnalysisCanvas().getCanvas(analTabs).cd(0);
        this.getAnalysisCanvas().getCanvas(analTabs)
                .draw(this.TDCHis.get(new Coordinate(i, j)));

    }
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas(analTabs).setGridX(false);
        this.getAnalysisCanvas().getCanvas(analTabs).setGridY(false);
        this.getAnalysisCanvas().getCanvas(analTabs).update();
        
        
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
       for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                this.Plot(i, j);
            }
        }
       this.getCalib().fireTableDataChanged();  
    }
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        System.out.println(str_sector + " " + str_layer + " " );
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);

       if(group.hasItem(sector,layer,0)==true){
           this.Plot(sector-1, layer-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
   
    }
}
