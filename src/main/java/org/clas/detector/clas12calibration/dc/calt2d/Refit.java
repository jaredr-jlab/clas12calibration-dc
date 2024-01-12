/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import org.jlab.rec.dc.cluster.Cluster;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.hit.FittedHit;

/**
 *
 * @author ziegler
 */
public class Refit {
    
    public Refit() {
    }
    public List<FittedCluster> recomposeClusters(List<FittedHit> fhits, boolean useOT) {
        Map<Integer, ArrayList<FittedHit>> grpHits = new HashMap<Integer, ArrayList<FittedHit>>();
        List<FittedCluster> clusters = new ArrayList<FittedCluster>();
        for (FittedHit hit : fhits) { 
            //if(hit.get_OutOfTimeFlag()==true) System.out.println(hit.printInfo());
            if (hit.get_AssociatedClusterID() == -1 || hit.get_AssociatedHBTrackID() == -1 ) {
                continue;
            }
            if(useOT && hit.get_OutOfTimeFlag()) {
                continue;
            }
            if (hit.get_AssociatedClusterID() != -1 &&
                    hit.get_AssociatedHBTrackID() != -1) {
                int index = hit.get_AssociatedHBTrackID()*10000+hit.get_AssociatedClusterID();
               
                if(grpHits.get(index)==null) { // if the list not yet created make it
                    grpHits.put(index, new ArrayList<FittedHit>()); 
                    grpHits.get(index).add(hit); // append hit
                    //System.out.println("appended first hit "+hit.get_Sector()+","+hit.get_Superlayer()+", "+hit.get_Layer()+","+hit.get_TDC()+
                    //        " to cluster "+index);
                } else {
                    grpHits.get(index).add(hit); // append hit
                    //System.out.println("appended subs hit "+hit.get_Sector()+","+hit.get_Superlayer()+", "+hit.get_Layer()+","+hit.get_TDC()+
                    //        " to cluster "+index);
                }
            }
        }
        Iterator<Map.Entry<Integer, ArrayList<FittedHit>>> itr = grpHits.entrySet().iterator(); 
          
        while(itr.hasNext()) {
            Map.Entry<Integer, ArrayList<FittedHit>> entry = itr.next(); 
             
            if(entry.getValue().size()>Double.parseDouble(T2DViewer.nWires.getText())) {
                Cluster cluster = new Cluster(entry.getValue().get(0).get_Sector(), 
                        entry.getValue().get(0).get_Superlayer(), entry.getValue().get(0).get_AssociatedClusterID());
                FittedCluster fcluster = new FittedCluster(cluster);
                
                fcluster.addAll(entry.getValue());
                if(fcluster.size()>Double.parseDouble(T2DViewer.nWires.getText())) {
                    clusters.add(fcluster);
                }
            }
        }
        
        return clusters;
    }
    private ClusterFitter cf = new ClusterFitter();
    public void reFit(List<FittedHit> hits, boolean useOT) {
        List<FittedCluster> clusters = this.recomposeClusters(hits, useOT);
        for(FittedCluster clus : clusters) {
            FittedCluster LRresolvClus = this.LRAmbiguityResolver(clus, cf);
            clus = LRresolvClus;
            if (clus == null) {
                continue;
            }
            cf.SetFitArray(clus, "TSC");
            cf.Fit(clus, true);
            cf.SetResidualDerivedParams(clus, false, false, T2DViewer.dcDetector); //calcTimeResidual=false, resetLRAmbig=false 
            
            double trkAngle = clus.get_clusterLineFitSlope();
            
            for(FittedHit h : clus) {
                
                //local angle updated
                double theta0 = Math.toDegrees(Math.acos(1-0.02*h.getB()));
                double alpha = Math.toDegrees(Math.atan(trkAngle));
                // correct alpha with theta0, the angle corresponding to the isochrone lines twist due to the electric field
                alpha-=(double)T2DCalib.polarity*theta0;
                h.setAlpha(alpha);
                double cosTkAng = 1./Math.sqrt(trkAngle*trkAngle + 1.);
                h.set_X(h.get_XWire() + h.get_LeftRightAmb() * (h.get_Doca() / cosTkAng) );
                
            }
            //refit after updating alpha
            cf.SetFitArray(clus, "TSC");
            cf.Fit(clus, true);
            cf.SetResidualDerivedParams(clus, true, false, T2DViewer.dcDetector); //calcTimeResidual=false, resetLRAmbig=false 
            
        }
        clusters.clear();
    }
    
    public FittedCluster LRAmbiguityResolver(FittedCluster fClus, 
            ClusterFitter cf) {
        //	int[] notResolvedLR = {0,0,0,0,0,0};
        //	if(fClus.get_Status()[1]==notResolvedLR) {
        //		return fClus;
        //	}
        //
        cf.SetFitArray(fClus, "TSC");
        cf.Fit(fClus, true);
        cf.SetResidualDerivedParams(fClus, false, false, T2DViewer.dcDetector); //calcTimeResidual=false, resetLRAmbig=false 
            
        double trkAngle = fClus.get_clusterLineFitSlope();
        double cosTkAng = 1./Math.sqrt(trkAngle*trkAngle + 1.);
        int index = 0;
        for(FittedHit hit: fClus) {
            
            if (hit.get_Doca() < 0.4 * hit.get_CellSize()) {
                hit.set_LeftRightAmb(0);
            }
            if (hit.get_LeftRightAmb() == 0) {
                index++;
            }

        }
        if (index == 0) {
            return fClus; // cluster OK
        }
        if (index > 6) {
            return null; // unresolveable cluster...
        }
        int arraySize = (int) Math.pow(2, (double) index);
        ArrayList<FittedCluster> arrayOfClus = new ArrayList<>(arraySize);

        //pass all acceptable clusters
        FittedCluster okClus = new FittedCluster(fClus.getBaseCluster());
        for (FittedHit hit : fClus) {
            if (hit.get_LeftRightAmb() != 0) {
                okClus.add(hit);
            }
        }
        //filter all other clusters
        FittedCluster notLRClus = new FittedCluster(fClus.getBaseCluster());
        for (FittedHit hit : fClus) {
            if (hit.get_LeftRightAmb() == 0) {
                notLRClus.add(hit);
            }
        }
        //make combinatorials

        FittedCluster totNotLRClus = new FittedCluster(fClus.getBaseCluster());
        FittedCluster posNotLRClus = new FittedCluster(fClus.getBaseCluster());
        FittedCluster negNotLRClus = new FittedCluster(fClus.getBaseCluster());

        for (FittedHit hit : notLRClus) {

            FittedHit newhitPos = new FittedHit(hit.get_Sector(), hit.get_Superlayer(), hit.get_Layer(), hit.get_Wire(),
                    hit.get_TDC(), hit.getJitter(), hit.get_Id());
            newhitPos.set_Doca(hit.get_Doca());
            newhitPos.set_DocaErr(hit.get_DocaErr());
            newhitPos.setT0(hit.getT0()); 
            newhitPos.set_Beta(hit.get_Beta()); 
            newhitPos.setB(hit.getB()); 
            newhitPos.set_DeltaTimeBeta(hit.get_DeltaTimeBeta());
            newhitPos.set_DeltaDocaBeta(hit.get_DeltaDocaBeta());
            newhitPos.setTStart(hit.getTStart());
            newhitPos.setTProp(hit.getTProp());
            newhitPos.betaFlag= hit.betaFlag;
            newhitPos.setTFlight(hit.getTFlight());
            newhitPos.set_Time(hit.get_Time());
            newhitPos.set_Id(hit.get_Id());
            newhitPos.set_TrkgStatus(0);
            newhitPos.calc_CellSize(T2DViewer.dcDetector);
            newhitPos.set_XWire(hit.get_XWire());
            newhitPos.set_Z(hit.get_Z());
            newhitPos.set_WireLength(hit.get_WireLength());
            newhitPos.set_WireMaxSag(hit.get_WireMaxSag());
            newhitPos.set_WireLine(hit.get_WireLine());
            newhitPos.set_LeftRightAmb(1);
            newhitPos.set_X(newhitPos.get_XWire() + newhitPos.get_LeftRightAmb() * newhitPos.get_Doca()/cosTkAng ); // assume the track angle is // to the layer, so that cosTrkAng =1
            
            newhitPos.set_AssociatedClusterID(hit.get_AssociatedClusterID());
            newhitPos.set_AssociatedHBTrackID(hit.get_AssociatedHBTrackID());

            FittedHit newhitNeg = new FittedHit(hit.get_Sector(), hit.get_Superlayer(), hit.get_Layer(), hit.get_Wire(),
                    hit.get_TDC(), hit.getJitter(), hit.get_Id());
            newhitNeg.set_Doca(hit.get_Doca());
            newhitNeg.set_DocaErr(hit.get_DocaErr());
            newhitNeg.setT0(hit.getT0()); 
            newhitNeg.set_Beta(hit.get_Beta());  
            newhitNeg.setB(hit.getB());  
            newhitNeg.set_DeltaTimeBeta(hit.get_DeltaTimeBeta());
            newhitNeg.set_DeltaDocaBeta(hit.get_DeltaDocaBeta());
            newhitNeg.setTStart(hit.getTStart());
            newhitNeg.setTProp(hit.getTProp());
            newhitNeg.betaFlag= hit.betaFlag;
            newhitNeg.setTFlight(hit.getTFlight());
            newhitNeg.set_Time(hit.get_Time());
            newhitNeg.set_Id(hit.get_Id());
            newhitNeg.set_TrkgStatus(0);
            newhitNeg.calc_CellSize(T2DViewer.dcDetector);
            newhitNeg.set_XWire(hit.get_XWire());
            newhitNeg.set_Z(hit.get_Z());
            newhitNeg.set_WireLength(hit.get_WireLength());
            newhitNeg.set_WireMaxSag(hit.get_WireMaxSag());
            newhitNeg.set_WireLine(hit.get_WireLine());
            newhitNeg.set_LeftRightAmb(-1);
            newhitNeg.set_X(newhitNeg.get_XWire() + newhitNeg.get_LeftRightAmb() * newhitNeg.get_Doca()/cosTkAng );  // assume the track angle is // to the layer

            newhitNeg.set_AssociatedClusterID(hit.get_AssociatedClusterID());
            newhitNeg.set_AssociatedHBTrackID(hit.get_AssociatedHBTrackID());

            totNotLRClus.add(newhitNeg);
            totNotLRClus.add(newhitPos);

            posNotLRClus.add(newhitPos);
            negNotLRClus.add(newhitNeg);
        }

        Collections.sort(totNotLRClus);

        if (index == 1) {
            arrayOfClus.add(posNotLRClus);
            arrayOfClus.add(negNotLRClus);
            arrayOfClus.get(0).addAll(okClus);
            arrayOfClus.get(1).addAll(okClus);
        }
        if (index == 2) {
            for (int i1 = 0; i1 < totNotLRClus.size(); i1++) {
                for (int i2 = 2; i2 < totNotLRClus.size(); i2++) {
                    if (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i2).get_Id()) {
                        continue;
                    }
                    FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                    newClus.addAll(okClus);
                    newClus.add(totNotLRClus.get(i1));
                    newClus.add(totNotLRClus.get(i2));
                    arrayOfClus.add(newClus);
                }
            }
        }

        if (index == 3) {
            for (int i1 = 0; i1 < totNotLRClus.size(); i1++) {
                for (int i2 = 2; i2 < totNotLRClus.size(); i2++) {
                    for (int i3 = 4; i3 < totNotLRClus.size(); i3++) {
                        if ((totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i2).get_Id())
                                || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i3).get_Id())
                                || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i3).get_Id())) {
                            continue;
                        }
                        FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                        newClus.addAll(okClus);
                        newClus.add(totNotLRClus.get(i1));
                        newClus.add(totNotLRClus.get(i2));
                        newClus.add(totNotLRClus.get(i3));
                        arrayOfClus.add(newClus);
                    }
                }
            }
        }

        if (index == 4) {
            for (int i1 = 0; i1 < totNotLRClus.size(); i1++) {
                for (int i2 = 2; i2 < totNotLRClus.size(); i2++) {
                    for (int i3 = 4; i3 < totNotLRClus.size(); i3++) {
                        for (int i4 = 6; i4 < totNotLRClus.size(); i4++) {
                            if ((totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i2).get_Id())
                                    || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i3).get_Id())
                                    || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i4).get_Id())
                                    || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i3).get_Id())
                                    || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i4).get_Id())
                                    || (totNotLRClus.get(i3).get_Id() == totNotLRClus.get(i4).get_Id())) {
                                continue;
                            }
                            FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                            newClus.addAll(okClus);
                            newClus.add(totNotLRClus.get(i1));
                            newClus.add(totNotLRClus.get(i2));
                            newClus.add(totNotLRClus.get(i3));
                            newClus.add(totNotLRClus.get(i4));
                            arrayOfClus.add(newClus);
                        }
                    }
                }
            }
        }

        if (index == 5) {
            for (int i1 = 0; i1 < totNotLRClus.size(); i1++) {
                for (int i2 = 2; i2 < totNotLRClus.size(); i2++) {
                    for (int i3 = 4; i3 < totNotLRClus.size(); i3++) {
                        for (int i4 = 6; i4 < totNotLRClus.size(); i4++) {
                            for (int i5 = 8; i5 < totNotLRClus.size(); i5++) {
                                if ((totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i2).get_Id())
                                        || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i3).get_Id())
                                        || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i4).get_Id())
                                        || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i5).get_Id())
                                        || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i3).get_Id())
                                        || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i4).get_Id())
                                        || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i5).get_Id())
                                        || (totNotLRClus.get(i3).get_Id() == totNotLRClus.get(i4).get_Id())
                                        || (totNotLRClus.get(i3).get_Id() == totNotLRClus.get(i5).get_Id())
                                        || (totNotLRClus.get(i4).get_Id() == totNotLRClus.get(i5).get_Id())) {
                                    continue;
                                }
                                FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                                newClus.addAll(okClus);
                                newClus.add(totNotLRClus.get(i1));
                                newClus.add(totNotLRClus.get(i2));
                                newClus.add(totNotLRClus.get(i3));
                                newClus.add(totNotLRClus.get(i4));
                                newClus.add(totNotLRClus.get(i5));
                                arrayOfClus.add(newClus);
                            }
                        }
                    }
                }
            }
        }

        if (index == 6) {
            for (int i1 = 0; i1 < totNotLRClus.size(); i1++) {
                for (int i2 = 2; i2 < totNotLRClus.size(); i2++) {
                    for (int i3 = 4; i3 < totNotLRClus.size(); i3++) {
                        for (int i4 = 6; i4 < totNotLRClus.size(); i4++) {
                            for (int i5 = 8; i5 < totNotLRClus.size(); i5++) {
                                for (int i6 = 10; i6 < totNotLRClus.size(); i6++) {
                                    if ((totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i2).get_Id())
                                            || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i3).get_Id())
                                            || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i4).get_Id())
                                            || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i5).get_Id())
                                            || (totNotLRClus.get(i1).get_Id() == totNotLRClus.get(i6).get_Id())
                                            || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i3).get_Id())
                                            || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i4).get_Id())
                                            || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i5).get_Id())
                                            || (totNotLRClus.get(i2).get_Id() == totNotLRClus.get(i6).get_Id())
                                            || (totNotLRClus.get(i3).get_Id() == totNotLRClus.get(i4).get_Id())
                                            || (totNotLRClus.get(i3).get_Id() == totNotLRClus.get(i5).get_Id())
                                            || (totNotLRClus.get(i3).get_Id() == totNotLRClus.get(i6).get_Id())
                                            || (totNotLRClus.get(i4).get_Id() == totNotLRClus.get(i5).get_Id())
                                            || (totNotLRClus.get(i4).get_Id() == totNotLRClus.get(i6).get_Id())
                                            || (totNotLRClus.get(i5).get_Id() == totNotLRClus.get(i6).get_Id())) {
                                        continue;
                                    }
                                    FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                                    newClus.addAll(okClus);
                                    newClus.add(totNotLRClus.get(i1));
                                    newClus.add(totNotLRClus.get(i2));
                                    newClus.add(totNotLRClus.get(i3));
                                    newClus.add(totNotLRClus.get(i4));
                                    newClus.add(totNotLRClus.get(i5));
                                    newClus.add(totNotLRClus.get(i6));
                                    arrayOfClus.add(newClus);
                                }
                            }
                        }
                    }
                }
            }
        }

        return cf.BestClusterSelector(arrayOfClus, "TSC");

    }

    
    
    
}
