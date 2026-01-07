/**
 * 
 * here, the total KMC time has been recorded: KMCTime_equilibrium && KMCTime_record = KMTTime
 * 
 * only KMCTime_record wiil be used to calculate the specific activity 
 * 
 */
package new_managers;


import java.io.File;

/**
 * @author liang
 *
 *
 * the adsorbate is *OH instead of *O
 */

/*
 * Created on October, 25, 2019
 *
 */

import java.util.Random;

//import PtNi_OH.RandomDecorate_Surface_PtOH_NiOH;
import matsci.Element;
import matsci.Species;
import matsci.engine.monte.BasicRecorder;
import matsci.io.clusterexpansion.PRIM;
import matsci.io.vasp.POSCAR;
import matsci.location.Coordinates;
import matsci.location.basis.AbstractBasis;
import matsci.location.basis.CartesianBasis;
import matsci.structure.Structure;
import matsci.structure.StructureBuilder;
import matsci.structure.decorate.function.AppliedDecorationFunction;
import matsci.structure.decorate.function.ce.AbstractAppliedCE;
import matsci.structure.decorate.function.ce.FastAppliedCE;
import matsci.structure.decorate.function.ce.GeneralAppliedCE;
import matsci.util.MSMath;

public class NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN_KMCTime {

  
  protected double m_AverageCurrent = 0;
  protected double m_AverageCurrentOccupied = 0;
  protected double m_AverageCurrentUnOccupied = 0;
  
  protected double m_AverageAdsCurrent = 0;
  
  protected double[] m_AverageOHBindingEnergySites;

  
  protected double m_AverageOHOccupancy = 0;
  protected double m_AverageOHBindingEnergy = 0;
  protected double m_AverageOHBindingEnergySq = 0;


  protected double m_AverageGCNFaceOHOccupancy = 0;
  protected double m_AverageGCNFaceOHBindingEnergy = 0;
  protected double m_AverageGCNFaceOHBindingEnergySq = 0;

  protected double m_AverageEdgeOHOccupancy = 0;
  protected double m_AverageEdgeOHBindingEnergy = 0;
  protected double m_AverageEdgeOHBindingEnergySq = 0;

  
  //for nanoparticle under grand canonical ensemble, the surface sites are not fixed (for PtNi surface, the surface sites are fixed)
  //protected boolean[] m_LastOccupiedOHPattern;

  
  protected final AbstractAppliedCE m_AppliedCE;

  protected static final Random m_Generator = new Random();

  protected double m_KMCTime = 0;
  protected double m_KMCTime_equil = 0;
  protected double m_KMCTime_record = 0;

  protected final Species m_OHSpec;
  
  protected final double m_GCN_Cutoff_FaceSites;
  
  protected final double m_KMCPreFactor;

  //protected final double m_KMCPreFactor_LC;
  
  //protected int m_KMCTotalNumEvents = 0;
  
  protected final int m_NumKMCIterations;
  
  protected final double m_EquilibraRatio_KMC;

  protected long m_NumAcceptedAdsEvents;
  protected long m_NumAcceptedDesEvents;
  
  protected int[] m_NumAcceptedAdsEventsSites; 
  //protected int[] m_NumAcceptedDesEventsSites; 
  
  protected final double m_MuOH_Inner_KMC;
  protected final double m_MuOH_Outer;
  
  protected final int m_numOHSites;
  protected final int m_numOHGCNFaceSites;
  protected final int m_numOHEdgeSites;
  
  protected final int[] m_OHSites;
  protected final int[] m_OHSitesCN;
  protected final double[] m_OHSitesGCN;

  protected final int[] m_OHStates;
  //protected final int[] m_VacancyStates;
  protected final int[] m_PtStates;
  //protected final int[] m_NiStates;
  
  //volcano plot peak for adsorbed *OH
   protected static final double m_VolcanoPeak = 1.03547445 + 0.118; //0.1 from JPCL (2012, 3, 2948-2951)
   //protected static final double m_VolcanoPeak = 1.03547445 + 0.1; //0.1 from JPCL (2012, 3, 2948-2951)
   //protected static final double m_VolcanoPeak = 1.03547445 + 0.136; //0.1 from JPCL (2012, 3, 2948-2951)


  
  //public NanoGCCoverageDepend_KMCCurrent_Recorder(AbstractAppliedCE appliedCE, double muOH_Inner_KMC, int numKMCIterations, double equilibraRatio_KMC, boolean[] lastOccupiedOxygenPattern, Species adsorbedOH) {
  public NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN_KMCTime(AbstractAppliedCE appliedCE, double muOH_Inner_KMC, int numKMCIterations, double equilibraRatio_KMC, int[] OHSites, int[] OHSitesCN, double[] OHSitesGCN, Species OH_Spec, int numOHEdgeSites, double GCN_Cutoff_FaceSites) {

    m_OHSpec = OH_Spec;
    
    m_GCN_Cutoff_FaceSites = GCN_Cutoff_FaceSites;
    
    m_AppliedCE = appliedCE;

    m_MuOH_Outer = appliedCE.getChemPot(OH_Spec);
    
    m_MuOH_Inner_KMC = muOH_Inner_KMC;
    
    m_AppliedCE.setChemPot(OH_Spec, m_MuOH_Inner_KMC);
    
    m_NumKMCIterations = numKMCIterations;
    
    m_EquilibraRatio_KMC = equilibraRatio_KMC;
    
    m_KMCTime = 0;
    m_KMCTime_equil = 0;
    m_KMCTime_record = 0;

    m_NumAcceptedAdsEvents = 0;
    m_NumAcceptedDesEvents = 0;

    m_numOHSites = OHSites.length;
    m_numOHGCNFaceSites = OHSites.length - numOHEdgeSites;
    m_numOHEdgeSites = numOHEdgeSites;
    
    System.out.println("\n\nwithin 2nd order recorder: ");
    System.out.println("# of surface sites for *OH, face(GCN>=" +  m_GCN_Cutoff_FaceSites + "), vertex/edge, and total: " + m_numOHGCNFaceSites + ", " + m_numOHEdgeSites + ", " + m_numOHSites + "\n\n");

    m_OHSites = OHSites;
    m_OHSitesCN = OHSitesCN;
    m_OHSitesGCN = OHSitesGCN;

    //m_VacancyStates = new int[m_OHSites.length];
    m_OHStates = new int[m_OHSites.length];
    m_PtStates = new int[m_OHSites.length];
    //m_NiStates = new int[m_OHSites.length];
    
    m_AverageOHBindingEnergySites = new double[m_OHSites.length];
    
    m_NumAcceptedAdsEventsSites = new int[m_OHSites.length];
    //m_NumAcceptedDesEventsSites = new int[numOSites];

    //m_LastOccupiedOHPattern = new boolean[numOHSites];
    
    //m_LastOccupiedOHPattern = lastOccupiedOxygenPattern;

    for(int ohSiteNum=0; ohSiteNum < m_OHSites.length; ohSiteNum++) {
        int sigmaIndex = m_OHSites[ohSiteNum];
        
        //m_VacancyStates[ohSiteNum] = appliedCE.getStateForSpecies(sigmaIndex, Species.vacancy);
        m_OHStates[ohSiteNum] = appliedCE.getStateForSpecies(sigmaIndex, m_OHSpec);
        m_PtStates[ohSiteNum] = appliedCE.getStateForSpecies(sigmaIndex, Species.platinum);
        //m_NiStates[ohSiteNum] = appliedCE.getStateForSpecies(sigmaIndex, Species.nickel);
        
        //System.out.println("ohSiteNum=" + ohSiteNum + ", OHState=" + m_OHStates[ohSiteNum]+ ", PtState=" + m_PtStates[ohSiteNum]+ ", NiState=" + m_NiStates[ohSiteNum] );
    }
    
    
    m_KMCPreFactor = 1.0 / (m_OHSites.length);

    //m_KMCPreFactor_LC = 1.0 / (this.m_numOHEdgeVertexSites);

    this.calculateAverages(); // this method already contains the function of clearSurfaceAdsorbates() method
    
    m_AppliedCE.setChemPot(m_OHSpec, m_MuOH_Outer); // set mu(O) back to the value for outer loop (Pt/Ni events) to make sure there is NO oxygen on the slab for outer loop. 
    
    this.clearSurfaceAdsorbates();
    
    System.out.println("\n\n**********************************************************************");
    System.out.println("Initialize 2nd-level recorder of KMC-current: NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN; muOH_Outer=" + m_MuOH_Outer + ", muOH_Inner_KMC=" + m_MuOH_Inner_KMC);
    System.out.println("peak of Sabatier volcano=" + m_VolcanoPeak);
    System.out.println("m_GCN_Cutoff_FaceSites=" + m_GCN_Cutoff_FaceSites);
    System.out.println("**********************************************************************\n\n");

    
  }


  
  public void calculateAverages() {
    
      m_AverageCurrent = 0;
      m_AverageCurrentOccupied = 0;
      m_AverageCurrentUnOccupied = 0;
      
      m_AverageAdsCurrent = 0;
      
      m_AverageOHOccupancy = 0;
      m_AverageOHBindingEnergy = 0;
      m_AverageOHBindingEnergySq = 0;
      
      m_AverageGCNFaceOHOccupancy = 0;
      m_AverageGCNFaceOHBindingEnergy = 0;
      m_AverageGCNFaceOHBindingEnergySq = 0;
      
      m_AverageEdgeOHOccupancy = 0;
      m_AverageEdgeOHBindingEnergy = 0;
      m_AverageEdgeOHBindingEnergySq = 0;
      
      this.m_KMCTime = 0;
      this.m_KMCTime_equil = 0;
      this.m_KMCTime_record = 0;

      this.m_NumAcceptedAdsEvents = 0;
      this.m_NumAcceptedDesEvents = 0;
      
      /*
      System.out.println("before KMC: ");
      System.out.println("this.m_KMCTime=" + this.m_KMCTime);
      System.out.println("this.m_NumAcceptedAdsEvents=" + this.m_NumAcceptedAdsEvents);
      System.out.println("this.m_NumAcceptedDesEvents=" + this.m_NumAcceptedDesEvents);
      */
      
      //for nanoparticle under grand canonical ensemble, the surface sites are not fixed (for PtNi surface, the surface sites are fixed)
      //initialize the oxygen pattern
      //int numOccupedO =0;
      /*
      for (int ohSiteNum = 0; ohSiteNum < m_OHSites.length; ohSiteNum++) { //go through oxygen sites on one side

          if(this.m_LastOccupiedOxygenPattern[ohSiteNum]) {
              int sigmaIndex = m_OxygenSites[ohSiteNum];
              int oxygenState = m_OxygenStates[ohSiteNum];
              m_AppliedCE.setSigmaState(sigmaIndex, OH_Spec);
              //numOccupedO++;
          }
      }
      */
      
      
      /*
      System.out.println("before the KMC, the occuped oxygen pattern from last time KMC: numOccupedO=" + numOccupedO);
      PRIM outfile = new PRIM(m_AppliedCE.getStructure(), true);
      outfile.useVASP5Format(true);
      outfile.writeFile("ce/fitECI_NiPt_AllRatio_rpbe_check/_KMC_tempStruct/beforeKMC-muPt=" + this.m_AppliedCE.getChemPot(Species.platinum) + "-muO=" + this.m_MuO_Inner_KMC+ ".vasp");
      */
      
      
      reject_free_KMC();

      this.m_KMCTime = this.m_KMCTime_equil + this.m_KMCTime_record;
      
      int local_numOHSites = 0;
      int local_numOHGCNFaceSites = 0;
      int local_numOHEdgeSites = 0;
      
      
      // the oxygen binding energy only record once after all KMC steps
      for (int ohSiteNum = 0; ohSiteNum < m_OHSites.length; ohSiteNum++) { //go through oxygen sites on one side
          local_numOHSites++;
          int sigmaIndex = m_OHSites[ohSiteNum];
          int OHState = m_OHStates[ohSiteNum];
          //int vacancyState = m_VacancyStates[ohSiteNum];
          int PtState = m_PtStates[ohSiteNum];

          if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == OHState) { //occupied by Pt-OH 
              
              //m_LastOccupiedOxygenPattern[oSiteNum] = true;
              
              m_AverageOHOccupancy += 1;   
              
              
              double delta_des = m_AppliedCE.getDelta(sigmaIndex, PtState); //Pt-OH --> Pt
              //double ohBindingEnergy = -1*delta_des + OHEnergyCorrectAccurate + 0.5*m_H2EnergyAccurate - m_H2OEnergyAccurate;
              double ohBindingEnergy = -1*delta_des;

              //double oBindingEnergy = m_AverageOBindingEnergySites[oSiteNum];
              m_AverageOHBindingEnergy += ohBindingEnergy;
              m_AverageOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
              

              if(m_OHSitesGCN[ohSiteNum] >= m_GCN_Cutoff_FaceSites) {
              //if(Math.abs(m_OHSitesGCN[ohSiteNum]-7.5)<=0.000001) {
              //if(m_OHSitesCN[ohSiteNum]==9) {
                  //System.out.println("occuped by *OH, CN=9: " + ohSiteNum + ", CN=" + m_OHSitesCN[ohSiteNum]);
                  local_numOHGCNFaceSites++;
                  
                  m_AverageGCNFaceOHOccupancy += 1;
                  m_AverageGCNFaceOHBindingEnergy += ohBindingEnergy;
                  m_AverageGCNFaceOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
              }
              else {
                  //System.out.println("occuped by *OH, CN<9: " + ohSiteNum + ", CN=" + m_OHSitesCN[ohSiteNum]);
                  local_numOHEdgeSites++;

                  m_AverageEdgeOHOccupancy += 1;
                  m_AverageEdgeOHBindingEnergy += ohBindingEnergy;
                  m_AverageEdgeOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
              }
          }
          
          else { // occupied by Pt-Vac; NOT occupied by Pt-OH 
              //m_LastOccupiedOxygenPattern[oSiteNum] = false;
              //m_AverageOHOccupancy += 1;   

              double delta_ads = m_AppliedCE.getDelta(sigmaIndex, OHState); //Pt --> Pt-OH
              //double ohBindingEnergy = delta_ads + OHEnergyCorrectAccurate + 0.5*m_H2EnergyAccurate - m_H2OEnergyAccurate;
              double ohBindingEnergy = delta_ads;

              //double oBindingEnergy = m_AverageOBindingEnergySites[oSiteNum];
              m_AverageOHBindingEnergy += ohBindingEnergy;
              m_AverageOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
              
              if(m_OHSitesGCN[ohSiteNum] >= m_GCN_Cutoff_FaceSites) {
              //if(Math.abs(m_OHSitesGCN[ohSiteNum]-7.5)<=0.000001) {
              //if(m_OHSitesCN[ohSiteNum]==9) {
                  //System.out.println("Un-occuped by *OH, CN=9: " + ohSiteNum + ", CN=" + m_OHSitesCN[ohSiteNum]);
                  local_numOHGCNFaceSites++;

                  m_AverageGCNFaceOHBindingEnergy += ohBindingEnergy;
                  m_AverageGCNFaceOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
              }
              else {
                  //System.out.println("Un-occuped by *OH, CN<9: " + ohSiteNum + ", CN=" + m_OHSitesCN[ohSiteNum]);
                  local_numOHEdgeSites++;

                  m_AverageEdgeOHBindingEnergy += ohBindingEnergy;
                  m_AverageEdgeOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
              }
          }

          
      } 
      
    
    System.out.println("# of local accunted face(GCN>=" +  m_GCN_Cutoff_FaceSites + "), edge/vertex, & all surface sites: " + local_numOHGCNFaceSites + ",    " + local_numOHEdgeSites + ",    " + local_numOHSites + "\n\n");  
    System.out.println("# of surface sites for *OH, face(GCN>=" +  m_GCN_Cutoff_FaceSites + "), vertex/edge, and total: " + m_numOHGCNFaceSites + ", " + m_numOHEdgeSites + ", " + m_numOHSites + "\n\n");
  
    m_AverageCurrent = this.m_KMCPreFactor * this.m_NumAcceptedDesEvents / this.m_KMCTime_record;
    m_AverageAdsCurrent = this.m_KMCPreFactor * this.m_NumAcceptedAdsEvents / this.m_KMCTime_record;
    
    m_AverageOHBindingEnergy /= (m_OHSites.length);
    m_AverageOHBindingEnergySq /= (m_OHSites.length);
    m_AverageOHOccupancy /= (m_OHSites.length);   

    m_AverageGCNFaceOHBindingEnergy /= (m_numOHGCNFaceSites);
    m_AverageGCNFaceOHBindingEnergySq /= (m_numOHGCNFaceSites);
    m_AverageGCNFaceOHOccupancy /= (m_numOHGCNFaceSites);   

    m_AverageEdgeOHBindingEnergy /= (m_numOHEdgeSites);
    m_AverageEdgeOHBindingEnergySq /= (m_numOHEdgeSites);
    m_AverageEdgeOHOccupancy /= (m_numOHEdgeSites);   

    

    

    
    
    /*
    //count the contribution of current from occupied sites or un-occupied sites
    for (int oSiteNum = 0; oSiteNum < m_OxygenSites.length/2; oSiteNum++) { //go through oxygen sites on one side
        int sigmaIndex = m_OxygenSites[oSiteNum];
        int oxygenState = m_OxygenStates[oSiteNum];

        //System.out.println(oSiteNum+ ": m_NumAcceptedAdsEventsSites[oSiteNum]="+ m_NumAcceptedAdsEventsSites[oSiteNum] + ", m_NumAcceptedDesEventsSites[oSiteNum]="+ m_NumAcceptedDesEventsSites[oSiteNum]);
        
        if((this.m_NumAcceptedAdsEventsSites[oSiteNum] + this.m_NumAcceptedDesEventsSites[oSiteNum])<1) {
            continue;
        }
        
        
        //double occupiedRatio = ((double)this.m_NumAcceptedAdsEventsSites[oSiteNum]) / (this.m_NumAcceptedAdsEventsSites[oSiteNum] + this.m_NumAcceptedDesEventsSites[oSiteNum]);
        //m_AverageCurrentOccupied += (occupiedRatio*this.m_NumAcceptedDesEventsSites[oSiteNum]);    
        //m_AverageCurrentUnOccupied += ((1-occupiedRatio)*this.m_NumAcceptedDesEventsSites[oSiteNum]);    
        
        
        if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == oxygenState) { //occupied by oxygen 
            m_AverageCurrentOccupied += this.m_NumAcceptedDesEventsSites[oSiteNum];    
        }
        else {
            m_AverageCurrentUnOccupied += this.m_NumAcceptedDesEventsSites[oSiteNum];    
        }
    }
    
    m_AverageCurrentOccupied /= this.m_KMCTime;
    m_AverageCurrentUnOccupied /= this.m_KMCTime;
    */
    
    
    
    System.out.println("after KMC: ");
    System.out.println("this.m_KMCTime=" + this.m_KMCTime);
    System.out.println("this.m_KMCTime_equil=" + this.m_KMCTime_equil);
    System.out.println("this.m_KMCTime_record=" + this.m_KMCTime_record);
    System.out.println("this.m_NumAcceptedAdsEvents=" + this.m_NumAcceptedAdsEvents);
    System.out.println("this.m_NumAcceptedDesEvents=" + this.m_NumAcceptedDesEvents);
    System.out.println("this.m_AverageCurrent=" + m_AverageCurrent);
    System.out.println("this.m_AverageOHOccupancy=" +  m_AverageOHOccupancy);
    System.out.println("m_AverageOBindingEnergy=" +  m_AverageOHBindingEnergy);
    
    
  }
  
  
  
  
  
  public void calculateAverages_old() {
    
      m_AverageCurrent = 0;
      m_AverageCurrentOccupied = 0;
      m_AverageCurrentUnOccupied = 0;
      m_AverageAdsCurrent = 0;
      m_AverageOHOccupancy = 0;
      m_AverageOHBindingEnergy = 0;
      m_AverageOHBindingEnergySq = 0;
      
      this.m_KMCTime = 0;
      this.m_KMCTime_equil = 0;
      this.m_KMCTime_record = 0;

      this.m_NumAcceptedAdsEvents = 0;
      this.m_NumAcceptedDesEvents = 0;
      
      /*
      System.out.println("before KMC: ");
      System.out.println("this.m_KMCTime=" + this.m_KMCTime);
      System.out.println("this.m_NumAcceptedAdsEvents=" + this.m_NumAcceptedAdsEvents);
      System.out.println("this.m_NumAcceptedDesEvents=" + this.m_NumAcceptedDesEvents);
      */
      
      reject_free_KMC();

      this.m_KMCTime = this.m_KMCTime_equil + this.m_KMCTime_record;

      
    /*  
    // the oxygen binding energy only record once after all KMC steps
    for (int oSiteNum = 0; oSiteNum < m_OxygenSites.length/2; oSiteNum++) { //go through oxygen sites on one side
        int sigmaIndex = m_OxygenSites[oSiteNum];
        int oxygenState = m_OxygenStates[oSiteNum];
        if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == oxygenState) { //occupied by oxygen 
            m_AverageOOccupancy += 1;    
            double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_VacancyStates[oSiteNum]) + m_OChemPot;
            double oBindingEnergy = -1*delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;
            m_AverageOBindingEnergy += oBindingEnergy;
            m_AverageOBindingEnergySq += oBindingEnergy * oBindingEnergy;
            this.m_AverageOBindingEnergySites[oSiteNum] = oBindingEnergy;
        }
        else { //occupied by Vacancy
            double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_OxygenStates[oSiteNum]) - m_OChemPot;
            double oBindingEnergy = delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;
            m_AverageOBindingEnergy += oBindingEnergy;
            m_AverageOBindingEnergySq += oBindingEnergy * oBindingEnergy;
            this.m_AverageOBindingEnergySites[oSiteNum] = oBindingEnergy;
        }
    }     
    */
      
      
      // the oxygen binding energy recorded for each KMC step
      int numOHSitesAdsorbed = 0;
      for (int ohSiteNum = 0; ohSiteNum < m_OHSites.length; ohSiteNum++) { //go through oxygen sites on one side
          int sigmaIndex = m_OHSites[ohSiteNum];
          int OHState = m_OHStates[ohSiteNum];
          //int vacancyState = m_VacancyStates[ohSiteNum];
          int PtState = m_PtStates[ohSiteNum];

          if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == OHState) { //occupied by oxygen 
              m_AverageOHOccupancy += 1;   
              
              /*
              double delta_des = m_AppliedCE.getDelta(sigmaIndex, vacancyState);
              double oBindingEnergy = -1*delta_des + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;
              
              //double oBindingEnergy = m_AverageOBindingEnergySites[oSiteNum];
              m_AverageOBindingEnergy += oBindingEnergy;
              m_AverageOBindingEnergySq += oBindingEnergy * oBindingEnergy;
              */
          }
          
          
          if(m_NumAcceptedAdsEventsSites[ohSiteNum]<=0) {
              continue;
          }
          
          numOHSitesAdsorbed++;
                    
          double ohBindingEnergy = m_AverageOHBindingEnergySites[ohSiteNum];
          
          
          m_AverageOHBindingEnergy += ohBindingEnergy;
          m_AverageOHBindingEnergySq += ohBindingEnergy * ohBindingEnergy;
          

      } 
      
      
    m_AverageCurrent = this.m_KMCPreFactor * this.m_NumAcceptedDesEvents / this.m_KMCTime_record;
    m_AverageAdsCurrent = this.m_KMCPreFactor * this.m_NumAcceptedAdsEvents / this.m_KMCTime_record;

    //m_AverageOBindingEnergy /= m_AverageOOccupancy;
    //m_AverageOBindingEnergySq /= m_AverageOOccupancy;
    //m_AverageOBindingEnergy /= (m_OxygenSites.length/2);
    //m_AverageOBindingEnergySq /= (m_OxygenSites.length/2);
    
    m_AverageOHBindingEnergy /= numOHSitesAdsorbed;
    m_AverageOHBindingEnergySq /= numOHSitesAdsorbed;
    
    m_AverageOHOccupancy /= (m_OHSites.length);   
    

    

    
    
    /*
    //count the contribution of current from occupied sites or un-occupied sites
    for (int oSiteNum = 0; oSiteNum < m_OxygenSites.length/2; oSiteNum++) { //go through oxygen sites on one side
        int sigmaIndex = m_OxygenSites[oSiteNum];
        int oxygenState = m_OxygenStates[oSiteNum];

        //System.out.println(oSiteNum+ ": m_NumAcceptedAdsEventsSites[oSiteNum]="+ m_NumAcceptedAdsEventsSites[oSiteNum] + ", m_NumAcceptedDesEventsSites[oSiteNum]="+ m_NumAcceptedDesEventsSites[oSiteNum]);
        
        if((this.m_NumAcceptedAdsEventsSites[oSiteNum] + this.m_NumAcceptedDesEventsSites[oSiteNum])<1) {
            continue;
        }
        
        
        //double occupiedRatio = ((double)this.m_NumAcceptedAdsEventsSites[oSiteNum]) / (this.m_NumAcceptedAdsEventsSites[oSiteNum] + this.m_NumAcceptedDesEventsSites[oSiteNum]);
        //m_AverageCurrentOccupied += (occupiedRatio*this.m_NumAcceptedDesEventsSites[oSiteNum]);    
        //m_AverageCurrentUnOccupied += ((1-occupiedRatio)*this.m_NumAcceptedDesEventsSites[oSiteNum]);    
        
        
        if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == oxygenState) { //occupied by oxygen 
            m_AverageCurrentOccupied += this.m_NumAcceptedDesEventsSites[oSiteNum];    
        }
        else {
            m_AverageCurrentUnOccupied += this.m_NumAcceptedDesEventsSites[oSiteNum];    
        }
    }
    
    m_AverageCurrentOccupied /= this.m_KMCTime;
    m_AverageCurrentUnOccupied /= this.m_KMCTime;
    */
    
    
    
    System.out.println("after KMC: ");
    System.out.println("this.m_KMCTime=" + this.m_KMCTime);
    System.out.println("this.m_KMCTime_equil=" + this.m_KMCTime_equil);
    System.out.println("this.m_KMCTime_record=" + this.m_KMCTime_record);
    System.out.println("this.m_NumAcceptedAdsEvents=" + this.m_NumAcceptedAdsEvents);
    System.out.println("this.m_NumAcceptedDesEvents=" + this.m_NumAcceptedDesEvents);
    System.out.println("this.m_AverageCurrent=" + m_AverageCurrent);
    System.out.println("this.m_AverageOHOccupancy=" +  m_AverageOHOccupancy);
    System.out.println("numOHSitesAdsorbed=" +  numOHSitesAdsorbed);
    System.out.println("m_AverageOHBindingEnergy=" +  m_AverageOHBindingEnergy);
    
    
  }
  
  
  /**
   * this clear the surface, so that the slab used for Pt/Ni events will not have surface adsorbates 
   */
  public void clearSurfaceAdsorbates() { //Pt-OH --> Pt
      
      //int num_OccupiedO = 0;
      for (int ohSiteNum = 0; ohSiteNum < m_OHSites.length; ohSiteNum++) { //go through oxygen sites on one side
          int sigmaIndex = m_OHSites[ohSiteNum];
          int OHState = m_OHStates[ohSiteNum];
          if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == OHState) {
              //num_OccupiedO++;
              m_AppliedCE.setSigmaState(sigmaIndex, this.m_PtStates[ohSiteNum]);
          }
      }   
  };
  
  
  
  public void reject_free_KMC() {
      
      this.m_AverageOHBindingEnergySites = new double[this.m_OHSites.length];
      
      m_NumAcceptedAdsEventsSites = new int[this.m_OHSites.length];
      //m_NumAcceptedDesEventsSites = new int[this.m_OxygenSites.length];

      double[] rateSums = new double[this.m_OHSites.length];

      for (int iterNum = 0; iterNum < this.m_NumKMCIterations; iterNum++){
          
          double rateSum = 0;
          
          for (int ohSiteNum = 0; ohSiteNum < this.m_OHSites.length;ohSiteNum++) {
            int sigmaIndex = this.m_OHSites[ohSiteNum];
            double rate = KMC_getRate(sigmaIndex, ohSiteNum);
            rateSum += rate;
            rateSums[ohSiteNum] = rateSum;
          }

          double triggerSum = rateSum * this.m_Generator.nextDouble();
          
          int triggerIndex = 0;
          int triggerOHSiteNum = 0;
          
          // Find the largest triggerIndex for which triggerSum<= rateSums[triggerIndex].  You could do this with a binarysearch, but simply cycling through the array might be fast enough.
          for (int ohSiteNum = 0; ohSiteNum < this.m_OHSites.length;ohSiteNum++) {
                            
              if (Double.compare(triggerSum, rateSums[ohSiteNum])>0) {
              //if(triggerSum > rateSums[oSiteNum]) {
                  continue;
              }
              else {
                  triggerOHSiteNum = ohSiteNum;
                  triggerIndex = this.m_OHSites[ohSiteNum];
                  //System.out.println("selected trigger event: rateSum=" + rateSum + ", triggerSum=" + triggerSum + ", oSiteNum=" + oSiteNum + ", triggerSigmaIndex=" + triggerIndex);
                  break;
              }
          }

          if(iterNum<this.m_NumKMCIterations*this.m_EquilibraRatio_KMC) {
              
              KMC_trigger_equilibra(triggerIndex, triggerOHSiteNum);  
              
              this.m_KMCTime_equil += (1.0/rateSum) * Math.log(1.0/this.m_Generator.nextDouble()); //rejection-free KMC
              
              
              //TODO: each 25% of edge sites of iterNum, print the snapshot with *OH
              if(iterNum<4*m_numOHSites) {
                  //int printSnapshot_deltaT = (int)(m_numOHEdgeSites*0.25);
                  //int printSnapshot_deltaT = 50;
                  int printSnapshot_deltaT = 100;

                  if((iterNum%printSnapshot_deltaT)==0) {
                      PRIM outfile = new PRIM(m_AppliedCE.getStructure(), true);
                      outfile.useVASP5Format(true);
                      String fileDir = "ce/PtNiOH/ORR-prediction/_KMC_tempStruct/";
                      
                      File dirFile  = new  File(fileDir);
                      if ( ! (dirFile.exists())  &&   ! (dirFile.isDirectory())) {
                              boolean  creadok  =  dirFile.mkdirs();
                      }
                      int numS = m_AppliedCE.getStructure().numDefiningSitesWithElement(Element.sulfur);
                      int numPt = m_AppliedCE.getStructure().numDefiningSitesWithElement(Element.platinum);
                      int numNi = m_AppliedCE.getStructure().numDefiningSitesWithElement(Element.nickel);
                      
                      //String structName = "Pt" + (numPt+numS) + "Ni"+ numNi + "-OHPattern-muPt=" + this.m_AppliedCE.getChemPot(Species.platinum) + "-muOH=" + this.m_MuOH_Inner_KMC;
                      String structName = "Pt" + (numPt+numS) + "Ni"+ numNi + "-OHPattern-iterNum=" + iterNum + "-KMCTime=" + this.m_KMCTime;

                      outfile.writeFile(fileDir + structName + ".vasp");

                      //Structure structure = this.m_AppliedCE.getStructure();
                      
                      //RandomDecorate_Surface_PtOH_NiOH.replaceSulfurWithPtOH_2(fileDir, structName);
                      replaceSulfurWithPtOH_2_local(fileDir, structName);
                  }
              }// end of if(iterNum<3*m_numOHSites) 
              
              
              
              continue;
          }
          
          else {
              
              KMC_trigger(triggerIndex, triggerOHSiteNum);  
              
              // Update the elapsed KMC time 
              //this.m_KMCTime += (1.0/m_KMCTotalNumEvents) * Math.log(1.0/this.m_Generator.nextDouble()); //rejection KMC
              this.m_KMCTime_record += (1.0/rateSum) * Math.log(1.0/this.m_Generator.nextDouble()); //rejection-free KMC

          }

          /*
          //Ideally we should also average the oxygen binding energy over the steps of the KMC run
          for (int oSiteNum = 0; oSiteNum < m_OxygenSites.length/2; oSiteNum++) { //go through oxygen sites on one side
              int sigmaIndex = m_OxygenSites[oSiteNum];
              //int oxygenState = m_OxygenStates[oSiteNum];
              int vacancyState = m_VacancyStates[oSiteNum];
              if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == vacancyState) { //occupied by Vacancy 
                  //double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_OxygenStates[oSiteNum]) - m_MuO_Outer;
                  double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_OxygenStates[oSiteNum]);
                  //double oBindingEnergy = delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy + m_MuO_Inner_KMC;
                  double oBindingEnergy = delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy ;
                  this.m_AverageOBindingEnergySites[oSiteNum] += oBindingEnergy;
              }
              else { //occupied by Oxygen
                  //double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_VacancyStates[oSiteNum]) + m_MuO_Outer;
                  double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_VacancyStates[oSiteNum]);
                  //double oBindingEnergy = -1*delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy + m_MuO_Inner_KMC;
                  double oBindingEnergy = -1*delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;
                  this.m_AverageOBindingEnergySites[oSiteNum] += oBindingEnergy;
              }
          } 
          */
          
      }
      
      //use to record the ObingingEnergy for each site
      /*
      for (int oSiteNum = 0; oSiteNum < m_OxygenSites.length/2; oSiteNum++) { //go through oxygen sites on one side 
          
          //this.m_AverageOBindingEnergySites[oSiteNum] /= this.m_NumKMCIterations;

          if(m_NumAcceptedAdsEventsSites[oSiteNum]<=0) {
              continue;
          }
          else {
              this.m_AverageOBindingEnergySites[oSiteNum] /= m_NumAcceptedAdsEventsSites[oSiteNum];
          }
      }
      */
      
      
      //TODO
      PRIM outfile = new PRIM(m_AppliedCE.getStructure(), true);
      outfile.useVASP5Format(true);
      String fileDir = "ce/PtNiOH/ORR-prediction/_KMC_tempStruct/";
      
      File dirFile  = new  File(fileDir);
      if ( ! (dirFile.exists())  &&   ! (dirFile.isDirectory())) {
              boolean  creadok  =  dirFile.mkdirs();
      }
      
      int numS = m_AppliedCE.getStructure().numDefiningSitesWithElement(Element.sulfur);
      int numPt = m_AppliedCE.getStructure().numDefiningSitesWithElement(Element.platinum);
      int numNi = m_AppliedCE.getStructure().numDefiningSitesWithElement(Element.nickel);
      
      //String structName = "Pt" + (numPt+numS) + "Ni"+ numNi + "-OHPattern-muPt=" + this.m_AppliedCE.getChemPot(Species.platinum) + "-muOH=" + this.m_MuOH_Inner_KMC;
      String structName = "Pt" + (numPt+numS) + "Ni"+ numNi + "-OHPattern-muOH=" + this.m_MuOH_Inner_KMC;

      outfile.writeFile(fileDir + structName + ".vasp");

      //Structure structure = this.m_AppliedCE.getStructure();
      
      //RandomDecorate_Surface_PtOH_NiOH.replaceSulfurWithPtOH_2(fileDir, structName);
      replaceSulfurWithPtOH_2_local(fileDir, structName);

  }
  

  public double KMC_getRate(int sigmaIndex, int ohSiteNum){
      
      Species spec = m_AppliedCE.getSpecies(sigmaIndex);

      if(spec==Species.nickel) {
          //System.out.println("\nthe surface spec within the method of KMC_getRate: " + spec.toString() + ", sigmaIndex=" + sigmaIndex + ", ohSiteNum=" + ohSiteNum);
      }
      
      //if (spec==Species.vacancy) { // adsorption event: Vac --> Oxygen, for Pt/Ni-O/Vac slab model
      if (spec==Species.platinum) { // adsorption event: Pt --> Pt-OH
    
          //m_NextState is oxygen 
          //double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_OxygenStates[oSiteNum]) - m_MuO_Outer;
          double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_OHStates[ohSiteNum]);
          //double delta = m_AppliedCE.getDelta(sigmaIndex, 0);

          //double oBindingEnergy = delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy + m_ShiftAdsEnergy;
          //double oBindingEnergy = delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy + m_MuO_Inner_KMC;
          //double oBindingEnergy = delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;
          //double ohBindingEnergy = delta + OHEnergyCorrectAccurate + 0.5*m_H2EnergyAccurate - m_H2OEnergyAccurate;
          double ohBindingEnergy = delta;

          //double barrier_ads = 0.297 + 0.53*(oBindingEnergy - m_VolcanoPeak);
          //double barrier_ads = 0.297 + Math.max(0, 0.53*(oBindingEnergy - m_VolcanoPeak));
          
          //double barrier_ads = Math.max(0, 1.06*(ohBindingEnergy - m_VolcanoPeak));
          
          //TODO
          //double barrier_ads = 0.297 + Math.max(0, 1.06*(ohBindingEnergy - m_VolcanoPeak));
          double barrier_ads = 0.297 + Math.max(0, 1.0*(ohBindingEnergy - m_VolcanoPeak));

          
          //double barrier_ads = 0.297 + Math.max(1.0*(m_VolcanoPeak - ohBindingEnergy), 1.06*(ohBindingEnergy - m_VolcanoPeak));

          //double barrier_ads = 0.297 + Math.max(0.5*(m_VolcanoPeak - oBindingEnergy), 0.53*(oBindingEnergy - m_VolcanoPeak));
          //double barrier_ads = 0.53*(-oBindingEnergy + m_VolcanoPeak);
          //double barrier_ads = oBindingEnergy;
          
          //System.out.println(ohSiteNum + ": adsorption event: Pt --> Pt-OH, delta=" + delta+ ", Eads(*OH)=" + ohBindingEnergy +  ", barrier_ads=" + barrier_ads + ", getRate()=" + (Math.exp(-barrier_ads / 0.0258519972)));
          //System.out.println("site: this.m_OHStates[ohSiteNum]=" + this.m_OHStates[ohSiteNum] + ", delta=" + delta);
          //System.out.println("adsorption event: Pt --> Pt-OH; site: this.m_OHStates[ohSiteNum]=" + this.m_OHStates[ohSiteNum] + ", delta=" + delta+ ", Eads(*OH)=" + ohBindingEnergy +  ", barrier_ads=" + barrier_ads + ", getRate()=" + (Math.exp(-barrier_ads / 0.0258519972)));

          return Math.exp(-1 * barrier_ads / 0.0258519972);
      }
      
      //else if(spec == this.m_OHSpec){ // desorption  event: Pt-OH --> Pt
      else { // desorption  event: Pt-OH --> Pt
    
          //m_NextState is vacancy 
          //double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_VacancyStates[oSiteNum]) + m_MuO_Outer;
          //double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_VacancyStates[ohSiteNum]);// for Pt/Ni-O/Vac slab model
          double delta = m_AppliedCE.getDelta(sigmaIndex, this.m_PtStates[ohSiteNum]);
          //double delta = m_AppliedCE.getDelta(sigmaIndex, 2);

          //double oBindingEnergy = -1*delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy + m_ShiftAdsEnergy;
          //double oBindingEnergy = -1*delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy + m_MuO_Inner_KMC;
          //double oBindingEnergy = -1*delta + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;
          //double ohBindingEnergy = -1*delta + OHEnergyCorrectAccurate + 0.5*m_H2EnergyAccurate - m_H2OEnergyAccurate;
          double ohBindingEnergy = -1*delta;

          //double barrier_des = 0.297 + 0.5*(m_VolcanoPeak - oBindingEnergy);
          
          //double barrier_des = Math.max(0, 1.0*(m_VolcanoPeak - ohBindingEnergy));
          double barrier_des = 0.297 + Math.max(0, 1.0*(m_VolcanoPeak - ohBindingEnergy));
          //double barrier_des = 0.297 + Math.max(1.0*(m_VolcanoPeak - ohBindingEnergy), 1.06*(ohBindingEnergy - m_VolcanoPeak));

          //double barrier_des = 0.297 + Math.max(0.5*(m_VolcanoPeak - oBindingEnergy), 0.53*(oBindingEnergy - m_VolcanoPeak));
          //double barrier_des = 0.5*(-m_VolcanoPeak + oBindingEnergy);
          //double barrier_des =  -1 * oBindingEnergy;

          //System.out.println(ohSiteNum + ": desorption  event: Pt-OH --> Pt, delta=" + delta+ ", Eads(*OH)=" + ohBindingEnergy + ", barrier_des=" + barrier_des + ", getRate()=" + (Math.exp(-barrier_des / 0.0258519972)));
          //System.out.println("desorption  event: Pt-OH --> Pt; site: this.m_OHStates[ohSiteNum]=" + this.m_OHStates[ohSiteNum] + ", delta=" + delta+ ", Eads(*OH)=" + ohBindingEnergy + ", barrier_des=" + barrier_des + ", getRate()=" + (Math.exp(-barrier_des / 0.0258519972)));

          return Math.exp(-1 * barrier_des / 0.0258519972);

      }
      
      /*
      else {
          
          System.out.println("ERROR: the surface site is occupied by element other than Pt and Pt-OH: " + spec.toString());
          return 0;
      }
      */
      
  } // end of getRate() method
  
  
  
  public void KMC_trigger(int sigmaIndex, int ohSiteNum) {
      
      Species spec = m_AppliedCE.getSpecies(sigmaIndex);

      //if (spec==Species.vacancy) { // adsorption event: Vac --> Oxygen, for Pt/Ni-O/Vac slab model
      if (spec==Species.platinum) { // adsorption event: Pt --> Pt-OH, for Pt/Ni/S/Vac NP model (S stands for Pt-OH)

          m_AppliedCE.setSigmaState(sigmaIndex, this.m_OHStates[ohSiteNum]);
          

          /*
          double delta_Ads = m_AppliedCE.setSigmaStateGetDelta(sigmaIndex, this.m_OxygenStates[oSiteNum]);
          
          double oBindingEnergy = delta_Ads + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;

          this.m_AverageOBindingEnergySites[oSiteNum] += oBindingEnergy;
          this.m_NumAcceptedAdsEventsSites[oSiteNum]++;
          */
          
          this.m_NumAcceptedAdsEvents++;

      }
      //else { // desorption  event: Oxygen --> Vac, for Pt/Ni-O/Vac slab model
      else { // desorption  event: Pt-OH --> Pt, for Pt/Ni/S/Vac NP model (S stands for Pt-OH)
          m_AppliedCE.setSigmaState(sigmaIndex, this.m_PtStates[ohSiteNum]);
          this.m_NumAcceptedDesEvents++;
          //m_NumAcceptedDesEventsSites[oSiteNum]++;
      }
      
  }
  
  
  public void KMC_trigger_equilibra(int sigmaIndex, int ohSiteNum) {
      
      Species spec = m_AppliedCE.getSpecies(sigmaIndex);

      //if (spec==Species.vacancy) { // adsorption event: Vac --> Oxygen, for Pt/Ni-O/Vac slab model
      if (spec==Species.platinum) { // adsorption event: Pt --> Pt-OH, for Pt/Ni/S/Vac NP model (S stands for Pt-OH)
          m_AppliedCE.setSigmaState(sigmaIndex, this.m_OHStates[ohSiteNum]);
          

          /*
          double delta_Ads = m_AppliedCE.setSigmaStateGetDelta(sigmaIndex, this.m_OxygenStates[oSiteNum]);
          
          double oBindingEnergy = delta_Ads + m_O2Energy / 2 + m_H2Energy - m_H2OEnergy;

          this.m_AverageOBindingEnergySites[oSiteNum] += oBindingEnergy;
          this.m_NumAcceptedAdsEventsSites[oSiteNum]++;
          */
          
          //this.m_NumAcceptedAdsEvents++;

      }
      //else { // desorption  event: Oxygen --> Vac, for Pt/Ni-O/Vac slab model
      else { // desorption  event: Pt-OH --> Pt, for Pt/Ni/S/Vac NP model (S stands for Pt-OH)
          m_AppliedCE.setSigmaState(sigmaIndex, this.m_PtStates[ohSiteNum]);
          //this.m_NumAcceptedDesEvents++;
          //m_NumAcceptedDesEventsSites[oSiteNum]++;
      }
      
  }
  
  
  public double getKMCTime() {
      return this.m_KMCTime;
  }
  
  
  
  public double getKMCTime_equil() {
      return this.m_KMCTime_equil;
  }
  
  
  
  public double getKMCTime_record() {
      return this.m_KMCTime_record;
  }
  
  
  public double getPreFactor() {
    
    return this.m_KMCPreFactor;
  }
  
  
  public double getAverageCurrent() {
    
    return m_AverageCurrent;
  }
  
  
  public double getAverageCurrentOccupied() {
    
    return m_AverageCurrentOccupied;
  }
  
  
  
  public double getAverageCurrentUnOccupied() {
    
    return m_AverageCurrentUnOccupied;
  }
  
  
  public double getAverageAdsCurrent() {
    
    return m_AverageAdsCurrent;
  }
  
  public double getAverageFaceOHOccupancy() {
    
    return m_AverageGCNFaceOHOccupancy;
  }
  
  public double getAverageFaceOHBindingEnergy() {
    
    return m_AverageGCNFaceOHBindingEnergy;
  }
  
  public double getAverageFaceOHBindingEnergySq() {
    
    return m_AverageGCNFaceOHBindingEnergySq;
  }
  
  public double getAverageEdgeOHOccupancy() {
      
    return m_AverageEdgeOHOccupancy;
  }
    
  public double getAverageEdgeOHBindingEnergy() {
      
    return m_AverageEdgeOHBindingEnergy;
  }
    
  public double getAverageEdgeOHBindingEnergySq() {
      
    return m_AverageEdgeOHBindingEnergySq;
  }
    
  public double getAverageOHOccupancy() {
        
    return m_AverageOHOccupancy;
  }
      
  public double getAverageOHBindingEnergy() {
        
    return m_AverageOHBindingEnergy;
  }
      
  public double getAverageOHBindingEnergySq() {
        
    return m_AverageOHBindingEnergySq;
  }
  
  public double[] getAverageOHBindingEnergySites() {
    
    return m_AverageOHBindingEnergySites;
  }
  
  /*
  public boolean[] getLastOccupiedOxygenPattern() {
      return this.m_LastOccupiedOxygenPattern;
  }
  */
  

    
    /**
     * 
     * @param fileDir
     * @param structName
     */
    public static void replaceSulfurWithPtOH_2_local(String fileDir, String structName) {
        POSCAR initial_poscar = new POSCAR(fileDir + structName + ".vasp");
        Structure structure = new Structure(initial_poscar);
        AbstractBasis basis = CartesianBasis.getInstance();
        
        int numS = structure.numDefiningSitesWithElement(Element.sulfur);
        int numPt = structure.numDefiningSitesWithElement(Element.platinum);
        int numNi = structure.numDefiningSitesWithElement(Element.nickel);
        int inital_totalAtoms = numS + numPt + numNi;
        
        if(numS==0){
            //System.out.println("  ERROR: no S-atom needs to be replaced by Pt-OH");
        }
        
        //int[] neighbors_num = new int[structure.numDefiningSites()];
        int num_replacedS = 0;
        
        for (int siteNum = 0; siteNum < inital_totalAtoms; siteNum++) {
            
            Structure.Site site = structure.getDefiningSite(siteNum);
            
            //Structure.Site[] neighbors = structure.getNearbySites(site.getCoords(), 2.6, false);
            // neighbors_num[siteNum] = neighbors.length;
            
            
            if(site.getSpecies()==Species.sulfur){ // if the surface site is sulfur, which stands for Pt-OH, then replace the S-atom with Pt-OH
                //superStructure.set
                num_replacedS++;
                
                double avg[]=new double[3];
                int avg_count=0;
                
                Structure.Site[] neighbors = structure.getNearbySites(site.getCoords(), 2.9, false);
                //Structure.Site[] neighbors = structure.getNearbySites(site.getCoords(), 2.6, false);
                for (int neighborNum = 0; neighborNum < neighbors.length; neighborNum++) {
                    Coordinates coord_neighbor = neighbors[neighborNum].getCoords();
                    double[] coord_neighbor_arr = coord_neighbor.getCoordArray(basis);
                    
                    //System.out.println("print out the coords for neighbor "+ siteNum +":" + coord_neighbor_arr[0] + ", " + coord_neighbor_arr[1] + ", " + coord_neighbor_arr[2]) ;
                    
                    avg=add(avg, coord_neighbor_arr) ;
                    avg_count++;
                }
                
                /*
                for (int siteNum1 = 0; siteNum1 < inital_totalAtoms; siteNum1++) {
                    
                    Structure.Site site1 = structure.getDefiningSite(siteNum1);
                    
                    Structure.Site[] neighbors = structure.getNearbySites(site1.getCoords(), 2.6, false);
                    for (int neighborNum = 0; neighborNum < neighbors.length; neighborNum++) {
                        //if (neighbors[neighborNum].getSpecies() == Species.vacancy) {
                        //  break;
                        //  }
                        if(neighbors[neighborNum].getSpecies() == Species.sulfur) {
                            
                            Coordinates coord_Sn = site1.getCoords();
                            double[] coord_Sn_arr = coord_Sn.getCoordArray(basis);
                            //test
                            System.out.println("print out the coords for site "+ siteNum1 +":" + coord_Sn_arr[0] + ", " + coord_Sn_arr[1] + ", " + coord_Sn_arr[2]) ;
                            
                            avg=add(avg, coord_Sn_arr) ;
                            avg_count++;
                            
                            //System.out.println("print out the coords for avg1: " + avg[0] + ", " + avg[1] + ", " + avg[2]) ;
                            
                            break;
                            }
                        } //end of neighbour loop
                    
                }
                */
                
                for( int i=0; i<3; i++) {
                    avg[i]=avg[i]/avg_count;
                }
                
                //test
                
                //System.out.println("print out the coords for avg2: " + avg[0] + ", " + avg[1] + ", " + avg[2]+ ", " +avg_count) ;
                    
                site.setSpecies(Species.platinum); // replace Sulfur with Pt atom
                    
                Coordinates coord_Pt = site.getCoords();
                    
                double[] coord_Pt_arr = coord_Pt.getCoordArray(basis);
                    
                //System.out.println("print out the coords for replaced S: " + coord_Pt_arr[0] + ", " + coord_Pt_arr[1] + ", " + coord_Pt_arr[2]) ;
                    
                //calculate the vector from average(neighbor_2ndlayer[0],neighbor_2ndlayer[1],neighbor_2ndlayer[2]) to replaced S
                double[] translationvector = minus(coord_Pt_arr, avg);
                double det_translationvector= Math.sqrt(translationvector[0]*translationvector[0]+translationvector[1]*translationvector[1]+translationvector[2]*translationvector[2]);
                for (int i=0; i<3; i++) {
                    translationvector[i]=translationvector[i]/det_translationvector;
                }
                    
                //test
                //System.out.println("print out the translation vector " + translationvector[0] + ", " + translationvector[1] + ", " + translationvector[2]) ;
                double O_Pt_translationvector[] = new double[3];
                double H_O_translationvector[] = new double[3];
                  
                for (int i=0; i<3; i++) {
                    O_Pt_translationvector[i]=translationvector[i]*1.98;
                    H_O_translationvector[i]=translationvector[i]*2.96;
                }
                    
                double[] coord_Oxygen = add(coord_Pt_arr , O_Pt_translationvector);
                //test
                //System.out.println("print out the Oxygen coord " + coord_Oxygen[0] + ", " + coord_Oxygen[1] + ", " + coord_Oxygen[2]) ;
                    
                double[] coord_Hydrogen = add(coord_Pt_arr , H_O_translationvector);
                    
                    
                //TODO, find a way to calculate the coordinates for the Oxygen and Hydrogen atoms
                StructureBuilder builder = new StructureBuilder(structure); // builder used to add atoms to the "structure"
                    
                /*  
                Coordinates coord_O = coord_Pt;
                Coordinates coord_H = coord_Pt;
                    
                Coordinates new_coords = new Coordinates(MSMath.arrayAdd(coord_Pt.getArrayCopy(), new double[] {-0.02, 0.02, 0.00}), basis);

                // this is a way to update the coords for O and H. 

                builder.addSite(coord_O, Species.oxygen);
                builder.addSite(coord_H, Species.hydrogen);
                    */
                    
                    
                Coordinates new_coords = new Coordinates(coord_Oxygen, basis);
                Coordinates new_coords_H = new Coordinates(coord_Hydrogen, basis);
                //Coordinates new_coords = new Coordinates(MSMath.arrayAdd(coord_Pt_arr, new double[] {-0.05*28.8, 0.05*28.8, 0.0*28.8}), basis);

                //testA = new_coords.getArrayCopy();
                    
                double[] test_arr2 = new_coords.getCoordArray(basis);
                    
                double[] test_arr2_H = new_coords_H.getCoordArray(basis);
                    

                //System.out.println("print out the coords for Oxygen: " + test_arr2[0] + ", " + test_arr2[1] + ", " + test_arr2[2]) ;
                //System.out.println("print out the coords for Hydrogen: " + test_arr2_H[0] + ", " + test_arr2_H[1] + ", " + test_arr2_H[2]) ;
               
                    
                builder.addSite(new_coords, Species.oxygen);
                    
                builder.addSite(new_coords_H, Species.hydrogen);

                structure = new Structure(builder);
                structure=structure.getCompactStructure().useSitesInCellAsDefining(false);
                    
                    //System.out.println("num of replaced S: " + num_replacedS);
            }
                
         }// end of for loop
            
         POSCAR outfile = new POSCAR(structure);
         //String subStructname = structName.substring(0, structName.length()-2);
         //outfile.writeFile(fileDir + subStructname + "-OH-perpendi.vasp");
         outfile.writeFile(fileDir + structName + "-PtOH-decorates.vasp");

            /*
            Structure struct_shrink = new Structure(outfile);
            struct_shrink = struct_shrink.scaleLatticeConstant(28.8 / 31.2);
            
            PRIM outfile_shrink = new PRIM(struct_shrink, true);
            outfile_shrink.setDescription("struct: latticeConstant-rescale");
            outfile_shrink.writeFile(fileDir + structName + "-PtOH-decorates-shrink.vasp");
            */
            

            /*Structure struct_expand = new Structure(outfile);
            struct_expand = struct_expand.scaleLatticeConstant(31.2 / 28.8);
            
            PRIM outfile_expand = new PRIM(struct_expand, true);
            outfile_expand.setDescription("struct: latticeConstant-rescale");
            outfile_expand.writeFile(fileDir + structName + "-PtOH-decorates-expand.vasp");
            */
    }
    
    
    

    
    //add the coordination arrays
    public static double[] add(double[] first, double[] second) {
        int length=first.length;
        double result[]=new double[length];
        for (int i = 0; i < length; i++) {
            result[i] = first[i] + second[i];
        }
        return result;
        
    }
    
    public static double[] minus(double[] first, double[] second) {
        int length=first.length;
        double result[]=new double[length];
        for (int i = 0; i < length; i++) {
            result[i] = first[i] - second[i];
        }
        return result;
        
    }
    
  
  
}

