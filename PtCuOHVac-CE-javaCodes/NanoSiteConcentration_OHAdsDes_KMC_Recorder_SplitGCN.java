package new_managers;

/**
 * @author Liang
 * 
 * Created on October 27,2019, used to simulate the coverage-dependent ORR currnet on the surface of NP (the adsorbate: *OH)
 */

import java.util.Arrays;

import matsci.Element;
import matsci.Species;
import matsci.engine.monte.GroundStateRecorder;
import matsci.engine.monte.IAllowsSnapshot;
import matsci.engine.monte.metropolis.Metropolis;
import matsci.location.Coordinates;
import matsci.location.symmetry.operations.SymmetryOperation;
import matsci.structure.Structure;
import matsci.structure.decorate.*;
import matsci.structure.decorate.function.*;
import matsci.structure.decorate.function.ce.AbstractAppliedCE;
import matsci.structure.decorate.function.ce.GeneralAppliedCE;
import matsci.structure.decorate.function.monte.SiteConcentrationRecorder;
import matsci.structure.superstructure.SuperStructure;
import matsci.util.arrays.ArrayUtils;



public class NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN extends GroundStateRecorder {
    
  protected AbstractAppliedCE m_AppliedFunction;
  
  protected final int m_CN_Max;
  
  protected final double m_GCN_Cutoff_FaceSites;
  
  protected final Species m_OHSpec;
  
  protected final double m_MuOH_Outer;

  protected final double m_MuOH_Inner_KMC;
  
  protected final int m_NumKMCIterations;
  
  protected final double m_EquilibraRatio_KMC;
    
  protected double[][] m_AverageStateCounts;
  
  protected double m_AverageCurrent;
  protected double m_AverageCurrentOccupied;
  protected double m_AverageCurrentUnOccupied;
  
  protected double m_AverageAdsCurrent;
  
  protected double m_AverageOHOccupancy;
  protected double m_AverageOHBindingEnergy;
  protected double m_AverageOHBindingEnergySq;
  
  protected double m_AverageFaceOHOccupancy;
  protected double m_AverageFaceOHBindingEnergy;
  protected double m_AverageFaceOHBindingEnergySq;
  
  protected double m_AverageEdgeOHOccupancy;
  protected double m_AverageEdgeOHBindingEnergy;
  protected double m_AverageEdgeOHBindingEnergySq;
  
  protected double m_MaxAvgCurrent = 0;
  
  //protected boolean[] m_LastOccupiedOxygenPattern;

  protected int m_NumOHSites;
  
  protected int m_NumOHGCNFaceSites;
  
  protected int m_NumOHEdgeSites;

  protected int[] m_OHSites;
  
  protected int[] m_OHSitesCN;
  
  protected double[] m_OHSitesGCN;

  protected double[] m_AverageOHBindingEnergySites;

  public Structure m_StructureMaxAvgCurrent;
  
  public double m_OutKMCTIme;
  
  //private int[] m_PtNiSigmaIndices;
  //private int[] m_SigmaIndexMap;
  
  /*
  private Metropolis m_MetropolisFreezeSlab;
  //private DecorationMixSlabFreezeManager_new m_ManagerFreezeSlab;
  private DecorationCanonicalSlabFreezeManager m_ManagerFreezeSlab;
  */
  

  
  //protected double m_ShiftAdsEnergy;
  
  
  public NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN(IAllowsSnapshot system, AbstractAppliedCE appliedFunction, double muOH_KMC, double numKMCIterationsPerSite, double equilibraRatio_KMC, Species OH_Spec) {
    //super(system, appliedFunction);
    super(system);
    
    m_AppliedFunction = appliedFunction;
    
    m_CN_Max = 12;
    
    m_GCN_Cutoff_FaceSites = 20.0 /3 ;
    
    m_OHSpec = OH_Spec;
    
    m_MuOH_Outer = m_AppliedFunction.getChemPot(m_OHSpec);
    
    m_MuOH_Inner_KMC = muOH_KMC;
            
    m_EquilibraRatio_KMC = equilibraRatio_KMC;
    
    m_AverageStateCounts = new double[appliedFunction.numSigmaSites()][];
    for (int sigmaIndex = 0; sigmaIndex < m_AverageStateCounts.length; sigmaIndex++) {
      int numStates = appliedFunction.getAllowedSpecies(sigmaIndex).length;
      m_AverageStateCounts[sigmaIndex] = new double[numStates];
    }

    
    
    int numOHSites = 0;
    int numOHGCNFaceSites = 0;
    int numOHEdgeSites = 0;
    
    int[] OHBindingSites = new int[0];
    int[] OHBindingSitesCN = new int[0];
    double[] OHBindingSitesGCN = new double[0];

    Structure structure = this.m_AppliedFunction.getStructure();
    SuperStructure superStruct = this.m_AppliedFunction.getSuperStructure();
    
    int numNiAtoms = 0;
    
    double min_GCN_SurfSite = Double.MAX_VALUE;
    double max_GCN_SurfSite = Double.MIN_VALUE;
    
    double min_CN_SurfSite = Double.MAX_VALUE;
    double max_CN_SurfSite = Double.MIN_VALUE;
    
    for (int siteNum = 0; siteNum < structure.numDefiningSites(); siteNum++) {
        
        Species spec = structure.getSiteSpecies(siteNum);
        
        
        if (spec == Species.vacancy) {
            continue;
        }

        if (spec == Species.nickel) {
            numNiAtoms++;
            //System.out.println("ATTENTION: # of Ni atoms=" + numNiAtoms + ", siteNum=" + siteNum);
            continue;
        }

        /*
        if (spec != Species.platinum) {
            continue;
        }
        */
        
        Structure.Site site = structure.getDefiningSite(siteNum);
        Structure.Site[] neighbors = structure.getNearbySites(site.getCoords(), 3, false);

        //System.out.println("siteNum=" + siteNum + ", neighbors.length=" + neighbors.length); // the neighbors.length == 12
        //neighbors_num[siteNum] = neighbors.length;
        
        int numNonVac = 0; //CN
        int numNonVac_2ndOrder = 0; // GCN * CN_max
        for(int neighborNum = 0; neighborNum < neighbors.length; neighborNum++) {
            Species neighbor_spec = neighbors[neighborNum].getSpecies();
            if(neighbor_spec != Species.vacancy) {
                numNonVac++;
                
                Structure.Site[] neighbors_2ndOrder = structure.getNearbySites(neighbors[neighborNum].getCoords(), 3, false);
                for(int neighborNum_2ndOrder = 0; neighborNum_2ndOrder < neighbors_2ndOrder.length; neighborNum_2ndOrder++) {
                    Species neighbor_spec_2ndOrder = neighbors_2ndOrder[neighborNum_2ndOrder].getSpecies();
                    if(neighbor_spec_2ndOrder != Species.vacancy) {
                        numNonVac_2ndOrder++;
                        //System.out.println("              numNonVac=" + numNonVac + ", numNonVac_2ndOrder=" + numNonVac_2ndOrder);
                    }
                }    
            }
        }
        
        //TODO: the CN of nearest neighboring sites, to calculate GCN
        double local_GCN = ((double)numNonVac_2ndOrder) / m_CN_Max;
        
        //System.out.println("siteNum=" + siteNum + ", CN=" + numNonVac + ", GCN=" + local_GCN);
        
        //if((local_GCN<=7.5)&&(local_GCN>=3)) {
        if((numNonVac<=9)&&(numNonVac>=3)) {
        //if(numNonVac<=9) {
            numOHSites++;
            
            Coordinates coords = structure.getSiteCoords(siteNum);
            //Coordinates coords = superStruct.getSiteCoords(siteNum)
            int sigmaIndex = superStruct.getSiteIndex(coords);

            OHBindingSites = ArrayUtils.appendElement(OHBindingSites, sigmaIndex);
            OHBindingSitesCN = ArrayUtils.appendElement(OHBindingSitesCN, numNonVac);
            OHBindingSitesGCN = ArrayUtils.appendElement(OHBindingSitesGCN,local_GCN);
            
            
            if((min_GCN_SurfSite>local_GCN) && (numNonVac==9)) {
                min_GCN_SurfSite = local_GCN;
                min_CN_SurfSite = numNonVac;
                //System.out.println("local min_GCN_SurfSite: min_GCN_SurfSite=" + min_GCN_SurfSite +  ", CN=" + numNonVac); // the neighbors.length == 12
            }
            if((max_GCN_SurfSite<local_GCN) && (numNonVac==9)) {
                max_GCN_SurfSite = local_GCN;
                max_CN_SurfSite = numNonVac;
                //System.out.println("local max_GCN_SurfSite: max_GCN_SurfSite=" + max_GCN_SurfSite +  ", CN=" + numNonVac); // the neighbors.length == 12
            }
            
            /*
            if(min_GCN_SurfSite>local_GCN) {
                min_GCN_SurfSite = local_GCN;
                min_CN_SurfSite = numNonVac;
                System.out.println("local min_GCN_SurfSite: min_GCN_SurfSite=" + min_GCN_SurfSite +  ", CN=" + numNonVac); // the neighbors.length == 12
            }
            if(max_GCN_SurfSite<local_GCN) {
                max_GCN_SurfSite = local_GCN;
                max_CN_SurfSite = numNonVac;
                System.out.println("local max_GCN_SurfSite: max_GCN_SurfSite=" + max_GCN_SurfSite +  ", CN=" + numNonVac); // the neighbors.length == 12
            }
            */
            
            
            //System.out.println("surface site: siteNum=" + siteNum +  ", sigmaIndex=" + sigmaIndex + ", neighbors.length=" + neighbors.length + ", numNonVac=" + numNonVac); // the neighbors.length == 12

            if(local_GCN >= m_GCN_Cutoff_FaceSites ) {
            //if(Math.abs(local_GCN-7.5)<=0.000001) {
                numOHGCNFaceSites++;
            }
            else {
                numOHEdgeSites++;
            }
            
        }
        //System.out.println(siteNum + ": # of vac_neighbors=" + neighbors.length);
        
    }
    
    
    /*
    int numOHSites = 0;
    int[] OHBindingSites = new int[0];
    for (int sigmaIndex = 0; sigmaIndex < m_AppliedFunction.numSigmaSites(); sigmaIndex++) {
      //if (m_AppliedFunction.allowsSpecies(sigmaIndex, Species.oxygen)) {numOSites++;}
      //TODO  
        Species spec = m_AppliedFunction.getSpecies(sigmaIndex);
        if (spec != Species.vacancy) {
            
            Structure.Site site = this.m_AppliedFunction.getStructure().getDefiningSite(sigmaIndex);
            Structure.Site[] neighbors = this.m_AppliedFunction.getStructure().getNearbySites(site.getCoords(), 3, false);
            
            int numNonVac = 0;
            for(int neighborNum = 0; neighborNum < neighbors.length; neighborNum++) {
                Species neighbor_spec = neighbors[neighborNum].getSpecies();
                if(neighbor_spec != Species.vacancy) {
                    numNonVac++;
                }
            }
            
            //System.out.println("sigmaIndex=" + sigmaIndex + ", spec=" + spec + ", numNonVac=" + numNonVac);

            
            if(numNonVac<=9) {
                numOHSites++;
                OHBindingSites = ArrayUtils.appendElement(OHBindingSites, sigmaIndex);
            }
            
        }
    }
    */
    
    System.out.println("\n\n########################################");
    System.out.println("max_GCN_SurfSite: max_GCN_SurfSite=" + max_GCN_SurfSite +  ", max_CN_SurfSite=" + max_CN_SurfSite); 
    System.out.println("min_GCN_SurfSite: min_GCN_SurfSite=" + min_GCN_SurfSite +  ", min_CN_SurfSite=" + min_CN_SurfSite); 
    System.out.println("########################################\n\n"); 

    
    System.out.println("# of surface sites for *OH, face(GCN>=" +  m_GCN_Cutoff_FaceSites + "), vertex/edge, and total: " + numOHGCNFaceSites + ", " + numOHEdgeSites + ", " + numOHSites);
    
    m_NumOHSites = numOHSites;
    m_NumOHGCNFaceSites = numOHGCNFaceSites;
    m_NumOHEdgeSites = numOHEdgeSites;
    
    if((m_NumOHGCNFaceSites+ m_NumOHEdgeSites)!=m_NumOHSites) {
        System.out.println("ATTENTION. VERY BAD news: (m_NumOHFaceSites+ m_NumOHEdgeSites)!=m_NumOHSites");
    }
    
    m_NumKMCIterations = (int)(numKMCIterationsPerSite * m_NumOHSites);

    
    m_OHSites = OHBindingSites;
    m_OHSitesCN = OHBindingSitesCN;
    m_OHSitesGCN = OHBindingSitesGCN;
    
    m_AverageOHBindingEnergySites = new double[m_NumOHSites];

    
    //m_LastOccupiedOxygenPattern = new boolean[numOSites];
    
    this.m_StructureMaxAvgCurrent = appliedFunction.getStructure();
    
    System.out.println("\n\n**********************************************************************");
    System.out.println("Initialize 1st-level recorder of KMC-current: NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN");
    System.out.println("numKMCIterationsPerSite=" + numKMCIterationsPerSite + ", m_NumKMCIterations=" + m_NumKMCIterations + ", m_EquilibraRatio_KMC=" + m_EquilibraRatio_KMC);
    System.out.println("m_GCN_Cutoff_FaceSites=" + m_GCN_Cutoff_FaceSites);
    System.out.println("**********************************************************************\n\n");


  }
  

  public void recordEndState(double value, double weight) {
    double oldTotalWeight = m_TotalWeight;
    double newTotalWeight = oldTotalWeight + weight;
    super.recordEndState(value, weight);
    double correctionRatio = oldTotalWeight / newTotalWeight;
    double incrementRatio = (1 - correctionRatio);
    for (int sigmaIndex = 0; sigmaIndex < m_AverageStateCounts.length; sigmaIndex++) {
      double[] stateCounts = m_AverageStateCounts[sigmaIndex];
      for (int state = 0; state < stateCounts.length; state++) {
        stateCounts[state] *= correctionRatio;
      }
    }
    
    //TODO
    //NanoGCCoverageDepend_KMCCurrent_Recorder currentRecorder = this.get_OHAdsDes_KMC_Current();
    NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN currentRecorder = this.get_OHAdsDes_KMC_Current();
    //NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN_KMCTime currentRecorder = this.get_OHAdsDes_KMC_Current();

    
    double tempAverageCurrent = currentRecorder.getAverageCurrent();
    
    m_AverageCurrent *= correctionRatio;
    m_AverageCurrent += incrementRatio * tempAverageCurrent;
    
    if(this.m_MaxAvgCurrent < tempAverageCurrent){
        this.m_MaxAvgCurrent = tempAverageCurrent;
        this.m_StructureMaxAvgCurrent = currentRecorder.m_AppliedCE.getStructure();
    }
    
    //m_LastOccupiedOxygenPattern = currentRecorder.getLastOccupiedOxygenPattern();
    
    m_AverageAdsCurrent *= correctionRatio;
    m_AverageAdsCurrent += incrementRatio * currentRecorder.getAverageAdsCurrent();
    
    /*
    m_AverageCurrentOccupied *= correctionRatio;
    m_AverageCurrentOccupied += incrementRatio * currentRecorder.getAverageCurrentOccupied();

    m_AverageCurrentUnOccupied *= correctionRatio;
    m_AverageCurrentUnOccupied += incrementRatio * currentRecorder.getAverageCurrentUnOccupied();
    */
    
    //TODO, clean the surface
    currentRecorder.clearSurfaceAdsorbates();
    
    m_AverageOHOccupancy *= correctionRatio;
    m_AverageOHOccupancy += incrementRatio * currentRecorder.getAverageOHOccupancy();
    
    m_AverageOHBindingEnergy *= correctionRatio;
    m_AverageOHBindingEnergy += incrementRatio * currentRecorder.getAverageOHBindingEnergy();
    
    m_AverageOHBindingEnergySq *= correctionRatio;
    m_AverageOHBindingEnergySq += incrementRatio * currentRecorder.getAverageOHBindingEnergySq();
    
    
    //for Face sites with CN=9
    m_AverageFaceOHOccupancy *= correctionRatio;
    m_AverageFaceOHOccupancy += incrementRatio * currentRecorder.getAverageFaceOHOccupancy();
    
    m_AverageFaceOHBindingEnergy *= correctionRatio;
    m_AverageFaceOHBindingEnergy += incrementRatio * currentRecorder.getAverageFaceOHBindingEnergy();
    
    m_AverageFaceOHBindingEnergySq *= correctionRatio;
    m_AverageFaceOHBindingEnergySq += incrementRatio * currentRecorder.getAverageFaceOHBindingEnergySq();
    
    
    //for edge (including vertex) sites with CN<9 && CN>=3
    m_AverageEdgeOHOccupancy *= correctionRatio;
    m_AverageEdgeOHOccupancy += incrementRatio * currentRecorder.getAverageEdgeOHOccupancy();
    
    m_AverageEdgeOHBindingEnergy *= correctionRatio;
    m_AverageEdgeOHBindingEnergy += incrementRatio * currentRecorder.getAverageEdgeOHBindingEnergy();
    
    m_AverageEdgeOHBindingEnergySq *= correctionRatio;
    m_AverageEdgeOHBindingEnergySq += incrementRatio * currentRecorder.getAverageEdgeOHBindingEnergySq();
    
    m_OutKMCTIme *= correctionRatio;
    m_OutKMCTIme += incrementRatio * currentRecorder.getKMCTime();
    
    /*
    double[] temp_avgOBindingEnergySites = currentRecorder.getAverageOBindingEnergySites();
    for(int oSiteNum=0; oSiteNum < this.m_AverageOBindingEnergySites.length/2; oSiteNum++) {
       m_AverageOBindingEnergySites[oSiteNum] *= correctionRatio; 
       m_AverageOBindingEnergySites[oSiteNum] += incrementRatio * temp_avgOBindingEnergySites[oSiteNum];
    }
    */
    
    m_AppliedFunction.setChemPot(this.m_OHSpec, m_MuOH_Outer);
    
    for (int sigmaIndex = 0; sigmaIndex < m_AppliedFunction.numSigmaSites(); sigmaIndex++) {
      int state = m_AppliedFunction.getQuickSigmaState(sigmaIndex);
      m_AverageStateCounts[sigmaIndex][state] += incrementRatio;
    }

  }
  
  public double getAverageStateCount(int sigmaIndex, int state) {
    return m_AverageStateCounts[sigmaIndex][state];
  }
  
  public double getOutKMCTime() {
      return m_OutKMCTIme;    
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
  
  public double getAverageOHOccupancy() {
    return m_AverageOHOccupancy;
  }
  
  public double getAverageFaceOHOccupancy() {
      return m_AverageFaceOHOccupancy;
  }
  
  public double getAverageEdgeOHOccupancy() {
      return m_AverageEdgeOHOccupancy;
  }

  
  public double getAverageOHBindingEnergy() {
    return m_AverageOHBindingEnergy;
  }
  
  public double getAverageFaceOHBindingEnergy() {
    return m_AverageFaceOHBindingEnergy;
  }

  
  public double getAverageEdgeOHBindingEnergy() {
    return m_AverageEdgeOHBindingEnergy;
  }

  
  public double getAverageOHBindingEnergySq() {
    return m_AverageOHBindingEnergySq;
  }
  
  public double getAverageFaceOHBindingEnergySq() {
    return m_AverageFaceOHBindingEnergySq;
  }
  
  public double getAverageEdgeOHBindingEnergySq() {
    return m_AverageEdgeOHBindingEnergySq;
  }
  
  
  public double[] getAverageOHBindingEnergySites() {
    return m_AverageOHBindingEnergySites;
  }
  
  
  public double getMaxAvgCurrent() {
      return this.m_MaxAvgCurrent;
    }
  
  public Structure getStructureMaxTempAvgCurrent() {
      return this.m_StructureMaxAvgCurrent;
  }
  
  public int getNumOHSites() {
    return this.m_NumOHSites;
  }
  
  public int getNumOHFaceSites() {
      return this.m_NumOHGCNFaceSites;
    }
  
  public int getNumOHEdgeSites() {
      return this.m_NumOHEdgeSites;
  }

  /*
  public NanoGCCoverageDepend_KMCCurrent_Recorder_SplitCN get_OHAdsDes_KMC_Current(){
      
    //return new NanoGCCoverageDepend_KMCCurrent_Recorder(m_AppliedFunction, this.m_MuOH_Inner_KMC, this.m_NumKMCIterations, this.m_EquilibraRatio_KMC, this.m_LastOccupiedOxygenPattern);
    //return new NanoGCCoverageDepend_KMCCurrent_Recorder(m_AppliedFunction, this.m_MuOH_Inner_KMC, this.m_NumKMCIterations, this.m_EquilibraRatio_KMC, m_OHSites, m_OHSpec, this.m_NumOHFaceSites, this.m_NumOHEdgeSites);
    
    return new NanoGCCoverageDepend_KMCCurrent_Recorder(m_AppliedFunction, this.m_MuOH_Inner_KMC, this.m_NumKMCIterations, this.m_EquilibraRatio_KMC, m_OHSites, m_OHSpec, this.m_NumOHEdgeSites);
    
  }
  */
  
  

  public NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN get_OHAdsDes_KMC_Current(){

      return new NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN(m_AppliedFunction, this.m_MuOH_Inner_KMC, this.m_NumKMCIterations, this.m_EquilibraRatio_KMC, m_OHSites, m_OHSitesCN, m_OHSitesGCN, m_OHSpec, this.m_NumOHEdgeSites, m_GCN_Cutoff_FaceSites);

    }

  
  /*
  public NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN_KMCTime get_OHAdsDes_KMC_Current(){

      return new NanoGCCoverageDepend_KMCCurrent_Recorder_SplitGCN_KMCTime(m_AppliedFunction, this.m_MuOH_Inner_KMC, this.m_NumKMCIterations, this.m_EquilibraRatio_KMC, m_OHSites, m_OHSitesCN, m_OHSitesGCN, m_OHSpec, this.m_NumOHEdgeSites, m_GCN_Cutoff_FaceSites);

    }
  */
  
}





