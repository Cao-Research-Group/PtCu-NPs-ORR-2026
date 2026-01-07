/**
 * 
 */
package PtCuOH;


/**
 * 
 */


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Arrays;
import java.util.StringTokenizer;

//import binary.NanoFitClusterExpansion;
import matsci.Element;
import matsci.Species;
import matsci.engine.monte.metropolis.Metropolis;
import matsci.io.app.log.Status;
import matsci.io.clusterexpansion.PRIM;
import matsci.io.clusterexpansion.StructureListFile;
import matsci.io.vasp.POSCAR;
import matsci.location.Coordinates;
import matsci.location.Vector;
import matsci.location.basis.CartesianBasis;
import matsci.location.basis.LinearBasis;
import matsci.model.linear.L2FullRLSFitter;
import matsci.model.reg.RegularizerSelector;
import matsci.structure.BravaisLattice;
import matsci.structure.Structure;
import matsci.structure.StructureBuilder;
import matsci.structure.decorate.function.ce.AbstractAppliedCE;
import matsci.structure.decorate.function.ce.ClusterExpansion;
import matsci.structure.decorate.function.ce.FastAppliedCE;
import matsci.structure.decorate.function.ce.clusters.ClusterGroup;
import matsci.structure.decorate.function.ce.reg.ILinearCERegGenerator;
import matsci.structure.decorate.function.ce.reg.PowExpRegGenerator;
import matsci.structure.decorate.function.ce.structures.StructureList;
import matsci.structure.function.ce.structures.ClustStructureData;
import matsci.structure.superstructure.SuperStructure;
import matsci.structure.symmetry.StructureMapper;
import matsci.util.MSMath;
import matsci.util.arrays.ArrayUtils;


/**
 * @author Liang
 *
 */
public class NanoFitClusterExpansion_PtCuOH {

    
    public static String ROOT = "ce/PtCuOH/";
    public static String STRUCT_DIR = ROOT + "/structures/";
    public static String CONTCAR_DIR = ROOT + "/contcars/";
    public static String CLUST_DIR = ROOT + "/clusters/";
    public static String CLUST_DIR_PTCUOH = ROOT + "/clusters_PtCuOH/";
    public static String GROUND_STATE_DIR = ROOT + "/groundStates/";
    
    public static String DFT_CE_ERROR_DIR = ROOT + "/DFT_CE_error/";
    
    //public static String  = ROOT + "/groundStates/";


    /*
     * RPBE-PREC=Accurate, getKPoints
     */
    public static double PtEnergyBulkAccurate = -5.4101415;
    public static double CuEnergyBulkAccurate = -3.2773616;
    public static double OHEnergyCorrectAccurate = 10.66793045;
    
    
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        
        Status.setLogLevelDetail();
        //Status.setLogLevelBasic();
        
        
        //fit CE for PtCuSVac (particles with *OH)
        ClusterExpansion ce = buildCE_PtCuOH();
        //addKnownStructures_PtCuOH_CalculateOHCoverage(ce, 1.0);      
        fitCE(ce,"energyList-trainingSet-PtCuSVac");
        
        
        
        //fit CE for PtCuVac (clean particles without *OH)
//        ClusterExpansion ce = buildCE_PtCuVac();
//        addKnownStructures_PtCuVac_cleanNP(ce);      
//        fitCE(ce);
//        
        
        //PREC=Accurate
        //checkNanoparticleRelaxation_driver(0.75, 1.0, ROOT + "/Pt-Ni-nanoparticles-CONTCAR-finish/gs-finished-PREC=Normal-CONTCAR-files/", ROOT + "/Pt-Ni-nanoparticles-CONTCAR-finish/gs-finished-PREC=Accurate-CONTCAR-files/", ROOT + "/Pt-Ni-nanoparticles-CONTCAR-finish/", "structList-gs-finished-PREC=Normal.txt");

        
    }
   
    /**
     * 
     * buildCE method
     */
    public static ClusterExpansion buildCE_PtCuOH() {
      
      PRIM prim = new PRIM(ROOT + "/PRIM-CuPtSVac.vasp", true);
      ClusterExpansion ce = new ClusterExpansion(prim);
      
      /*
       * (10, 6, 5, 4, 4, 0) are the parameters for Pt3Ni-Mo Nano Lett paper
       */
      double pairCutoff = 7; 
      double tripleCutoff = 6;
      double quadCutoff = 5;
      double quintCutoff = 0;// include
      double hexCutoff = 0; //include
      double seCuCutoff = 0;
 
      /*
       * 2nd set: (6,6,5,4,4,0)
       */
      /*
      pairCutoff = 6;
      tripleCutoff = 6;
      quadCutoff = 5;
      quintCutoff = 4;
      hexCutoff = 4;
      seCuCutoff = 0;
      */
      
      /*
       * 3rd set: (6,6,4,4,4,0)
       */
      /*
      pairCutoff = 6;
      tripleCutoff = 6;
      quadCutoff = 4;
      quintCutoff = 4;
      hexCutoff = 4;
      seCuCutoff = 0;
      */
      
      
      /*
       * 4th set: (5,4,0,0,0,0)
       */
      /*
      pairCutoff = 5;
      tripleCutoff = 4;
      quadCutoff = 0;
      quintCutoff = 0;
      hexCutoff = 0;
      seCuCutoff = 0;
      */
      
      System.out.println("Pair Cutoff: " + pairCutoff);
      System.out.println("Triple Cutoff: " + tripleCutoff);
      System.out.println("Quad Cutoff: " + quadCutoff);
      System.out.println("Quint Cutoff: " + quintCutoff);
      System.out.println("Hex Cutoff: " + hexCutoff);
      System.out.println("SeCu Cutoff: " + seCuCutoff);
      ce.appendClusters(1, 1);
      ce.appendClusters(2, pairCutoff);
      ce.appendClusters(3, tripleCutoff);
      ce.appendClusters(4, quadCutoff);
      ce.appendClusters(5, quintCutoff);
      ce.appendClusters(6, hexCutoff);
      ce.appendClusters(7, seCuCutoff);
      
      Species gold = Species.get("Au");  // first element in PRIM file
      Species palladium = Species.get("Pd"); // second element in PRIM file
      Species[] speciesInOrder = new Species[] {gold, palladium};
      int numFunctionGroups = 0;
      for (int clustNum = 0; clustNum < ce.numGroups(); clustNum++) {
        ClusterGroup group = ce.getClusterGroup(clustNum);
        POSCAR outfile = new POSCAR(group.getSampleStructure(gold, palladium), speciesInOrder);
        outfile.writeFile(CLUST_DIR_PTCUOH + clustNum + ".vasp");
        //outfile.writeVICSOutFile(CLUST_DIR_PTNIO + clustNum + ".out");
        numFunctionGroups += group.numFunctionGroups();
      }
      System.out.println(ce.numGroups() + " cluster groups generated");
      System.out.println(numFunctionGroups + " function groups generated");
      return ce;

    }

    
    

      
      
      
    
    /**
     * 
     * buildCE method
     */
    public static ClusterExpansion buildCE_PtCuVac() {
      
      PRIM prim = new PRIM(ROOT + "/PRIM-PtCu.vasp", true);
      ClusterExpansion ce = new ClusterExpansion(prim);
      
      /*
       * (10, 6, 5, 4, 4, 0) are the parameters for Pt3Ni-Mo Nano Lett paper
       */
      double pairCutoff = 7; 
      double tripleCutoff = 6;
      double quadCutoff = 5;
      double quintCutoff = 4;// include
      double hexCutoff = 4; //include
      double seCuCutoff = 0;
 
      /*
       * 2nd set: (6,6,5,4,4,0)
       */
      /*
      pairCutoff = 6;
      tripleCutoff = 6;
      quadCutoff = 5;
      quintCutoff = 4;
      hexCutoff = 4;
      seCuCutoff = 0;
      */
      
      /*
       * 3rd set: (6,6,4,4,4,0)
       */
      /*
      pairCutoff = 6;
      tripleCutoff = 6;
      quadCutoff = 4;
      quintCutoff = 4;
      hexCutoff = 4;
      seCuCutoff = 0;
      */
      
      
      /*
       * 4th set: (5,4,0,0,0,0)
       */
      /*
      pairCutoff = 5;
      tripleCutoff = 4;
      quadCutoff = 0;
      quintCutoff = 0;
      hexCutoff = 0;
      seCuCutoff = 0;
      */
      
      System.out.println("Pair Cutoff: " + pairCutoff);
      System.out.println("Triple Cutoff: " + tripleCutoff);
      System.out.println("Quad Cutoff: " + quadCutoff);
      System.out.println("Quint Cutoff: " + quintCutoff);
      System.out.println("Hex Cutoff: " + hexCutoff);
      System.out.println("SeCu Cutoff: " + seCuCutoff);
      ce.appendClusters(1, 1);
      ce.appendClusters(2, pairCutoff);
      ce.appendClusters(3, tripleCutoff);
      ce.appendClusters(4, quadCutoff);
      ce.appendClusters(5, quintCutoff);
      ce.appendClusters(6, hexCutoff);
      ce.appendClusters(7, seCuCutoff);
      
      Species gold = Species.get("Au");  // first element in PRIM file
      Species palladium = Species.get("Pd"); // second element in PRIM file
      Species[] speciesInOrder = new Species[] {gold, palladium};
      int numFunctionGroups = 0;
      for (int clustNum = 0; clustNum < ce.numGroups(); clustNum++) {
        ClusterGroup group = ce.getClusterGroup(clustNum);
        POSCAR outfile = new POSCAR(group.getSampleStructure(gold, palladium), speciesInOrder);
        outfile.writeFile(CLUST_DIR + clustNum + ".vasp");
        //outfile.writeVICSOutFile(CLUST_DIR + clustNum + ".out");
        numFunctionGroups += group.numFunctionGroups();
      }
      System.out.println(ce.numGroups() + " cluster groups generated");
      System.out.println(numFunctionGroups + " function groups generated");
      return ce;

    }

    
    /**
     * 
     * buildCE method
     */
    public static ClusterExpansion buildCE_dummy() {
      
      PRIM prim = new PRIM(ROOT + "/PRIM-PtCu.vasp", true);
      ClusterExpansion ce = new ClusterExpansion(prim);
      
      double pairCutoff = 4; 
      double tripleCutoff = 4;
      
      ce.appendClusters(1, 1);
      ce.appendClusters(2, pairCutoff);
      //ce.appendClusters(3, tripleCutoff);
      
      return ce;

    }
    
    
    /**
     * addKnownStructures() method
     * the txt file contains free energies, while, the formation energies are calculated and added to CE
     * @param ce
     * @return
     */
    public static ClusterExpansion addKnownStructures_PtCuOH_CalculateOHCoverage(ClusterExpansion ce, double cutoffOHML,String StructList) {
      
      //StructureListFile enerIn = new StructureListFile(ROOT + "energyList-trainingSet-PtNiVac.txt"); //chosen
      StructureListFile enerIn = new StructureListFile(ROOT + StructList+".txt"); //chosen

  
      /*
       * different x,y,z supercels, use different appliedCE object.   
       */
      
      /*
       * may have problem when bulk struct added to training set. 
       */
      int[][] superToDirect = new int[][] {
          {-8, 8, 8},
          {8, -8, 8},
          {8, 8, -8}
      };

      FastAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);
      
      int num_struct_withOH = 0;
      int num_struct_withSingleOH = 0;
      
      int num_struct_nonSurfaceOH =0;
      
      double max_OHCoverage = Double.NEGATIVE_INFINITY;
      double min_OHCoverage = Double.POSITIVE_INFINITY;
      
      int numStructLess85 = 0;
      int numStructLess100 = 0;
      
      int OHML1stBin = 0;
      int OHML2ndBin = 0;
      int OHML3rdBin = 0;
      int OHML4thBin = 0;
      int OHML5thBin = 0;
      int OHML6thBin = 0;
      int OHML7thBin = 0;
      int OHML8thBin = 0;
      int OHML9thBin = 0;
      int OHML10thBin = 0;
      int OHMLAllBin = 0;

      
      try{
          
          FileWriter fw = new FileWriter(ROOT + "/_added-energyList-PtCu-OH.txt", false);
          BufferedWriter bw = new BufferedWriter(fw); 
          
          bw.write("structIndex structName numPt numCu numOH numsurfacePt energy FE OHcoverage OH_CN=9 OH_CN=8 OH_CN=7 OH_CN=6 OH_CN=5 OH_CN<5 GCN>=6.67" + "\r\n");
          bw.flush();
      
          //TODO
          for (int entryNum = 0; entryNum < enerIn.numEntries(); entryNum++) {
              String fileName = enerIn.getFileName(entryNum);
              //POSCAR infile = new POSCAR(STRUCT_DIR + fileName, true);
              POSCAR infile = new POSCAR(ROOT + "/structures/" + fileName, true);

              //System.out.println(fileName);

              Structure structure = new Structure(infile);

              int numPt = structure.numDefiningSitesWithElement(Element.platinum);
              int numCu = structure.numDefiningSitesWithElement(Element.copper);
              int numNi = structure.numDefiningSitesWithElement(Element.nickel);
              int numP = structure.numDefiningSitesWithElement(Element.phosphorus);
              int numS = structure.numDefiningSitesWithElement(Element.sulfur);

              //int totalAtoms = numPt + numNi + numMo + numCu;
              int totalAtoms = numPt + numCu + numP + numS;
              
              if(numP!=0){
                  //System.out.println("the structure contains Ni-*OH: " + fileName);
                  continue;
              }
              
              
              if(numS!=0){
                  //System.out.println("the structure contains Ni-*OH: " + fileName);
                  //continue; //exlude NPs with *OH
              }
              
              //if (((totalAtoms) > 50) && ((totalAtoms) < 85)) {continue;}
              if (((totalAtoms) > 50) && ((totalAtoms) < 100)) {continue;}
              if (numNi>0) { continue;}
              //if (((totalAtoms) > 1) && ((totalAtoms) < 100)) {continue;}
              
              
              //TODO
              /*
               * calculate the *OH coverage for each training structure
               */
              int numSurfaceS = 0; //only CN<=9, considered as surface site
              int numSurfaceSites = 0;
              int numSurfacePt = 0;
              int numSurfaceS_CN9=0;
              int numSurfaceS_CN8=0;
              int numSurfaceS_CN7=0;
              int numSurfaceS_CN6=0;
              int numSurfaceS_CN5=0;
              int numSurfaceS_CNlessthan5=0;
              int numGCNmorethan6=0;
              for (int siteNum = 0; siteNum < structure.numDefiningSites(); siteNum++) {
                  Structure.Site site = structure.getDefiningSite(siteNum);
                  Structure.Site[] neighbors = structure.getNearbySites(site.getCoords(), 3, false);
                              
                  //if((neighbors.length < 12) && (site.getSpecies()==Species.nickel)){
                  if(neighbors.length <= 9){
                      numSurfaceSites++;
                  }
                  if((neighbors.length <= 9) && (site.getSpecies()==Species.platinum)){
                	  numSurfacePt++;
                  }
                  if((neighbors.length <= 9) && (site.getSpecies()==Species.sulfur)){
                	  numSurfacePt++;
                  }
                  if((neighbors.length <= 9) && (site.getSpecies()==Species.sulfur)){
                      numSurfaceS++;
                  }
                  if((neighbors.length == 9) && (site.getSpecies()==Species.sulfur)){
                	  numSurfaceS_CN9++;
                  }
                  if((neighbors.length == 8) && (site.getSpecies()==Species.sulfur)){
                	  numSurfaceS_CN8++;
                  }
                  if((neighbors.length == 7) && (site.getSpecies()==Species.sulfur)){
                	  numSurfaceS_CN7++;
                  }
                  if((neighbors.length == 6) && (site.getSpecies()==Species.sulfur)){
                	  numSurfaceS_CN6++;
                  }
                  if((neighbors.length == 5) && (site.getSpecies()==Species.sulfur)){
                	  numSurfaceS_CN5++;
                  }
                  if((neighbors.length < 5) && (site.getSpecies()==Species.sulfur)){
                	  numSurfaceS_CNlessthan5++;
                  }
                  if(site.getSpecies()==Species.sulfur) {
                	  double GCN=0;
                	  for(int neighborSiteNum=0;neighborSiteNum<neighbors.length;neighborSiteNum++) {
                		  Structure.Site site2=neighbors[neighborSiteNum];
                		  Structure.Site[] neighbors2 = structure.getNearbySites(site2.getCoords(), 3, false);
    
                		  GCN+=neighbors2.length;
                	  }
                	  
                	  if(GCN>=80) {
                		  numGCNmorethan6++;  
                	  }
                  }
              }
              
              double SCoverage_temp = ((double)numSurfaceS) / numSurfaceSites;

              if(SCoverage_temp > cutoffOHML) {
                  continue;
              }
              
              if(max_OHCoverage < SCoverage_temp) {
                  max_OHCoverage = SCoverage_temp; 
              }
              if(min_OHCoverage > SCoverage_temp) {
                  min_OHCoverage = SCoverage_temp; 
              }
              
              
              if((SCoverage_temp>0.0)&&(SCoverage_temp<=0.1)) {
                  OHML1stBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.1)&&(SCoverage_temp<=0.2)) {
                  OHML2ndBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.2)&&(SCoverage_temp<=0.3)) {
                  OHML3rdBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.3)&&(SCoverage_temp<=0.4)) {
                  OHML4thBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.4)&&(SCoverage_temp<=0.5)) {
                  OHML5thBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.5)&&(SCoverage_temp<=0.6)) {
                  OHML6thBin++;
                  OHMLAllBin++;
                  
              }
              if((SCoverage_temp>0.6)&&(SCoverage_temp<=0.7)) {
                  OHML7thBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.7)&&(max_OHCoverage<=0.8)) {
                  OHML8thBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.8)&&(SCoverage_temp<=0.9)) {
                  OHML9thBin++;
                  OHMLAllBin++;
              }
              if((SCoverage_temp>0.9)&&(SCoverage_temp<=1.0)) {
                  OHML10thBin++;
                  OHMLAllBin++;
              }
              
              
              
              if(numSurfaceS != numS){
                  num_struct_nonSurfaceOH++;
                  System.out.println("struct=" + fileName + ", has surface Pt-OH sites with CN>9");
                  continue;
              }
              
              if(numS != 0){
                  num_struct_withOH++;
              }
              
              //exclude NPs with only one *OH
//              if(numS == 1){
//                  num_struct_withSingleOH++;
//                  continue;
//              }

              double delta = Math.abs(structure.getDefiningVolume() - appliedCE.getSuperStructure().getDefiningVolume());
              
              //the supercell of this structure is the same as the implicit one in the appliedCE
              if (delta < 0.1) {
                  //appliedNanoparticleCE.decorateFromStructure(structure);
                  //appliedCE = appliedNanoparticleCE;
                  
                  appliedCE.decorateFromStructure(structure);

              } 
              //the supercell of this structure is different from the implicit one in the appliedCE
              else {
                  System.out.println("!!! Attention: " + fileName + " has different shape of super cell.");
                  appliedCE = new FastAppliedCE(ce, ce.getBaseStructure().getSuperStructure(infile));
              }
              
              //SuperStructure superStructure = appliedCE.getSuperStructure();
              
              
              //appliedCE.decorateFromStructure(new Structure(infile));
              appliedCE.getSuperStructure().setDescription(fileName);
              
              
              double value = enerIn.getEnergy(entryNum); 
              
              if(value==Double.NaN){
                  continue;
              }

              
              // formation energy: formation E = E(PtxNiy)-xE(Pt)-yE(Ni)
              value -= numPt * PtEnergyBulkAccurate; // From bulk calCulations
              value -= numCu * CuEnergyBulkAccurate; // From bulk calCulations
              
              //value -= numCu * CuEnergyBulkAccurate; // From bulk calCulations
              //value -= numMo * MoEnergyBulkAccurate; // From bulk calCulations

              //TODO
              // Sulfur stands for Pt-OH: "value += numS * OHEnergyCorrect;" used to correct the *OH from the group of "Pt-OH", which is referenced to 0.5*E(H2)-E(H2O)
              //                          "value -= numS * PtEnergyBulk;"    used to account for the Pt from the group of "Pt-OH", , which is referenced to E(Pt-bulk)
              value -= numS * PtEnergyBulkAccurate; // From bulk calCulations
              value -= numP * CuEnergyBulkAccurate; // From bulk calCulations
              
              //TODO
              //correction to *OH term, referenced to 0.5*E(H2)-E(H2O)
              value += numP * OHEnergyCorrectAccurate;
              value += numS * OHEnergyCorrectAccurate;

              
              value /= appliedCE.numPrimCells();
              
              //Status.detail(fileName + "      " + enerIn.getEnergy(entryNum));

              Status.detail("Adding structure " + (entryNum + 1) + "/" + enerIn.numEntries() + ", totalAtoms=" + totalAtoms + ", " + fileName + ". numPt: " +  numPt +", numCu: " + numCu + ", numS: " + numS + ", value / site: " + value + ", appliedCE.numPrimCells(): " + appliedCE.numPrimCells()+ ", OHCoverage: " + SCoverage_temp);
              //Status.detail("Adding structure " + (entryNum + 1) + "/" + enerIn.numEntries() + ", " + fileName + ". numPt: " +  numPt +", numNi: " + numNi + ", numMo: " + numMo + ", value: " + value);
              
              if(totalAtoms < 85) {
                  System.out.println("added structure with totalAtoms<85: " + fileName + ", totalAtoms=" + totalAtoms);
                  numStructLess85++;
              }
              
              if(totalAtoms < 100) {
                  //System.out.println("added structure with totalAtoms<85: " + fileName + ", totalAtoms=" + totalAtoms);
                  numStructLess100++;
              }
              
              bw.write((entryNum + 1) + "   " + fileName + "   " + numPt + "   " + numCu + "   " + numS + "   " +numSurfacePt+ "   " + enerIn.getEnergy(entryNum) + "   " + value + "   " + SCoverage_temp + " "+numSurfaceS_CN9+" "+numSurfaceS_CN8+" "+numSurfaceS_CN7+" "+numSurfaceS_CN6+" "+numSurfaceS_CN5+" "+numSurfaceS_CNlessthan5+" "+numGCNmorethan6+"\r\n");
              bw.flush();
              
              appliedCE.activateAllGroups();
              ce.getStructureList().addCorrelations(appliedCE, value);
          }
          
          bw.flush();
          bw.close();
          fw.close();

      }
      
      catch (IOException e) {
      e.printStackTrace();
      }

      
      Status.basic("Added " + ce.getStructureList().numKnownStructures() + " structures.\n\n");

      Status.basic("The # structures with non-surface Pt-OH (CN>9) on the surface is: " + num_struct_nonSurfaceOH + "\n\n");

      Status.basic("The # structures with Pt-OH and Pt-SingleOH on the surface is: " + num_struct_withOH + ", " + num_struct_withSingleOH + "\n\n");

      Status.basic("The, cutoff, max and min OHCoverage for added training set: cutoff=" + cutoffOHML + ", max=" + max_OHCoverage + ", min=" + min_OHCoverage + "\n\n");

      System.out.println("# of structs with atoms < 85: " +  numStructLess85);
      System.out.println("# of structs with atoms < 100: " +  numStructLess100);

      System.out.println("\n\n******************************************************");
      System.out.println("the analysis of *OH coverages:");

      System.out.println("# of structs with *OH coverage < 0.1 ML: " +  OHML1stBin);
      System.out.println("# of structs with *OH coverage < 0.2 ML: " +  OHML2ndBin);
      System.out.println("# of structs with *OH coverage < 0.3 ML: " +  OHML3rdBin);
      System.out.println("# of structs with *OH coverage < 0.4 ML: " +  OHML4thBin);
      System.out.println("# of structs with *OH coverage < 0.5 ML: " +  OHML5thBin);
      System.out.println("# of structs with *OH coverage < 0.6 ML: " +  OHML6thBin);
      System.out.println("# of structs with *OH coverage < 0.7 ML: " +  OHML7thBin);
      System.out.println("# of structs with *OH coverage < 0.8 ML: " +  OHML8thBin);
      System.out.println("# of structs with *OH coverage < 0.9 ML: " +  OHML9thBin);
      System.out.println("# of structs with *OH coverage < 1.0 ML: " +  OHML10thBin);
      System.out.println("\n# of structs with *OH coverage > 0 ML: " +  OHMLAllBin);
      System.out.println("******************************************************\n\n");
      
      return ce;
    }
    


    
    
    /**
     * addKnownStructures() method
     * the txt file contains free energies, while, the formation energies are calculated and added to CE
     * @param ce
     * @return
     */
    public static ClusterExpansion addKnownStructures_PtCuVac_cleanNP(ClusterExpansion ce) {
      
      StructureListFile enerIn = new StructureListFile(ROOT + "energyList-trainingSet-PtCuSVac-20240727.txt"); //chosen
      //StructureListFile enerIn = new StructureListFile(ROOT + "energyList-trainingSet-PtNiSVac.txt"); //chosen

  
      /*
       * different x,y,z supercels, use different appliedCE object.   
       */
      
      /*
       * may have problem when bulk struct added to training set. 
       */
      int[][] superToDirect = new int[][] {
          {-8, 8, 8},
          {8, -8, 8},
          {8, 8, -8}
      };

      FastAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);

      
      try{
          
          FileWriter fw = new FileWriter(ROOT + "/_added-energyList-PtCu-cleanNP.txt", false);
          BufferedWriter bw = new BufferedWriter(fw); 
          
          bw.write("structIndex structName numPt numNi numOH energy FE" + "\r\n");
          bw.flush();
      
          //TODO
          for (int entryNum = 0; entryNum < enerIn.numEntries(); entryNum++) {
              String fileName = enerIn.getFileName(entryNum);
              //POSCAR infile = new POSCAR(STRUCT_DIR + fileName, true);
              POSCAR infile = new POSCAR(ROOT + "/structures/" + fileName, true);

              //System.out.println(fileName);

              Structure structure = new Structure(infile);

              int numPt = structure.numDefiningSitesWithElement(Element.platinum);
              int numCu = structure.numDefiningSitesWithElement(Element.copper);
              
              int numP = structure.numDefiningSitesWithElement(Element.phosphorus);
              int numS = structure.numDefiningSitesWithElement(Element.sulfur);

              //int totalAtoms = numPt + numNi + numMo + numCu;
              int totalAtoms = numPt + numCu + numP + numS;

              
              if(numP!=0){
                  //System.out.println("the structure contains Ni-*OH: " + fileName);
                  continue;
              }
              
              
              if(numS!=0){
                  //System.out.println("the structure contains Ni-*OH: " + fileName);
                  continue; //exlude NPs with *OH
              }
              
              //if (((totalAtoms) > 50) && ((totalAtoms) < 85)) {continue;}
              //if (((totalAtoms) > 50) && ((totalAtoms) < 100)) {continue;}
              if (((totalAtoms) > 1) && ((totalAtoms) < 100)) {continue;}
              

              double delta = Math.abs(structure.getDefiningVolume() - appliedCE.getSuperStructure().getDefiningVolume());
              
              //the supercell of this structure is the same as the implicit one in the appliedCE
              if (delta < 0.1) {
                  //appliedNanoparticleCE.decorateFromStructure(structure);
                  //appliedCE = appliedNanoparticleCE;
                  
                  appliedCE.decorateFromStructure(structure);

              } 
              //the supercell of this structure is different from the implicit one in the appliedCE
              else {
                  System.out.println("!!! Attention: " + fileName + " has different shape of super cell.");
                  appliedCE = new FastAppliedCE(ce, ce.getBaseStructure().getSuperStructure(infile));
              }
              
              //SuperStructure superStructure = appliedCE.getSuperStructure();
              
              
              //appliedCE.decorateFromStructure(new Structure(infile));
              appliedCE.getSuperStructure().setDescription(fileName);
              
              
              double value = enerIn.getEnergy(entryNum); 
              
              if(value==Double.NaN){
                  continue;
              }

              
              // formation energy: formation E = E(PtxNiy)-xE(Pt)-yE(Ni)
              value -= numPt * PtEnergyBulkAccurate; // From bulk calCulations
              value -= numCu * CuEnergyBulkAccurate; // From bulk calCulations             

              //TODO
              // Sulfur stands for Pt-OH: "value += numS * OHEnergyCorrect;" used to correct the *OH from the group of "Pt-OH", which is referenced to 0.5*E(H2)-E(H2O)
              //                          "value -= numS * PtEnergyBulk;"    used to account for the Pt from the group of "Pt-OH", , which is referenced to E(Pt-bulk)
              value -= numS * PtEnergyBulkAccurate; // From bulk calCulations
              value -= numP * CuEnergyBulkAccurate; // From bulk calCulations
              
              //TODO
              //correction to *OH term, referenced to 0.5*E(H2)-E(H2O)
              value += numP * OHEnergyCorrectAccurate;
              value += numS * OHEnergyCorrectAccurate;

              
              value /= appliedCE.numPrimCells();
              
              //Status.detail(fileName + "      " + enerIn.getEnergy(entryNum));

              Status.detail("Adding structure " + (entryNum + 1) + "/" + enerIn.numEntries() + ", totalAtoms=" + totalAtoms + ", " + fileName + ". numPt: " +  numPt +", numNi: " + numCu + ", numS: " + numS + ", value / site: " + value + ", appliedCE.numPrimCells(): " + appliedCE.numPrimCells());
              //Status.detail("Adding structure " + (entryNum + 1) + "/" + enerIn.numEntries() + ", " + fileName + ". numPt: " +  numPt +", numNi: " + numNi + ", numMo: " + numMo + ", value: " + value);
              
              bw.write((entryNum + 1) + "   " + fileName + "   " + numPt + "   " + numCu + "   " + enerIn.getEnergy(entryNum) + "   " + value + "\r\n");
              bw.flush();
              
              appliedCE.activateAllGroups();
              ce.getStructureList().addCorrelations(appliedCE, value);
          }
          
          bw.flush();
          bw.close();
          fw.close();

      }
      
      catch (IOException e) {
      e.printStackTrace();
      }

      
      Status.basic("Added " + ce.getStructureList().numKnownStructures() + " structures.\n\n");
      
      return ce;
    }
    



    public static ClusterExpansion fitCE(ClusterExpansion ce,String StructList) {
        
        /**
         * Build the cluster expansion and generate NN clusters with up to 4 sites.
         */ 
        //ClusterExpansion ce = buildCE();

        /**
         * Read in correlations from a file called "energies.txt"  
         * This should just be a list of fileNames followed by energies (per unit cell)
         * Each fileName should correspond to a POSCAR in the "nanowire/2nm/ce/poscars" folder
         */    
        // Order the species in this array the same way they are ordered in the POSCAR files.
        //addKnownStructures(ce);
        //addKnownStructuresWithoutMo(ce);
        //addKnownStructuresRuleOutStructWithHighFormationE(ce, FECutoff);
    	addKnownStructures_PtCuOH_CalculateOHCoverage(ce, 1.0,StructList);
        
        ClusterGroup[] activeGroups = ce.getAllGroups(true);
        
        System.out.println(activeGroups.length + " cluster groups.");
        
        //ILinearCERegGenerator generator = new ExpExpRegGenerator(activeGroups);
        
        //ILinearCERegGenerator generator = new ExpPowRegGenerator(activeGroups);
        
        ILinearCERegGenerator generator = new PowExpRegGenerator(activeGroups);
        //ILinearCERegGenerator generator = new PowPowRegGenerator(activeGroups);
        generator.setMultiplicityExponent(0);

        double[] initState = new double[generator.numParameters()];
        Arrays.fill(initState, 1);
        //initState = new double[] {0.8823124822899626, 0.16174455365217855, 0.11868777414268145, 0.0922548990460885};
        //initState = new double[] {0.6580633972197518, 0.14685427063104411, 0.11907807350546687, 0.1319838515075602};
        //initState = new double[] {1E-8, 1E-8, 10, 10};
        
        //20191121
    // initState = new double[] {1E-8, 1E-8, 10, 10};
     initState = new double[] {1E-8, 1E-8, 10, 10};
//    initState = new double[] { 
//    		9.999896287184006E-9, 
//    		9.47389118410421E-9, 
//    		4.835384023847213, 
//    		3.3922112504813082};
        //initState = new double[] {1.0000017820487327E-8, 9.394658097219738E-9, 4.227456360502615, 2.9162053123814786};
        
        //20200122, 20200123
        // the set for trainingSet including Pt(111) slabs
        //initState = new double[] {1.0138695063838821E-8, 9.417254064670763E-9, 0.6590492647968678, 6.077005051205868};

        //CE: 20200124_n=354
        // the set for trainingSet including Pt(111) slabs
        //initState = new double[] {9.978191734636064E-9, 9.407220917090453E-9, 4.084263541427147, 2.974329932793443};

        //CE: 20200124_n=363
        // trainingSet excluding Pt(111) slabs
        //initState = new double[] {1.0138695063838821E-8, 9.417254064670763E-9, 0.6590492647968678, 6.077005051205868};

        //CE: 20200213_n=383
        // trainingSet excluding Pt(111) slabs
        //initState = new double[] {1.0138496562698203E-8, 9.318735467371137E-9, 0.6516398871070076, 5.036628490864908};
        //CE: 20200213_n=425-singleOH
        // trainingSet excluding Pt(111) slabs
        //initState = new double[] {1.0138690077370717E-8, 9.317896093307552E-9, 0.6512046578222509, 5.033281982901187};
        
        
        System.out.println();
        System.out.println("Initial Parameters: ");
        for (int paramNum = 0; paramNum < initState.length; paramNum++) {
          System.out.println(initState[paramNum] + ", ");
        }
        System.out.println();
        System.out.println();
        
        double[] values = ce.getStructureList().getValues(true);
        double[][] correlations = ce.getStructureList().getMultCorrelations(activeGroups, true);
        double[] regularizer = generator.getRegularizer(initState, null);
        
        L2FullRLSFitter fitter = new L2FullRLSFitter(correlations, values);
        fitter = (L2FullRLSFitter) fitter.setRegularizer(generator, initState);
            
        System.out.println("************ Before oCuimizing score *********");
        System.out.println("GCV: " + fitter.getGeneralizedCVScore() + ", ");
        System.out.println("LOO CV: " + fitter.getLKOCVScore(1) + ", ");
        System.out.println();
        
        fitter.useLKOCVScore(1);
        //fitter.useGCVScore();
        //fitter.useUnweightedLOOCVScore();
        fitter.setNoSingularGuarantee(false);
        fitter.allowIncrementalUpdates(false);
        
        RegularizerSelector selector = new RegularizerSelector(fitter, generator, initState);
        //selector.setMaxAllowedGrad(5E-3);
        //selector.setMaxAllowedGrad(1E-5);
        selector.setMaxAllowedGrad(1E-8);
        selector.findMinPolakRibiere();
        fitter = (L2FullRLSFitter) selector.getFitter();
        
        System.out.println("************ After oCuimizing score *********");
        System.out.println("Score: " + fitter.scoreFit() + ", ");
        System.out.println("GCV: " + fitter.getGeneralizedCVScore() + ", ");
        System.out.println("LOO CV: " + fitter.getLKOCVScore(1) + ", ");
        System.out.println("Unweighted LOO CV: " + fitter.getUnweightedLOOCVScore() + ", ");
        System.out.println("RMS Error: " + fitter.getRMSError() + ", ");
        System.out.println();
        double[] finalState = selector.getBoundedParameters(null);
        System.out.println("Final Parameters: ");
        for (int paramNum = 0; paramNum < finalState.length; paramNum++) {
          System.out.println(finalState[paramNum] + ", ");
        }
        System.out.println();
        double[] eci = fitter.getCoefficients();
        System.out.println("\nECI: ");
        for (int varNum = 0; varNum < eci.length; varNum++) {
          System.out.print(eci[varNum] + ", ");
        }
        System.out.println();
        System.out.println();

        ce.setECI(eci, activeGroups, false);
        
        int numFunctionGroups = 0;
        for (int clustNum = 0; clustNum < ce.numGroups(); clustNum++) {
          ClusterGroup group = ce.getClusterGroup(clustNum);
          for (int functNum = 0; functNum < group.numFunctionGroups(); functNum++) {
            String description = "Group number " + clustNum + ", Funct group number: " + functNum + ", Sites per cluster: " + group.numSitesPerCluster() + ", Multiplicity: " + group.getMultiplicity() + ", Max distance: " + group.getMaxDistance() + ", ECI: " + group.getECI(functNum) + ", Correlation coefficient: " + group.getCorrelationCoefficient(functNum);        
            Status.detail(description);
            numFunctionGroups++;
          }
        }
        Status.detail("Num total function groups: " + numFunctionGroups);
        
        Status.detail("");
        Status.detail("Structure information: ");

        double[] cvErrors = fitter.getCVErrors(false);
        StructureList list = ce.getStructureList();
        Status.detail("the # of structures: " + list.numKnownStructures());
        double totalSqCVErrorPerAtom = 0;
        double totalAbsCVErrorPerAtom = 0;
        int numNonEmpty = 0;
        for (int structNum = 0; structNum < list.numKnownStructures(); structNum++) {
          ClustStructureData structData = list.getKnownStructureData(structNum);
          double value = structData.getValue();
          double predValue = list.predictValue(activeGroups, structNum);
          SuperStructure structure = structData.newSuperStructure(ce);
          double totalatom = (structure.numDefiningSites() - structure.numDefiningSitesWithSpecies(Species.vacancy));
          double ratio = (double) structure.numDefiningSites() / totalatom;
          double cvErrorPerAtom = cvErrors[structNum] * ratio;
          if (!Double.isInfinite(ratio)) {
            totalSqCVErrorPerAtom += cvErrorPerAtom * cvErrorPerAtom;
            totalAbsCVErrorPerAtom += Math.abs(cvErrorPerAtom);
            numNonEmpty++;
          }
          double weight = fitter.getWeight(structNum);
          int numCu=structure.numDefiningSitesWithSpecies(Species.copper);
          int numPt=structure.numDefiningSitesWithSpecies(Species.platinum);
          int numS=structure.numDefiningSitesWithSpecies(Species.sulfur);
          double KnowntotalEnergy=value*structure.numDefiningSites()+numCu * CuEnergyBulkAccurate + numPt * PtEnergyBulkAccurate+numS*PtEnergyBulkAccurate-numS*OHEnergyCorrectAccurate;
          double PretotalEnergy=predValue*structure.numDefiningSites()+numCu * CuEnergyBulkAccurate + numPt * PtEnergyBulkAccurate+numS*PtEnergyBulkAccurate-numS*OHEnergyCorrectAccurate;
          double CVtotalEnergy=(value - cvErrors[structNum])*structure.numDefiningSites()+numCu * CuEnergyBulkAccurate + numPt * PtEnergyBulkAccurate+numS*PtEnergyBulkAccurate-numS*OHEnergyCorrectAccurate;
          double KnowntotalFE=value*structure.numDefiningSites();
          double PretotalFE=predValue*structure.numDefiningSites();
          double CVtotalFE=(value - cvErrors[structNum])*structure.numDefiningSites();
          
          //Status.detail(structData.getDescription() + ", Known Value: " + value + ", Predicted value: " + predValue + ", CV value: " + (value - cvErrors[structNum]) + ", CV Error: " + cvErrors[structNum]*structure.numDefiningSites() + ", CV Error per atom: " + cvErrorPerAtom + ", Weight: " + weight);
          Status.detail(structData.getDescription() + ", Known Value energy: " + KnowntotalEnergy + ", Predicted value energy: " + PretotalEnergy + ", CV value energy: " + CVtotalEnergy + ", Known Value: " + value*ratio + ", Predicted value: " + predValue*ratio + ", CV value: " + (value - cvErrors[structNum])*ratio+", CV Error: " + cvErrors[structNum]*structure.numDefiningSites() + ", CV Error per atom: " + cvErrorPerAtom + ", Weight: " + weight);
        }
        
        
        System.out.println("\n\n\nlist out structures with absolute value of CV Error per atom larger than 0.005 eV.\n");
        int numStructLargeError = 0;
        for (int structNum = 0; structNum < list.numKnownStructures(); structNum++) {
            ClustStructureData structData = list.getKnownStructureData(structNum);
            double value = structData.getValue();
            double predValue = list.predictValue(activeGroups, structNum);
            SuperStructure structure = structData.newSuperStructure(ce);
            double ratio = (double) structure.numDefiningSites() / (structure.numDefiningSites() - structure.numDefiningSitesWithSpecies(Species.vacancy));
            double cvErrorPerAtom = cvErrors[structNum] * ratio;
            if (!Double.isInfinite(ratio)) {
              //totalSqCVErrorPerAtom += cvErrorPerAtom * cvErrorPerAtom;
              //totalAbsCVErrorPerAtom += Math.abs(cvErrorPerAtom);
              //numNonEmpty++;
            }
            double weight = fitter.getWeight(structNum);
            if(Math.abs(cvErrorPerAtom) >= 0.005){
                numStructLargeError++;
                Status.detail(structData.getDescription() + ", Known Value: " + value + ", Predicted value: " + predValue + ", CV value: " + (value - cvErrors[structNum]) + ", CV Error: " + cvErrors[structNum] + ", CV Error per atom: " + cvErrorPerAtom + ", Weight: " + weight);
            }
        }
        
        System.out.println("\nThe # of structures with CV Error per atom larger than 0.005 eV is: " + numStructLargeError );
        
        
        System.out.println("\n\n\nlist out structures with absolute value of CV Error per atom larger than 0.010 eV.\n");
        int numStructLargeError_10 = 0;
        for (int structNum = 0; structNum < list.numKnownStructures(); structNum++) {
            ClustStructureData structData = list.getKnownStructureData(structNum);
            double value = structData.getValue();
            double predValue = list.predictValue(activeGroups, structNum);
            SuperStructure structure = structData.newSuperStructure(ce);
            double ratio = (double) structure.numDefiningSites() / (structure.numDefiningSites() - structure.numDefiningSitesWithSpecies(Species.vacancy));
            double cvErrorPerAtom = cvErrors[structNum] * ratio;
            if (!Double.isInfinite(ratio)) {
              //totalSqCVErrorPerAtom += cvErrorPerAtom * cvErrorPerAtom;
              //totalAbsCVErrorPerAtom += Math.abs(cvErrorPerAtom);
              //numNonEmpty++;
            }
            double weight = fitter.getWeight(structNum);
            if(Math.abs(cvErrorPerAtom) >= 0.010){
                numStructLargeError_10++;
                Status.detail(structData.getDescription() + ", Known Value: " + value + ", Predicted value: " + predValue + ", CV value: " + (value - cvErrors[structNum]) + ", CV Error: " + cvErrors[structNum] + ", CV Error per atom: " + cvErrorPerAtom + ", Weight: " + weight);
            }
        }
        
        System.out.println("\nThe # of structures with CV Error per atom larger than 0.010 eV is: " + numStructLargeError_10);
        
        
        System.out.println("\n\n\nlist out structures with absolute value of CV Error per atom larger than 0.020 eV.\n");
        int numStructLargeError_20 = 0;
        for (int structNum = 0; structNum < list.numKnownStructures(); structNum++) {
            ClustStructureData structData = list.getKnownStructureData(structNum);
            double value = structData.getValue();
            double predValue = list.predictValue(activeGroups, structNum);
            SuperStructure structure = structData.newSuperStructure(ce);
            double ratio = (double) structure.numDefiningSites() / (structure.numDefiningSites() - structure.numDefiningSitesWithSpecies(Species.vacancy));
            double cvErrorPerAtom = cvErrors[structNum] * ratio;
            if (!Double.isInfinite(ratio)) {
              //totalSqCVErrorPerAtom += cvErrorPerAtom * cvErrorPerAtom;
              //totalAbsCVErrorPerAtom += Math.abs(cvErrorPerAtom);
              //numNonEmpty++;
            }
            double weight = fitter.getWeight(structNum);
            if(Math.abs(cvErrorPerAtom) >= 0.020){
                numStructLargeError_20++;
                Status.detail(structData.getDescription() + ", Known Value: " + value + ", Predicted value: " + predValue + ", CV value: " + (value - cvErrors[structNum]) + ", CV Error: " + cvErrors[structNum] + ", CV Error per atom: " + cvErrorPerAtom + ", Weight: " + weight);
            }
        }
        
        System.out.println("\nThe # of structures with CV Error per atom larger than 0.020 eV is: " + numStructLargeError_20);
        
        
        System.out.println("\n\n");
        Status.basic("RMS CV score per atom: " + Math.sqrt(totalSqCVErrorPerAtom / numNonEmpty));
        Status.basic("Mean Abs CV score per atom: " + totalAbsCVErrorPerAtom / numNonEmpty);
        Status.basic("LOOCV score: " + fitter.getLOOCVScore());
        Status.basic("Final score: " + fitter.scoreFit());
        Status.basic("RMS Error: " + fitter.getRMSError());
        
        Status.basic("");
        
        return ce;
      }
       

    
    public static ClusterExpansion getPreFittedCE20250428(){

        ClusterExpansion ce = buildCE_PtCuOH();
        ClusterGroup[] activeGroups = ce.getAllGroups(true);

        
        double[] eci = new double[]{
        		0.42439570447820074, -0.41267771618458937, -0.26284628208416194, 0.21941955203431007, -6.330960420207754E-4, -5.146630434421651E-5, 4.497776297135436E-4, 6.935542071022857E-4, -2.935919966405437E-4, -1.482059043204087E-4, 2.2026480509579338E-4, 4.260087382476213E-4, 8.020761917942043E-4, 3.1745194743049154E-4, 1.849234307946812E-4, 1.0754361093606724E-4, -2.889935528735954E-4, 1.384137257338658E-4, 5.651158871111256E-4, 9.19500258533823E-5, 1.5321159509401605E-4, -2.1882142397944282E-4, -3.838600081731841E-4, 0.0016179046830853416, -9.300697979297556E-5, -2.950427381047671E-4, -7.018070137929655E-4, 4.118749202678962E-4, -0.0013086306685434757, 7.22570916878329E-4, 0.0010648497171956706, 0.002001934163167851, -4.196953081991861E-5, -3.384855987973445E-4, -1.8170980395687542E-4, -0.0013921913377827888, 5.788218616009644E-5, 0.0021203977899296456, 3.45780957911046E-4, -6.328736702245763E-4, -0.008185237724536165, 0.01398826345921389, 0.0017948332436903414, 0.005776206184853638, -0.006656512992508745, 0.0012884671180817264, 1.4181465615523397E-4, -8.92522974604848E-6, 8.572410978626224E-5, -6.968540424901609E-5, 4.569947350109051E-5, 4.021466637626689E-5, 5.1735667301835986E-5, -3.285524538182427E-6, 3.94519158657472E-5, 3.4697373164291334E-5, -1.064184772946869E-5, -6.315248706660614E-5, 2.0702128674110476E-6, 8.104275940354344E-5, 5.033295049011732E-5, 2.5218083947967096E-5, -7.060373506659334E-5, -3.355850547115402E-6, -4.807704386899231E-5, -1.2236067206205029E-4, -7.758089669585369E-5, 1.2950636473223198E-4, 1.6553371558564178E-5, -7.299911738169532E-5, 9.835217041045605E-5, -6.102505236165205E-5, -7.705081036386515E-5, -8.384406452865631E-5, -8.021106226182434E-5, -3.6608179799460516E-5, -5.9068121134223106E-5, -1.8778101270814848E-4, -1.0277503615757633E-4, -4.355532258871336E-5, -8.143247002575984E-5, -4.199584197451395E-5, 7.376348493667861E-5, -1.9126238220991756E-4, -3.2416439151247006E-5, 1.0181492546606358E-4, 3.811666563352653E-5, -4.88670291548185E-5, 1.2336174399791948E-4, 5.9320968113549255E-5, -3.4153629220184945E-5, 6.504163945945246E-5, 2.084993798643492E-5, -7.593518530454868E-6, -2.4132241847071162E-4, -7.736344358154346E-6, -8.84320764481597E-5, -7.576269848003775E-6, -2.1368950533112718E-4, -7.571116208091033E-5, 9.010300985152694E-5, -1.0156780874369387E-5, -4.458918751688249E-5, -1.4110252019040176E-4, 3.605118424232501E-5, -3.5668070308608746E-5, -9.980802001503634E-5, -3.322812752321047E-5, -1.1479254214413412E-5, -4.459519343497781E-5, -1.4590625635815466E-4, -2.1166061200929287E-5, 9.15013900711012E-5, 4.629217247855994E-5, 4.2167808587748046E-6, 6.089157276883278E-5, -3.5823506479279454E-5, -1.8578538541376421E-6, 1.337279096216143E-5, 7.173536418902601E-5, 3.982793469636814E-6, -1.2479659597425252E-4, 1.347109777820963E-4, -4.985571164607683E-6, -8.746799543771222E-5, 7.724112417882567E-5, 2.8082187238334666E-5, 1.7264762137722198E-5, -3.615212416818725E-6, -6.656644922738153E-5, -4.323821145907347E-5, 7.240737716401185E-5, 4.029781586070008E-5, -7.196145131729836E-5, 3.0819826714378884E-5, -2.7689325961893774E-5, -5.222674038648026E-5, 3.905436145894975E-5, 5.607368234776533E-5, -5.659685056337701E-5, 6.572592145059897E-5, -1.8300919762508083E-4, 1.2641362855513733E-5, 1.3466016545976567E-5, -2.9419211628058106E-5, 5.981174411579665E-5, -1.1852378218116244E-4, -2.274641111624648E-4, 1.4839724259442175E-4, 1.0192540115847373E-4, 2.5101409272620945E-4, 1.3612084455950219E-4, -1.0731468901455099E-4, 1.2703106273778686E-4, 8.537457196129694E-5, -1.0475260648526441E-5, -6.47267935973625E-5, 1.0068884260520244E-4, -1.7056924037445187E-5, 1.2367527661236166E-4, 3.235944879673852E-5, -6.866763138969889E-6, -8.454878892077966E-5, 1.739059638828505E-5, 2.0353924660466982E-6, 9.101816071425961E-5, 8.33165885584259E-5, -1.6707549389250592E-5, -1.94735851173633E-5, 1.524654780054313E-4, 2.2198416488738645E-5, -2.015416822187012E-5, -5.9435798323972944E-5, 9.356259170881231E-5, 5.3161773992602135E-5, -1.192267752448406E-4, -3.748404041079359E-5, -2.36928458376362E-4, -1.4251531771952708E-5, -1.3629406989519588E-4, -5.1126470915902056E-5, 2.22600305412584E-5, 5.384677479774308E-5, 2.1133798086617195E-5, -5.209143760849522E-5, 4.3262280386868746E-5, -2.3378773147050835E-5, -1.5605380404404637E-4, 5.472603045318658E-5, -1.2123931229054974E-5, -2.8034828240954685E-5, -4.869727018565889E-5, -5.554836874197E-5, 1.273926902990944E-4, -1.0006442057639703E-4, -2.5174491682676547E-5, 6.96297577994498E-5, -2.1152413017862942E-6, -2.8995549853007285E-5, 9.207394853803785E-5, 9.360375929977349E-5, 2.0158843288493877E-5, 3.794674722051514E-5, -2.835248892451344E-5, 2.4906257877328752E-5, -5.1143840666539875E-6, 6.4802155332818745E-6, 4.4686848921365874E-5, -2.079290908056034E-5, -3.278725080108416E-5, -1.2691912605840294E-5, 1.303077247576806E-4, -3.406114765515736E-5, 4.799665993085456E-5, 1.0461417244748278E-4, -2.5690686131403517E-5, 1.8865164331992312E-5, -7.689230809304615E-5, -7.98895981992398E-5, 1.806674285666915E-5, 1.3246453671206943E-4, 1.6026181034282244E-4, -1.8339540879618033E-5, 4.939538963137646E-5, 2.0627355659649793E-5, 7.54506531747065E-5, -7.932751952183237E-5, -6.087126241630797E-5, -1.0944874733193055E-4, 1.61488852951093E-4, 1.2926643424239492E-4, 1.2686562662755918E-4, 7.554193021493127E-5, 4.360780459046476E-5, 4.2732966423410235E-5, 3.5020744879000337E-4, -3.1031527211958926E-4, 8.389395126771892E-5, 1.085281300426677E-5, -9.861453831527811E-6, -6.117485840795994E-5, 1.584886569592752E-4, -1.482841202597729E-4, -1.0412926291917743E-4, -3.851711019147719E-4, 2.768022485959698E-4, 9.610496365800416E-5, 2.918250760136354E-4, 1.7991066942782408E-4, 1.8645910306314128E-4, 1.0134167852389464E-4, 5.505388301249581E-5, -6.187457270798863E-5, -1.603760261526671E-4, -1.4450231079195242E-4, -1.2436736772002133E-4, -8.862734741126033E-5, 2.9475034039969095E-4, 1.4480379102616708E-4, -1.6432700740211942E-4, 8.296901199045248E-5, 2.8545943711100197E-5, 8.258896613070818E-5, -1.7691679214952571E-4, 6.625059424855695E-5, 2.9205323768590266E-4, 1.7791775801855036E-4, 1.026893295244316E-4, -5.1732932036294144E-5, 2.2173616252223194E-5, -6.281580601116832E-5, 2.336580859958748E-5, 8.816046076298711E-5, -1.2516854489816667E-4, 1.476844677336317E-5, -9.30847513766625E-5, -1.9659755234549547E-4, -3.3487970128081614E-5, -7.129024678677841E-5, 2.231426167114553E-5, -8.284907813385645E-5, 2.2982112423367007E-5, -2.723088878350821E-5, -5.660480943882329E-5, -1.1554646393241441E-5, -1.5239051003285153E-5, -4.5952417215042174E-5, -6.434034105337605E-5, -8.073376614561867E-5, -1.1670958777473817E-4, -6.387652466376126E-6, 9.31907989788032E-5, 6.891567304322804E-5, -7.453741786329557E-5, -1.3688297113872983E-5, -4.321890743782392E-5, 5.537329831431688E-5, 2.6621873855131595E-5, 1.2808696115906848E-4, 3.452980005629741E-6, -1.4878867466541074E-4, -1.6242333665142036E-5, 2.1405563197657296E-5, -1.7468708617508072E-5, 7.523095374832589E-5, -3.1475409484197844E-5, -2.626126022742259E-5, 1.4908092470095058E-6, 4.804919735922758E-5, -2.023079435781866E-5, 3.7245291582994715E-5, 2.130633533873119E-5, 2.4350187511831822E-5, 1.8073005627245592E-5, -6.801286597782775E-5, -4.214998947421264E-5, -2.404759845366783E-5, 1.1189527848529113E-4, 9.57639954218964E-5, 1.3606899047086118E-5, -8.559188110763942E-6, -1.5539562475894164E-4, -1.2343284810446137E-4, 1.6865586080633266E-4, 1.2907841328710368E-4, -1.7140540620842054E-5, 4.526203149597674E-5, -3.640271775944004E-5, 2.3347647852382455E-4, -4.199217181932724E-5, -1.9232793631485112E-4, -4.135309972124329E-4, 8.204084754591737E-5, 6.776603814494936E-5, -1.907719161777649E-4, 2.8112581148820315E-5, 2.5869105396778E-4, 1.7557042128840683E-4, 2.4119886396278822E-4, -2.4480188802333824E-4, 1.759649319411083E-4, -6.476252080066247E-5, -4.188539628210801E-5, 6.356096753680832E-5, 6.238932663774066E-5, 8.455088354190551E-5, 4.290726707670278E-6, 1.7145212234360497E-4, -3.1576697147214783E-4, 3.540768307885231E-4, 2.789713714745762E-4, -2.112828248462817E-4, -5.754919111172763E-6, 3.6561778993443594E-5, 1.9637256068626152E-4, -3.086596477336334E-4, 4.7351257164984784E-5, 7.577768818893549E-5, -3.340365247220256E-5, 4.401443189194776E-6, 9.043532221965356E-5, 4.7567680343696406E-5, 1.2288431346656886E-4, -7.739628297058487E-5, 2.7863352052112655E-5, -2.5483006236452083E-5, -7.803605832398808E-6, 6.332039764836428E-5, 3.0704181230665345E-5, 1.562784406586304E-5, 5.788427821333542E-6, -4.253721695626489E-4, 3.2678331374869085E-4, 5.67568872502644E-5, -1.8269069057053056E-5, 2.2658611358505887E-4, 2.9764548510184154E-5, -1.1970644122163976E-4, 2.622229978509541E-5, 4.594554370301873E-5, -3.100957838948759E-4, 2.2222253920069733E-4, -6.471073896563428E-5, -5.366157671060634E-5, 1.7517761318619485E-4, 1.7539826805135217E-4, 2.0887897748678164E-6, -2.501421795710057E-5, 5.5883175162307896E-5, 1.14288000476713E-4, -1.2016807769339027E-4, 5.723449697504421E-5, 1.3851509425199388E-4, 2.637305724081246E-5, 6.921541487150962E-5, 1.069070031765873E-5, 8.839731870076417E-5, 8.659621704061593E-5, -6.306947760253978E-5, 1.9975045932579816E-4, -4.040317898554302E-5, 2.8084218532096333E-4, -4.6894766407401534E-5, 2.4951284468915924E-5, 6.884501536914833E-5, -4.7601899930713176E-4, 2.3176806061311894E-4, -5.4169316611341056E-5, -6.680006823109319E-5, -2.2301502131447798E-4, -1.1677230923629045E-5, -2.4188228776109155E-4, 9.857609826885438E-6, -9.63471246810399E-5, 6.729988861767927E-5, -1.433974000935763E-4, -3.863628545482112E-4, -1.277196766234709E-4, 5.1302634666960945E-5, 2.1285256606230069E-4, 4.822701044321365E-5, -1.6971210788052637E-4, -1.1671649301764689E-4, -4.3897183033805853E-4, 4.985052525427469E-5, -1.6747828701990652E-4, 3.428464831294886E-4, 1.545397954219041E-4, 2.011107531037365E-4, -3.6863549770345394E-4, 3.2715246530345965E-4, -2.6864705732679295E-4, -2.6642939249529873E-4, -3.1850901432109446E-5, 2.2637765562327772E-5, -7.194091586060214E-5, -0.0024567436792979347, -0.0010487689934656624, 0.0021478757255084524, 0.0019171273256768373, 5.35447767842534E-4, 3.337381099824399E-6, 5.28256855234055E-4, 1.4368647625095932E-4, 5.765823218428778E-5, -1.6829962688750654E-5, 8.33645760092938E-6, -9.95815223404353E-7, 2.94590081444078E-5, -3.335214699038986E-5, 4.3301944880041326E-8, 7.421795302459816E-6, -1.909883228455242E-5, -1.1987354666894573E-5, -2.1501211291997946E-5, 1.961717199052554E-5, -4.877932785508909E-6, -1.885090235348441E-5, 1.5281324798395446E-5, -7.660933660858307E-7, -1.397482924689851E-5, -3.294284830705936E-6, -3.946031247368739E-6, -8.40293669058384E-6, 1.1984135490053678E-6, -2.168938692704821E-6, -2.0978665032523464E-5, -2.369199861677884E-6, -1.6529660686551845E-5, -1.9996204916017737E-5, 5.502994162185845E-6, -6.772397446598577E-6, -2.997645880152555E-6, 1.4865188448535613E-6, -2.667129903153806E-6, 1.1783844425050302E-5, -9.338261469411238E-6, -8.560960395822236E-6, 1.0073993813158703E-5, 4.421968056986247E-6, 3.049572920988102E-5, -3.672347446552135E-6, -1.7016644861889488E-6, 2.8483953550612756E-5, 1.6394352584787625E-5, 2.1640816927273078E-5, 1.9382852110256567E-5, -2.0752127162119993E-5, -1.1637519110680112E-5, -1.6062745583668436E-5, -8.622719846541546E-6, -6.55960151377255E-7, -9.58463211898214E-7, 1.4044237240338E-5, -3.166090327327842E-6, 6.888323024418011E-6, -2.0028006628945716E-5, 1.015671099173672E-5, -6.379446362580511E-6, -1.4888521459924513E-5, 1.1471181238602995E-5, -1.1395711043294537E-5, 1.2375547666413755E-5, 1.72985615393903E-5, 1.1624436869268201E-5, 8.955190903128102E-6, -2.3302842793273756E-5, -1.7836333452371858E-5, 1.4342518847726589E-5, -1.3373480299247106E-5, -1.721656912353044E-5, -2.725228018127246E-6, -1.9448847452935986E-5, -1.4538357149913514E-6, 1.2331261230161608E-5, 2.590458603438442E-5, 4.373489646621216E-6, 2.8350811970574014E-6, 7.5404284748548794E-6, -1.125536171880042E-5, 1.4899397589186042E-6, 9.3257331369688E-6, -2.3914621365764277E-6, 1.1715669650738322E-5, -3.5844430028939694E-6, 1.3697390528867072E-5, 1.1946599390598011E-5, -1.177991086942845E-5, 7.715090245590667E-6, 2.165825456061296E-5, 9.029411163424296E-6, -2.6757070652012535E-6, 2.8457815449862978E-6, -7.056497370545867E-6, 7.273559539267699E-7, 5.725320911340948E-6, -6.653994890248871E-6, 1.158623752580204E-5, -5.466971818439718E-6, 6.637671547289224E-6, 1.296875941602521E-5, 1.1684907532988733E-5, -4.3755745510647495E-6, 8.958755981055199E-7, 2.9864814032462203E-6, 1.6221119608760687E-6, -1.0261098489984521E-5, -3.504193566985993E-6, 1.1056039061690254E-5, 1.9896117966329563E-5, 1.1886484295731137E-6, 9.333941068806473E-7, 2.3678428989655498E-6, 2.668714737082828E-6, -1.5859998942144293E-6, 1.8027876102963067E-5, 3.2792424770373843E-6, -1.4326997696041503E-5, -1.8305607259579137E-5, 2.6378930803677703E-5, -2.3592136623084388E-5, 1.8919173308704247E-5, 1.1516413618208737E-5, -2.10837936570063E-5, 1.1155608617783447E-5, 7.661797295190759E-6, 1.6620481466083324E-5, 5.490185594512554E-6, -2.6195672407055073E-5, 5.599192904232747E-6, -6.341529023267886E-6, -2.7242878727792158E-5, 1.9829159401532626E-5, 6.252973827255157E-6, 8.198587099016164E-6, -8.6224991947036E-6, -5.5600883264433215E-6, -2.5164555078349E-5, 7.298885949594269E-6, -7.37362990207123E-6, -2.6416356265344942E-5, 8.646416215896124E-6, -3.313086227238471E-6, 3.376947968175906E-5, 7.147137024796774E-6, 1.1494334472138268E-5, -4.137221005436412E-5, 2.938427629494301E-5, 8.365324321753242E-7, -3.650799238209823E-5, 2.9652401863982E-5, 1.1145977293781223E-5, -3.175875228570227E-5, -3.7522153889086272E-6, -9.850723356333649E-7, 1.8610015615953368E-5, 1.3105515031592262E-5, 1.368259211967024E-5, 1.823054004856592E-5, -1.232507427901923E-5, -7.315254347574937E-6, 8.783463944159187E-6, -1.1466575796446887E-5, -6.1236969595322645E-6, 4.675655001736001E-6, 1.9797017913585763E-6, 1.1344729136192333E-5, 1.6196522948708214E-5, 4.156456856080695E-6, 9.254786241943069E-6, 2.185514046689277E-5, 2.1109808285907014E-6, -5.875486532309904E-6, -4.0053011545535325E-5, 3.0068212210811705E-5, 1.40589444298303E-5, -2.4677005923331996E-5, 1.5351728894080763E-5, 6.964818746027501E-6, -1.6989856452846383E-5, 4.715315586508936E-6, 4.997749790638939E-6, 1.831496408405582E-5, -1.7839891309448085E-5, -6.842787676104856E-6, -4.381642996802854E-6, 5.4559894162847095E-6, -4.823462499331521E-6, 1.3463688346622774E-5, -1.4074179782519805E-5, -4.972508225507474E-6, -8.205785808786809E-6, 2.6513708273230574E-6, 3.599390129037166E-6, -1.4734676234834376E-6, 9.801898267530057E-6, 3.321422477366859E-5, -3.1000897098146414E-5, 1.0234423012244614E-6, -6.8154499733835E-6, 4.378927126128321E-5, -2.3930738138642355E-5, 1.3685401771247607E-5, 1.3227660417895054E-5, -3.098703552216128E-5, -4.469688972598591E-6, -1.3785603143546491E-5, 1.537460452684554E-5, -2.2282470239498586E-6, -1.8701171664189248E-5, 3.012107277801032E-5, -7.999683043619782E-6, -1.2177481733027101E-5, 6.127765717222136E-6, -8.135020161038897E-6, -8.659676779074893E-6, -2.5156067289811013E-6, -3.773575962127268E-6, -4.001641856443694E-5, 1.1085104242334253E-5, -1.2488353715192288E-5, -1.0721893400865434E-5, -4.538024966232318E-6, -1.3726479505524192E-5, -1.343111006211345E-5, 1.5560246444083423E-5, -4.717700957682742E-6, -4.0401193068872326E-5, 3.289611647248197E-5, -2.2780008298445287E-6, -6.519985923188545E-6, 1.0521332373390627E-5, 2.793752457107973E-6, -1.6147498333719118E-5, -1.9126434075908416E-6, 5.820079847615521E-6, -9.161472898200351E-6, -8.46068924422079E-6, -3.244580092917268E-5, -3.446726614951545E-5, 1.90421761669321E-5, -9.829947339298414E-6, -1.7990555953500964E-5, 1.4589944494628956E-5, 5.575240883265442E-6, -1.748190654692061E-5, 3.3893229910327005E-6, 1.1188486886866411E-5, -5.59766396256503E-6, 5.836892163047815E-6, 2.5898918055031195E-5, -4.517618030529225E-5, 2.600790403374867E-6, -1.1531722863626984E-5, 6.805273616241822E-5, -4.885933342797799E-5, 2.29537206864229E-5, 3.263568633475947E-5, -3.554708933262152E-5, 6.087193803616707E-6, 2.9355806572053944E-5, -1.2396235035945155E-5, 9.103219886646974E-7, -1.3492291064183537E-5, 1.6790984967395006E-6, -5.7044956896868185E-6, -1.4547349195054987E-5, 1.207220492655198E-5, -5.6773870674171195E-6, -3.779292442333132E-6, -1.0877212004389153E-5, -1.5172616235705622E-5, 1.236441315144313E-6, -1.418337198383357E-8, 1.427859197103247E-6, 1.777640433435132E-6, -1.756883453489077E-5, -3.170246784866236E-6, 2.718594163063446E-5, -2.738639601081866E-5, -2.4369228588856824E-5, -4.925008859096858E-5, 2.230020436927785E-5, 4.895805402644303E-6, -8.199780728258485E-6, 7.955388875790755E-6, -1.1195187066704443E-7, -3.228566660572717E-5, 3.211618355330407E-5, 1.3482333567229818E-5, 4.031718813439153E-6, -3.3763455301195205E-6, -2.196869178469141E-5, -1.1582282366759747E-5, -8.977301451615696E-6, 4.481152687219395E-6, 1.2292093628408097E-5, 9.137001350034875E-6, 1.0876054940524293E-5, -1.2265341298115656E-6, -6.806289222758416E-6, -3.6424006618681476E-6, -1.0676314156458921E-5, 1.1427496030801121E-5, 6.786196980596709E-6, -1.3140521321394066E-6, -1.604310299703135E-5, -1.8726169204093604E-5, -9.249075306635296E-6, 1.2042614361469119E-6, -3.5231035282329293E-6, -3.483826771132764E-6, 1.269047937003942E-5, 7.464234184170748E-6, -9.796547420061821E-6, 2.3794754594611348E-5, 1.650063018989248E-5, 1.3984465255869076E-5, -2.5007689498403318E-5, -1.0165527401830021E-5, -3.546427539320024E-6, -1.4466636092336227E-5, 3.036209102608156E-7, 4.117295205802174E-6, 3.001980649085245E-6, -5.024171499173544E-7, -5.508161335661399E-6, -2.964366156734538E-6, -2.722221755298293E-6, 3.4755665311372863E-7, 3.32307318035973E-6, 2.824735153765523E-5, -1.757833109754508E-5, 1.2805475287213785E-5, -8.606491619441749E-6, 1.5163277235082009E-5, -3.3373746435820875E-5, -6.660610172191765E-6, 2.8929344312624006E-7, -2.846580559149811E-5, -2.4879839374949975E-5, -5.379158287046465E-7, -4.69605166548287E-7, -4.0635758324089957E-7, -9.626871651537833E-6, 1.1531244961757864E-5, -1.3774531863524828E-6, -2.083545349019361E-6, -1.6049102396513762E-5, -7.2559558267204635E-6, -4.519066241680898E-6, 1.29261439530369E-6, 2.084526416342824E-5, 1.7023069543666523E-5, 1.1347520238813927E-5, 7.178318471862642E-6, 4.1589773367564295E-6, -8.348197514468575E-7, -1.1344146806750721E-6, 2.224563606254004E-6, 9.774057306948837E-6, 4.579712976864178E-6, -4.751947261397681E-6, 9.484885960798356E-6, 3.680368467272148E-7, -4.001972579072838E-6, -2.2382628730808346E-6, -9.74878750290356E-7, -1.4104027217471359E-5, 1.3991617571608092E-5, 1.4327038321378799E-5, 1.3830847990944223E-5, -1.6620064624611394E-5, -2.1431041200815903E-5, -7.049307096118791E-6, 2.52336037651589E-6, -7.251045911841822E-6, 1.171461448606654E-6, -7.283399649080612E-6, 8.734413900080701E-7, 7.128475364694131E-6, -1.087996877577953E-5, -5.756373210366247E-6, 1.5153193665610847E-6, -6.9256624535755755E-6, -5.114608384803113E-6, -2.8415673644151696E-6, 1.2615175690147869E-6, 4.035673062821078E-6, -2.2074110694226753E-5, 9.070606955501829E-6, 9.676669948416317E-6, 3.859578095609731E-6, -6.469080107372391E-6, 2.7500102291258544E-6, -2.3275619100079914E-5, 2.7825855318857826E-6, -1.9072243773614328E-5, -4.118318572396244E-5, -8.362705620055047E-6, -2.370531726089177E-5, -8.438927293701395E-6, -1.1062824312969063E-5, -1.0460126917797782E-5, -7.810168562251852E-6, 4.2286075404201E-6, 1.1715882262753122E-5, -1.593223657274552E-5, 4.7348915830185134E-6, -1.8273793878610493E-5, -4.688624372514716E-6, -1.4653478584295346E-6, 1.0139686934130393E-6, -1.2379044805083984E-5, -7.657528835714692E-6, 1.335667718529109E-6, -1.8889841516999545E-5, -2.7218891332010664E-5, -2.150134824357693E-5, -1.3778543201101168E-5, -3.842055789382223E-6, -3.814886915975112E-6, -1.2212924253454313E-5, -7.259251779794407E-6, 2.1615245173432773E-6, -1.0598733747266623E-5, 1.931362437009509E-5, -7.962843665092895E-6, -2.1263343311416286E-5, -2.207412479273991E-6, -1.1525844588626544E-5, 2.107445646846312E-5, -1.7156792034930422E-6, 1.3367938200201795E-5, 1.6775967897803695E-5, 2.076798074528846E-5, 1.141677524583648E-5, 1.0008460260679114E-6, 1.4643019458230114E-6, 1.4880933199646132E-6, -2.79244610771872E-5, 9.62570478788699E-6, 5.091325584249268E-6, 9.383901382609139E-7, -2.333347806745146E-5, -2.046524744613032E-5, -5.738295251892303E-6, -2.4261136484749965E-6, -3.803860136097903E-6, -1.5801763202689598E-5, 4.908620597425268E-5, 2.4085036671605067E-5, -2.6546563356273645E-5, 5.4427789990686704E-5, 2.7639486547340053E-5, -6.894080657026541E-6, 1.7376737070222067E-5, 1.1624241077603419E-5, -2.684358161672768E-5, 7.950405018209566E-6, -4.005304698188189E-6, -2.322258746332434E-6, -2.4488560158134748E-5, -1.2125825926766198E-5, 1.5666105331046953E-6, -7.749498176740782E-6, 9.407167793107005E-6, -6.460825834685158E-6, 3.585705679184665E-5, 2.2941818933639932E-5, -2.4250176441817286E-5, 3.4812852974658304E-5, 3.1842905963862765E-5, -7.04366680847994E-6, 1.8570335088573947E-5, 2.9271863961768214E-5, 4.3181614914681015E-6, -6.204986323922381E-6, 3.976523800043857E-6, 1.0520110095100978E-5, -5.1481993209373175E-6, -4.093058701984398E-6, 5.2174748399382714E-6, -3.6002967280115636E-5, -1.4599612411038236E-5, -2.5836916607378108E-6, -2.8130416213235226E-5, -1.731723519970106E-5, -2.1022388187099533E-6, -1.7937553487081676E-5, -2.896848700064715E-6, 4.5736495122425044E-6, -1.3785919110510456E-5, -2.1187956585829414E-5, -2.8793936656735294E-6, -2.4505614504311338E-5, -5.286341311692196E-6, 1.6740214276887364E-5, -1.8307115003238228E-5, 4.092696243612054E-6, 2.310187821241358E-5, -1.8040112781224013E-5, 9.965215884543549E-6, 2.5116167327816907E-5, -1.6160362153898322E-5, 5.4997247923306005E-6, -1.334694463711334E-5, 1.5437063154486938E-5, 1.5018795967661897E-5, 3.478622912611801E-6, 1.7510336917415936E-6, -5.9156244946803604E-6, -1.590726571635767E-5, 9.822478236783116E-6, -6.181439981090692E-6, -1.787266004408053E-5, 1.460318634390577E-6, 5.113305442686067E-6, -3.202576073208993E-6, 1.445447546721823E-5, -1.1837857333517332E-5, 3.5811215623633708E-6, -2.938663838780676E-6, -1.209876668635361E-5, -4.377505620220025E-5, 2.598253698855163E-5, 6.967197634131367E-6, -1.686784999631319E-5, -1.5800972204237096E-5, -2.8386011812455564E-5, -1.563858607025667E-5, 4.518046785802893E-6, -1.5980275248668983E-5, -9.604486167774304E-6, -4.318407847424006E-5, -4.4110302736859E-5, 1.597024597276009E-5, 3.0538684781269976E-5, 2.1470679112891786E-5, 1.524579326295239E-5, -1.836643118063863E-5, -3.3867082276581306E-6, 2.4782735773606207E-6, -3.351783465790556E-5, -2.2850825522290757E-5, 8.674423409141263E-6, 9.749004505940628E-6, -7.710877004477035E-7, 1.041802491785239E-5, 3.789445638870126E-6, -1.0723312883625935E-6, -4.066747237274279E-5, 1.4199396299464928E-5, -4.821233196694065E-6, -1.1287450991128573E-5, -1.2087352782356314E-5, -1.1028974942994124E-5, 2.1644926147960323E-5, -3.661447637598805E-5, -7.372215127921834E-6, 2.0112207595289998E-5, -4.3667381242182626E-5, -2.273183351741287E-5, 1.3377364804834091E-5, 7.221630302867155E-6, 3.788877193464796E-7, -8.431271385517941E-7, 1.703704492973107E-5, 9.979349803367066E-6, 3.594372625268655E-5, -4.363813103066519E-5, -1.7955794933934704E-5, 6.1709138727203405E-6, -3.5346032974739594E-6, 8.20541798311758E-6, -1.7486756150977487E-5, 3.4337017017778415E-5, 3.9711336869568264E-5, 4.081027519935079E-6, 1.4757155238171067E-5, 1.6709656829407005E-5, -1.4586383046940113E-5, -1.3073993745107634E-6, 1.1027653557010852E-5, -3.6858684871305357E-6, -1.6572918649510049E-6, 2.0776104220627915E-6, 1.8034696458143E-5, 1.581569025890274E-5, 1.9564312086269714E-6, 9.55222280397093E-6, 8.280500325864339E-6, -3.063144659583121E-7, 1.2374915427480239E-5, 1.4729358271381096E-5, 2.580747212842238E-6, 1.091819666507697E-5, 1.0588661245551113E-5, 1.0901586786313342E-5, 2.6264851289488356E-5, 8.699494151030024E-6, 1.0633044082915613E-5, 1.348120922552154E-5, 2.011576090127972E-5, 5.5379147688161545E-6, -7.100050384373224E-6, -5.356362988114156E-6, -1.2381245886258947E-5, 1.4238563902237156E-5, 8.643356961268722E-6, 1.6314092213427778E-6, 5.277384464320769E-6, 2.1678488138460877E-5, 1.054457541492873E-5, 5.344729053641251E-7, -7.222136094943779E-6, -2.894358660626402E-7, 7.946559401077538E-6, -1.470979003270756E-6, -3.0065677619865727E-6, 1.0400700113789678E-6, 2.3362420648936986E-5, 4.391097689744605E-5, 3.5840915776319884E-5, 1.3721533261852603E-5, -1.0194234487486548E-5, 1.4478647850677678E-5, 1.7339236967938224E-5, -9.127462015283286E-6, -1.7897857530313477E-5, -3.7565044131958325E-6, 2.0881421995692234E-5, 4.630135226630933E-5, 4.1467562947692695E-5, 2.0551418393157522E-5, -8.02045935297934E-6, 2.549381930354224E-5, 1.3584377100062626E-5, -1.7174500576102496E-5, 2.4493333968889433E-7, 7.23263102002554E-6, 1.7122972810563932E-6, -3.0751951967234945E-5, -8.985520307952521E-6, 1.3260975938815187E-5, -2.460013717271378E-6, -1.0264088348712307E-5, -1.1515418446447767E-5, -1.1543014425546838E-7, 1.1866205642842298E-6, 2.697933474885353E-6, -1.1015986863493338E-5, -1.879143524671084E-5, -8.572584161840088E-6, 9.859184620361183E-6, -4.919658897147711E-6, 7.593251633260852E-7, 1.850770175298392E-6, 1.515161479802296E-6, 9.687979266027186E-6, 1.0079529358506E-5, 4.339634076028567E-6, -1.6657076039074126E-6, 2.371391214936458E-6, 1.1900577007076252E-5, 9.44858029826301E-6, -4.226154620866861E-6, -1.253759898425192E-5, 5.464222360724615E-6, 1.247304598158902E-5, 1.4650061728207348E-5, -6.101839613196103E-5, 7.561030475377976E-5, 1.871148190839972E-5, -4.572787150584033E-5, -2.245174131432215E-5, -4.28367032485194E-6, 3.9470062353580415E-5, -9.486505167735637E-5, -1.7660886841741608E-5, 4.341081285313443E-5, 1.5067227084749127E-5, -1.571365607440367E-5, 5.0834095019561265E-6, -4.6618690497254867E-5, -4.330227277268839E-5, 2.374642573299557E-6, -2.8902998248092615E-6, -5.400918741342575E-6, -4.2545504436008325E-5, 2.1911298408186508E-5, -2.7557812594697362E-5, 2.0840082014857847E-5, 2.707906854088754E-5, -2.625539361662127E-6, -2.288555326777908E-5, -3.384587344827661E-6, -7.744110679101801E-6, 3.0642076707380894E-5, 1.501391997932311E-5, 1.7525672994168212E-5, -1.6295126759791403E-5, 5.544116823144596E-6, 1.2144805415958276E-5, 1.649776940952901E-5, 7.66125521590999E-6, 1.528987332765925E-5, -2.055707015044322E-6, 2.8964802637440265E-5, 2.1126062647039575E-5, -2.360457256013991E-5, -6.815887523468206E-6, 3.1439115029023926E-6, -2.4821992186894667E-5, -2.823761394949451E-5, 1.289693291501269E-5, 1.2164599525434777E-5, -1.1230428452562585E-5, -6.3918620767204615E-6, -3.1765725912106905E-6, -1.4152695485818156E-5, -2.0634874757030977E-6, -1.0475040643547796E-5, -4.5115684709925996E-5, -4.961424689243349E-6, -2.9650608578637994E-5, -2.0775794807914583E-5, -6.072034948064603E-6, 3.0769400951189936E-5, -1.0345969021122864E-4, -4.7847305687076513E-5, 9.37313764056277E-5, 8.67065014705853E-5, 8.986150044047042E-5, -3.7634003104553004E-5, -3.454951401454798E-5, -1.044138464939004E-4, -4.5483235867865845E-5, 1.5578649457531624E-6, 1.2365248380774574E-8, -5.692423599157337E-8, 3.84389468868207E-5, 1.7010821542484678E-6, 
      };
        
     
        System.out.println("the # of ECI are: " + eci.length);
        
        ce.setECI(eci, activeGroups, false);
        return ce;
    }
    public static void checkNanoparticleRelaxation_driver(double relaxRatio, double cellRatio, String poscarDir, String contcarDir, String structListDir, String structList){
        
        
        StructureListFile enerIn = new StructureListFile(structListDir + structList);

        for (int entryNum = 0; entryNum < enerIn.numEntries(); entryNum++) {
            String fileName = enerIn.getFileName(entryNum);
            
            //POSCAR infile = new POSCAR(STRUCT_DIR + fileName, true);

            //Structure structure = new Structure(infile);

            /*
            int numPt = structure.numDefiningSitesWithElement(Element.platinum);
            int numNi = structure.numDefiningSitesWithElement(Element.nickel);
            int numMo = structure.numDefiningSitesWithElement(Element.molybdenum);

            int numH = structure.numDefiningSitesWithElement(Element.hydrogen);
            */
            
            //simpleRelaxationCleanNP_new(poscarDir, contcarDir, fileName);
            
            simpleRelaxationOHDecorated_new(poscarDir, contcarDir, fileName, relaxRatio, cellRatio);

        }
        
    }
    
    
    
    public static boolean simpleRelaxationCleanNP_new(String poscarDir, String contcarDir, String fileName) {
        
        
        String contcarPath = contcarDir + "/" + fileName.substring(0, fileName.length()-5) + "/" + "CONTCAR";
        
        if (! new File(contcarPath).exists()) {
          Status.detail("No CONTCAR found for " + fileName);
          return true;
        }
        POSCAR poscar = new POSCAR(poscarDir + "/" + fileName);
        POSCAR contcar = new POSCAR(contcarPath);
        poscar.setVectorPeriodicity(0, false);
        poscar.setVectorPeriodicity(1, false);
        poscar.setVectorPeriodicity(2, false);
        contcar.setVectorPeriodicity(0, false);
        contcar.setVectorPeriodicity(1, false);
        contcar.setVectorPeriodicity(2, false);
        Structure poscarStruct = new Structure(poscar);
        Structure contcarStruct = new Structure(contcar);

        double oldStatus = Status.getLogLevel();
        Status.setLogLevelError();
        StructureMapper mapper = new StructureMapper(poscarStruct, contcarStruct);
        mapper.setMaxAllowedRMSError(1);
        mapper.setMaxAllowedOffsetDifferenceScale(0.75);
        boolean returnValue = mapper.mapExists();
        Status.setLogLevel(oldStatus);
        if (!returnValue) {
            System.out.println("Not Accepted: " + fileName);
        }
        else{
            System.out.println("Accepted: " + fileName);
        }
        return returnValue;
    }
    
    
    public static boolean simpleRelaxationOHDecorated_new(String poscarDir, String contcarDir, String fileName, double relaxRatio, double cellRatio) {
        
        //System.out.println("struct: " + fileName);
        
        //String contcarPath = contcarDir + "/" + fileName.substring(0, fileName.length()-5) + "/" + "CONTCAR";
        String contcarPath = contcarDir + "/" + fileName + "/" + "CONTCAR";
        
        //String poscarPath = poscarDir + "/" + fileName + "/" + "POSCAR";
        String poscarPath = poscarDir + "/" + fileName + "/" + "CONTCAR-decorates.vasp";

        
        if (! new File(contcarPath).exists()) {
          Status.detail("No CONTCAR found for " + fileName);
          return true;
        }
        //POSCAR poscar = new POSCAR(poscarDir + "/" + fileName);
        POSCAR poscar = new POSCAR(poscarPath);
        //poscar.setScaleFactor(cellRatio);
        
        POSCAR contcar = new POSCAR(contcarPath);
        poscar.setVectorPeriodicity(0, false);
        poscar.setVectorPeriodicity(1, false);
        poscar.setVectorPeriodicity(2, false);
        contcar.setVectorPeriodicity(0, false);
        contcar.setVectorPeriodicity(1, false);
        contcar.setVectorPeriodicity(2, false);
        Structure poscarStruct = new Structure(poscar);
        Structure contcarStruct = new Structure(contcar);

        int numPt = poscarStruct.numDefiningSitesWithElement(Element.platinum);
        int numNi = poscarStruct.numDefiningSitesWithElement(Element.nickel);
        int numMo = poscarStruct.numDefiningSitesWithElement(Element.molybdenum);

        int numO = poscarStruct.numDefiningSitesWithElement(Element.oxygen);
        int numH = poscarStruct.numDefiningSitesWithElement(Element.hydrogen);
        int initial_totalAtoms = poscarStruct.numDefiningSites();
        
        //System.out.println("numPt, numNi, numO, numH: " + numPt + ", " + numNi + ", " + numO + ", " + numH + ", " + poscarStruct.numDefiningSites());
        
        if(numH!=0){
            //TODO
            // remove all hydorgen atoms for both CONTCAR & POSCAR structures.
            StructureBuilder poscarBuilder = new StructureBuilder(poscarStruct); // builder used to add atoms to the "structure"
            StructureBuilder contcarBuilder = new StructureBuilder(contcarStruct); // builder used to add atoms to the "structure"
            
            //System.out.println("# of DefiningSites for poscar, contcar: " + poscarStruct.numDefiningSites() + ",    " + poscarStruct.numDefiningSites());
            int num_removedH = 0;
            
            poscarBuilder.removeAllSites();
            contcarBuilder.removeAllSites();
            
            for (int siteNum = 0; siteNum < initial_totalAtoms; siteNum++) {
                
                Structure.Site site = poscarStruct.getDefiningSite(siteNum);
                Coordinates coords = site.getCoords();
                Species spec = site.getSpecies();

                Structure.Site site_contcar = contcarStruct.getDefiningSite(siteNum);
                Coordinates coords_contcar = site_contcar.getCoords();
                
                if(spec!=Species.hydrogen){
                    poscarBuilder.addSite(coords, spec);
                    contcarBuilder.addSite(coords_contcar, spec);
                }
                else{
                    num_removedH++;
                }
            }
            
            poscarStruct = new Structure(poscarBuilder);
            contcarStruct = new Structure(contcarBuilder);
            
            POSCAR outfilePOS= new POSCAR(poscarStruct);
            outfilePOS.writeFile(poscarDir + "/" + fileName  + "/POSCAR-H-Removed.vasp");
            
            POSCAR outfileCONT= new POSCAR(contcarStruct);
            outfileCONT.writeFile(contcarDir + "/" + fileName  + "/CONTCAR-H-Removed.vasp");

            //System.out.println("after using builder class to remove H atoms: ");
            //System.out.println("# of DefiningSites for poscar, contcar: " + poscarStruct.numDefiningSites() + ",    " + poscarStruct.numDefiningSites());

            /*
            for (int siteNum = 0; siteNum < initial_totalAtoms; siteNum++) {
                
                Structure.Site site = poscarStruct.getDefiningSite(siteNum);

                if(site.getSpecies()==Species.hydrogen){ // if the surface site is sulfur, which stands for Pt-OH, then replace the S-atom with Pt-OH
                    Structure.Site siteCONTCAR = contcarStruct.getDefiningSite(siteNum);

                    poscarStruct.removeDefiningSite(siteNum);
                    contcarStruct.removeDefiningSite(siteNum);
                    siteNum--;
                    
                    initial_totalAtoms--;
                    
                    
                    //site.setSpecies(Species.vacancy);
                    //siteCONTCAR.setSpecies(Species.vacancy);
                    
                    
                    num_removedH++;
                }
            }
            */
            
            
            
            if(num_removedH!=numH){
                System.out.println("ERROR: the numH != num_removedH");
                return false;
            }
            
        }
        
        //System.out.println("struct: " + fileName + ", numH after removing: " + poscarStruct.numDefiningSitesWithElement(Element.hydrogen) + ", " + poscarStruct.numDefiningSites());
        
        double oldStatus = Status.getLogLevel();
        Status.setLogLevelError();
        //StructureMapper mapper = new StructureMapper(poscarStruct, contcarStruct);
        StructureMapper mapper = new StructureMapper(contcarStruct, poscarStruct);
        mapper.setMaxAllowedRMSError(1);
        mapper.setMaxAllowedOffsetDifferenceScale(relaxRatio);
        //System.out.println("\n      m_MaxAllowedOffsetDifference=" + mapper.getMaxAllowedOffsetDifference());
        boolean returnValue = mapper.mapExists();
        Status.setLogLevel(oldStatus);
        if (!returnValue) {
            System.out.println("Not Accepted: " + fileName);
        }
        else{
            //System.out.println("Accepted: " + fileName);
        }
        return returnValue;
    }
    
    
    
    
    public static boolean simpleRelaxation(String fileName) {
        
        String contcarPath = CONTCAR_DIR + "/" + "CONTCAR-" + fileName;
        
        if (! new File(contcarPath).exists()) {
          Status.detail("No CONTCAR found for " + fileName);
          return true;
        }
        POSCAR poscar = new POSCAR(STRUCT_DIR + "/" + fileName);
        POSCAR contcar = new POSCAR(contcarPath);
        poscar.setVectorPeriodicity(0, false);
        poscar.setVectorPeriodicity(1, false);
        poscar.setVectorPeriodicity(2, false);
        contcar.setVectorPeriodicity(0, false);
        contcar.setVectorPeriodicity(1, false);
        contcar.setVectorPeriodicity(2, false);
        Structure poscarStruct = new Structure(poscar);
        Structure contcarStruct = new Structure(contcar);

        double oldStatus = Status.getLogLevel();
        Status.setLogLevelError();
        //StructureMapper mapper = new StructureMapper(poscarStruct, contcarStruct);
        StructureMapper mapper = new StructureMapper(contcarStruct, poscarStruct);
        mapper.setMaxAllowedRMSError(1);
        mapper.setMaxAllowedOffsetDifferenceScale(0.75);
        boolean returnValue = mapper.mapExists();
        Status.setLogLevel(oldStatus);
        if (!returnValue) {
            System.out.println("Not Accepted: " + fileName);
        }
        return returnValue;
    }
    
    
    
    
    
}

