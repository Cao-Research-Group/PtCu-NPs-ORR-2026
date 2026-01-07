/**
 * 
 */
package PtCuOH;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.text.NumberFormat;
import java.util.Random;
import java.util.StringTokenizer;

import PtCuOH.NanoFitClusterExpansion_PtCuOH;
//import a_surface.DecorationGCHalfSlabFreezeManager;
//import a_surface.DecorationGCHalfSlabFreeze_OAdsEds_KMC_Manager;
//import a_surface.DecorationGCHalfSlabFreeze_OAdsEds_Manager;
//import a_surface.SiteConcentrationInnerLoopRecorder;
//import a_surface.SiteConcentration_OAdsDes_KMC_Recorder;
//import a_surface.SiteConcentration_OHAdsDes_KMC_Recorder;
//import binary.NanoFitClusterExpansion;
//import managers.NanoManagerLowMemVacancyHop_AllowGrowth_LT;
import managers_KMC_OHCurrent.NanoSiteConcentration_OHAdsDes_KMC_Recorder;
import managers_KMC_OHCurrent.NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN;
import managers_KMC_OHCurrent.NanoSiteConcentration_OHAdsDes_KMC_Recorder_allSites;
//import managers_KMC_OHCurrent.SiteConcentration_OHAdsDes_KMC_Recorder_Split;
import matsci.Element;
import matsci.Species;
import matsci.engine.monte.metropolis.Metropolis;
import matsci.io.app.log.Status;
import matsci.io.clusterexpansion.PRIM;
import matsci.io.clusterexpansion.StructureListFile;
import matsci.io.vasp.POSCAR;
import matsci.location.Coordinates;
import matsci.location.symmetry.operations.SpaceGroup;
import matsci.location.symmetry.operations.SymmetryOperation;
import matsci.structure.Structure;
import matsci.structure.decorate.DecorationTemplate;
import matsci.structure.decorate.function.ce.AbstractAppliedCE;
import matsci.structure.decorate.function.ce.ClusterExpansion;
import matsci.structure.decorate.function.ce.FastAppliedCE;
import matsci.structure.decorate.function.ce.clusters.ClusterGroup;
import matsci.structure.decorate.function.monte.DecorationCanonicalManager;
import matsci.structure.decorate.function.monte.NanoCanonicalManagerLowMemFixedShape;
import matsci.structure.decorate.function.monte.SiteConcentrationRecorder;
import matsci.structure.superstructure.SuperStructure;
//import purePt_111_OH.FitCE_purePt_111_OH;

/**
 * @author liang
 *
 */
public class PtCuOH_NanoFixedShape_Snapshot_KMC_Current {



    public static String ROOT = "ce/PtCu/";
    public static String STRUCT_DIR = ROOT + "/structures/";
    public static String CLUST_DIR = ROOT + "/clusters/";
    public static String CLUST_DIR_PTNIO = ROOT + "/clusters_PtNiO/";
    public static String GROUND_STATE_DIR = ROOT + "/groundStates/";
    public static String FIXEDSHAPE_DIR = ROOT + "/fixedShape/";
    public static String LAYERRECORDER_DIR = ROOT + "/fixedShape/layer-composition/";
    public static String SITERECORDER_DIR = ROOT + "/fixedShape/site-composition/";
    public static String SITERECORDER_DIR_NEW = ROOT + "/fixedShape/site-composition-new/";
    public static String VACANCYSITE_DIR = ROOT + "/vacancySite/";
    public static String DIFFUSIONENNI_DIR = ROOT + "/diffusionEnergyNi/";
    public static String TEMPSTRUCT_DIR = ROOT + "/tempStruct/";
    
    public static String ORR_PREDICTION_DIR = ROOT + "/ORR-prediction/";
    
    public static String OXYGENADSORPTION_DIR = ROOT + "/oxygenAdsorption/";
    
    public static String KMC_DIR = ROOT + "/dummy_diffusion/";
    
    
    public static NumberFormat numberFormat = NumberFormat.getNumberInstance();
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub

        
        String[] abc = new String[args.length];
        for(int i = 0; i < args.length & i < 100; i++){
               abc[i]=args[i];
        }
        
        //Status.setLogLevelBasic();
        Status.setLogLevelDetail(); //print out the generated expectation matrix
        
        
        
        //TODO fit CE immediately
        //fit CE for PtCuSVac (particles with *OH)
//        ClusterExpansion ce = NanoFitClusterExpansion_PtCuOH.buildCE_PtCuOH();
//        ce = NanoFitClusterExpansion_PtCuOH.fitCE(ce,"energyList-PREC_Accu-PtNiCuSVac-20250912_test");
        
        
        
        
        //TODO Pt-Ni-Cu-OH-Vac CE
        //ClusterExpansion ce = NanoFitClusterExpansion_PtNiCuOH.getPreFittedCE20240912_PRIM_PtNiCuOH_PtNi111_761_maxOH056ML(); 
        //ClusterExpansion ce = NanoFitClusterExpansion_PtNiCuOH.getPreFittedCE20241002_PRIM_PtNiCuOH_PtNi111_818_maxOH056ML(); 
        ClusterExpansion ce = NanoFitClusterExpansion_PtCuOH.getPreFittedCE20250428(); 

        
        
        //TODO nanoparticles
        int cellSize = 9;
        //cellSize = 21;
        
        double muOH_Inner_KMC = 0.0;
        
        double numKMCIterationsPerSite = 30;
        
        double equilibraRatio_KMC = 0.5;//?TODO
        
        Species OH_Spec = Species.sulfur;
        
        String fileDir_outer = "/483-cellSize=9-KMC-vary-compPt/";
        //fileDir_outer = "/885-cellSize=11-KMC-vary-compPt/";

        
        //2249-atom, Pt-Ni-Cu
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/2249-cellSize=15-KMC-vary-compPt/", "LBL.random.cellSize=15-Pt1349Cu450Ni449", 1, 1, 1, 0, 10, 20, 15, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
       //2249-atom, Pt-Ni
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/2249-cellSize=15-KMC-vary-compPt/", "LBL.random.cellSize=15-oct-1-truncated-Pt1349Ni900", 1, 1, 1, 0, 10, 20, 15, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
       //2249-atom, Pt-Cu
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/2249-cellSize=15-KMC-vary-compPt/", "LBL.random.cellSize=15-dummy-Pt1349Cu900", 1, 1, 1, 0, 10, 20, 11, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
  
        
        //4573-atom, Pt-Ni
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/4573-cellSize=19-KMC-vary-compPt/", "LBL.random.cellSize=19-oct-1-truncated-Pt2744Ni1829", 1, 1, 1, 0, 10, 20, 19, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
  

        //6175-atom, Pt-Ni-Cu
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/6175-cellSize=21-KMC-vary-comp-20240923/", "cellsize=21-oct-2.2.6", 1, 1, 1, 0, 10, 20, 21, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
      //6175-atom, Pt-Ni
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/6175-cellSize=21-KMC-vary-comp-20240923/", "cellsize=21-oct-0.4.6", 1, 1, 1, 0, 10, 20, 21, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
          //6175-atom, Pt-Cu
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/6175-cellSize=21-KMC-vary-comp-20240923/", "cellsize=21-oct-4.0.6", 1, 1, 1, 0, 10, 20, 21, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
     
        //muOH_Inner_KMC is a variable
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, abc[0], abc[1], abc[2], Integer.parseInt(abc[3]), Integer.parseInt(abc[4]), 1, 0, 10, 20, Integer.parseInt(abc[5]), Double.parseDouble(abc[6]), Double.parseDouble(abc[7]), equilibraRatio_KMC, OH_Spec) ;

        
        //benchmark
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/purePt-octahedral-NP-vary-size/", "cellSize=8-oct-1-truncated-Pt225-singleS", 1, 1, 1, 0, 10, 20, 8, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
        
        // groundstate.184.19.22-1 gs,
        //muOH_Inner_KMC = 0.0 ;
        //runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-SplitGCN", "/purePt-octahedral-NP-vary-size/", "groundstate.184.19.22-1", 1, 1, 1, 0, 10, 20, 8, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
        
        //before 2022-05-20
        //re-generated jar file on 2024-11-13 for Xiao Qi
        runGCMonteCarlo_KMC_Current_driver_eachSite(ce, abc[0], abc[1], abc[2], Integer.parseInt(abc[3]), Integer.parseInt(abc[4]), 1, 0, 10, 20, Integer.parseInt(abc[5]), muOH_Inner_KMC, Integer.parseInt(abc[6]), equilibraRatio_KMC, OH_Spec) ;
      //runGCMonteCarlo_KMC_Current_driver_eachSite(ce, "CE=20250403_Pt29Ni7_noPtCugs_Pt_skin", "//", "LBL.random.cellSize=21-GC-thermo-muPt=0.5", 2, 2, 1, 0, 10, 20, 21, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
        // runGCMonteCarlo_KMC_Current_driver_splitGCN(ce, "CE=20201101_355_Pt111=1.035-peakS=0.118-eachSite", "3275-cellSize=17-KMC-vary-compPt", "LBL.random.cellSize=17-oct-1-truncated-Pt2784Cu491", 1, 1, 1, 0, 10, 20, 17, muOH_Inner_KMC, 30, equilibraRatio_KMC, OH_Spec) ;
       // String ceVersion, String fileDir_outer, String fileName, int paraRunStartIndex, int paraRunEndIndex, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterationsPerSite, double equilibraRatio_KMC, Species OH_Spec)
        //2022-05-20, muOH as a variable
        //runGCMonteCarlo_KMC_Current_driver_eachSite(ce, abc[0], abc[1], abc[2], Integer.parseInt(abc[3]), Integer.parseInt(abc[4]), 1, 0, 10, 20, Integer.parseInt(abc[5]), Double.parseDouble(abc[6]), Integer.parseInt(abc[7]), equilibraRatio_KMC, OH_Spec) ;
          


    }
    
    
    
    
    public static void runGCMonteCarlo_KMC_Current_driver_splitGCN(ClusterExpansion ce, String ceVersion, String fileDir_outer, String fileName, int paraRunStartIndex, int paraRunEndIndex, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterationsPerSite, double equilibraRatio_KMC, Species OH_Spec) {
        
        try{

            String predictedSA_path = ORR_PREDICTION_DIR + "/" + ceVersion + "/";
            
            File dirPredictedSA  = new  File(predictedSA_path);
            if ( ! (dirPredictedSA.exists())  &&   ! (dirPredictedSA.isDirectory())) {
                    boolean  creadok  =  dirPredictedSA.mkdirs();
            }
            
            String sub_fileName = fileName.substring(20);
            
            //FileWriter fw0 = new FileWriter(ORR_PREDICTION_DIR + fileDir_outer + "/" + ceVersion + "-KMC_SA_" + fileName + "-size=" + cellSize + ".txt", false);
            //FileWriter fw0 = new FileWriter(ORR_PREDICTION_DIR + fileDir_outer + "/" + ceVersion + "-IteN=" + numKMCIterations + "-" + sub_fileName + ".txt", false);
            //FileWriter fw0 = new FileWriter(predictedSA_path + "/" + ceVersion + "-IteN=" + numKMCIterations + "-" + sub_fileName + ".txt", false);
            FileWriter fw0 = new FileWriter(predictedSA_path + "/IteN=" + numKMCIterationsPerSite + "-" + fileName + "-s=" + paraRunStartIndex + "-e=" + paraRunEndIndex +  ".txt", false);

            BufferedWriter bw0 = new BufferedWriter(fw0); 
          
            //bw0.write("structName   temp numPt numNi numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgFaceOHOccupancy avgEdgeOHOccupancy avgEads avgFaceEads avgEdgeEads avgEads2 stdDev " + "\r\n");
            bw0.write("structName   temp numPt numCu numPtSurface numPt_6.67<=GCN numPt_EdgeVertex ratio_6.67<=GCN ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgGCNFaceOHOccupancy avgEdgeOHOccupancy avgEads avgGCNFaceEads avgEdgeEads avgEads2 stdDev " + "\r\n");
            bw0.flush();
            
            //StructureListFile enerIn = new StructureListFile(ORR_PREDICTION_DIR + fileDir_outer + structList);

            //for (int entryNum = 0; entryNum < paraRunNum; entryNum++) {
            for (int entryNum = paraRunStartIndex; entryNum <= paraRunEndIndex; entryNum++) {
    
                System.out.println("\n\n**************************************************************");
                System.out.println("path and structure name used for KMC SAMA: ");
                System.out.println("fileDir_outer: " + fileDir_outer);
                System.out.println("fileName: " + fileName + "-" + entryNum);
                System.out.println("**************************************************************\n\n");

                //String structName = enerIn.getFileName(entryNum);
                //int size1 = (int)enerIn.getEnergy(entryNum);
                
                String structName = "snapShot-afterKMC.vasp";
                int size1 = cellSize;
                
                //POSCAR infile = new POSCAR(STRUCT_DIR + fileName, true);

                //Structure structure = new Structure(infile);

                /*
                int numPt = structure.numDefiningSitesWithElement(Element.platinum);
                int numNi = structure.numDefiningSitesWithElement(Element.nickel);
                int numMo = structure.numDefiningSitesWithElement(Element.molybdenum);

                int numH = structure.numDefiningSitesWithElement(Element.hydrogen);
                */
                
                //simpleRelaxationCleanNP_new(poscarDir, contcarDir, fileName);
                
                /*
                int [] returnArr2 = SplitNPLayerbyLayer1stLayer(ORR_PREDICTION_DIR + fileDir_outer + "/" + fileName + "-" + entryNum + "/", structName) ;
                double ratio_surface = ((double)returnArr2[1]) / returnArr2[0];
                */
                
                //double [] returnArr = runGCMonteCarlo_KMC_Current_test(ce, "/" + fileDir_outer + "/", fileName + "-" + entryNum + "/", structName, numPasses2, startTemp2, endTemp2, tempIncre, cellSize, muOH_Inner_KMC, numKMCIterationsPerSite, equilibraRatio_KMC, OH_Spec) ;
                double [] returnArr = runGCMonteCarlo_KMC_Current_splitGCN(ce, "/" + fileDir_outer + "/", fileName + "-" + entryNum + "/", structName, numPasses2, startTemp2, endTemp2, tempIncre, cellSize, muOH_Inner_KMC, numKMCIterationsPerSite, equilibraRatio_KMC, OH_Spec) ;

                //double [] returnArr = calculateAvgCur_listActivityEachSite(ce, fileDir, structName, size1);
                //double [] returnArr = calculateAvgCur(ce, fileDir, structName, size1);

                //System.out.println("Temp: " + temp + "\tConcentration: " + concentration + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                Status.detail("structName: " + structName + "\tKMC_avgCur: " + returnArr[7] + "\tOH_coverage: " + returnArr[10]);

                //Status.detail("Temp: " + temp + "\ttrigRatio: " + recorder.getTriggerRatio() + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                //bw0.write(structName + "    " + returnArr2[0] + "  " + returnArr2[1] + "    " + ratio_surface + "  " + returnArr[0] + "  " + returnArr[1] + "  " + returnArr[2] + "  " + returnArr[3] + "  " + returnArr[4] + "  " + returnArr[5] + "    " + returnArr[6] + "  " + returnArr[7] + "  " + returnArr[8] + "  " + returnArr[9] + "  " + returnArr[10] + "  " + returnArr[11] + "  " + returnArr[12] + "    " + returnArr[13] + "\r\n");
                bw0.write(structName +  "  " + returnArr[0] + "  " + returnArr[1] + "  " + returnArr[3] + "  " + returnArr[19] + "  " + returnArr[20] + "  " + returnArr[21] + "  " + returnArr[22] + "  " + returnArr[23] +
                
                        "  " + returnArr[4] + "  " + returnArr[5] + "  " + returnArr[6] + "    " + returnArr[7] + "  " + returnArr[8] + "  " + returnArr[9] + "  " + returnArr[10] + "  " + returnArr[11] + "  " + returnArr[12] + "  " + returnArr[13] + "    " + returnArr[14] + "    " + returnArr[15] + "    " + returnArr[16] + "    " + returnArr[17] + "    " + returnArr[18] + "\r\n");

                //bw0.write("structName   temp numPt numNi numCu numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgEads avgEads2 stdDev " + "\r\n");

                
                bw0.flush();

            }
            
            
            
            bw0.flush();
            bw0.close();
            fw0.close();
            

        }
            
        catch (IOException e) {
            e.printStackTrace();
        }

    }
    
    
    
    
    public static void runGCMonteCarlo_KMC_Current_driver_eachSite(ClusterExpansion ce, String ceVersion, String fileDir_outer, String fileName, int paraRunStartIndex, int paraRunEndIndex, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterationsPerSite, double equilibraRatio_KMC, Species OH_Spec) {
        
        try{

            String predictedSA_path = ORR_PREDICTION_DIR + "/" + ceVersion + "/";
            
            File dirPredictedSA  = new  File(predictedSA_path);
            if ( ! (dirPredictedSA.exists())  &&   ! (dirPredictedSA.isDirectory())) {
                    boolean  creadok  =  dirPredictedSA.mkdirs();
            }
            
            String sub_fileName = fileName.substring(20);
            
            //FileWriter fw0 = new FileWriter(ORR_PREDICTION_DIR + fileDir_outer + "/" + ceVersion + "-KMC_SA_" + fileName + "-size=" + cellSize + ".txt", false);
            //FileWriter fw0 = new FileWriter(ORR_PREDICTION_DIR + fileDir_outer + "/" + ceVersion + "-IteN=" + numKMCIterations + "-" + sub_fileName + ".txt", false);
            //FileWriter fw0 = new FileWriter(predictedSA_path + "/" + ceVersion + "-IteN=" + numKMCIterations + "-" + sub_fileName + ".txt", false);
            //FileWriter fw0 = new FileWriter(predictedSA_path + "/IteN=" + numKMCIterationsPerSite + "-" + sub_fileName + "-s=" + paraRunStartIndex + "-e=" + paraRunEndIndex +  ".txt", false);
            FileWriter fw0 = new FileWriter(predictedSA_path + "/IteN=" + numKMCIterationsPerSite + "-" + fileName + "-muOH=" + muOH_Inner_KMC + "-s=" + paraRunStartIndex + "-e=" + paraRunEndIndex +  ".txt", false);

            BufferedWriter bw0 = new BufferedWriter(fw0); 
          
            //bw0.write("structName   temp numPt numNi numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgFaceOHOccupancy avgEdgeOHOccupancy avgEads avgFaceEads avgEdgeEads avgEads2 stdDev " + "\r\n");
            bw0.write("structName   temp numPt numNi numCu numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgFaceOHOccupancy avgEdgeOHOccupancy avgEads avgFaceEads avgEdgeEads avgEads2 stdDev " + "\r\n");
            bw0.flush();
            
            //StructureListFile enerIn = new StructureListFile(ORR_PREDICTION_DIR + fileDir_outer + structList);

            //for (int entryNum = 0; entryNum < paraRunNum; entryNum++) {
            for (int entryNum = paraRunStartIndex; entryNum <= paraRunEndIndex; entryNum++) {
    
                System.out.println("\n\n**************************************************************");
                System.out.println("path and structure name used for KMC SAMA: ");
                System.out.println("fileDir_outer: " + fileDir_outer);
                System.out.println("fileName: " + fileName + "-" + entryNum);
                System.out.println("**************************************************************\n\n");

                //String structName = enerIn.getFileName(entryNum);
                //int size1 = (int)enerIn.getEnergy(entryNum);
                
                String structName = "snapShot-afterKMC.vasp";
                int size1 = cellSize;
                
                //POSCAR infile = new POSCAR(STRUCT_DIR + fileName, true);

                //Structure structure = new Structure(infile);

                /*
                int numPt = structure.numDefiningSitesWithElement(Element.platinum);
                int numNi = structure.numDefiningSitesWithElement(Element.nickel);
                int numMo = structure.numDefiningSitesWithElement(Element.molybdenum);

                int numH = structure.numDefiningSitesWithElement(Element.hydrogen);
                */
                
                //simpleRelaxationCleanNP_new(poscarDir, contcarDir, fileName);
                
                /*
                int [] returnArr2 = SplitNPLayerbyLayer1stLayer(ORR_PREDICTION_DIR + fileDir_outer + "/" + fileName + "-" + entryNum + "/", structName) ;
                double ratio_surface = ((double)returnArr2[1]) / returnArr2[0];
                */
                
                //double [] returnArr = runGCMonteCarlo_KMC_Current_test(ce, "/" + fileDir_outer + "/", fileName + "-" + entryNum + "/", structName, numPasses2, startTemp2, endTemp2, tempIncre, cellSize, muOH_Inner_KMC, numKMCIterationsPerSite, equilibraRatio_KMC, OH_Spec) ;
                double [] returnArr = runGCMonteCarlo_KMC_Current_eachSite(ce, "/" + fileDir_outer + "/", fileName + "-" + entryNum + "/", structName, numPasses2, startTemp2, endTemp2, tempIncre, cellSize, muOH_Inner_KMC, numKMCIterationsPerSite, equilibraRatio_KMC, OH_Spec) ;

                //double [] returnArr = calculateAvgCur_listActivityEachSite(ce, fileDir, structName, size1);
                //double [] returnArr = calculateAvgCur(ce, fileDir, structName, size1);

                //System.out.println("Temp: " + temp + "\tConcentration: " + concentration + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                Status.detail("structName: " + structName + "\tKMC_avgCur: " + returnArr[7] + "\tOH_coverage: " + returnArr[10]);

                //Status.detail("Temp: " + temp + "\ttrigRatio: " + recorder.getTriggerRatio() + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                //bw0.write(structName + "    " + returnArr2[0] + "  " + returnArr2[1] + "    " + ratio_surface + "  " + returnArr[0] + "  " + returnArr[1] + "  " + returnArr[2] + "  " + returnArr[3] + "  " + returnArr[4] + "  " + returnArr[5] + "    " + returnArr[6] + "  " + returnArr[7] + "  " + returnArr[8] + "  " + returnArr[9] + "  " + returnArr[10] + "  " + returnArr[11] + "  " + returnArr[12] + "    " + returnArr[13] + "\r\n");
                bw0.write(structName +  "  " + returnArr[0] + "  " + returnArr[1] + "  " + returnArr[2] + "  " + returnArr[3] + "  " + returnArr[19] + "  " + returnArr[20] + "  " + returnArr[21] + "  " + returnArr[22] + "  " + returnArr[23] +
                
                        "  " + returnArr[4] + "  " + returnArr[5] + "  " + returnArr[6] + "    " + returnArr[7] + "  " + returnArr[8] + "  " + returnArr[9] + "  " + returnArr[10] + "  " + returnArr[11] + "  " + returnArr[12] + "  " + returnArr[13] + "    " + returnArr[14] + "    " + returnArr[15] + "    " + returnArr[16] + "    " + returnArr[17] + "    " + returnArr[18] + "\r\n");

                //bw0.write("structName   temp numPt numNi numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgEads avgEads2 stdDev " + "\r\n");

                
                bw0.flush();

            }
            
            
            
            bw0.flush();
            bw0.close();
            fw0.close();
            

        }
            
        catch (IOException e) {
            e.printStackTrace();
        }

    }
    
    
    
    
    
    
    
    
    
    
    public static double[] runGCMonteCarlo_KMC_Current_splitGCN(ClusterExpansion ce, String fileDir_outer, String fileDir_inner, String initialStruct, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterations, double equilibraRatio_KMC, Species OH_Spec) {
        
        double[] returnArr = new double[24];

        
        double kb = 8.6173423E-5;
        //double stepsFor4 = 10;
        
        //String fileDir = ORR_PREDICTION_DIR + "_variedOCoverage/" + initialStruct + "/";
        String fileDir = ORR_PREDICTION_DIR + fileDir_outer + fileDir_inner + "/";

        
        File dirFile  = new  File(fileDir);
        if ( ! (dirFile.exists())  &&   ! (dirFile.isDirectory())) {
                boolean  creadok  =  dirFile.mkdirs();
        }
        
        
        String fileDirInner = fileDir + "muOH_KMC=" + muOH_Inner_KMC + "-numKMCItes=" + numKMCIterations + "-KMCEquRatio=" + equilibraRatio_KMC + "/";
        File dirFileInner  = new  File(fileDirInner);
        if ( ! (dirFileInner.exists())  &&   ! (dirFileInner.isDirectory())) {
                boolean  creadok  =  dirFileInner.mkdirs();
        }
        
        ClusterGroup[] activeGroups = ce.getAllGroups(true);
        
        //double CutoffPerSite = 0;
        double CutoffPerSite = 1E-5;
        System.out.println("Using Cutoff per site of " + CutoffPerSite + " to trim insignificant groups...");
        activeGroups = ce.getSignificantGroups(activeGroups, ce.numSigmaSites() * CutoffPerSite);
        System.out.println(activeGroups.length + " significant groups remaining.\n\n");
        
        double startTemp = startTemp2;
        double endTemp = endTemp2;
        double tempIncrement = tempIncre;
        
        int numPasses = numPasses2;
        System.out.println("the number of numPasses : " + numPasses);

        /*
         * chemical potential of all species
         */
        
        double chemPotPt = 0.0;
        double chemPotNi = 0.0; 
       
        
        int[][] superToDirect = new int[][] {
                {-cellSize, cellSize, cellSize},
                {cellSize, -cellSize, cellSize},
                {cellSize, cellSize, -cellSize}
        };
        
        AbstractAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        //AbstractAppliedCE appliedCE = new MakeMonteAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        System.out.println("Supercell has " + appliedCE.numPrimCells() + " primitive cells.");
        appliedCE.activateGroups(activeGroups);
        System.out.println("the number of sites is: " + appliedCE.numSigmaSites());
        
        SuperStructure groundState2 = appliedCE.getSuperStructure();
        
        POSCAR outfile2= new POSCAR(groundState2);
        outfile2.writeFile(fileDirInner + "groundState_0.vasp");
        //outfile2.writeVICSOutFile(fileDirInner + "groundState_0.out");
        System.out.println("write out the initial structure before decoration!!!!!!");
        
        Metropolis metropolis = new Metropolis();
        metropolis.setOptimistic(true);
        
        //appliedCE.decorateRandomly(new int[] {0, 1, 2}, new int[] {numNi, numPt, numVac});

        
        //POSCAR infile = new POSCAR(fileDir + initialStruct + ".vasp", true);
        POSCAR infile = new POSCAR(fileDir + initialStruct, true);
        Structure structure = new Structure(infile);
        appliedCE.decorateFromStructure(structure);
        
        SuperStructure groundState3 = appliedCE.getSuperStructure();
        
        POSCAR outfile3= new POSCAR(groundState3);
        outfile3.writeFile(fileDirInner + "groundState_1.vasp");
        //outfile3.writeVICSOutFile(fileDirInner + "groundState_1.out");
        System.out.println("write out the structure after decoration!!!!!!");
        


        
        int numPt = structure.numDefiningSitesWithSpecies(Species.platinum);
        //int numNi = structure.numDefiningSitesWithSpecies(Species.nickel);
        int numCu = structure.numDefiningSitesWithSpecies(Species.copper);
        int numVac = appliedCE.numSigmaSites() - numCu - numPt;
        int numAtoms = numPt  + numCu;
        
        //DecorationCanonicalManager manager = new DecorationCanonicalManager(appliedCE); // Use this line if you want a canonical run

        //DecorationCrossCanonicalManager manager = new DecorationCrossCanonicalManager(appliedCE); // Use this line if you want a canonical run
        NanoCanonicalManagerLowMemFixedShape manager = new NanoCanonicalManagerLowMemFixedShape(appliedCE);
        
        metropolis.setTemp(startTemp * kb); // Whatever you want
        
        //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
        metropolis.setNumIterations(1);
        Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps...");
        metropolis.runBasic(manager);

        //System.out.println("Initial FBeta: " + totalFBeta);
        //System.out.println("Running for real...");
        
        
        if ((endTemp - startTemp) / tempIncrement < 0) {
            tempIncrement *= -1;
        }
        int numSteps = (int) Math.floor((endTemp - startTemp) / tempIncrement);
        System.out.println("the # of step of Monte Carlo simulation is: " + numSteps + "\n\n");
        
        Structure groundState = null;
        
        
        try{
            
            //FileWriter fw = new FileWriter(fileDirInner + "compositionProfile.txt", false);
            FileWriter fw = new FileWriter(fileDir + "/muOH_KMC=" + muOH_Inner_KMC + "-nKMCItes=" + numKMCIterations + "-avgCurrent.txt", false);
            BufferedWriter bw = new BufferedWriter(fw); 
            
            bw.write("temp numPt numCu numPtSurface numPt_6.67<=GCN numPt_EdgeVertex ratio_6.67<=GCN ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgGCNFaceOHOccupancy avgEdgeOHOccupancy avgEads avgGCNFaceEads avgEdgeEads avgEads2 stdDev " + "\r\n");
            bw.flush();
            
            for (int stepNum = 0; stepNum <= numSteps; stepNum++) {
                
                double temp = startTemp + stepNum * tempIncrement;
                metropolis.setTemp(kb * temp);
                
                /*
                 * If want to run grand canonical ensemble, in which the compositions are allowed to vary.
                 * thus, we think we can find some interesting phenomena.
                 * I will create a mirror symmetrical grand cannonical ensemble manager, named as DecorationMSymGCManager.
                 */
                //DecorationMSymGCManager manager = new DecorationMSymGCManager(appliedCE, mirrorOp);
                
                //SiteConcentrationRecorder recorder = new SiteConcentrationRecorder(manager, appliedCE); 
                
                //TODO
                NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 
                //NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN_KMCTime recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder_SplitGCN_KMCTime(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 
                //SiteConcentration_OHAdsDes_KMC_Recorder_Split recorder = new SiteConcentration_OHAdsDes_KMC_Recorder_Split(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 

                // just for symmetrical Monte Carlo simulations
               // Equilibration            
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
                //structMetropolis.runBasic(manager);

                // Run for real
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());         
                metropolis.setNumIterations(1);         
                Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps!!!");
                metropolis.run(manager, recorder);
                        
                Status.detail("Temperature: " + metropolis.getTemp());
                Status.detail("Trigger ratio: " + recorder.getTriggerRatio());
                Status.detail("Average energy: " + recorder.getAverageValue());
                

             
                /*
                 * print out the ground state for every temperature step
                 * 
                 */
                SuperStructure groundState4 = appliedCE.getSuperStructure();            
                POSCAR outfile4= new POSCAR(groundState4);                        
                outfile4.writeFile(fileDirInner + "snapShot." + stepNum + "-temp=" + temp + ".vasp");
                //outfile4.writeVICSOutFile(fileDirInner + "snapShot." + stepNum  + "-temp=" + temp + ".out"); 

                System.out.println("\nOuter metropolis run to get the cleanslab under ChemPot(Pt-Ni): ");
                System.out.println("Temperature: " + metropolis.getTemp());
                System.out.println("# of iterations: " + metropolis.getNumIterations());
                System.out.println("Trigger ratio: " + recorder.getTriggerRatio());
                System.out.println("Average energy: " + recorder.getAverageValue());
                
                System.out.println("Average current: " + recorder.getAverageCurrent());
                System.out.println("Average current from occupied sites: " + recorder.getAverageCurrentOccupied());
                System.out.println("Average current from unoccupied sites: " + recorder.getAverageCurrentUnOccupied());
                System.out.println("Average Ads current: " + recorder.getAverageAdsCurrent());
                System.out.println("Max current: " + recorder.getMaxAvgCurrent());
                System.out.println("Average OH Occupancy: " + recorder.getAverageOHOccupancy());
                System.out.println("Average OH Binding Energy: " + recorder.getAverageOHBindingEnergy());
                System.out.println("Average OH Binding Energy Sq: " + recorder.getAverageOHBindingEnergySq());
                double stdDev = Math.sqrt(recorder.getAverageOHBindingEnergySq() - (recorder.getAverageOHBindingEnergy() *  recorder.getAverageOHBindingEnergy()));
                System.out.println("Average OH Binding Energy Std Dev: " + stdDev);
                
                System.out.println("Average Face(GCN=7.5) OH Occupancy: " + recorder.getAverageFaceOHOccupancy());
                System.out.println("Average Face(GCN=7.5) OH Binding Energy: " + recorder.getAverageFaceOHBindingEnergy());
                
                System.out.println("Average Edge OH Occupancy: " + recorder.getAverageEdgeOHOccupancy());
                System.out.println("Average Edge OH Binding Energy: " + recorder.getAverageEdgeOHBindingEnergy());
                
                
                /*
                // the manager: DecorationGCHalfSlabFreeze_OAdsEds_KMC_Manager, which contains KMC (reject-free)
                long totalNumGeneratedAdsEvents = manager.getNumGeneratedAdsEvents();
                long totalNumGeneratedDesEvents = manager.getNumGeneratedDesEvents();

                long totalNumAcceptedAdsEvents = manager.getNumAcceptedAdsEvents();
                long totalNumAcceptedDesEvents = manager.getNumAcceptedDesEvents();

                double totalKMCTime = manager.getKMCTime();
                
                double avgAdsBarrier = manager.getAvgAdsBarrierEvents();
                double avgEdsBarrier = manager.getAvgDesBarrierEvents();

                System.out.println("totalNumGeneratedAdsEvents: " + totalNumGeneratedAdsEvents);
                System.out.println("totalNumGeneratedDesEvents: " + totalNumGeneratedDesEvents);
                System.out.println("totalNumAcceptedAdsEvents: " + totalNumAcceptedAdsEvents);
                System.out.println("totalNumAcceptedDesEvents: " + totalNumAcceptedDesEvents);
                System.out.println("avgAdsBarrier: " + avgAdsBarrier);     
                System.out.println("avgEdsBarrier: " + avgEdsBarrier);    
                System.out.println("totalKMCTime: " + totalKMCTime);    
                System.out.println("Average current: " + totalNumAcceptedDesEvents / totalKMCTime);
                */
                
                double avgE = recorder.getAverageValue();
                double avgE2 = recorder.getAverageSquaredValue();
                double sigSq = avgE2 - avgE * avgE;
                double cV = (avgE2 - avgE*avgE) / (kb * temp * temp);
                cV /= appliedCE.numSigmaSites();

                double maxCurrent = recorder.getMaxAvgCurrent();
                double avgCurrent = recorder.getAverageCurrent();
                double avgCurrentOccupied = recorder.getAverageCurrentOccupied();
                double avgCurrentUnOccupied = recorder.getAverageCurrentUnOccupied();                
                double avgAdsCurrent = recorder.getAverageAdsCurrent();
                double avgOHOccupancy = recorder.getAverageOHOccupancy();
                double avgFaceOHOccupancy = recorder.getAverageFaceOHOccupancy();
                double avgEdgeOHOccupancy = recorder.getAverageEdgeOHOccupancy();
                double avgEads =  recorder.getAverageOHBindingEnergy();
                double avgFaceEads =  recorder.getAverageFaceOHBindingEnergy();
                double avgEdgeEads =  recorder.getAverageEdgeOHBindingEnergy();
                double avgEads2 = recorder.getAverageOHBindingEnergySq();
                
                double avgKMCTime = recorder.getOutKMCTime();
                
                int totalSurfacePt = recorder.getNumOHSites();
                int totalSurfacePt_111 = recorder.getNumOHFaceSites();
                int totalSurfacePt_EdgeVertex = recorder.getNumOHEdgeSites();
                
                Structure sturctureMaxAvgCurrent = recorder.getStructureMaxTempAvgCurrent();
                
                POSCAR outfileMaxAvgCurrent= new POSCAR(sturctureMaxAvgCurrent);                                        outfileMaxAvgCurrent.writeFile(fileDirInner + "snapshot." + stepNum + "-temp=" + temp + "-maxAvgCurrent.vasp");
                
                
                System.out.println("\n\n****************************************************\n");
                System.out.println("the KMC time: " +  avgKMCTime + "\n");
                System.out.println("****************************************************\n\n");
                
                
                //System.out.println("Temp: " + temp + "\tConcentration: " + concentration + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                Status.detail("Temp: " + temp + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                
                Structure snapShot = appliedCE.getStructure(); 
              
                //System.out.println(stepNum);
                System.out.println("the num of steps is: " + stepNum + "\n");
                
                
                //String stringNum = temp + "   " + numPt + "   " + numNi + "   " + numOTem + "   " + compositionPtTemp_bulk + "   " + compositionPtTemp + "   " + latticeP_Temp_bulk + "   " + latticeP_Temp + "   " + strain_Temp_bulk + "   " + strain_Temp + "   " + shiftAdsEnergy + "  " + avgE + "   " + avgE2 + "   " + Cv_Slab + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOOccupancy + "    " + avgEads + "    " + avgEads2 + "   " + stdDev + "   " + compositionNiSplitTemp[0] + "   " + compositionNiSplitTemp[1] + "   " + compositionNiSplitTemp[2] + "   " + compositionNiSplitTemp[3] + "   " + NiComp_567_Freeze[0] + "   " + NiComp_567_Freeze[1] + "   " + NiComp_567_Freeze[2] + "\r\n";
                
                String stringNum = temp + "   " + numPt  + "   " + numCu + "   " + totalSurfacePt + "   " + totalSurfacePt_111 + "   " + totalSurfacePt_EdgeVertex + "   " + (((double)totalSurfacePt_111) / totalSurfacePt) + "   " + (((double)totalSurfacePt) / numPt) + "   " + avgE + "   " + avgE2 + "   " + cV + "   " + muOH_Inner_KMC + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOHOccupancy + "    " + avgFaceOHOccupancy + "    " + avgEdgeOHOccupancy + "    " + avgEads + "    " + avgFaceEads + "    " + avgEdgeEads + "    " + avgEads2 + "   " + stdDev + "\r\n";
                
                System.out.println("Temp and Cv are: " + stringNum);
                bw.write(stringNum);    
                bw.flush();

                
                returnArr[0] = temp;
                returnArr[1] = numPt; 
               // returnArr[2] = numNi; 
                
                returnArr[3] = numCu; 

                returnArr[4] = avgE;
                returnArr[5] = avgE2;
                returnArr[6] = cV;
                returnArr[7] = muOH_Inner_KMC;
                returnArr[8] = avgCurrent;
                returnArr[9] = avgAdsCurrent; 
                returnArr[10] = maxCurrent ;
                returnArr[11] = avgOHOccupancy;
                returnArr[12] = avgFaceOHOccupancy;
                returnArr[13] = avgEdgeOHOccupancy;
                returnArr[14] = avgEads;
                returnArr[15] = avgFaceEads;
                returnArr[16] = avgEdgeEads;
                returnArr[17] = avgEads2;
                returnArr[18] = stdDev;

                returnArr[19] = totalSurfacePt;
                returnArr[20] = totalSurfacePt_111;
                returnArr[22] = totalSurfacePt_EdgeVertex;
                returnArr[22] = ((double)totalSurfacePt_111) / totalSurfacePt;
                returnArr[23] = ((double)totalSurfacePt) / numPt;

            } 
            
            bw.flush();
            bw.close();
            fw.close();

        }
            
        catch (IOException e) {
            e.printStackTrace();
        }
            

        
        groundState = appliedCE.getSuperStructure().getCompactStructure();
        int numPtF = groundState.numDefiningSitesWithSpecies(Species.platinum);
       // int numNiF = groundState.numDefiningSitesWithSpecies(Species.nickel);
        int numCuF = groundState.numDefiningSitesWithSpecies(Species.copper);
        POSCAR outfile5= new POSCAR(groundState);
        outfile5.setDescription("fixedShape");
        outfile5.writeFile(fileDirInner  + "groundState." + numPtF  + "."+ numCuF + ".vasp");
        //outfile5.writeVICSOutFile(fileDirInner  + "groundState." + numPtF +"."+ numNiF + "."+ numCuF + ".out");

        
        return returnArr;
        
      } 
      
    
    
    
    
    
    public static double[] runGCMonteCarlo_KMC_Current_eachSite(ClusterExpansion ce, String fileDir_outer, String fileDir_inner, String initialStruct, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterations, double equilibraRatio_KMC, Species OH_Spec) {
        
        double[] returnArr = new double[24];

        
        double kb = 8.6173423E-5;
        //double stepsFor4 = 10;
        
        //String fileDir = ORR_PREDICTION_DIR + "_variedOCoverage/" + initialStruct + "/";
        String fileDir = ORR_PREDICTION_DIR + fileDir_outer + fileDir_inner + "/";

        
        File dirFile  = new  File(fileDir);
        if ( ! (dirFile.exists())  &&   ! (dirFile.isDirectory())) {
                boolean  creadok  =  dirFile.mkdirs();
        }
        
        
        String fileDirInner = fileDir + "muOH_KMC=" + muOH_Inner_KMC + "-numKMCItes=" + numKMCIterations + "-KMCEquRatio=" + equilibraRatio_KMC + "/";
        File dirFileInner  = new  File(fileDirInner);
        if ( ! (dirFileInner.exists())  &&   ! (dirFileInner.isDirectory())) {
                boolean  creadok  =  dirFileInner.mkdirs();
        }
        
        ClusterGroup[] activeGroups = ce.getAllGroups(true);
        
        //double CutoffPerSite = 0;
        double CutoffPerSite = 1E-5;
        System.out.println("Using Cutoff per site of " + CutoffPerSite + " to trim insignificant groups...");
        activeGroups = ce.getSignificantGroups(activeGroups, ce.numSigmaSites() * CutoffPerSite);
        System.out.println(activeGroups.length + " significant groups remaining.\n\n");
        
        double startTemp = startTemp2;
        double endTemp = endTemp2;
        double tempIncrement = tempIncre;
        
        int numPasses = numPasses2;
        System.out.println("the number of numPasses : " + numPasses);

        /*
         * chemical potential of all species
         */
        
        double chemPotPt = 0.0;
        double chemPotNi = 0.0; 
       
        
        int[][] superToDirect = new int[][] {
                {-cellSize, cellSize, cellSize},
                {cellSize, -cellSize, cellSize},
                {cellSize, cellSize, -cellSize}
        };
        
        AbstractAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        //AbstractAppliedCE appliedCE = new MakeMonteAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        System.out.println("Supercell has " + appliedCE.numPrimCells() + " primitive cells.");
        appliedCE.activateGroups(activeGroups);
        System.out.println("the number of sites is: " + appliedCE.numSigmaSites());
        
        SuperStructure groundState2 = appliedCE.getSuperStructure();
        
        POSCAR outfile2= new POSCAR(groundState2);
        outfile2.writeFile(fileDirInner + "groundState_0.vasp");
        //outfile2.writeVICSOutFile(fileDirInner + "groundState_0.out");
        System.out.println("write out the initial structure before decoration!!!!!!");
        
        Metropolis metropolis = new Metropolis();
        metropolis.setOptimistic(true);
        
        //appliedCE.decorateRandomly(new int[] {0, 1, 2}, new int[] {numNi, numPt, numVac});

        
        //POSCAR infile = new POSCAR(fileDir + initialStruct + ".vasp", true);
        POSCAR infile = new POSCAR(fileDir + initialStruct, true);
        Structure structure = new Structure(infile);
        appliedCE.decorateFromStructure(structure);
        
        SuperStructure groundState3 = appliedCE.getSuperStructure();
        
        POSCAR outfile3= new POSCAR(groundState3);
        outfile3.writeFile(fileDirInner + "groundState_1.vasp");
        //outfile3.writeVICSOutFile(fileDirInner + "groundState_1.out");
        System.out.println("write out the structure after decoration!!!!!!");
        


        
        int numPt = structure.numDefiningSitesWithSpecies(Species.platinum);
        int numNi = structure.numDefiningSitesWithSpecies(Species.nickel);
        int numCu = structure.numDefiningSitesWithSpecies(Species.copper);
        int numVac = appliedCE.numSigmaSites() - numCu - numNi - numPt;
        int numAtoms = numPt + numNi + numCu;
        
        //DecorationCanonicalManager manager = new DecorationCanonicalManager(appliedCE); // Use this line if you want a canonical run

        //DecorationCrossCanonicalManager manager = new DecorationCrossCanonicalManager(appliedCE); // Use this line if you want a canonical run
        NanoCanonicalManagerLowMemFixedShape manager = new NanoCanonicalManagerLowMemFixedShape(appliedCE);
        
        metropolis.setTemp(startTemp * kb); // Whatever you want
        
        //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
        metropolis.setNumIterations(1);
        Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps...");
        metropolis.runBasic(manager);

        //System.out.println("Initial FBeta: " + totalFBeta);
        //System.out.println("Running for real...");
        
        
        if ((endTemp - startTemp) / tempIncrement < 0) {
            tempIncrement *= -1;
        }
        int numSteps = (int) Math.floor((endTemp - startTemp) / tempIncrement);
        System.out.println("the # of step of Monte Carlo simulation is: " + numSteps + "\n\n");
        
        Structure groundState = null;
        
        
        try{
            
            //FileWriter fw = new FileWriter(fileDirInner + "compositionProfile.txt", false);
            FileWriter fw = new FileWriter(fileDir + "/muOH_KMC=" + muOH_Inner_KMC + "-nKMCItes=" + numKMCIterations + "-avgCurrent.txt", false);
            BufferedWriter bw = new BufferedWriter(fw); 
            
            bw.write("temp numPt numNi numCu numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgFaceOHOccupancy avgEdgeOHOccupancy avgEads avgFaceEads avgEdgeEads avgEads2 stdDev " + "\r\n");
            bw.flush();
            
            
            
            FileWriter fw_eachSite = new FileWriter(fileDir + "/muOH_KMC=" + muOH_Inner_KMC + "-nKMCItes=" + numKMCIterations + "-avgEads-eachSite.txt", false);
            BufferedWriter bw_eachSite = new BufferedWriter(fw_eachSite); 
            
            bw_eachSite.write("siteNum surfaceSiteIndex surfaceSiteSigmaIndex avgEads CN avgEads2 stdDev " + "\r\n");
            bw_eachSite.flush();
            
            for (int stepNum = 0; stepNum <= numSteps; stepNum++) {
                
                double temp = startTemp + stepNum * tempIncrement;
                metropolis.setTemp(kb * temp);
                
                /*
                 * If want to run grand canonical ensemble, in which the compositions are allowed to vary.
                 * thus, we think we can find some interesting phenomena.
                 * I will create a mirror symmetrical grand cannonical ensemble manager, named as DecorationMSymGCManager.
                 */
                //DecorationMSymGCManager manager = new DecorationMSymGCManager(appliedCE, mirrorOp);
                
                //SiteConcentrationRecorder recorder = new SiteConcentrationRecorder(manager, appliedCE); 
                
                //TODO
                //NanoSiteConcentration_OHAdsDes_KMC_Recorder recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 
                //SiteConcentration_OHAdsDes_KMC_Recorder_Split recorder = new SiteConcentration_OHAdsDes_KMC_Recorder_Split(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 
                NanoSiteConcentration_OHAdsDes_KMC_Recorder_allSites recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder_allSites(manager, appliedCE, structure, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec,fileDirInner); 

                // just for symmetrical Monte Carlo simulations
               // Equilibration            
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
                //structMetropolis.runBasic(manager);

                // Run for real
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());         
                metropolis.setNumIterations(1);         
                Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps!!!");
                metropolis.run(manager, recorder);
                        
                Status.detail("Temperature: " + metropolis.getTemp());
                Status.detail("Trigger ratio: " + recorder.getTriggerRatio());
                Status.detail("Average energy: " + recorder.getAverageValue());
                

             
                /*
                 * print out the ground state for every temperature step
                 * 
                 */
                SuperStructure groundState4 = appliedCE.getSuperStructure();            
                POSCAR outfile4= new POSCAR(groundState4);                        
                //outfile4.writeFile(fileDirInner + "snapShot." + stepNum + "-temp=" + temp + ".vasp");
                //outfile4.writeVICSOutFile(fileDirInner + "snapShot." + stepNum  + "-temp=" + temp + ".out"); 

                System.out.println("\nOuter metropolis run to get the cleanslab under ChemPot(Pt-Ni): ");
                System.out.println("Temperature: " + metropolis.getTemp());
                System.out.println("# of iterations: " + metropolis.getNumIterations());
                System.out.println("Trigger ratio: " + recorder.getTriggerRatio());
                System.out.println("Average energy: " + recorder.getAverageValue());
                
                System.out.println("Average current: " + recorder.getAverageCurrent());
                System.out.println("Average current from occupied sites: " + recorder.getAverageCurrentOccupied());
                System.out.println("Average current from unoccupied sites: " + recorder.getAverageCurrentUnOccupied());
                System.out.println("Average Ads current: " + recorder.getAverageAdsCurrent());
                System.out.println("Max current: " + recorder.getMaxAvgCurrent());
                System.out.println("Average OH Occupancy: " + recorder.getAverageOHOccupancy());
                System.out.println("Average OH Binding Energy: " + recorder.getAverageOHBindingEnergy());
                System.out.println("Average OH Binding Energy Sq: " + recorder.getAverageOHBindingEnergySq());
                double stdDev = Math.sqrt(recorder.getAverageOHBindingEnergySq() - (recorder.getAverageOHBindingEnergy() *  recorder.getAverageOHBindingEnergy()));
                System.out.println("Average OH Binding Energy Std Dev: " + stdDev);
                
                System.out.println("Average Face OH Occupancy: " + recorder.getAverageFaceOHOccupancy());
                System.out.println("Average Face OH Binding Energy: " + recorder.getAverageFaceOHBindingEnergy());
                
                System.out.println("Average Edge OH Occupancy: " + recorder.getAverageEdgeOHOccupancy());
                System.out.println("Average Edge OH Binding Energy: " + recorder.getAverageEdgeOHBindingEnergy());
                
                
                /*
                // the manager: DecorationGCHalfSlabFreeze_OAdsEds_KMC_Manager, which contains KMC (reject-free)
                long totalNumGeneratedAdsEvents = manager.getNumGeneratedAdsEvents();
                long totalNumGeneratedDesEvents = manager.getNumGeneratedDesEvents();

                long totalNumAcceptedAdsEvents = manager.getNumAcceptedAdsEvents();
                long totalNumAcceptedDesEvents = manager.getNumAcceptedDesEvents();

                double totalKMCTime = manager.getKMCTime();
                
                double avgAdsBarrier = manager.getAvgAdsBarrierEvents();
                double avgEdsBarrier = manager.getAvgDesBarrierEvents();

                System.out.println("totalNumGeneratedAdsEvents: " + totalNumGeneratedAdsEvents);
                System.out.println("totalNumGeneratedDesEvents: " + totalNumGeneratedDesEvents);
                System.out.println("totalNumAcceptedAdsEvents: " + totalNumAcceptedAdsEvents);
                System.out.println("totalNumAcceptedDesEvents: " + totalNumAcceptedDesEvents);
                System.out.println("avgAdsBarrier: " + avgAdsBarrier);     
                System.out.println("avgEdsBarrier: " + avgEdsBarrier);    
                System.out.println("totalKMCTime: " + totalKMCTime);    
                System.out.println("Average current: " + totalNumAcceptedDesEvents / totalKMCTime);
                */
                
                double avgE = recorder.getAverageValue();
                double avgE2 = recorder.getAverageSquaredValue();
                double sigSq = avgE2 - avgE * avgE;
                double cV = (avgE2 - avgE*avgE) / (kb * temp * temp);
                cV /= appliedCE.numSigmaSites();

                double maxCurrent = recorder.getMaxAvgCurrent();
                double avgCurrent = recorder.getAverageCurrent();
                double avgCurrentOccupied = recorder.getAverageCurrentOccupied();
                double avgCurrentUnOccupied = recorder.getAverageCurrentUnOccupied();                
                double avgAdsCurrent = recorder.getAverageAdsCurrent();
                double avgOHOccupancy = recorder.getAverageOHOccupancy();
                double avgFaceOHOccupancy = recorder.getAverageFaceOHOccupancy();
                double avgEdgeOHOccupancy = recorder.getAverageEdgeOHOccupancy();
                double avgEads =  recorder.getAverageOHBindingEnergy();
                double avgFaceEads =  recorder.getAverageFaceOHBindingEnergy();
                double avgEdgeEads =  recorder.getAverageEdgeOHBindingEnergy();
                double avgEads2 = recorder.getAverageOHBindingEnergySq();
                
                double[] avgEadsAllSites = recorder.getAverageOHBindingEnergyEachSite();
                int[] allSurfaceSitesIndices = recorder.getSurfaceOHSitesIndices();
                int[] allSurfaceSitesSigmaIndices = recorder.getSurfaceOHSitesSigmaIndices();

                int[] allSurfaceSitesCN = recorder.getSurfaceOHSitesCN();

                for(int siteNum=0; siteNum< avgEadsAllSites.length; siteNum++) {
                    String inputStr = siteNum + "   " + allSurfaceSitesIndices[siteNum]  + "   " +   allSurfaceSitesSigmaIndices[siteNum]  + "   " +   avgEadsAllSites[siteNum] + "   " + allSurfaceSitesCN[siteNum] + "\r\n";
                    bw_eachSite.write(inputStr);    
                    bw_eachSite.flush();
                }
                
                int totalSurfacePt = recorder.getNumOHSites();
                int totalSurfacePt_111 = recorder.getNumOHFaceSites();
                int totalSurfacePt_EdgeVertex = recorder.getNumOHEdgeSites();
                
                Structure sturctureMaxAvgCurrent = recorder.getStructureMaxTempAvgCurrent();
                
                POSCAR outfileMaxAvgCurrent= new POSCAR(sturctureMaxAvgCurrent);                                        
                //outfileMaxAvgCurrent.writeFile(fileDirInner + "snapshot." + stepNum + "-temp=" + temp + "-maxAvgCurrent.vasp");
                
                
                //System.out.println("Temp: " + temp + "\tConcentration: " + concentration + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                Status.detail("Temp: " + temp + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                
                Structure snapShot = appliedCE.getStructure(); 
              
                //System.out.println(stepNum);
                System.out.println("the num of steps is: " + stepNum + "\n");
                
                
                //String stringNum = temp + "   " + numPt + "   " + numNi + "   " + numCu + "   " + numOTem + "   " + compositionPtTemp_bulk + "   " + compositionPtTemp + "   " + latticeP_Temp_bulk + "   " + latticeP_Temp + "   " + strain_Temp_bulk + "   " + strain_Temp + "   " + shiftAdsEnergy + "  " + avgE + "   " + avgE2 + "   " + Cv_Slab + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOOccupancy + "    " + avgEads + "    " + avgEads2 + "   " + stdDev + "   " + compositionNiSplitTemp[0] + "   " + compositionNiSplitTemp[1] + "   " + compositionNiSplitTemp[2] + "   " + compositionNiSplitTemp[3] + "   " + NiComp_567_Freeze[0] + "   " + NiComp_567_Freeze[1] + "   " + NiComp_567_Freeze[2] + "\r\n";
                
                String stringNum = temp + "   " + numPt + "   " + numNi + "   " + numCu + "   " + totalSurfacePt + "   " + totalSurfacePt_111 + "   " + totalSurfacePt_EdgeVertex + "   " + (((double)totalSurfacePt_111) / totalSurfacePt) + "   " + (((double)totalSurfacePt) / numPt) + "   " + avgE + "   " + avgE2 + "   " + cV + "   " + muOH_Inner_KMC + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOHOccupancy + "    " + avgFaceOHOccupancy + "    " + avgEdgeOHOccupancy + "    " + avgEads + "    " + avgFaceEads + "    " + avgEdgeEads + "    " + avgEads2 + "   " + stdDev + "\r\n";
                
                System.out.println("Temp and Cv are: " + stringNum);
                bw.write(stringNum);    
                bw.flush();

                
                returnArr[0] = temp;
                returnArr[1] = numPt; 
                returnArr[2] = numNi; 
                
                returnArr[3] = numCu; 
                
                returnArr[4] = avgE;
                returnArr[5] = avgE2;
                returnArr[6] = cV;
                returnArr[7] = muOH_Inner_KMC;
                returnArr[8] = avgCurrent;
                returnArr[9] = avgAdsCurrent; 
                returnArr[10] = maxCurrent ;
                returnArr[11] = avgOHOccupancy;
                returnArr[12] = avgFaceOHOccupancy;
                returnArr[13] = avgEdgeOHOccupancy;
                returnArr[14] = avgEads;
                returnArr[15] = avgFaceEads;
                returnArr[16] = avgEdgeEads;
                returnArr[17] = avgEads2;
                returnArr[18] = stdDev;

                returnArr[19] = totalSurfacePt;
                returnArr[20] = totalSurfacePt_111;
                returnArr[21] = totalSurfacePt_EdgeVertex;
                returnArr[22] = ((double)totalSurfacePt_111) / totalSurfacePt;
                returnArr[23] = ((double)totalSurfacePt) / numPt;

            } 
            
            bw.flush();
            bw.close();
            fw.close();

        }
            
        catch (IOException e) {
            e.printStackTrace();
        }
            

        
        groundState = appliedCE.getSuperStructure().getCompactStructure();
        int numPtF = groundState.numDefiningSitesWithSpecies(Species.platinum);
        int numNiF = groundState.numDefiningSitesWithSpecies(Species.nickel);
        POSCAR outfile5= new POSCAR(groundState);
        outfile5.setDescription("fixedShape");
       // outfile5.writeFile(fileDirInner  + "groundState." + numPtF +"."+ numNiF + ".vasp");
        //outfile5.writeVICSOutFile(fileDirInner  + "groundState." + numPtF +"."+ numNiF + ".out");

        
        return returnArr;
        
      } 
      
 public static double[] runGCMonteCarlo_KMC_Current_slab_test_new(ClusterExpansion ce, String ceVersion, String fileDir_outer, String fileDir_inner, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterations, double equilibraRatio_KMC, Species OH_Spec) {
        
        double[] returnArr = new double[20];

        
        double kb = 8.6173423E-5;
        //double stepsFor4 = 10;
        
        //String fileDir = ORR_PREDICTION_DIR + "_variedOCoverage/" + initialStruct + "/";
        String fileDir = ORR_PREDICTION_DIR + fileDir_outer + "/" + fileDir_inner + "/";

        
        File dirFile  = new  File(fileDir);
        if ( ! (dirFile.exists())  &&   ! (dirFile.isDirectory())) {
                boolean  creadok  =  dirFile.mkdirs();
        }
        
        
        String fileDirInner = fileDir + "muOH_KMC=" + muOH_Inner_KMC + "-numKMCItes=" + numKMCIterations + "-KMCEquRatio=" + equilibraRatio_KMC + "/";
        File dirFileInner  = new  File(fileDirInner);
        if ( ! (dirFileInner.exists())  &&   ! (dirFileInner.isDirectory())) {
                boolean  creadok  =  dirFileInner.mkdirs();
        }
        
        ClusterGroup[] activeGroups = ce.getAllGroups(true);
        
        //double CutoffPerSite = 0;
        double CutoffPerSite = 1E-5;
        System.out.println("Using Cutoff per site of " + CutoffPerSite + " to trim insignificant groups...");
        activeGroups = ce.getSignificantGroups(activeGroups, ce.numSigmaSites() * CutoffPerSite);
        System.out.println(activeGroups.length + " significant groups remaining.\n\n");
        
        double startTemp = startTemp2;
        double endTemp = endTemp2;
        double tempIncrement = tempIncre;
        
        int numPasses = numPasses2;
        System.out.println("the number of numPasses : " + numPasses);

        /*
         * chemical potential of all species
         */
        
        double chemPotPt = 0.0;
        double chemPotNi = 0.0; 
       
        int x_slab = 2;
        int y_slab = 2;
        int z_slab = 6;
        
        x_slab = cellSize;
        y_slab = cellSize;
        
        /*
        //bad for "9-layer-slab-25.11-cleanslab-prim=144.vasp"
        //bad for "9-layer-slab-27.9-cleanslab-300K-prim=144.vasp"
        int[][] superToDirect = new int[][] {
                {x_slab, 0, -x_slab},
                {0, -y_slab, y_slab},
                {z_slab, z_slab, z_slab}
        }; 
        */
        
        
        
        int[][] superToDirect = new int[][] {
            {0, -y_slab, y_slab},
            {x_slab, 0, -x_slab},
            {-z_slab, -z_slab, -z_slab}
        }; 
        
    
        /*
        int[][] superToDirect = new int[][] {
                {-cellSize, cellSize, cellSize},
                {cellSize, -cellSize, cellSize},
                {cellSize, cellSize, -cellSize}
        };
        */
        
        
        AbstractAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        //AbstractAppliedCE appliedCE = new MakeMonteAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        System.out.println("Supercell has " + appliedCE.numPrimCells() + " primitive cells.");
        appliedCE.activateGroups(activeGroups);
        System.out.println("the number of sites is: " + appliedCE.numSigmaSites());
        
        SuperStructure groundState2 = appliedCE.getSuperStructure();
        
        POSCAR outfile2= new POSCAR(groundState2);
        outfile2.writeFile(fileDirInner + "groundState_0.vasp");
        //outfile2.writeVICSOutFile(fileDirInner + "groundState_0.out");
        System.out.println("write out the initial structure before decoration!!!!!!");
        
        Metropolis metropolis = new Metropolis();
        metropolis.setOptimistic(true);
        
        //appliedCE.decorateRandomly(new int[] {0, 1, 2}, new int[] {numNi, numPt, numVac});

        
        //POSCAR infile = new POSCAR(fileDir + initialStruct + ".vasp", true);
        //POSCAR infile = new POSCAR(fileDir + initialStruct, true);
        POSCAR infile = new POSCAR(fileDir + "snapShot-afterKMC.vasp", true);
        Structure structure = new Structure(infile);
        appliedCE.decorateFromStructure(structure);
        
        SuperStructure groundState3 = appliedCE.getSuperStructure();
        
        POSCAR outfile3= new POSCAR(groundState3);
        outfile3.writeFile(fileDirInner + "groundState_1.vasp");
        //outfile3.writeVICSOutFile(fileDirInner + "groundState_1.out");
        System.out.println("write out the structure after decoration!!!!!!");
        


        
        int numPt = structure.numDefiningSitesWithSpecies(Species.platinum);
       // int numNi = structure.numDefiningSitesWithSpecies(Species.nickel);
        int numCu = structure.numDefiningSitesWithSpecies(Species.copper);
        int numVac = appliedCE.numSigmaSites()  - numCu/ - numPt;
        int numAtoms = numPt + numCu;
        
        //DecorationCanonicalManager manager = new DecorationCanonicalManager(appliedCE); // Use this line if you want a canonical run

        //DecorationCrossCanonicalManager manager = new DecorationCrossCanonicalManager(appliedCE); // Use this line if you want a canonical run
        NanoCanonicalManagerLowMemFixedShape manager = new NanoCanonicalManagerLowMemFixedShape(appliedCE);
        
        metropolis.setTemp(startTemp * kb); // Whatever you want
        
        //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
        metropolis.setNumIterations(1);
        Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps...");
        metropolis.runBasic(manager);

        //System.out.println("Initial FBeta: " + totalFBeta);
        //System.out.println("Running for real...");
        
        
        if ((endTemp - startTemp) / tempIncrement < 0) {
            tempIncrement *= -1;
        }
        int numSteps = (int) Math.floor((endTemp - startTemp) / tempIncrement);
        System.out.println("the # of step of Monte Carlo simulation is: " + numSteps + "\n\n");
        
        Structure groundState = null;
        
        
        try{
            
            //FileWriter fw = new FileWriter(fileDirInner + "compositionProfile.txt", false);
            FileWriter fw = new FileWriter(fileDir + "/muOH_KMC=" + muOH_Inner_KMC + "-nKMCItes=" + numKMCIterations + "-avgCurrent.txt", false);
            BufferedWriter bw = new BufferedWriter(fw); 
            
            bw.write("temp numPt  numCu numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgEads avgEads2 stdDev " + "\r\n");
            bw.flush();
            
            for (int stepNum = 0; stepNum <= numSteps; stepNum++) {
                
                double temp = startTemp + stepNum * tempIncrement;
                metropolis.setTemp(kb * temp);
                
                /*
                 * If want to run grand canonical ensemble, in which the compositions are allowed to vary.
                 * thus, we think we can find some interesting phenomena.
                 * I will create a mirror symmetrical grand cannonical ensemble manager, named as DecorationMSymGCManager.
                 */
                //DecorationMSymGCManager manager = new DecorationMSymGCManager(appliedCE, mirrorOp);
                
                //SiteConcentrationRecorder recorder = new SiteConcentrationRecorder(manager, appliedCE); 
                
                //TODO
                NanoSiteConcentration_OHAdsDes_KMC_Recorder recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 
                //SiteConcentration_OHAdsDes_KMC_Recorder_Split recorder = new SiteConcentration_OHAdsDes_KMC_Recorder_Split(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 

                // just for symmetrical Monte Carlo simulations
               // Equilibration            
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
                //structMetropolis.runBasic(manager);

                // Run for real
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());         
                metropolis.setNumIterations(1);         
                Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps!!!");
                metropolis.run(manager, recorder);
                        
                Status.detail("Temperature: " + metropolis.getTemp());
                Status.detail("Trigger ratio: " + recorder.getTriggerRatio());
                Status.detail("Average energy: " + recorder.getAverageValue());
                

             
                /*
                 * print out the ground state for every temperature step
                 * 
                 */
                SuperStructure groundState4 = appliedCE.getSuperStructure();            
                POSCAR outfile4= new POSCAR(groundState4);                        
                outfile4.writeFile(fileDirInner + "snapShot." + stepNum + "-temp=" + temp + ".vasp");
                //outfile4.writeVICSOutFile(fileDirInner + "snapShot." + stepNum  + "-temp=" + temp + ".out"); 

                System.out.println("\nOuter metropolis run to get the cleanslab under ChemPot(Pt-Ni): ");
                System.out.println("Temperature: " + metropolis.getTemp());
                System.out.println("# of iterations: " + metropolis.getNumIterations());
                System.out.println("Trigger ratio: " + recorder.getTriggerRatio());
                System.out.println("Average energy: " + recorder.getAverageValue());
                
                System.out.println("Average current: " + recorder.getAverageCurrent());
                System.out.println("Average current from occupied sites: " + recorder.getAverageCurrentOccupied());
                System.out.println("Average current from unoccupied sites: " + recorder.getAverageCurrentUnOccupied());
                System.out.println("Average Ads current: " + recorder.getAverageAdsCurrent());
                System.out.println("Max current: " + recorder.getMaxAvgCurrent());
                System.out.println("Average OH Occupancy: " + recorder.getAverageOHOccupancy());
                System.out.println("Average OH Binding Energy: " + recorder.getAverageOHBindingEnergy());
                System.out.println("Average OH Binding Energy Sq: " + recorder.getAverageOHBindingEnergySq());
                double stdDev = Math.sqrt(recorder.getAverageOHBindingEnergySq() - (recorder.getAverageOHBindingEnergy() *  recorder.getAverageOHBindingEnergy()));
                System.out.println("Average OH Binding Energy Std Dev: " + stdDev);
                
                
                /*
                // the manager: DecorationGCHalfSlabFreeze_OAdsEds_KMC_Manager, which contains KMC (reject-free)
                long totalNumGeneratedAdsEvents = manager.getNumGeneratedAdsEvents();
                long totalNumGeneratedDesEvents = manager.getNumGeneratedDesEvents();

                long totalNumAcceptedAdsEvents = manager.getNumAcceptedAdsEvents();
                long totalNumAcceptedDesEvents = manager.getNumAcceptedDesEvents();

                double totalKMCTime = manager.getKMCTime();
                
                double avgAdsBarrier = manager.getAvgAdsBarrierEvents();
                double avgEdsBarrier = manager.getAvgDesBarrierEvents();

                System.out.println("totalNumGeneratedAdsEvents: " + totalNumGeneratedAdsEvents);
                System.out.println("totalNumGeneratedDesEvents: " + totalNumGeneratedDesEvents);
                System.out.println("totalNumAcceptedAdsEvents: " + totalNumAcceptedAdsEvents);
                System.out.println("totalNumAcceptedDesEvents: " + totalNumAcceptedDesEvents);
                System.out.println("avgAdsBarrier: " + avgAdsBarrier);     
                System.out.println("avgEdsBarrier: " + avgEdsBarrier);    
                System.out.println("totalKMCTime: " + totalKMCTime);    
                System.out.println("Average current: " + totalNumAcceptedDesEvents / totalKMCTime);
                */
                
                double avgE = recorder.getAverageValue();
                double avgE2 = recorder.getAverageSquaredValue();
                double sigSq = avgE2 - avgE * avgE;
                double cV = (avgE2 - avgE*avgE) / (kb * temp * temp);
                cV /= appliedCE.numSigmaSites();

                double maxCurrent = recorder.getMaxAvgCurrent();
                double avgCurrent = recorder.getAverageCurrent();
                double avgCurrentOccupied = recorder.getAverageCurrentOccupied();
                double avgCurrentUnOccupied = recorder.getAverageCurrentUnOccupied();                
                double avgAdsCurrent = recorder.getAverageAdsCurrent();
                double avgOHOccupancy = recorder.getAverageOHOccupancy();
                double avgEads =  recorder.getAverageOHBindingEnergy();
                double avgEads2 = recorder.getAverageOHBindingEnergySq();
                
                
                
                
                
                int totalSurfacePt = recorder.getNumOHSites();
                int totalSurfacePt_111 = recorder.getNumOHFaceSites();
                int totalSurfacePt_EdgeVertex = recorder.getNumOHEdgeSites();
                
                Structure sturctureMaxAvgCurrent = recorder.getStructureMaxTempAvgCurrent();
                
                POSCAR outfileMaxAvgCurrent= new POSCAR(sturctureMaxAvgCurrent);                                        outfileMaxAvgCurrent.writeFile(fileDirInner + "snapshot." + stepNum + "-temp=" + temp + "-maxAvgCurrent.vasp");
                
                
                //System.out.println("Temp: " + temp + "\tConcentration: " + concentration + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                Status.detail("Temp: " + temp + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                
                Structure snapShot = appliedCE.getStructure(); 
              
                //System.out.println(stepNum);
                System.out.println("the num of steps is: " + stepNum + "\n");
                
                
                //String stringNum = temp + "   " + numPt + "   " + numNi + "   " + numCu + "   " + numOTem + "   " + compositionPtTemp_bulk + "   " + compositionPtTemp + "   " + latticeP_Temp_bulk + "   " + latticeP_Temp + "   " + strain_Temp_bulk + "   " + strain_Temp + "   " + shiftAdsEnergy + "  " + avgE + "   " + avgE2 + "   " + Cv_Slab + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOOccupancy + "    " + avgEads + "    " + avgEads2 + "   " + stdDev + "   " + compositionNiSplitTemp[0] + "   " + compositionNiSplitTemp[1] + "   " + compositionNiSplitTemp[2] + "   " + compositionNiSplitTemp[3] + "   " + NiComp_567_Freeze[0] + "   " + NiComp_567_Freeze[1] + "   " + NiComp_567_Freeze[2] + "\r\n";
                
                String stringNum = temp + "   " + numPt  + "   " + numCu + "   " + totalSurfacePt + "   " + totalSurfacePt_111 + "   " + totalSurfacePt_EdgeVertex + "   " + (((double)totalSurfacePt_111) / totalSurfacePt) + "   " + (((double)totalSurfacePt) / numPt) + "   " + avgE + "   " + avgE2 + "   " + cV + "   " + muOH_Inner_KMC + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOHOccupancy + "    " + avgEads + "    " + avgEads2 + "   " + stdDev + "\r\n";
                
                System.out.println("Temp and Cv are: " + stringNum);
                bw.write(stringNum);    
                bw.flush();

                
                returnArr[0] = temp;
                returnArr[1] = numPt; 
            //    returnArr[2] = numNi; 
                returnArr[3] = numCu; 
                returnArr[4] = avgE;
                returnArr[5] = avgE2;
                returnArr[6] = cV;
                returnArr[7] = muOH_Inner_KMC;
                returnArr[8] = avgCurrent;
                returnArr[9] = avgAdsCurrent; 
                returnArr[10] = maxCurrent ;
                returnArr[11] = avgOHOccupancy;
                returnArr[12] = avgEads;
                returnArr[13] = avgEads2;
                returnArr[14] = stdDev;

                returnArr[15] = totalSurfacePt;
                returnArr[16] = totalSurfacePt_111;
                returnArr[17] = totalSurfacePt_EdgeVertex;
                returnArr[18] = ((double)totalSurfacePt_111) / totalSurfacePt;
                returnArr[19] = ((double)totalSurfacePt) / numPt;

            } 
            
            bw.flush();
            bw.close();
            fw.close();

        }
            
        catch (IOException e) {
            e.printStackTrace();
        }
            

        
        groundState = appliedCE.getSuperStructure().getCompactStructure();
        int numPtF = groundState.numDefiningSitesWithSpecies(Species.platinum);
       // int numNiF = groundState.numDefiningSitesWithSpecies(Species.nickel);
        int numCuF = groundState.numDefiningSitesWithSpecies(Species.copper);
        POSCAR outfile5= new POSCAR(groundState);
        outfile5.setDescription("fixedShape");
        outfile5.writeFile(fileDirInner  + "groundState." + numPtF  + "."+ numCuF +  ".vasp");
        //outfile5.writeVICSOutFile(fileDirInner  + "groundState." + numPtF +"."+ numNiF + "."+ numCuF + ".out");

        
        return returnArr;
        
      } 
      
public static double[] runGCMonteCarlo_KMC_Current_slab_test_eachSite_new(ClusterExpansion ce, String ceVersion, String fileDir_outer, String fileDir_inner, int numPasses2, int startTemp2, int endTemp2, int tempIncre, int cellSize, double muOH_Inner_KMC, double numKMCIterations, double equilibraRatio_KMC, Species OH_Spec) {
        
        double[] returnArr = new double[20];

        
        double kb = 8.6173423E-5;
        //double stepsFor4 = 10;
        
        //String fileDir = ORR_PREDICTION_DIR + "_variedOCoverage/" + initialStruct + "/";
        String fileDir = ORR_PREDICTION_DIR + fileDir_outer + "/" + fileDir_inner + "/";

        
        File dirFile  = new  File(fileDir);
        if ( ! (dirFile.exists())  &&   ! (dirFile.isDirectory())) {
                boolean  creadok  =  dirFile.mkdirs();
        }
        
        
        String fileDirInner = fileDir + "muOH_KMC=" + muOH_Inner_KMC + "-numKMCItes=" + numKMCIterations + "-KMCEquRatio=" + equilibraRatio_KMC + "/";
        File dirFileInner  = new  File(fileDirInner);
        if ( ! (dirFileInner.exists())  &&   ! (dirFileInner.isDirectory())) {
                boolean  creadok  =  dirFileInner.mkdirs();
        }
        
        ClusterGroup[] activeGroups = ce.getAllGroups(true);
        
        //double CutoffPerSite = 0;
        double CutoffPerSite = 1E-5;
        System.out.println("Using Cutoff per site of " + CutoffPerSite + " to trim insignificant groups...");
        activeGroups = ce.getSignificantGroups(activeGroups, ce.numSigmaSites() * CutoffPerSite);
        System.out.println(activeGroups.length + " significant groups remaining.\n\n");
        
        double startTemp = startTemp2;
        double endTemp = endTemp2;
        double tempIncrement = tempIncre;
        
        int numPasses = numPasses2;
        System.out.println("the number of numPasses : " + numPasses);

        /*
         * chemical potential of all species
         */
        
        double chemPotPt = 0.0;
        double chemPotNi = 0.0; 
       
        int x_slab = 2;
        int y_slab = 2;
        int z_slab = 6;
        
        x_slab = cellSize;
        y_slab = cellSize;
        
        /*
        //bad for "9-layer-slab-25.11-cleanslab-prim=144.vasp"
        //bad for "9-layer-slab-27.9-cleanslab-300K-prim=144.vasp"
        int[][] superToDirect = new int[][] {
                {x_slab, 0, -x_slab},
                {0, -y_slab, y_slab},
                {z_slab, z_slab, z_slab}
        }; 
        */
        
        
        
        int[][] superToDirect = new int[][] {
            {0, -y_slab, y_slab},
            {x_slab, 0, -x_slab},
            {-z_slab, -z_slab, -z_slab}
        }; 
        
    
        /*
        int[][] superToDirect = new int[][] {
                {-cellSize, cellSize, cellSize},
                {cellSize, -cellSize, cellSize},
                {cellSize, cellSize, -cellSize}
        };
        */
        
        
        AbstractAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        //AbstractAppliedCE appliedCE = new MakeMonteAppliedCE(ce, superToDirect);  // Applies the cluster expansion to a superstructure
        System.out.println("Supercell has " + appliedCE.numPrimCells() + " primitive cells.");
        appliedCE.activateGroups(activeGroups);
        System.out.println("the number of sites is: " + appliedCE.numSigmaSites());
        
        SuperStructure groundState2 = appliedCE.getSuperStructure();
        
        POSCAR outfile2= new POSCAR(groundState2);
        outfile2.writeFile(fileDirInner + "groundState_0.vasp");
        //outfile2.writeVICSOutFile(fileDirInner + "groundState_0.out");
        System.out.println("write out the initial structure before decoration!!!!!!");
        
        Metropolis metropolis = new Metropolis();
        metropolis.setOptimistic(true);
        
        //appliedCE.decorateRandomly(new int[] {0, 1, 2}, new int[] {numNi, numPt, numVac});

        
        //POSCAR infile = new POSCAR(fileDir + initialStruct + ".vasp", true);
        //POSCAR infile = new POSCAR(fileDir + initialStruct, true);
        POSCAR infile = new POSCAR(fileDir + "snapShot-afterKMC.vasp", true);
        Structure structure = new Structure(infile);
        appliedCE.decorateFromStructure(structure);
        
        SuperStructure groundState3 = appliedCE.getSuperStructure();
        
        POSCAR outfile3= new POSCAR(groundState3);
        outfile3.writeFile(fileDirInner + "groundState_1.vasp");
        //outfile3.writeVICSOutFile(fileDirInner + "groundState_1.out");
        System.out.println("write out the structure after decoration!!!!!!");
        


        
        int numPt = structure.numDefiningSitesWithSpecies(Species.platinum);
        int numNi = structure.numDefiningSitesWithSpecies(Species.nickel);
        int numCu = structure.numDefiningSitesWithSpecies(Species.copper);
        int numVac = appliedCE.numSigmaSites()  - numCu - numNi - numPt;
        int numAtoms = numPt + numNi + numCu;
        
        //DecorationCanonicalManager manager = new DecorationCanonicalManager(appliedCE); // Use this line if you want a canonical run

        //DecorationCrossCanonicalManager manager = new DecorationCrossCanonicalManager(appliedCE); // Use this line if you want a canonical run
        NanoCanonicalManagerLowMemFixedShape manager = new NanoCanonicalManagerLowMemFixedShape(appliedCE);
        
        metropolis.setTemp(startTemp * kb); // Whatever you want
        
        //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
        metropolis.setNumIterations(1);
        Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps...");
        metropolis.runBasic(manager);

        //System.out.println("Initial FBeta: " + totalFBeta);
        //System.out.println("Running for real...");
        
        
        if ((endTemp - startTemp) / tempIncrement < 0) {
            tempIncrement *= -1;
        }
        int numSteps = (int) Math.floor((endTemp - startTemp) / tempIncrement);
        System.out.println("the # of step of Monte Carlo simulation is: " + numSteps + "\n\n");
        
        Structure groundState = null;
        
        
        try{
            
            //FileWriter fw = new FileWriter(fileDirInner + "compositionProfile.txt", false);
            FileWriter fw = new FileWriter(fileDir + "/muOH_KMC=" + muOH_Inner_KMC + "-nKMCItes=" + numKMCIterations + "-avgCurrent.txt", false);
            BufferedWriter bw = new BufferedWriter(fw); 
            
            bw.write("temp numPt numNi  numCu numPtSurface numPt_111 numPt_EdgeVertex ratio_111 ratio_MA avgE avgE2 cV muOH_Inner_KMC avgCurrent avgAdsCurrent maxCurrent avgOHOccupancy avgEads avgEads2 stdDev " + "\r\n");
            bw.flush();
            FileWriter fw_eachSite = new FileWriter(fileDir + "/muOH_KMC=" + muOH_Inner_KMC + "-nKMCItes=" + numKMCIterations + "-avgEads-eachSite.txt", false);
            BufferedWriter bw_eachSite = new BufferedWriter(fw_eachSite); 
            
            bw_eachSite.write("siteNum surfaceSiteIndex surfaceSiteSigmaIndex avgEads CN avgEads2 stdDev " + "\r\n");
            bw_eachSite.flush();
            for (int stepNum = 0; stepNum <= numSteps; stepNum++) {
                
                double temp = startTemp + stepNum * tempIncrement;
                metropolis.setTemp(kb * temp);
                
                /*
                 * If want to run grand canonical ensemble, in which the compositions are allowed to vary.
                 * thus, we think we can find some interesting phenomena.
                 * I will create a mirror symmetrical grand cannonical ensemble manager, named as DecorationMSymGCManager.
                 */
                //DecorationMSymGCManager manager = new DecorationMSymGCManager(appliedCE, mirrorOp);
                
                //SiteConcentrationRecorder recorder = new SiteConcentrationRecorder(manager, appliedCE); 
                
                //TODO
               // NanoSiteConcentration_OHAdsDes_KMC_Recorder recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec);
                //NanoSiteConcentration_OHAdsDes_KMC_Recorder recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 
                NanoSiteConcentration_OHAdsDes_KMC_Recorder_allSites recorder = new NanoSiteConcentration_OHAdsDes_KMC_Recorder_allSites(manager, appliedCE, structure, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec,fileDirInner); 
                 
                //SiteConcentration_OHAdsDes_KMC_Recorder_Split recorder = new SiteConcentration_OHAdsDes_KMC_Recorder_Split(manager, appliedCE, muOH_Inner_KMC, numKMCIterations, equilibraRatio_KMC, OH_Spec); 

                // just for symmetrical Monte Carlo simulations
               // Equilibration            
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());
                //structMetropolis.runBasic(manager);

                // Run for real
                //structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites());         
                metropolis.setNumIterations(1);         
                Status.detail("Running equilibration for " + metropolis.getNumIterations() + " steps!!!");
                metropolis.run(manager, recorder);
                        
                Status.detail("Temperature: " + metropolis.getTemp());
                Status.detail("Trigger ratio: " + recorder.getTriggerRatio());
                Status.detail("Average energy: " + recorder.getAverageValue());
                

             
                /*
                 * print out the ground state for every temperature step
                 * 
                 */
                SuperStructure groundState4 = appliedCE.getSuperStructure();            
                POSCAR outfile4= new POSCAR(groundState4);                        
               // outfile4.writeFile(fileDirInner + "snapShot." + stepNum + "-temp=" + temp + ".vasp");
                //outfile4.writeVICSOutFile(fileDirInner + "snapShot." + stepNum  + "-temp=" + temp + ".out"); 

                System.out.println("\nOuter metropolis run to get the cleanslab under ChemPot(Pt-Ni): ");
                System.out.println("Temperature: " + metropolis.getTemp());
                System.out.println("# of iterations: " + metropolis.getNumIterations());
                System.out.println("Trigger ratio: " + recorder.getTriggerRatio());
                System.out.println("Average energy: " + recorder.getAverageValue());
                
                System.out.println("Average current: " + recorder.getAverageCurrent());
                System.out.println("Average current from occupied sites: " + recorder.getAverageCurrentOccupied());
                System.out.println("Average current from unoccupied sites: " + recorder.getAverageCurrentUnOccupied());
                System.out.println("Average Ads current: " + recorder.getAverageAdsCurrent());
                System.out.println("Max current: " + recorder.getMaxAvgCurrent());
                System.out.println("Average OH Occupancy: " + recorder.getAverageOHOccupancy());
                System.out.println("Average OH Binding Energy: " + recorder.getAverageOHBindingEnergy());
                System.out.println("Average OH Binding Energy Sq: " + recorder.getAverageOHBindingEnergySq());
                double stdDev = Math.sqrt(recorder.getAverageOHBindingEnergySq() - (recorder.getAverageOHBindingEnergy() *  recorder.getAverageOHBindingEnergy()));
                System.out.println("Average OH Binding Energy Std Dev: " + stdDev);
                
                
                /*
                // the manager: DecorationGCHalfSlabFreeze_OAdsEds_KMC_Manager, which contains KMC (reject-free)
                long totalNumGeneratedAdsEvents = manager.getNumGeneratedAdsEvents();
                long totalNumGeneratedDesEvents = manager.getNumGeneratedDesEvents();

                long totalNumAcceptedAdsEvents = manager.getNumAcceptedAdsEvents();
                long totalNumAcceptedDesEvents = manager.getNumAcceptedDesEvents();

                double totalKMCTime = manager.getKMCTime();
                
                double avgAdsBarrier = manager.getAvgAdsBarrierEvents();
                double avgEdsBarrier = manager.getAvgDesBarrierEvents();

                System.out.println("totalNumGeneratedAdsEvents: " + totalNumGeneratedAdsEvents);
                System.out.println("totalNumGeneratedDesEvents: " + totalNumGeneratedDesEvents);
                System.out.println("totalNumAcceptedAdsEvents: " + totalNumAcceptedAdsEvents);
                System.out.println("totalNumAcceptedDesEvents: " + totalNumAcceptedDesEvents);
                System.out.println("avgAdsBarrier: " + avgAdsBarrier);     
                System.out.println("avgEdsBarrier: " + avgEdsBarrier);    
                System.out.println("totalKMCTime: " + totalKMCTime);    
                System.out.println("Average current: " + totalNumAcceptedDesEvents / totalKMCTime);
                */
                
                double avgE = recorder.getAverageValue();
                double avgE2 = recorder.getAverageSquaredValue();
                double sigSq = avgE2 - avgE * avgE;
                double cV = (avgE2 - avgE*avgE) / (kb * temp * temp);
                cV /= appliedCE.numSigmaSites();

                double maxCurrent = recorder.getMaxAvgCurrent();
                double avgCurrent = recorder.getAverageCurrent();
                double avgCurrentOccupied = recorder.getAverageCurrentOccupied();
                double avgCurrentUnOccupied = recorder.getAverageCurrentUnOccupied();                
                double avgAdsCurrent = recorder.getAverageAdsCurrent();
                double avgOHOccupancy = recorder.getAverageOHOccupancy();
                double avgEads =  recorder.getAverageOHBindingEnergy();
                double avgEads2 = recorder.getAverageOHBindingEnergySq();
                
                double[] avgEadsAllSites = recorder.getAverageOHBindingEnergyEachSite();
                int[] allSurfaceSitesIndices = recorder.getSurfaceOHSitesIndices();
                int[] allSurfaceSitesSigmaIndices = recorder.getSurfaceOHSitesSigmaIndices();
                int[] allSurfaceSitesCN = recorder.getSurfaceOHSitesCN();
                
                for(int siteNum=0; siteNum< avgEadsAllSites.length; siteNum++) {
                    String inputStr = siteNum + "   " + allSurfaceSitesIndices[siteNum]  + "   " +   allSurfaceSitesSigmaIndices[siteNum]  + "   " +   avgEadsAllSites[siteNum] + "   " + allSurfaceSitesCN[siteNum] + "\r\n";
                    bw_eachSite.write(inputStr);    
                    bw_eachSite.flush();
                }
                
                
                int totalSurfacePt = recorder.getNumOHSites();
                int totalSurfacePt_111 = recorder.getNumOHFaceSites();
                int totalSurfacePt_EdgeVertex = recorder.getNumOHEdgeSites();
                
                
                
                
                
                Structure sturctureMaxAvgCurrent = recorder.getStructureMaxTempAvgCurrent();
                
                POSCAR outfileMaxAvgCurrent= new POSCAR(sturctureMaxAvgCurrent);                                        
               // outfileMaxAvgCurrent.writeFile(fileDirInner + "snapshot." + stepNum + "-temp=" + temp + "-maxAvgCurrent.vasp");
                
                
                //System.out.println("Temp: " + temp + "\tConcentration: " + concentration + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                Status.detail("Temp: " + temp + "\tE: " + avgE + "\tcv: " + cV+ "\tsigSq: " + sigSq);
                
                Structure snapShot = appliedCE.getStructure(); 
              
                //System.out.println(stepNum);
                System.out.println("the num of steps is: " + stepNum + "\n");
                
                
                //String stringNum = temp + "   " + numPt + "   " + numNi + "   " + numCu + "   " + numOTem + "   " + compositionPtTemp_bulk + "   " + compositionPtTemp + "   " + latticeP_Temp_bulk + "   " + latticeP_Temp + "   " + strain_Temp_bulk + "   " + strain_Temp + "   " + shiftAdsEnergy + "  " + avgE + "   " + avgE2 + "   " + Cv_Slab + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOOccupancy + "    " + avgEads + "    " + avgEads2 + "   " + stdDev + "   " + compositionNiSplitTemp[0] + "   " + compositionNiSplitTemp[1] + "   " + compositionNiSplitTemp[2] + "   " + compositionNiSplitTemp[3] + "   " + NiComp_567_Freeze[0] + "   " + NiComp_567_Freeze[1] + "   " + NiComp_567_Freeze[2] + "\r\n";
                
                String stringNum = temp + "   " + numPt + "   " + numNi + "   " + numCu + "   " + totalSurfacePt + "   " + totalSurfacePt_111 + "   " + totalSurfacePt_EdgeVertex + "   " + (((double)totalSurfacePt_111) / totalSurfacePt) + "   " + (((double)totalSurfacePt) / numPt) + "   " + avgE + "   " + avgE2 + "   " + cV + "   " + muOH_Inner_KMC + "   " + avgCurrent + "   " + avgAdsCurrent + "   " + maxCurrent + "   " + avgOHOccupancy + "    " + avgEads + "    " + avgEads2 + "   " + stdDev + "\r\n";
                
                System.out.println("Temp and Cv are: " + stringNum);
                bw.write(stringNum);    
                bw.flush();

                
                returnArr[0] = temp;
                returnArr[1] = numPt; 
                returnArr[2] = numNi; 
                returnArr[3] = numCu; 
                returnArr[4] = avgE;
                returnArr[5] = avgE2;
                returnArr[6] = cV;
                returnArr[7] = muOH_Inner_KMC;
                returnArr[8] = avgCurrent;
                returnArr[9] = avgAdsCurrent; 
                returnArr[10] = maxCurrent ;
                returnArr[11] = avgOHOccupancy;
                returnArr[12] = avgEads;
                returnArr[13] = avgEads2;
                returnArr[14] = stdDev;

                returnArr[15] = totalSurfacePt;
                returnArr[16] = totalSurfacePt_111;
                returnArr[17] = totalSurfacePt_EdgeVertex;
                returnArr[18] = ((double)totalSurfacePt_111) / totalSurfacePt;
                returnArr[19] = ((double)totalSurfacePt) / numPt;

            } 
            
            bw.flush();
            bw.close();
            fw.close();

        }
            
        catch (IOException e) {
            e.printStackTrace();
        }
            

        
        groundState = appliedCE.getSuperStructure().getCompactStructure();
        int numPtF = groundState.numDefiningSitesWithSpecies(Species.platinum);
        int numNiF = groundState.numDefiningSitesWithSpecies(Species.nickel);
        int numCuF = groundState.numDefiningSitesWithSpecies(Species.copper);
        POSCAR outfile5= new POSCAR(groundState);
        //outfile5.setDescription("fixedShape");
        //outfile5.writeFile(fileDirInner  + "groundState." + numPtF + "."+ numNiF + "."+ numCuF +  ".vasp");
        //outfile5.writeVICSOutFile(fileDirInner  + "groundState." + numPtF +"."+ numNiF + "."+ numCuF + ".out");

        
        return returnArr;
        
      } 
   

    
    
    
    

}
