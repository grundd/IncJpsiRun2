// AnalysisConfig.h
// David Grund, Feb 26, 2022
// Configure values of the parameters

#include "TSystem.h"
#include "TString.h"

TString str_subfolder = "";
TString str_in_DT_fldr = "";
TString str_in_DT_tree = "";
TString str_in_MC_fldr_rec = "";
TString str_in_MC_tree_rec = "";
TString str_in_MC_fldr_gen = "";
TString str_in_MC_tree_gen = "";

void InitAnalysis(Int_t iAnalysis){
    // pass1, 4bins, original MC PID
    if(iAnalysis == 0)
    {
        str_subfolder = "pass1_4bins/";
        nPtBins = 4;
        isPass3 = kFALSE;
        isPIDCalibrated = kFALSE;
        isNParInDSCBFixed = kFALSE;
    }
    // pass1, 5bins, original MC PID
    if(iAnalysis == 1)
    {
        str_subfolder = "pass1_5bins/";
        nPtBins = 5;
        isPass3 = kFALSE;
        isPIDCalibrated = kFALSE;
        isNParInDSCBFixed = kFALSE;
    }
    // pass3, 4bins, original MC PID
    if(iAnalysis == 10){
        str_subfolder = "pass3_4bins/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isPIDCalibrated = kFALSE;
        isNParInDSCBFixed = kFALSE;
    }
    // pass3, 5bins, original MC PID
    if(iAnalysis == 11)
    {
        str_subfolder = "pass3_5bins/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isPIDCalibrated = kFALSE;
        isNParInDSCBFixed = kFALSE;
    }
    // pass3, 4bins, calibrated MC PID
    if(iAnalysis == 12)
    {
        str_subfolder = "pass3_4bins_calibPID/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isPIDCalibrated = kTRUE;
        isNParInDSCBFixed = kFALSE;
    }
    // pass3, 5bins, calibrated MC PID
    if(iAnalysis == 13)
    {
        str_subfolder = "pass3_5bins_calibPID/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isPIDCalibrated = kTRUE;
        isNParInDSCBFixed = kFALSE;
    }
    // pass3, 4bins, calibrated MC PID, N params in DSCB fixed to 10.
    if(iAnalysis == 14)
    {
        str_subfolder = "pass3_4bins_calibPID_Nfix/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isPIDCalibrated = kTRUE;
        isNParInDSCBFixed = kTRUE;
    }
    // pass3, 5bins, calibrated MC PID, N params in DSCB fixed to 10.
    if(iAnalysis == 15)
    {
        str_subfolder = "pass3_5bins_calibPID_Nfix/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isPIDCalibrated = kTRUE;
        isNParInDSCBFixed = kTRUE;
    }
    // pass3, 5bins, calibrated MC PID, cut_fVertexZ set to 10 cm
    if(iAnalysis == 100)
    {
        str_subfolder = "pass3_5bins_calibPID_Nfix_Zcut10.0/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isPIDCalibrated = kTRUE;
        isNParInDSCBFixed = kTRUE;
        cut_fVertexZ = 10.0;
    }
    if(iAnalysis == 101)
    {
        str_subfolder = "pass3_5bins_calibPID_Nfix_Zcut12.5/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isPIDCalibrated = kTRUE;
        isNParInDSCBFixed = kTRUE;
        cut_fVertexZ = 12.5;
    }
    // set reduced run lists for the given pass
    SetReducedRunList(isPass3);
    // set the path to the input data
    if(!isPass3)
    {
        str_in_DT_fldr = "Trees/AnalysisData_pass1/";
        str_in_DT_tree = "AnalysisOutput/fTreeJPsi";
        // whether to use MC files with calibrated NSigmas
        if(!isPIDCalibrated) str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass1/";
        else                 str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass1/PIDCalibrated/";
        str_in_MC_tree_rec = "AnalysisOutput/fTreeJPsiMCRec";
        str_in_MC_fldr_gen = "Trees/AnalysisDataMC_pass1/";
        str_in_MC_tree_gen = "AnalysisOutput/fTreeJPsiMCGen";
    } 
    else 
    {
        str_in_DT_fldr = "Trees/AnalysisData_pass3/";
        str_in_DT_tree = "AnalysisOutput/fTreeJpsi";
        // whether to use MC files with calibrated NSigmas
        if(!isPIDCalibrated) str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass3/";
        else                 str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass3/PIDCalibrated/";        
        str_in_MC_tree_rec = "AnalysisOutput/fTreeJpsi";
        str_in_MC_fldr_gen = "Trees/AnalysisDataMC_pass3/";
        str_in_MC_tree_gen = "AnalysisOutput/fTreeJpsiMCGen";
    }
    // create a new one of the same name (nothing happens if it already exists)
    gSystem->Exec("mkdir -p Results/" + str_subfolder);
    gSystem->Exec("mkdir -p Trees/" + str_subfolder);

    Printf("*** Analysis initiated successfully. ***");
    Printf("*** Results will be stored in the folder: Results/%s ***", str_subfolder.Data());

    return;
}