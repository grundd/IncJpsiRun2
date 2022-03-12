// AnalysisConfig.h
// David Grund, Feb 26, 2022
// Configure values of the parameters

// ******************************
Bool_t START_FROM_CLEAN = kFALSE;
Int_t iWhichAnalysisToRun = 0;
// ******************************

TString str_subfolder = "";
TString str_in_DT_fldr = "";
TString str_in_DT_tree = "";
TString str_in_MC_fldr_rec = "";
TString str_in_MC_tree_rec = "";
TString str_in_MC_fldr_gen = "";
TString str_in_MC_tree_gen = "";

void SetReducedRunList(Bool_t pass3){

    if(!pass3){
        nRuns_18q = sizeof(DPG_LHC18q_pass1_reduced) / sizeof(DPG_LHC18q_pass1_reduced[0]); // 123 runs
        for(Int_t i = 0; i < nRuns_18q; i++) runList_18q.push_back(DPG_LHC18q_pass1_reduced[i]);
        nRuns_18r = sizeof(DPG_LHC18r_pass1_reduced) / sizeof(DPG_LHC18r_pass1_reduced[0]); // 96 runs
        for(Int_t i = 0; i < nRuns_18r; i++) runList_18r.push_back(DPG_LHC18r_pass1_reduced[i]);
    } else {
        nRuns_18q = sizeof(DPG_LHC18q_pass3_reduced) / sizeof(DPG_LHC18q_pass3_reduced[0]); // 122 runs
        for(Int_t i = 0; i < nRuns_18q; i++) runList_18q.push_back(DPG_LHC18q_pass3_reduced[i]);
        nRuns_18r = sizeof(DPG_LHC18r_pass3_reduced) / sizeof(DPG_LHC18r_pass3_reduced[0]); // 96 runs
        for(Int_t i = 0; i < nRuns_18r; i++) runList_18r.push_back(DPG_LHC18r_pass3_reduced[i]);
    }
    Printf("Number of runs in LHC18q run list: %i (%i)", nRuns_18q, (Int_t)runList_18q.size());
    Printf("Number of runs in LHC18r run list: %i (%i)", nRuns_18r, (Int_t)runList_18r.size());

    return;
}

void InitAnalysis(Int_t iAnalysis){
    // pass1, 4bins, original MC PID
    if(iAnalysis == 0){
        str_subfolder = "pass1_4bins/";
        nPtBins = 4;
        isPass3 = kFALSE;
        isZNcut = kFALSE;
        isPIDCalibrated = kFALSE;
    }
    // pass1, 5bins, original MC PID
    if(iAnalysis == 1){
        str_subfolder = "pass1_5bins/";
        nPtBins = 5;
        isPass3 = kFALSE;
        isZNcut = kFALSE;
        isPIDCalibrated = kFALSE;
    }
    // pass3, 4bins, original MC PID
    if(iAnalysis == 10){
        str_subfolder = "pass3_4bins/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isZNcut = kFALSE;
        isPIDCalibrated = kFALSE;
    }
    // pass3, 5bins, original MC PID
    if(iAnalysis == 11){
        str_subfolder = "pass3_5bins/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isZNcut = kFALSE;
        isPIDCalibrated = kFALSE;
    }
    // pass3, 4bins, calibrated MC PID
    if(iAnalysis == 12){
        str_subfolder = "pass3_4bins_calibPID/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isZNcut = kFALSE;
        isPIDCalibrated = kTRUE;
    }
    // pass3, 5bins, calibrated MC PID
    if(iAnalysis == 13){
        str_subfolder = "pass3_5bins_calibPID/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isZNcut = kFALSE;
        isPIDCalibrated = kTRUE;
    }
    // set reduced run lists for the given pass
    SetReducedRunList(isPass3);
    // set the path to the input data
    if(!isPass3){
        str_in_DT_fldr = "Trees/AnalysisData_pass1/";
        str_in_DT_tree = "AnalysisOutput/fTreeJPsi";
        // whether to use MC files with calibrated NSigmas
        if(!isPIDCalibrated) str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass1/";
        else                 str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass1/PIDCalibrated/";
        str_in_MC_tree_rec = "AnalysisOutput/fTreeJPsiMCRec";
        str_in_MC_fldr_gen = "Trees/AnalysisDataMC_pass1/";
        str_in_MC_tree_gen = "AnalysisOutput/fTreeJPsiMCGen";
    } else {
        str_in_DT_fldr = "Trees/AnalysisData_pass3/";
        str_in_DT_tree = "AnalysisOutput/fTreeJpsi";
        // whether to use MC files with calibrated NSigmas
        if(!isPIDCalibrated) str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass3/";
        else                 str_in_MC_fldr_rec = "Trees/AnalysisDataMC_pass3/PIDCalibrated/";        
        str_in_MC_tree_rec = "AnalysisOutput/fTreeJpsi";
        str_in_MC_fldr_gen = "Trees/AnalysisDataMC_pass3/";
        str_in_MC_tree_gen = "AnalysisOutput/fTreeJpsiMCGen";
    }
    // if we want to start with a clean folder
    if(START_FROM_CLEAN){
        // delete previous subfolder inside the Results/ and Trees/ folder (if exists)
        gSystem->Exec("rm -r Results/" + str_subfolder);
        gSystem->Exec("rm -r Trees/" + str_subfolder);
    }
    // create a new one of the same name (nothing happens if it already exists)
    gSystem->Exec("mkdir -p Results/" + str_subfolder);
    gSystem->Exec("mkdir -p Trees/" + str_subfolder);

    Printf("*** Analysis initiated successfully. ***");
    Printf("*** Results will be stored in the folder: Results/%s ***", str_subfolder.Data());

    return;
}