// RunAnalysis.C
// David Grund, Feb 26, 2022
// Run the whole analysis

// my headers
#include "AnalysisHeaders.h"
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "CountEvents.h"
#include "CountEvents_MC.h"
#include "RunListCheck.h"
#include "GetTriggerCounters.h"
#include "IntegratedLuminosity.h"
#include "InvMassFit_MC.h"
#include "InvMassFit.h"
#include "BinsThroughMassFit.h"
#include "SetPtBinning.h"
#include "AxE_PtBins.h"
#include "PtFit_SubtractBkg.h"
#include "PtFit_PrepareMCTemplates.h"
#include "PtFit_NoBkg.h"
#include "STARlight_OptimalRA.h"

const Int_t nSteps = 15;
Bool_t AnalysisStepsDone[nSteps+1] = { kFALSE };
Int_t iProgress = 0;
ofstream *progress_file = NULL;

void CheckProgressFile(){

    ifstream ifs;
    ifs.open(("Results/" + str_subfolder + "progress_file.txt").Data());
    // if the file already exists, then check the progress
    if(!ifs.fail() && !START_FROM_CLEAN){
        // file exists and progress was loaded
        for(Int_t i = 0; i < nSteps+1; i++){
            Int_t iStep;
            Int_t isDone = 0;
            ifs >> iStep >> isDone;
            if(isDone == 1){
                AnalysisStepsDone[i] = kTRUE;
                Printf("Step %i of the analysis: old results were found, will be skipped.", i);
            } else {
                Printf("Step %i of the analysis: calculations will be performed.", i);
            }
        }
    } else {
        // file doesn't exist
    }
    ifs.close(); 
    return;
}

void UpdateProgressFile(){

    *progress_file << iProgress << "\t" << "1" << "\n";
    iProgress++;

    return;
};

void RunAnalysis(){

    // Initialize the analysis
    InitAnalysis(iWhichAnalysisToRun);

    // check if the progress file already exists, if yes, load progress
    CheckProgressFile();

    // create a new progress file
    progress_file = new ofstream(("Results/" + str_subfolder + "progress_file.txt").Data());
    UpdateProgressFile();

    // 1) Count events (data)
    if(!AnalysisStepsDone[1] || START_FROM_CLEAN) CountEvents_main();
    UpdateProgressFile();

    // 2) Count events (MC)
    if(!AnalysisStepsDone[2] || START_FROM_CLEAN) CountEvents_MC_main();
    UpdateProgressFile();

    // 3) Run list check
    if(!AnalysisStepsDone[3] || START_FROM_CLEAN) RunListCheck_main();
    UpdateProgressFile();

    // 4) Get trigger counters for both periods
    if(!AnalysisStepsDone[4] || START_FROM_CLEAN) GetTriggerCounters_main();
    UpdateProgressFile();

    // 5) Calculate the integrated luminosity
    if(!AnalysisStepsDone[5] || START_FROM_CLEAN) IntegratedLuminosity_main();
    UpdateProgressFile();

    // 6) MC invariant mass fits of coh, inc, all and allbins
    if(!AnalysisStepsDone[6] || START_FROM_CLEAN) InvMassFit_MC_main(0);
    UpdateProgressFile();

    // 7) Invariant mass fits of coh, inc, all and allbins
    if(!AnalysisStepsDone[7] || START_FROM_CLEAN) InvMassFit_main(0);
    UpdateProgressFile();

    // 8) Set pT binning via the invariant mass fitting
    if(!AnalysisStepsDone[8] || START_FROM_CLEAN) BinsThroughMassFit_main();
    UpdateProgressFile();

    // Set pT binning (must be always done as it is required in the following steps!)
    SetPtBinning_main();

    // 9) MC invariant mass fits in pT bins
    if(!AnalysisStepsDone[9] || START_FROM_CLEAN) InvMassFit_MC_main(1);
    UpdateProgressFile();

    // 10) Invariant mass fits in pT bins
    if(!AnalysisStepsDone[10] || START_FROM_CLEAN) InvMassFit_main(1);
    UpdateProgressFile();

    // 11) AxE in pT bins
    if(!AnalysisStepsDone[11] || START_FROM_CLEAN) AxE_PtBins_main();
    UpdateProgressFile();

    // 12) Create pT bins for the fit of the pT distribution and subtract background in these bins
    if(!AnalysisStepsDone[12] || START_FROM_CLEAN) PtFit_SubtractBkg_main();
    UpdateProgressFile();    

    // 13) Prepare MC templates (PDFs) to be used in pT fits
    if(!AnalysisStepsDone[13] || START_FROM_CLEAN) PtFit_PrepareMCTemplates_main();
    UpdateProgressFile();    

    // 14) Various pT fits with background subtracted
    if(!AnalysisStepsDone[14] || START_FROM_CLEAN) PtFit_NoBkg_main();
    UpdateProgressFile();  

    // 15) Find the optimal RA (for which there is a minimum in chi2 of the pT fit)
    if(!AnalysisStepsDone[15] || START_FROM_CLEAN) STARlight_OptimalRA_main();
    UpdateProgressFile();  

    // Save the progress to the file
    progress_file->close();
    Printf("*** Analysis status printed to %s.***", ("Results/" + str_subfolder + "progress_file.txt").Data());

    Printf("\n******************");
    Printf("*** FINISHED. ****");
    Printf("******************");

    return;
}