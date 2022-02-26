// RunAnalysis.C
// David Grund, Feb 26, 2022
// Run the whole analysis

// my headers
#include "AnalysisConfig.h"
#include "AnalysisManager.h"
#include "CountEvents.h"
#include "CountEvents_MC.h"
#include "GetTriggerCounters.h"
#include "IntegratedLuminosity.h"

const Int_t nSteps = 7;
Bool_t AnalysisStepsDone[nSteps+1] = { kFALSE };
Int_t iProgress = 0;
ofstream progress_file;

void CheckProgressFile(){

    ifstream ifs;
    ifs.open(("Results/" + str_subfolder + "progress_file.txt").Data());
    // if the file already exists, then check the progress
    if(!ifs.fail() && !START_FROM_CLEAN){
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
        ifs.close(); 
    } else {
        // file doesn't exist
        return;
    }
    // file exists and progress was loaded
    return;
}

void UpdateProgressFile(){

    progress_file << iProgress << "\t" << "1" << "\n";
    iProgress++;

    return;
};

void RunAnalysis(){

    // Initialize the analysis
    InitAnalysis(iWhichAnalysisToRun);

    // check if the progress file already exists, if yes, load progress
    CheckProgressFile();

    // create a new progress file
    progress_file.open(("Results/" + str_subfolder + "progress_file.txt").Data());
    UpdateProgressFile();

    // 1) Count events (data)
    if(!AnalysisStepsDone[1] || START_FROM_CLEAN) CountEvents_main();
    UpdateProgressFile();

    // 2) Count events (MC)
    if(!AnalysisStepsDone[2] || START_FROM_CLEAN) CountEvents_MC_main();
    UpdateProgressFile();

    // 3) Ge trigger counters for both periods
    if(!AnalysisStepsDone[3] || START_FROM_CLEAN) GetTriggerCounters_main();
    UpdateProgressFile();

    // 4) Calculate the integrated luminosity
    if(!AnalysisStepsDone[4] || START_FROM_CLEAN) IntegratedLuminosity_main();
    UpdateProgressFile();

    /*
    // 5) Invariant mass fits of coh, inc, all and allbins
    if(!AnalysisStepsDone[5] || START_FROM_CLEAN);
    UpdateProgressFile();

    // 6) Calculating widths of pT bins
    if(!AnalysisStepsDone[6] || START_FROM_CLEAN);
    UpdateProgressFile();

    // 7) Invariant mass fits in pT bins
    if(!AnalysisStepsDone[7] || START_FROM_CLEAN);
    UpdateProgressFile();
    */

    // Save the progress to the file
    progress_file.close();
    Printf("*** Analysis status printed to %s.***", ("Results/" + str_subfolder + "progress_file.txt").Data());

    Printf("\n******************");
    Printf("*** FINISHED. ****");
    Printf("******************");

    return;
}