// PtFit_SystUncertainties.C
// David Grund, June 21, 2022

#include "PtFit_Utilities.h"

// #############################################################################################

void PtFit_SystUncertainties(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning_PtFit();
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_SystUncertainties/");

    // Try various values of R
    for(Int_t ifD = -3; ifD <= 3; ifD++) PtFit_NoBkg_DoFit(4,ifD);
    
    return;
}