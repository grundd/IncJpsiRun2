// PtFit_NoBkg.C
// David Grund, Mar 20, 2022

#include "PtFit_Utilities.h"

// #############################################################################################

void PtFit_NoBkg(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning_PtFit();
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_NoBkg/");

    // CohJ created from official MC sample kCohJpsiToMu
    PtFit_NoBkg_DoFit(0);

    // CohJ: fit using a "Gaussian shape" (free parameter b)
    PtFit_NoBkg_DoFit(2);

    // CohJ: fit using a pure STARlight formfactor (free parameter R_A)
    PtFit_NoBkg_DoFit(3);

    // CohJ: to find the optimal value of R_A
    for(Int_t i = 1001; i < 1014; i++) PtFit_NoBkg_DoFit(i);
    
    // Fit using the optimal value of R_A
    PtFit_NoBkg_DoFit(4);

    return;
}