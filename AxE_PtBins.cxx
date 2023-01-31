// AxE_PtBins.cxx
// David Grund, Mar 20, 2022
// To calculate the acceptance x efficiency from MC data in defined pt bins

// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "AxE_Utilities.h"

void AxE_PtBins(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_PtBins/");

    // original AxE
    AxE_PtBins_Calculate(kFALSE,cut_fVertexZ);

    // re-weighted AxE
    AxE_PtBins_Calculate(kTRUE,cut_fVertexZ);

    return;
}