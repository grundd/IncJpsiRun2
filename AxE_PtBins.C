// AxE_PtBins.C
// David Grund, Mar 20, 2022
// To calculate the acceptance x efficiency from MC data in defined pt bins

// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "AxE_PtBins_Utilities.h"

void AxE_PtBins(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_PtBins/");

    AxE_PtBins_Calculate(cut_fVertexZ);

    return;
}