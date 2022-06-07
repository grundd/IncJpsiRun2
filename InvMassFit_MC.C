// InvMassFit_MC.C
// David Grund, Jun 07, 2022
// To perform fit of the invariant mass distribution of MC data

// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "InvMassFit_MC_Utilities.h"

// Main functions
void InvMassFit_MC_SetFit(Int_t opt);
// Support function
// see InvMassFit_MC_Utilities.h

void InvMassFit_MC(Int_t iAnalysis, Int_t optMain)
{
    InitAnalysis(iAnalysis);

    if(optMain == 0)
    {
        InvMassFit_MC_PrepareData();
        // inc
        InvMassFit_MC_SetFit(0);
        // coh
        InvMassFit_MC_SetFit(1);
        // all
        InvMassFit_MC_SetFit(2);
        // allbins
        InvMassFit_MC_SetFit(3);
    }

    if(optMain == 1)
    {
        SetPtBinning();
        // in bins
        InvMassFit_MC_SetFit(4);
        InvMassFit_MC_SetFit(5);
        InvMassFit_MC_SetFit(6);
        InvMassFit_MC_SetFit(7);
        if(nPtBins == 5) InvMassFit_MC_SetFit(8);
    }

    Printf("Done.");

    return;
}

void InvMassFit_MC_SetFit(Int_t opt)
{
    // Prepare path
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/");
    TString str = "Results/" + str_subfolder + "InvMassFit_MC/";

    switch(opt){
        case 0:
            gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/inc/");
            str = str + "inc/inc";
            break;
        case 1:
            gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/coh/");
            str = str + "coh/coh";
            break;
        case 2:
            gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/all/");
            str = str + "all/all";
            break;
        case 3:
            gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/allbins/");
            str = str + "allbins/allbins";
            break;
        case 4:
            if(nPtBins == 4) str = str + "4bins/bin1"; 
            if(nPtBins == 5) str = str + "5bins/bin1"; 
            break;
        case 5:
            if(nPtBins == 4) str = str + "4bins/bin2"; 
            if(nPtBins == 5) str = str + "5bins/bin2"; 
            break;
        case 6:
            if(nPtBins == 4) str = str + "4bins/bin3"; 
            if(nPtBins == 5) str = str + "5bins/bin3"; 
            break;
        case 7:
            if(nPtBins == 4) str = str + "4bins/bin4"; 
            if(nPtBins == 5) str = str + "5bins/bin4"; 
            break;
        case 8:
            if(nPtBins == 5) str = str + "5bins/bin5"; 
            break;
    }

    if(opt > 3 && nPtBins == 4) gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/4bins/");
    if(opt > 3 && nPtBins == 5) gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit_MC/5bins/");

    InvMassFit_MC_DoFit(opt, str);

    return;
}