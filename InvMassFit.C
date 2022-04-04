// InvMassFit.C
// David Grund, Mar 20, 2022
// To perform fit of the invariant mass distribution of measured data

// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "InvMassFit_Utilities.h"

// Main function
void InvMassFit_SetFit(Int_t opt);
// Support function
// see InvMassFit_Utilities.h

void InvMassFit(Int_t iAnalysis, Int_t optMain)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "InvMassFit/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit/");

    if(optMain == 0)
    {
        InvMassFit_PrepareData(0);
        // inc
        InvMassFit_SetFit(0);
        // coh
        InvMassFit_SetFit(1);
        // all
        InvMassFit_SetFit(2);
        // allbins
        InvMassFit_SetFit(3);
    }

    if(optMain == 1)
    {
        SetPtBinning();
        // in bins
        InvMassFit_SetFit(4);
        InvMassFit_SetFit(5);
        InvMassFit_SetFit(6);
        InvMassFit_SetFit(7);
        if(nPtBins == 5) InvMassFit_SetFit(8);
    }

    Printf("Done.");

    return;
}

void InvMassFit_SetFit(Int_t opt)
{
    // Inv mass cuts
    Double_t fMCutLow   = 2.2;
    Double_t fMCutUpp   = 4.5;

    // CB tail parameters: load MC values
    Double_t fAlpha_L;
    Double_t fAlpha_R;
    Double_t fN_L;
    Double_t fN_R;

    char name[20];
    Double_t values[4];
    Double_t errors[4];

    TString* path = new TString("Results/" + str_subfolder + "InvMassFit_MC/");
    switch(opt){
        case 0: // 'inc': incoherent-enriched sample
            path->Append("inc/inc.txt");
            break;
        case 1: // 'coh': coherent-enriched sample
            path->Append("coh/coh.txt");
            break;
        case 2: // 'all': total sample (pT < 2.0 GeV/c)
            path->Append("all/all.txt");
            break;
        case 3: // 'allbins': sample with pT from 0.2 to 1 GeV/c 
            path->Append("allbins/allbins.txt");
            break;
        case 4: // pT bin 1
            if(nPtBins == 4) path->Append("4bins/bin1.txt");
            if(nPtBins == 5) path->Append("5bins/bin1.txt");
            break;
        case 5: // pT bin 2
            if(nPtBins == 4) path->Append("4bins/bin2.txt");
            if(nPtBins == 5) path->Append("5bins/bin2.txt");
            break;
        case 6: // pT bin 3
            if(nPtBins == 4) path->Append("4bins/bin3.txt");
            if(nPtBins == 5) path->Append("5bins/bin3.txt");
            break;
        case 7: // pT bin 4
            if(nPtBins == 4) path->Append("4bins/bin4.txt");
            if(nPtBins == 5) path->Append("5bins/bin4.txt");
            break;
        case 8: // pT bin 5
            if(nPtBins == 5) path->Append("5bins/bin5.txt");
            break;
    }

    ifstream ifs;
    ifs.open(path->Data());
    if(ifs.fail()){
        Printf("\n");
        Printf("*** Warning! ***");
        Printf("*** MC values for tail parameters not found. Terminating... *** \n");
        return;
    } else {
        Int_t i_line = 0;
        while(!ifs.eof()){
            ifs >> name >> values[i_line] >> errors[i_line];
            i_line++;
        }
    }
    ifs.close();

    fAlpha_L = values[0];
    fAlpha_R = values[1];
    fN_L = values[2];
    fN_R = values[3];

    // Prepare the output path
    TString names[4] = {"inc", "coh", "all", "allbins"};
    TString str_out = "Results/" + str_subfolder + "InvMassFit/";

    if(opt == 0 || // inc
       opt == 1 || // coh
       opt == 2 || // all
       opt == 3)   // allbins
    {
        gSystem->Exec("mkdir -p Results/" + str_subfolder + "InvMassFit/" + names[opt]);
        str_out = str_out + names[opt] + "/" + names[opt];
    }
    if(opt > 3) // in pT bins
    {
        gSystem->Exec("mkdir -p Results/" + str_subfolder + Form("InvMassFit/%ibins/", nPtBins));
        str_out = str_out + Form("%ibins/bin%i", nPtBins, opt-3); 
    }   

    InvMassFit_DoFit(opt, fMCutLow, fMCutUpp, fAlpha_L, fAlpha_R, fN_L, fN_R, str_out);

    return;
}