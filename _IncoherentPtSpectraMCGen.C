// IncoherentPtSpectraMCGen.C
// David Grund, Jul 1, 2022

// root headers
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
// my headers
#include "AnalysisManager.h"

void _IncoherentPtSpectraMCGen()
{
    TH1D *hJ = new TH1D("hJ","hJ",140,0.,1.4);
    TH1D *hP = new TH1D("hP","hP",140,0.,1.4);

    // get IncJ MC gen tree
    TFile *fJ = TFile::Open("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    if(fJ) Printf("MC gen file loaded.");
    TTree *tJ = dynamic_cast<TTree*> (fJ->Get("AnalysisOutput/fTreeJpsiMCGen"));
    if(tJ) Printf("MC gen tree loaded.");
    ConnectTreeVariablesMCGen(tJ);
    for(Int_t iEntry = 0; iEntry < tJ->GetEntries(); iEntry++)
    {
        tJ->GetEntry(iEntry);
        hJ->Fill(fPtGen);
    }

    // get IncP MC gen tree
    TFile *fP = TFile::Open("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
    if(fP) Printf("MC gen file loaded.");
    TTree *tP = dynamic_cast<TTree*> (fP->Get("AnalysisOutput/fTreeJpsiMCGen"));
    if(tP) Printf("MC gen tree loaded.");
    ConnectTreeVariablesMCGen(tP);
    for(Int_t iEntry = 0; iEntry < tP->GetEntries(); iEntry++)
    {
        tP->GetEntry(iEntry);
        hP->Fill(fPtGen);
    }

    // normalize and plot histograms
    hJ->Scale(1./hJ->Integral());
    hP->Scale(1./hP->Integral());
    // calculate the ratio
    TH1D *hR = (TH1D*)hP->Clone("hR");
    hR->SetTitle("hR");
    hR->Sumw2();
    hR->Divide(hJ);

    TCanvas *c = new TCanvas("c","c",900,600);
    c->cd();
    hR->Draw();
    gSystem->Exec("mkdir -p Results/_IncoherentPtSpectraMCGen/");
    c->Print("Results/_IncoherentPtSpectraMCGen/RatioPtPsi2sToJpsi.pdf");

    return;
}