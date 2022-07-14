// _Guillermo_BinsPtFit.C
// To create histograms with inv mass distributions in pT bins used in the pT fit
// July 14, 2022

#include <fstream>
#include <stdio.h> // printf
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TH1.h"

const Int_t nBins = 42;
Double_t PtBins[nBins+1] = { 0 };
TH1D *h[nBins] = { NULL };
Double_t fPt, fM, fY;

void _Guillermo_BinsPtFit()
{
    // load the boundaries
    ifstream ifs;
    ifs.open("Results/5bins_pass3/PtFit_SubtractBkg/bins_defined.txt");
    Printf("Loaded boundaries:");
    for(Int_t i = 0; i <= nBins; i++)
    {
        ifs >> PtBins[i];
        Printf("%.3f", PtBins[i]);
    }
    ifs.close();
    // open the file with invariant mass tree
    TFile *f = TFile::Open("Trees/5bins_pass3/PtFit/PtFit.root","read");
    TTree *t = NULL;
    f->GetObject("tPtFit",t);
    if(t) Printf("Tree %s loaded.", t->GetName());
    t->SetBranchAddress("fM", &fM);
    t->SetBranchAddress("fPt", &fPt);
    // create histograms 
    TList *l = new TList();
    for(Int_t i = 0; i < nBins; i++)
    {
        h[i] = new TH1D(Form("hBin%i",i),Form("hBin%i: pT %.3f to %.3f GeV",i,PtBins[i],PtBins[i+1]), 115, 2.2, 4.5);
        l->Add(h[i]);
    } 
    // fill them
    Printf("%lli entries found in the tree.", t->GetEntries());
    Int_t nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++)
    {
        t->GetEntry(iEntry);
        for(Int_t i = 0; i < nBins; i++)
        {
            if(fPt > PtBins[i] && fPt <= PtBins[i+1]) h[i]->Fill(fM);
        }
        if((iEntry+1) % 1000 == 0)
        {
            nEntriesAnalysed += 1000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    // save the histograms to a root file
    gSystem->Exec("mkdir -p Trees/_Guillermo_BinsPtFit/");
    TFile *f_out = new TFile("Trees/_Guillermo_BinsPtFit/InvMassDistribution_PtBins.root","RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f_out->ls();
    f_out->Close();

    return;
}