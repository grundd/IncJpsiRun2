// ResolutionPt.C
// David Grund, May 23, 2022

// root headers
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

Double_t res;   // (pt_rec - pt_gen) / pt_gen [-]
Double_t diff;  // (pt_rec - pt_gen) [GeV/c]
Double_t resLow = -0.5;
Double_t resUpp = 0.5;
Int_t nBins = 50;
Int_t nBins2 = 100;
Double_t BinSize = (resUpp - resLow) / (Double_t)nBins;

void CalculateFWHMAndFWTM(TH1D *h, Int_t iBin);
void CalculateResPerBin();

void ResolutionPt(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "ResolutionPt/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "ResolutionPt/");

    CalculateResPerBin();

    // Load data file
    TFile *file = TFile::Open(("Trees/" + str_subfolder + "ResolutionPt/histograms.root").Data(), "read");
    if(file) Printf("File %s loaded.", file->GetName());    
    else {
        Printf("Input file missing. Terminating...");
        return;
    }

    TList *list = (TList*) file->Get("Output");
    if(list) Printf("List %s loaded.", list->GetName()); 
    Printf("\n***\n");

    TTree *tResBins[5] = { NULL };
    TH1D *hResBins[5] = { NULL };
    TH1D *hDiffBins[5] = { NULL };

    // Loop over bins
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        tResBins[iBin] = (TTree*)list->FindObject(Form("tResBin%i", iBin+1));
        hResBins[iBin] = (TH1D*)list->FindObject(Form("hResBin%i", iBin+1));
        hDiffBins[iBin] = (TH1D*)list->FindObject(Form("hDiffBins%i", iBin+1));
        if(hResBins[iBin]){
            Printf("Tree %s found.", tResBins[iBin]->GetName());
            Printf("Tree %s contains %lli entries.", tResBins[iBin]->GetName(), tResBins[iBin]->GetEntries());
            Printf("Histogram %s found.", hResBins[iBin]->GetName());
            // Fit the resolution in bin using double-sided CB function
            //FitResInBins(hResBins[iBin], iBin+1, 2);
            Printf("\n***\n");
        } 
        if(hDiffBins[iBin]){
            // Calculate the FWHM and FWTM for each bin
            CalculateFWHMAndFWTM(hDiffBins[iBin], iBin+1);
        }  
    }

    return;
}

void CalculateFWHMAndFWTM(TH1D *h, Int_t iBin)
{
    TCanvas *c = new TCanvas("c","c",900,600);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.10);

    // Calculate FWHM
    Int_t bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
    Int_t bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
    Double_t FWHM = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
    // Calculate FWTM
    bin1 = h->FindFirstBinAbove(h->GetMaximum()/10);
    bin2 = h->FindLastBinAbove(h->GetMaximum()/10);
    Double_t FWTM = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
    // Draw histogram
    // Vertical axis
    h->GetYaxis()->SetTitle("Counts per 3 MeV");
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.95);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    h->GetXaxis()->SetTitle("#it{p}_{T}^{rec} #minus #it{p}_{T}^{gen} [GeV/#it{c}]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetLabelOffset(0.01); 
    h->GetXaxis()->SetDecimals(1);
    h->Draw();
    // Style
    h->SetLineWidth(3.);
    h->SetLineColor(215);
    // Legend
    TLegend *l1 = new TLegend(0.02,0.44,0.5,0.92);
    l1->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    l1->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l1->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.3f, %.3f) GeV/#it{c}", ptBoundaries[iBin-1], ptBoundaries[iBin]),"");
    l1->AddEntry((TObject*)0,Form("FWHM = %.0f MeV/#it{c}", FWHM * 1000),"");
    l1->AddEntry((TObject*)0,Form("FWTM = %.0f MeV/#it{c}", FWTM * 1000),"");
    l1->AddEntry((TObject*)0,Form("FWTM/FWHM = %.2f", FWTM / FWHM),"");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();

    TString str = "Results/" + str_subfolder + "ResolutionPt/" + Form("FWHM_bin%i", iBin);
    c->Print((str + ".pdf").Data());
    c->Print((str + ".png").Data());
}

void CalculateResPerBin()
{
    TFile *file = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tRec = dynamic_cast<TTree*> (file->Get(str_in_MC_tree_rec.Data()));
    if(tRec) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(tRec);

    TFile *f = new TFile("Trees/" + str_subfolder + "ResolutionPt/histograms.root","RECREATE");
    TList *l = new TList();

    TH1D *hResBins[5] = { NULL };
    TH1D *hDiffBins[5] = { NULL };
    TTree *tResBins[5] = { NULL };

    for(Int_t i = 0; i < nPtBins; i++)
    {
        hResBins[i] = new TH1D(Form("hResBin%i", i+1), Form("hResBin%i", i+1), nBins, resLow, resUpp);
        hDiffBins[i] = new TH1D(Form("hDiffBins%i", i+1), Form("hDiffBins%i", i+1), nBins2, -0.15, 0.15);
        tResBins[i] = new TTree(Form("tResBin%i", i+1), Form("tResBin%i", i+1));
        tResBins[i]->Branch("res", &res, "res/D");
        tResBins[i]->Branch("diff", &diff, "diff/D");
        l->Add(hResBins[i]);
        l->Add(hDiffBins[i]);
        l->Add(tResBins[i]);
    }

    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++)
    {
        tRec->GetEntry(iEntry);
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            if(EventPassedMCRec(0,4,iBin+1) && fPtGen > ptBoundaries[iBin] && fPtGen < ptBoundaries[iBin+1]){
                res = (fPt - fPtGen) / fPtGen;
                diff = fPt - fPtGen;
                hResBins[iBin]->Fill(res);
                hDiffBins[iBin]->Fill(diff);
                tResBins[iBin]->Fill();
            } 
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");

    l->Write("Output", TObject::kSingleKey);
    f->ls();

    return;
}