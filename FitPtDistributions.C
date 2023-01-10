// FitPtDistributions.C
// David Grund, Jan 05, 2023

// cpp headers
#include <fstream>
// root headers
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TSystem.h"
#include "TStyle.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

Int_t iCount(0);

void TH1_SetStyle(TH1F* h, Color_t c, Int_t style = 1)
{
    h->SetLineColor(c);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
}

void TF1_SetStyle(TF1 *f, Color_t c, Int_t style = 1)
{
    f->SetLineColor(c);
    f->SetLineWidth(2);
    f->SetLineStyle(style);
}

TF1* PtShape(TString name, Float_t par0, Float_t par1)
{
    TF1* f = new TF1(name.Data(),"[0] * abs(x) * exp(-[1] * pow(x,2))",0.,2.); // x = pT
    f->SetParameter(0,par0);
    f->SetParameter(1,par1);
    return f;
}

TF1* Pt2Shape(TString name, Float_t par0, Float_t par1)
{
    TF1* f = new TF1(name.Data(),"[0] * exp(-[1] * x)",0.,4.); // x = pT^2
    f->SetParameter(0,par0);
    f->SetParameter(1,par1);
    return f;
}

void PlotHistos(TString name, TString opt, Bool_t log, TH1F* h1, TF1* f1 = NULL, TH1F* h2 = NULL, TF1* f2 = NULL)
{
    Bool_t vsPt = kFALSE;
    if(name.Contains("_pT")) vsPt = kTRUE;
    TCanvas c("c","c",700,600);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.03);
    c.SetTopMargin(0.085);
    c.SetBottomMargin(0.12);
    if(log) c.SetLogy();
    TH1_SetStyle(h1,kBlue);
    // x-axis
    if(vsPt) h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    else     h1->GetXaxis()->SetTitle("#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})");
    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetDecimals(1);
    // y-axis
    if(log) {
        Float_t minimum = h1->GetMinimum();
        if(minimum < 0.1) minimum = 0.1;
        h1->GetYaxis()->SetRangeUser(0.5*minimum,1.5*h1->GetMaximum());
    }
    h1->GetYaxis()->SetTitle("value");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetMaxDigits(3);
    // draw it
    h1->Draw(Form("%s",opt.Data()));
    if(f1) {
        TF1_SetStyle(f1,kGreen+1,9);
        f1->Draw("SAME"); 
        TLegend* l = new TLegend(0.50,0.78,0.80,0.92);
        if(vsPt) l->AddEntry(f1,"fit: #it{N}#it{p}_{T} exp(#minus#it{b}#it{p}_{T}^{2})","L");
        else     l->AddEntry(f1,"fit: #it{N} exp(#minus#it{b}#it{p}_{T}^{2})","L");
        l->AddEntry((TObject*)0,Form("#it{N} = %.0f #pm %.0f",f1->GetParameter(0),f1->GetParError(0)),"");
        l->AddEntry((TObject*)0,Form("#it{b} = %.3f #pm %.3f",f1->GetParameter(1),f1->GetParError(1)),"");
        l->SetTextSize(0.03);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetMargin(0.25);
        l->Draw();
    }
    if(h2) { 
        TH1_SetStyle(h2,kRed);
        h2->Draw(Form("%s SAME",opt.Data()));
    }
    if(f2) {
        TF1_SetStyle(f2,kViolet,9);
        f1->Draw("SAME"); 
    }
    c.Print(("Results/" + str_subfolder + "FitPtDistributions/" + name + ".pdf").Data());
    return;
}

void FitMC()
{
    TFile *fRec = TFile::Open("Trees/AnalysisDataMC_pass3/PIDCalibrated/AnalysisResults_MC_kIncohJpsiToMu.root","read");
    if(fRec) Printf("MC rec file loaded.");
    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJpsi"));
    if(tRec) Printf("MC rec tree loaded.");
    ConnectTreeVariablesMCRec(tRec);

    TFile *fGen = TFile::Open("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root","read");
    if(fGen) Printf("MC gen file loaded.");
    TTree *tGen = dynamic_cast<TTree*> (fGen->Get("AnalysisOutput/fTreeJpsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    ConnectTreeVariablesMCGen(tGen);

    TH1F* hRec_pT = new TH1F("hRec_pT","#it{N}_{rec} vs #it{p}_{T}",100,0.,2.); // vs pT
    TH1F* hGen_pT = new TH1F("hGen_pT","#it{N}_{gen} vs #it{p}_{T}",100,0.,2.);
    TH1F* hRec = new TH1F("hRec","#it{N}_{rec} vs #it{p}_{T}^{2}",100,0.,4.); // vs pT2 (default)
    TH1F* hGen = new TH1F("hGen","#it{N}_{gen} vs #it{p}_{T}^{2}",100,0.,4.);

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) {
        tRec->GetEntry(iEntry);
        // m between 2.2 and 4.5 GeV/c^2 and any pT
        if(EventPassedMCRec(0, -1)) { 
            hRec_pT->Fill(fPt);
            hRec->Fill(fPt*fPt);
        }
    }
    TF1* fRecPt = PtShape("fRecPt",2.0e4,5.33);
    hRec_pT->Fit(fRecPt);
    PlotHistos("hRec_pT","HIST",kTRUE,hRec_pT,fRecPt);
    TF1* fRecPt2 = Pt2Shape("fRecPt2",2.0e4,5.33);
    hRec->Fit(fRecPt2);
    PlotHistos("hRec","HIST",kTRUE,hRec,fRecPt2);

    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++) {
        tGen->GetEntry(iEntry);
        // m between 2.2 and 4.5 GeV/c^2 and any pT
        if(EventPassedMCGen(-1)) { 
            hGen_pT->Fill(fPtGen);
            hGen->Fill(fPtGen*fPtGen);
        }
    }
    TF1* fGenPt = PtShape("fGenPt",6.1e5,3.99);
    hGen_pT->Fit(fGenPt);
    PlotHistos("hGen_pT","HIST",kTRUE,hGen_pT,fGenPt);
    TF1* fGenPt2 = Pt2Shape("fGenPt2",6.1e5,3.99);
    hGen->Fit(fGenPt2);
    PlotHistos("hGen","HIST",kTRUE,hGen,fGenPt2);
    return;
}

void FitData()
{
    TFile *f = TFile::Open("Trees/AnalysisData_pass3/AnalysisResults.root","read");
    if(f) Printf("Data file loaded.");
    TTree *t = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJpsi"));
    if(t) Printf("Data tree loaded.");
    ConnectTreeVariables(t);

    TH1F* hData_pT = new TH1F("hData_pT","#it{N}_{data} vs #it{p}_{T}",100,0.2,1.0);
    TH1F* hData = new TH1F("hData","#it{N}_{data} vs #it{p}_{T}^{2}",100,0.04,1.0);

    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++) {
        t->GetEntry(iEntry);
        // m between 3.0 and 3.2 GeV/c^2 and any pT
        if(EventPassed(1, -1)) {
            hData_pT->Fill(fPt);
            hData->Fill(fPt*fPt);
        }
    }
    TF1* fDataPt = PtShape("fDataPt",500.,5.33);
    hData_pT->Fit(fDataPt);
    PlotHistos("hData_pT","HIST",kTRUE,hData_pT,fDataPt);
    TF1* fDataPt2 = Pt2Shape("fDataPt2",500.,5.33);
    hData->Fit(fDataPt2);
    PlotHistos("hData","HIST",kFALSE,hData,fDataPt2);
    return;
}

void FitPtDistributions(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "FitPtDistributions/");

    // fit MC
    if(kTRUE) FitMC();

    // fit data
    if(kTRUE) FitData();

    return;
}