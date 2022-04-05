// Plot2DHistogramPtAndT.c
// David Grund, Oct 19, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"

Double_t fPt2Gm, fPt2VM, fPt2Pm;
Double_t fPtGm, fPtVM, fPtPm;
Int_t nGenEv = 500000;
Int_t nBins = 1000;
Double_t fPt2Low, fPt2Upp;

TString TxtPtGamma = "DependenceOnT/SL_simulations_10-19-2021/PtGamma.txt";
TString TxtPtVMPom = "DependenceOnT/SL_simulations_10-19-2021/PtVMpomeron.txt";

void PlotResults(Int_t opt);
void CalculateAvgTPerBin();
void CorrectionPt2ToT();
void PrepareTree();

void Plot2DHistogramPtAndT(){

    //PrepareTree();

    //PlotResults(0);

    //PlotResults(1);

    //CalculateAvgTPerBin();

    CorrectionPt2ToT();

    return;
}

void CalculateAvgTPerBin(){

    TFile *f = TFile::Open("DependenceOnT/SL_simulations_10-19-2021/tree.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TList *l = (TList*) f->Get("TreeList");
    if(l) Printf("List %s loaded.", l->GetName()); 

    TTree *tPtVMPom = (TTree*)l->FindObject("tPtVMPom");
    if(tPtVMPom) Printf("Tree %s loaded.", tPtVMPom->GetName());

    tPtVMPom->SetBranchAddress("fPt2VM", &fPt2VM);
    tPtVMPom->SetBranchAddress("fPt2Pm", &fPt2Pm);
    tPtVMPom->SetBranchAddress("fPtVM", &fPtVM);
    tPtVMPom->SetBranchAddress("fPtPm", &fPtPm);    

    SetPtBinning();

    Double_t nPt2VMPerBin[nPtBins] = { 0 };
    Double_t SumOfTPerBin[nPtBins] = { 0 };
    Double_t AvgOfTPerBin[nPtBins] = { 0 };

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++){
        tPtVMPom->GetEntry(iEntry);
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            if(fPtVM > ptBoundaries[iBin] && fPtVM <= ptBoundaries[iBin + 1]){
                nPt2VMPerBin[iBin]++;
                SumOfTPerBin[iBin] += fPt2Pm;
            }
        }
    }

    TString str = Form("DependenceOnT/output_%ibins.txt", nPtBins);
    ofstream outfile(str.Data());
    outfile << std::fixed << std::setprecision(6);

    // Calculate the average values
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        AvgOfTPerBin[iBin] = SumOfTPerBin[iBin] / nPt2VMPerBin[iBin];
        Printf("Bin %i: avg t value = %.6f", iBin+1, AvgOfTPerBin[iBin]);
        outfile << iBin+1 << "\t" << AvgOfTPerBin[iBin] << "\n";
    }

    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());

    return;
}

void PlotResults(Int_t opt){

    if(opt == 0){
        fPt2Low = 0.00; // GeV^2
        fPt2Upp = 2.56; // GeV^2
    } else if(opt == 1){
        fPt2Low = 0.04; // GeV^2
        fPt2Upp = 1.00; // GeV^2
    }

    TFile *f = TFile::Open("DependenceOnT/SL_simulations_10-19-2021/tree.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TList *l = (TList*) f->Get("TreeList");
    if(l) Printf("List %s loaded.", l->GetName()); 

    TTree *tPtGamma = (TTree*)l->FindObject("tPtGamma");
    if(tPtGamma) Printf("Tree %s loaded.", tPtGamma->GetName());

    TTree *tPtVMPom = (TTree*)l->FindObject("tPtVMPom");
    if(tPtVMPom) Printf("Tree %s loaded.", tPtVMPom->GetName());

    tPtGamma->SetBranchAddress("fPt2Gm", &fPt2Gm);
    tPtGamma->SetBranchAddress("fPtGm", &fPtGm);

    tPtVMPom->SetBranchAddress("fPt2VM", &fPt2VM);
    tPtVMPom->SetBranchAddress("fPt2Pm", &fPt2Pm);
    tPtVMPom->SetBranchAddress("fPtVM", &fPtVM);
    tPtVMPom->SetBranchAddress("fPtPm", &fPtPm);

    TH2D *Hist = new TH2D("Hist", "#it{t} vs #it{p}_{T}^{2} of J/#psi", nBins, fPt2Low, fPt2Upp, nBins, fPt2Low, fPt2Upp);
    // on a horizontal axis: J/psi transverse momentum squared, i.e. fPt2VM
    // on a vertical axis: Mandelstam t, i.e. fPt2Pm

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++){
        tPtGamma->GetEntry(iEntry);
        tPtVMPom->GetEntry(iEntry);
        //Printf("fPt2VM: %.4f, fPt2Gm+fPt2Pm: %.4f, fPt2Gm: %.4f, fPt2Pm: %.4f", fPt2VM, fPt2Gm+fPt2Pm, fPt2Gm, fPt2Pm);

        Hist->Fill(fPt2VM, fPt2Pm);
    }

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetLogz();
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.145);
    c->SetRightMargin(0.11);
    c->SetLeftMargin(0.1);

    // a vertical axis
    Hist->GetYaxis()->SetTitle("|#it{t}| or #it{p}_{T, pom}^{2} (GeV^{2}/#it{c}^{2})");
    Hist->GetYaxis()->SetTitleSize(0.05);
    Hist->GetYaxis()->SetLabelSize(0.05);
    Hist->GetYaxis()->SetTitleOffset(0.915);
    Hist->GetYaxis()->SetDecimals(1);
    // a horizontal axis
    Hist->GetXaxis()->SetTitle("#it{p}_{T, J/#psi}^{2} (GeV^{2}/#it{c}^{2})");
    Hist->GetXaxis()->SetTitleSize(0.05);
    Hist->GetXaxis()->SetTitleOffset(1.3);
    Hist->GetXaxis()->SetLabelSize(0.05);
    Hist->GetXaxis()->SetLabelOffset(0.02);
    Hist->GetXaxis()->SetDecimals(1);
    // draw the histogram
    Hist->GetZaxis()->SetLabelSize(0.05);
    Hist->Draw("COLZ");

    if(opt == 0){
        c->Print("DependenceOnT/HistPtVsT_0.png");
        c->Print("DependenceOnT/HistPtVsT_0.pdf");
    } else if(opt == 1){
        c->Print("DependenceOnT/HistPtVsT_1.png");
        c->Print("DependenceOnT/HistPtVsT_1.pdf");
    }

    return;
}

void CorrectionPt2ToT(){

    TFile *f = TFile::Open("DependenceOnT/SL_simulations_10-19-2021/tree.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TList *l = (TList*) f->Get("TreeList");
    if(l) Printf("List %s loaded.", l->GetName()); 

    TTree *tPtGamma = (TTree*)l->FindObject("tPtGamma");
    if(tPtGamma) Printf("Tree %s loaded.", tPtGamma->GetName());

    TTree *tPtVMPom = (TTree*)l->FindObject("tPtVMPom");
    if(tPtVMPom) Printf("Tree %s loaded.", tPtVMPom->GetName());

    tPtGamma->SetBranchAddress("fPt2Gm", &fPt2Gm);
    tPtGamma->SetBranchAddress("fPtGm", &fPtGm);

    tPtVMPom->SetBranchAddress("fPt2VM", &fPt2VM);
    tPtVMPom->SetBranchAddress("fPt2Pm", &fPt2Pm);
    tPtVMPom->SetBranchAddress("fPtVM", &fPtVM);
    tPtVMPom->SetBranchAddress("fPtPm", &fPtPm);

    SetPtBinning();

    Double_t tBoundaries[nPtBins + 1] = { 0 };
    for(Int_t i = 0; i < nPtBins + 1; i++) tBoundaries[i] = ptBoundaries[i]*ptBoundaries[i];

    TH1D *hEventsInT = new TH1D("hEventsInT", "hEventsInT", nPtBins, tBoundaries);
    TH1D *hEventsInPt2 = new TH1D("hEventsInPt2", "hEventsInPt2", nPtBins, tBoundaries);
    TH1D *hCorrection = NULL;

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++){
        tPtVMPom->GetEntry(iEntry);
        hEventsInT->Fill(fPt2Pm);
        hEventsInPt2->Fill(fPt2VM);
    } 

    hCorrection = (TH1D*)hEventsInPt2->Clone("hCorrection");
    hCorrection->SetTitle("hCorrection");
    hCorrection->Sumw2();
    hCorrection->Divide(hEventsInT);   

    hCorrection->SetMarkerStyle(21);
    hCorrection->SetMarkerColor(kBlue);
    hCorrection->SetMarkerSize(1.0);
    hCorrection->SetLineColor(kBlue);
    hCorrection->SetLineWidth(2.0);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetLogz();
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.145);
    c->SetRightMargin(0.02);
    c->SetLeftMargin(0.14);

    // a vertical axis
    hCorrection->GetYaxis()->SetTitle("#it{N}[#it{p}_{T,J/#psi}^{2} #in (|#it{t}|_{min}, |#it{t}|_{max})]/#it{N}[#it{p}_{T,pom}^{2} #in (|#it{t}|_{min}, |#it{t}|_{max})]");
    hCorrection->GetYaxis()->SetTitleSize(0.048);
    hCorrection->GetYaxis()->SetLabelSize(0.05);
    hCorrection->GetYaxis()->SetTitleOffset(1.4);
    hCorrection->GetYaxis()->SetDecimals(3);
    // a horizontal axis
    hCorrection->GetXaxis()->SetTitle("|#it{t}| or #it{p}_{T, pom}^{2} (GeV^{2}/#it{c}^{2})");
    hCorrection->GetXaxis()->SetTitleSize(0.05);
    hCorrection->GetXaxis()->SetTitleOffset(1.3);
    hCorrection->GetXaxis()->SetLabelSize(0.05);
    hCorrection->GetXaxis()->SetDecimals(1);
    // draw the histogram
    hCorrection->Draw("P E1");
    // Legend
    TLegend *leg = new TLegend(0.32,0.76,0.65,0.96);
    leg->AddEntry((TObject*)0,Form("STARlight Simulation"),""); 
    leg->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    leg->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    leg->SetTextSize(0.056);
    leg->SetBorderSize(0); // no border
    leg->SetFillStyle(0);  // legend is transparent
    leg->Draw();

    TString str = Form("DependenceOnT/CorrectionPt2ToT_%ibins", nPtBins);

    c->Print((str + ".pdf").Data());
    c->Print((str + ".png").Data());

    return;
}

void PrepareTree(){

    TTree *tPtGamma = new TTree("tPtGamma", "tPtGamma");
    tPtGamma->Branch("fPt2Gm", &fPt2Gm, "fPt2Gm/D");
    tPtGamma->Branch("fPtGm", &fPtGm, "fPtGm/D");

    TTree *tPtVMPom = new TTree("tPtVMPom", "tPtVMPom");
    tPtVMPom->Branch("fPt2VM", &fPt2VM, "fPt2VM/D");
    tPtVMPom->Branch("fPt2Pm", &fPt2Pm, "fPt2Pm/D");
    tPtVMPom->Branch("fPtVM", &fPtVM, "fPtVM/D");
    tPtVMPom->Branch("fPtPm", &fPtPm, "fPtPm/D");

    ifstream ifs;
    ifs.open(TxtPtGamma.Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < nGenEv; i++){
            ifs >> fPt2Gm;
            fPtGm = TMath::Sqrt(fPt2Gm);
            tPtGamma->Fill();
        }
        ifs.close();
    } else {
        Printf("File %s missing. Terminating.", TxtPtGamma.Data());
        return;
    }

    ifs.open(TxtPtVMPom.Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < nGenEv; i++){
            // see the src/eventfilewriter.cpp in the STARlight source code
            ifs >> fPt2VM >> fPt2Pm;
            fPtVM = TMath::Sqrt(fPt2VM);
            fPtPm = TMath::Sqrt(fPt2Pm);
            tPtVMPom->Fill();
        }
        ifs.close();
    } else {
        Printf("File %s missing. Terminating.", TxtPtVMPom.Data());
        return;
    }

    TList *l = new TList();
    l->Add(tPtGamma);
    l->Add(tPtVMPom);

    TFile *f = new TFile("DependenceOnT/SL_simulations_10-19-2021/tree.root","RECREATE");
    l->Write("TreeList", TObject::kSingleKey);
    f->ls();

    return;
}