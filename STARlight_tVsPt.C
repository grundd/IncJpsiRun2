// STARlight_tVsPt.C
// David Grund, May 24, 2022

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
#include "_STARlight_Utilities.h"
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

Int_t nGenEv = 6e6;
Int_t nBins = 1000;

void PlotResults(Double_t pT2_min, Double_t pT2_max); // pT2 in [GeV^2]
void CalculateAvgTPerBin();
void CorrectionPt2ToT(Int_t opt_err_bars);

void STARlight_tVsPt(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "STARlight_tVsPt/");

    PlotResults(0.00, 2.56);

    PlotResults(0.04, 1.00);

    CalculateAvgTPerBin();

    CorrectionPt2ToT(0);

    CorrectionPt2ToT(1);

    return;
}

void PlotResults(Double_t pT2_min, Double_t pT2_max)
{
    TFile *f = TFile::Open("Trees/STARlight/IncJ_tVsPt/tree_tPtGammaVMPom.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TTree *tPtGammaVMPom = dynamic_cast<TTree*> (f->Get("tPtGammaVMPom"));
    if(tPtGammaVMPom) Printf("Input tree loaded.");

    ConnectTreeVariables_tPtGammaVMPom(tPtGammaVMPom);

    TH2D *Hist = new TH2D("Hist", "#it{t} vs #it{p}_{T}^{2} of J/#psi", nBins, pT2_min, pT2_max, nBins, pT2_min, pT2_max);
    // on a horizontal axis: J/psi transverse momentum squared
    // on a vertical axis: Mandelstam |t|

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++)
    {
        tPtGammaVMPom->GetEntry(iEntry);
        Hist->Fill(fPtVM*fPtVM, fPtPm*fPtPm);
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

    TString path_out = "Results/" + str_subfolder + "STARlight_tVsPt/" + Form("2dhist_%.2f-%.2f", pT2_min, pT2_max);

    c->Print((path_out + ".pdf").Data());
    c->Print((path_out + ".png").Data());

    return;
}

void CalculateAvgTPerBin()
// calculate the average value of |t| (p_T,pom^2) and the average value of p_T,J/psi^2 in each bin as predicted by STARlight
{
    TFile *f = TFile::Open("Trees/STARlight/IncJ_tVsPt/tree_tPtGammaVMPom.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TTree *tPtGammaVMPom = dynamic_cast<TTree*> (f->Get("tPtGammaVMPom"));
    if(tPtGammaVMPom) Printf("Input tree loaded.");

    ConnectTreeVariables_tPtGammaVMPom(tPtGammaVMPom);

    Double_t nPt2VMPerBin[5] = { 0 };
    // to calculate average |t|:
    Double_t SumOfTPerBin[5] = { 0 };
    Double_t AvgOfTPerBin[5] = { 0 };
    // to calculate average p_T,J/psi^2
    Double_t SumOfPt2VMPerBin[5] = { 0 };
    Double_t AvgOfPt2VMPerBin[5] = { 0 };

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++)
    {
        tPtGammaVMPom->GetEntry(iEntry);
        for(Int_t iBin = 0; iBin < nPtBins; iBin++)
        {
            if(fPtVM > ptBoundaries[iBin] && fPtVM <= ptBoundaries[iBin + 1])
            {
                nPt2VMPerBin[iBin]++;
                SumOfTPerBin[iBin] += fPtPm * fPtPm;
                SumOfPt2VMPerBin[iBin] += fPtVM * fPtVM;
            }
        }
    }

    TString str_1 = Form("Results/%sSTARlight_tVsPt/AvgTPerBin.txt", str_subfolder.Data());
    ofstream outfile_1(str_1.Data());
    outfile_1 << std::fixed << std::setprecision(6);
    TString str_2 = Form("Results/%sSTARlight_tVsPt/AvgPt2VMPerBin.txt", str_subfolder.Data());
    ofstream outfile_2(str_2.Data());
    outfile_2 << std::fixed << std::setprecision(6);

    // Calculate the average values
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        // to calculate average |t|:
        AvgOfTPerBin[iBin] = SumOfTPerBin[iBin] / nPt2VMPerBin[iBin];
        // to calculate average p_T,J/psi^2
        AvgOfPt2VMPerBin[iBin] = SumOfPt2VMPerBin[iBin] / nPt2VMPerBin[iBin];
        // print the results
        Printf("Bin %i: avg |t| value = %.6f, avg p_T,Jpsi^2 value = %.6f", iBin+1, AvgOfTPerBin[iBin], AvgOfPt2VMPerBin[iBin]);
        outfile_1 << iBin+1 << "\t" << AvgOfTPerBin[iBin] << "\n";
        outfile_2 << iBin+1 << "\t" << AvgOfPt2VMPerBin[iBin] << "\n";
    }

    outfile_1.close();
    outfile_2.close();
    Printf("*** Results printed to %s. ***", str_1.Data());
    Printf("*** Results printed to %s. ***", str_2.Data());

    return;
}

void CorrectionPt2ToT(Int_t opt_err_bars)
{
    TFile *f = TFile::Open("Trees/STARlight/IncJ_tVsPt/tree_tPtGammaVMPom.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TTree *tPtGammaVMPom = dynamic_cast<TTree*> (f->Get("tPtGammaVMPom"));
    if(tPtGammaVMPom) Printf("Input tree loaded.");

    ConnectTreeVariables_tPtGammaVMPom(tPtGammaVMPom);

    Double_t* tBoundaries = NULL;
    Double_t tBoundaries_4bins[5] = { 0 };
    Double_t tBoundaries_5bins[6] = { 0 };
    
    if(nPtBins == 4)      tBoundaries = &tBoundaries_4bins[0];
    else if(nPtBins == 5) tBoundaries = &tBoundaries_5bins[0];

    for(Int_t i = 0; i < nPtBins + 1; i++) tBoundaries[i] = ptBoundaries[i]*ptBoundaries[i];

    TH1D *hEventsInT = new TH1D("hEventsInT", "hEventsInT", nPtBins, tBoundaries);
    TH1D *hEventsInPt2 = new TH1D("hEventsInPt2", "hEventsInPt2", nPtBins, tBoundaries);
    TH1D *hCorrection = NULL;

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++)
    {
        tPtGammaVMPom->GetEntry(iEntry);
        hEventsInT->Fill(fPtPm * fPtPm);
        hEventsInPt2->Fill(fPtVM * fPtVM);
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
    TString path_out = "";
    if(opt_err_bars == 0)
    {
        path_out = "Results/" + str_subfolder + "STARlight_tVsPt/CorrectionPt2ToT_errbars0";
        hCorrection->Draw("P");
    } 
    else if(opt_err_bars == 1)
    {
        path_out = "Results/" + str_subfolder + "STARlight_tVsPt/CorrectionPt2ToT_errbars1";
        hCorrection->SetFillColor(kBlack);
        hCorrection->SetFillStyle(3001);
        hCorrection->Draw("P E2");
    }
    // Legend
    TLegend *leg = new TLegend(0.12,0.76,0.40,0.96);
    leg->AddEntry((TObject*)0,Form("STARlight Simulation"),""); 
    leg->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    leg->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0); // no border
    leg->SetFillStyle(0);  // legend is transparent
    leg->Draw();
 
    c->Print((path_out + ".pdf").Data());
    c->Print((path_out + ".png").Data());

    return;
}