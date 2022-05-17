// MigrationPtRecGen.C
// David Grund, May 17, 2022

// root headers
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

// vs pT or |t|, print absolute numbers or percentages
void PlotHistBins(Bool_t vsPt, Bool_t inPercent);

void MigrationPtRecGen(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "MigrationPtRecGen/");

    PlotHistBins(kTRUE, kTRUE);

    //PlotHistBins(kFALSE, kTRUE);

    return;
}

void PlotHistBins(Bool_t vsPt, Bool_t inPercent)
{
    // Load data
    TFile *f = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
    if(f) Printf("Input data loaded.");

    TTree *t = dynamic_cast<TTree*> (f->Get(str_in_MC_tree_rec.Data()));
    if(t) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(t);

    Printf("%lli entries found in the tree.", t->GetEntries());
    Int_t nEntriesAnalysed = 0;

    // Create 2D histogram
    // pT gen on the horizontal axis, pT rec on the vertical axis
    Double_t* boundaries_pT;
    Double_t* boundaries_pT2;
    Double_t boundaries_pT_4[7] = { 0. };
    Double_t boundaries_pT_5[8] = { 0. };
    Double_t boundaries_pT2_4[7] = { 0. };
    Double_t boundaries_pT2_5[8] = { 0. };
    if(nPtBins == 4)
    {
        boundaries_pT = &boundaries_pT_4[0];
        boundaries_pT2 = &boundaries_pT2_4[0];
    } 
    else if(nPtBins == 5) 
    {
        boundaries_pT = &boundaries_pT_5[0];
        boundaries_pT2 = &boundaries_pT2_5[0];
    }
    for(Int_t iBin = 1; iBin <= nPtBins+1; iBin++)
    {
        boundaries_pT[iBin] = ptBoundaries[iBin-1];
        boundaries_pT2[iBin] = ptBoundaries[iBin-1]*ptBoundaries[iBin-1];
    }
    boundaries_pT[nPtBins+2] = 1.2;
    boundaries_pT2[nPtBins+2] = 1.2;

    // 2d histograms with pT bins
    TH2D *hMigration = NULL;
    TH1D *hScaleByTotalNRec = NULL;
    if(vsPt){
        hMigration = new TH2D("hMigration", "pt gen vs pt rec", nPtBins+2, boundaries_pT, nPtBins+2, boundaries_pT);
        hScaleByTotalNRec = new TH1D("hScaleByTotalNRec", "Total NRec per bin in pt gen", nPtBins+2, boundaries_pT);

    } else {
        hMigration = new TH2D("hMigration", "pt^{2} gen vs pt^{2} rec", nPtBins+2, boundaries_pT2, nPtBins+2, boundaries_pT2);
        hScaleByTotalNRec = new TH1D("hScaleByTotalNRec", "Total NRec per bin in pt gen", nPtBins+2, boundaries_pT2);
    } 
    // 2d histogram showing only pT gen vs. pT gen without using the main binning
    TH2D *hPtGenVsRec = new TH2D("hPtGenVsRec","hPtGenVsRec",360,0.0,1.2,240,0.0,1.2);
    // 1d histograms for pT distributions
    TH1D *hDistPtRec = new TH1D("hDistPtRec","hDistPtRec",240,0.0,1.2);
    TH1D *hDistPtGen = new TH1D("hDistPtGen","hDistPtGen",240,0.0,1.2);

    // Fill the histograms
    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++)
    {
        t->GetEntry(iEntry);

        // mass between 3.0 and 3.2, no pT cut
        if(EventPassedMCRec(0,-1))
        {
            if(vsPt){
                hMigration->Fill(fPtGen,fPt);
                hScaleByTotalNRec->Fill(fPtGen);
            } else {
                hMigration->Fill(fPtGen*fPtGen,fPt*fPt);
                hScaleByTotalNRec->Fill(fPtGen*fPtGen);
            }
            hPtGenVsRec->Fill(fPtGen,fPt);
            hDistPtRec->Fill(fPt);
            hDistPtGen->Fill(fPtGen);
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    // Scale the histogram
    for(Int_t iBinX = 1; iBinX <= nPtBins+2; iBinX++){
        Double_t NRec_tot = hScaleByTotalNRec->GetBinContent(iBinX);
        Printf("BinX %i: %.0f", iBinX, NRec_tot);
        for(Int_t iBinY = 1; iBinY <= nPtBins+2; iBinY++){
            Double_t ValueScaled;
            if(inPercent) ValueScaled = hMigration->GetBinContent(iBinX,iBinY) / NRec_tot * 100;
            else ValueScaled = hMigration->GetBinContent(iBinX,iBinY) / NRec_tot;
            hMigration->SetBinContent(iBinX,iBinY,ValueScaled);
        }
    }

    Printf("Number of entries in the histogram: %.0f", hMigration->GetEntries());

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    if(inPercent) gStyle->SetPaintTextFormat("4.2f");
    else gStyle->SetPaintTextFormat("4.3f");

    //****************************************************************************************
    // plot 2d histograms with pT bins
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    if(!vsPt){
        c1->SetLogx();
        c1->SetLogy();
    }
    c1->SetTopMargin(0.03);
    c1->SetBottomMargin(0.145);
    c1->SetRightMargin(0.1);
    c1->SetLeftMargin(0.075);

    hMigration->SetMarkerSize(1.8);
    // horizontal axis
    if(vsPt){
        hMigration->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
        hMigration->GetXaxis()->SetLabelOffset(0.015);
    } else {
        hMigration->GetXaxis()->SetTitle("(#it{p}_{T}^{gen})^{2} (GeV^{2}/#it{c}^{2})");
    }   
    hMigration->GetXaxis()->SetTitleSize(0.05);
    hMigration->GetXaxis()->SetTitleOffset(1.3);
    hMigration->GetXaxis()->SetLabelSize(0.05);
    hMigration->GetXaxis()->SetDecimals(1);
    // vertical axis
    if(vsPt) hMigration->GetYaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    else     hMigration->GetYaxis()->SetTitle("(#it{p}_{T}^{rec})^{2} (GeV^{2}/#it{c}^{2})");
    hMigration->GetYaxis()->SetTitleSize(0.05);
    hMigration->GetYaxis()->SetLabelSize(0.05);
    hMigration->GetYaxis()->SetTitleOffset(0.65);
    hMigration->GetYaxis()->SetDecimals(1);
    // Set ranges
    hMigration->GetXaxis()->SetRangeUser(0.0,1.2);
    hMigration->GetYaxis()->SetRangeUser(0.0,1.2);
    // Z-axis
    hMigration->GetZaxis()->SetLabelSize(0.05);
    hMigration->Draw("COLZ TEXT");


    TString path_out1 = "Results/" + str_subfolder + "MigrationPtRecGen/";
    if(vsPt) path_out1 += "migr_bins_pT";
    else     path_out1 += "migr_bins_pT2";
    c1->Print((path_out1 + ".pdf").Data());
    c1->Print((path_out1 + ".png").Data());

    //****************************************************************************************
    // plot 2d histogram showing only pT gen vs. pT gen without using the main binning
    TCanvas *c2 = new TCanvas("c2","c2",900,600);
    c2->SetTopMargin(0.03);
    c2->SetBottomMargin(0.145);
    c2->SetRightMargin(0.1);
    c2->SetLeftMargin(0.075);

    // horizontal axis
    hPtGenVsRec->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hPtGenVsRec->GetXaxis()->SetLabelOffset(0.015);
    hPtGenVsRec->GetXaxis()->SetTitleSize(0.05);
    hPtGenVsRec->GetXaxis()->SetTitleOffset(1.3);
    hPtGenVsRec->GetXaxis()->SetLabelSize(0.05);
    hPtGenVsRec->GetXaxis()->SetDecimals(1);
    // vertical axis
    hPtGenVsRec->GetYaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    hPtGenVsRec->GetYaxis()->SetTitleSize(0.05);
    hPtGenVsRec->GetYaxis()->SetLabelSize(0.05);
    hPtGenVsRec->GetYaxis()->SetTitleOffset(0.65);
    hPtGenVsRec->GetYaxis()->SetDecimals(1);
    // Set ranges
    hPtGenVsRec->GetXaxis()->SetRangeUser(0.0,1.2);
    hPtGenVsRec->GetYaxis()->SetRangeUser(0.0,1.2);
    // Z-axis
    hPtGenVsRec->GetZaxis()->SetLabelSize(0.05);
    hPtGenVsRec->Draw("COLZ");

    TString path_out2 = "Results/" + str_subfolder + "MigrationPtRecGen/migr_hist";
    c2->Print((path_out2 + ".pdf").Data());
    c2->Print((path_out2 + ".png").Data());

    //****************************************************************************************
    // plot 1d distributions
    TCanvas *c3 = new TCanvas("c3","c3",900,600);
    c3->SetTopMargin(0.03);
    c3->SetBottomMargin(0.145);
    c3->SetRightMargin(0.03);
    c3->SetLeftMargin(0.13);

    // scale the histograms
    hDistPtRec->Scale(1.0/hDistPtRec->Integral());
    hDistPtGen->Scale(1.0/hDistPtGen->Integral());
    // horizontal axis
    hDistPtRec->GetXaxis()->SetTitle("#it{p}_{T}^{gen} or #it{p}_{T}^{rec} (GeV/#it{c})");
    hDistPtRec->GetXaxis()->SetLabelOffset(0.015);
    hDistPtRec->GetXaxis()->SetTitleSize(0.05);
    hDistPtRec->GetXaxis()->SetTitleOffset(1.35);
    hDistPtRec->GetXaxis()->SetLabelSize(0.05);
    hDistPtRec->GetXaxis()->SetDecimals(1);
    // vertical axis
    hDistPtRec->GetYaxis()->SetTitle("Counts (normalized to 1.0)");
    hDistPtRec->GetYaxis()->SetTitleSize(0.05);
    hDistPtRec->GetYaxis()->SetLabelSize(0.05);
    hDistPtRec->GetYaxis()->SetTitleOffset(1.2);
    hDistPtRec->GetYaxis()->SetDecimals(1);
    // set colors
    hDistPtRec->SetLineWidth(2.0);
    hDistPtRec->SetLineColor(kRed);
    hDistPtGen->SetLineWidth(2.0);
    hDistPtGen->SetLineColor(kBlue);
    // draw
    hDistPtRec->Draw("HIST");
    hDistPtGen->Draw("HIST SAME");

    // legend
    TLegend* leg = new TLegend(0.65,0.80,0.95,0.95);
    leg->AddEntry(hDistPtRec,"distribution in #it{p}_{T}^{rec}","L");
    leg->AddEntry(hDistPtGen,"distribution in #it{p}_{T}^{gen}","L");
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    TString path_out3 = "Results/" + str_subfolder + "MigrationPtRecGen/migr_dist";
    c3->Print((path_out3 + ".pdf").Data());
    c3->Print((path_out3 + ".png").Data());

    return;
}