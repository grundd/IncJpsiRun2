// STARlight_tVsPt2.cxx
// David Grund, March 2, 2023

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
Int_t nBins = 200;

void PlotResults(Double_t pT2_min, Double_t pT2_max); // pT2 in [GeV^2]
void CalculateAvgTPerBin();
void CorrectionPt2ToT();

template <typename TH> // TH2D or TProfile
void SetHisto(TH* h, TString xTitle, TString yTitle, float yOffset = 1.25)
{
    // x-axis
    h->GetXaxis()->SetTitle(xTitle.Data());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetLabelOffset(0.012);
    h->GetXaxis()->SetDecimals(1);
    // y-axis
    h->GetYaxis()->SetTitle(yTitle.Data());
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(yOffset);
    h->GetYaxis()->SetDecimals(1);
    // z-axis    
    h->GetZaxis()->SetLabelSize(0.05);
    return;
}

void STARlight_tVsPt2(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "STARlight_tVsPt2/");

    PlotResults(0.00, 2.56);

    PlotResults(0.04, 1.00);

    CalculateAvgTPerBin();

    CorrectionPt2ToT();

    return;
}

void PlotResults(Double_t pT2_min, Double_t pT2_max)
{
    TFile *f = TFile::Open("Trees/STARlight/IncJ_tVsPt/tree_tPtGammaVMPom.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TTree *t = dynamic_cast<TTree*> (f->Get("tPtGammaVMPom"));
    if(t) Printf("Input tree loaded.");

    ConnectTreeVariables_tPtGammaVMPom(t);

    TH2D* h = new TH2D("h","|#it{t}| vs #it{p}_{T}^{2} of J/#psi",nBins,pT2_min,pT2_max,nBins,pT2_min,pT2_max);
    // on a horizontal axis: J/psi transverse momentum squared
    // on a vertical axis: Mandelstam |t|
    // new: (March 2, 2023)
    Double_t maxDiff = +0.1;
    Double_t minDiff = -0.1;
    TH2D* hDisp = new TH2D("hDisp","(#it{p}_{T}^{2}#minus|#it{t}|)/|#it{t}| vs |#it{t}|",nBins,pT2_min,pT2_max,nBins,minDiff,maxDiff);
    TProfile* hMean = new TProfile("hMean","(#it{p}_{T}^{2}#minus|#it{t}|)/|#it{t}| vs |#it{t}|",nBins,pT2_min,pT2_max,minDiff,maxDiff);

    Double_t* boundaries_pT;
    Double_t boundaries_pT_4[7] = { 0. };
    Double_t boundaries_pT_5[8] = { 0. };
    if(nPtBins == 4) boundaries_pT = &boundaries_pT_4[0];
    else if(nPtBins == 5) boundaries_pT = &boundaries_pT_5[0];
    for(Int_t iBin = 1; iBin <= nPtBins+1; iBin++) boundaries_pT[iBin] = ptBoundaries[iBin-1];
    boundaries_pT[nPtBins+2] = 1.2;
    TH2D *hBins = new TH2D("hBins", "pt gen vs pt rec", nPtBins+2, boundaries_pT, nPtBins+2, boundaries_pT);
    TH1D *hScaleByTotal = new TH1D("hScaleByTotal", "", nPtBins+2, boundaries_pT);

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++)
    {
        t->GetEntry(iEntry);
        Double_t pt2 = fPtVM*fPtVM;
        Double_t abst = fPtPm*fPtPm;
        Double_t relDiff = (pt2 - abst) / abst;
        h->Fill(pt2, abst);
        hDisp->Fill(abst, relDiff);
        hMean->Fill(abst, relDiff);
        hBins->Fill(TMath::Sqrt(abst), TMath::Sqrt(pt2));
        hScaleByTotal->Fill(TMath::Sqrt(abst));
    }

    // Scale the histogram
    for(Int_t iBinX = 1; iBinX <= nPtBins+2; iBinX++){
        Double_t N_tot = hScaleByTotal->GetBinContent(iBinX);
        Printf("BinX %i: %.0f", iBinX, N_tot);
        for(Int_t iBinY = 1; iBinY <= nPtBins+2; iBinY++) {
            Double_t ValueScaled = hBins->GetBinContent(iBinX,iBinY) / N_tot * 100;
            hBins->SetBinContent(iBinX,iBinY,ValueScaled);
        }
    }

    TCanvas *c1 = new TCanvas("c1","",900,800);
    c1->SetLogz();
    c1->SetTopMargin(0.03);
    c1->SetBottomMargin(0.14);
    c1->SetRightMargin(0.13);
    c1->SetLeftMargin(0.14);
    SetHisto(h,"#it{p}_{T,J/#psi}^{2} (GeV^{2}/#it{c}^{2})","|#it{t}|/#it{c}^{2} or #it{p}_{T,pom}^{2} (GeV^{2}/#it{c}^{2})");
    h->Draw("COLZ");
    TString sOut = "Results/" + str_subfolder + "STARlight_tVsPt2/" + Form("dist_%.2f_%.2f", pT2_min, pT2_max);
    c1->Print((sOut + ".pdf").Data());
    //c1->Print((sOut + ".png").Data());

    TLegend *ltw = new TLegend(0.20,0.88,0.35,0.94);
    ltw->AddEntry((TObject*)0,"#bf{This work}","");
    ltw->SetMargin(0.);
    ltw->SetTextSize(0.05);
    ltw->SetBorderSize(0);
    ltw->SetFillStyle(0);
    ltw->Draw();

    sOut = "Results/" + str_subfolder + "_rozprava/t_vs_pt_hist.pdf";
    c1->Print(sOut.Data());
    
    TCanvas *c2 = new TCanvas("c2","",900,800);
    c2->SetLogz();
    c2->SetTopMargin(0.03);
    c2->SetBottomMargin(0.14);
    c2->SetRightMargin(0.13);
    c2->SetLeftMargin(0.16);
    SetHisto(hDisp,"(#it{p}_{T}^{2} #minus |#it{t}|) / |#it{t}| (-)","|#it{t}| (GeV^{2})", 1.5);
    hDisp->Draw("COLZ");
    hMean->SetLineColor(kBlack);
    hMean->SetLineWidth(2);
    hMean->Draw("E0 SAME");
    TLegend l(0.61,0.91,0.86,0.96);
    l.AddEntry(hMean,"mean value","ELP");
    l.SetTextSize(0.045);
    l.SetBorderSize(0);
    //l.SetFillStyle(0);
    l.SetMargin(0.2);
    l.Draw();
    sOut = "Results/" + str_subfolder + "STARlight_tVsPt2/" + + Form("diff_%.2f_%.2f", pT2_min, pT2_max);
    c2->Print((sOut + ".pdf").Data());
    //c2->Print((sOut + ".png").Data());

    // matrix with bins as we have in migration
    TCanvas *c3 = new TCanvas("c3","",1200,600);
    c3->SetTopMargin(0.03);
    c3->SetBottomMargin(0.145);
    c3->SetRightMargin(0.1);
    c3->SetLeftMargin(0.085);

    hBins->SetMarkerSize(2.);
    // horizontal axis
    hBins->GetXaxis()->SetTitle("#sqrt{|#it{t}|}/#it{c} or #it{p}_{T,pom} (GeV/#it{c})");
    hBins->GetXaxis()->SetLabelOffset(0.015);
    hBins->GetXaxis()->SetTitleSize(0.055);
    hBins->GetXaxis()->SetTitleOffset(1.2);
    hBins->GetXaxis()->SetLabelSize(0.055);
    hBins->GetXaxis()->SetDecimals(1);
    // vertical axis
    hBins->GetYaxis()->SetTitle("#it{p}_{T,J/#psi} (GeV/#it{c})");
    hBins->GetYaxis()->SetTitleSize(0.055);
    hBins->GetYaxis()->SetLabelSize(0.055);
    hBins->GetYaxis()->SetTitleOffset(0.7);
    hBins->GetYaxis()->SetDecimals(1);
    // Set ranges
    hBins->GetXaxis()->SetRangeUser(0.0,1.2);
    hBins->GetYaxis()->SetRangeUser(0.0,1.2);
    // Z-axis
    hBins->GetZaxis()->SetLabelSize(0.055);
    hBins->Draw("COLZ TEXT");

    sOut = "Results/" + str_subfolder + "STARlight_tVsPt2/bins.pdf";
    c3->Print(sOut.Data());

    return;
}

void CalculateAvgTPerBin()
// calculate the average value of |t| (p_T,pom^2) and the average value of p_T,J/psi^2 in each bin as predicted by STARlight
{
    TFile *f = TFile::Open("Trees/STARlight/IncJ_tVsPt/tree_tPtGammaVMPom.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TTree *t = dynamic_cast<TTree*> (f->Get("tPtGammaVMPom"));
    if(t) Printf("Input tree loaded.");

    ConnectTreeVariables_tPtGammaVMPom(t);

    Double_t nPt2VMPerBin[5] = { 0 };
    // to calculate average |t|:
    Double_t SumOfTPerBin[5] = { 0 };
    Double_t AvgOfTPerBin[5] = { 0 };
    // to calculate average p_T,J/psi^2
    Double_t SumOfPt2VMPerBin[5] = { 0 };
    Double_t AvgOfPt2VMPerBin[5] = { 0 };

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++)
    {
        t->GetEntry(iEntry);
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

    TString str_1 = Form("Results/%sSTARlight_tVsPt2/AvgTPerBin.txt", str_subfolder.Data());
    ofstream outfile_1(str_1.Data());
    outfile_1 << std::fixed << std::setprecision(6);
    TString str_2 = Form("Results/%sSTARlight_tVsPt2/AvgPt2VMPerBin.txt", str_subfolder.Data());
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

void CorrectionPt2ToT()
{
    TFile *f = TFile::Open("Trees/STARlight/IncJ_tVsPt/tree_tPtGammaVMPom.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TTree *t = dynamic_cast<TTree*> (f->Get("tPtGammaVMPom"));
    if(t) Printf("Input tree loaded.");

    ConnectTreeVariables_tPtGammaVMPom(t);

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
        t->GetEntry(iEntry);
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

    TCanvas *c = new TCanvas("c","c",900,800);
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.04);
    c->SetLeftMargin(0.18);

    // a vertical axis
    hCorrection->GetYaxis()->SetTitle("#it{N}[#it{p}_{T,J/#psi}^{2} #in (|#it{t}|_{min}, |#it{t}|_{max})]/#it{N}[#it{p}_{T,pom}^{2} #in (|#it{t}|_{min}, |#it{t}|_{max})]");
    hCorrection->GetYaxis()->SetTitleSize(0.042);
    hCorrection->GetYaxis()->SetLabelSize(0.05);
    hCorrection->GetYaxis()->SetTitleOffset(2.05);
    hCorrection->GetYaxis()->SetDecimals(3);
    // a horizontal axis
    hCorrection->GetXaxis()->SetTitle("|#it{t}|/#it{c}^{2} or #it{p}_{T,pom}^{2} (GeV^{2}/#it{c}^{2})");
    hCorrection->GetXaxis()->SetTitleSize(0.05);
    hCorrection->GetXaxis()->SetTitleOffset(1.2);
    hCorrection->GetXaxis()->SetLabelSize(0.05);
    hCorrection->GetXaxis()->SetDecimals(1);
    // draw the hogram
    TString sOut = "Results/" + str_subfolder + "STARlight_tVsPt2/ratios_bins.pdf";
    hCorrection->Draw("P");
    // Legend
    /*
    TLegend *leg = new TLegend(0.18,0.78,0.40,0.96);
    leg->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    leg->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    leg->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0); // no border
    leg->SetFillStyle(0);  // legend is transparent
    leg->Draw();
    */
 
    c->Print((sOut + ".pdf").Data());

    TLegend *ltw = new TLegend(0.25,0.88,0.40,0.94);
    ltw->AddEntry((TObject*)0,"#bf{This work}","");
    ltw->SetMargin(0.);
    ltw->SetTextSize(0.05);
    ltw->SetBorderSize(0);
    ltw->SetFillStyle(0);
    ltw->Draw();

    sOut = "Results/" + str_subfolder + "_rozprava/t_vs_pt_ratio.pdf";
    c->Print((sOut + ".pdf").Data());

    return;
}