// plot.cxx
// David Grund, April 26, 2023

#include "TGraphAsymmErrors.h"
#include "TH1F.h"

double limits_coh[7] = { 0., 0.00072, 0.0016, 0.0026, 0.004, 0.0062, 0.0121 };
double abst_coh[6] = {};
double cs_coh_val[6] = { 8.15, 5.75, 4.23, 2.87, 1.48, 0.4 };// mb
double cs_coh_err_stat[6] = { 0.5, 0.27, 0.2, 0.15, 0.09, 0.04 };
double cs_coh_err_syst_uncr[6] = { 0.18, 0.06, 0.03, 0.04, 0.02, 0.01 };
double cs_coh_err_syst_corr[6] = { 0.46, 0.34, 0.25, 0.17, 0.09, 0.03 };
double cd_coh_err_syst_phot[6] = { 0.20, 0.16, 0.11, 0.08, 0.04, 0.03 };
double limits_inc[6] = { 0.04, 0.08, 0.152, 0.258, 0.477, 1. }; // mub
double abst_inc[5] = {};
double cs_inc_val[5] = { 21.9, 19.1, 13.1, 8.1, 4.6 };
double cs_inc_err_stat[5] = { 2.1, 1.9, 1.6, 1.1, 0.6 };
double cs_inc_err_syst_uncr[5] = { 0.3, 0.3, 0.4, 0.1, 0.1 };
double cs_inc_err_syst_corr[5] = { 2.1, 1.5, 0.9, 0.6, 0.3 };

void SetFrame(TH1* fr, Double_t textSize)
{
    // title and label sizes
    fr->GetXaxis()->SetTitleSize(textSize);
    fr->GetYaxis()->SetTitleSize(textSize);
    fr->GetXaxis()->SetLabelSize(textSize);
    fr->GetYaxis()->SetLabelSize(textSize);
    // font types
    fr->GetXaxis()->SetTitleFont(42);
    fr->GetYaxis()->SetTitleFont(42);
    fr->GetXaxis()->SetLabelFont(42);
    fr->GetYaxis()->SetLabelFont(42);
    // divisions on the x-axis
    fr->GetXaxis()->SetNdivisions(505);
    return;
}

void SetStyleUncr (TGraphAsymmErrors* g, Marker_t m)
{
    g->SetMarkerStyle(m);
    g->SetMarkerSize(0.75);
    g->SetLineColor(kBlack);
    g->SetLineWidth(1);
    g->SetMarkerColor(kBlack);
    return;
}

void SetupSysErrorBox(TGraph* g, Color_t color)
{
    g->SetLineWidth(0);
    g->SetMarkerSize(0);
    g->SetFillStyle(1001);
    g->SetFillColorAlpha(color,0.3);
    return;
}

void plot () 
{   
    // coh
    TGraphAsymmErrors* gr_coh_uncr = new TGraphAsymmErrors(6);
    TGraphAsymmErrors* gr_coh_corr = new TGraphAsymmErrors(6);
    for(int i = 0; i < 6; i++) {
        abst_coh[i] = (limits_coh[i+1] + limits_coh[i]) / 2.;
        double abst_err_low = abst_coh[i] - limits_coh[i];
        double abst_err_upp = limits_coh[i+1] - abst_coh[i];
        double cs_err_uncr = TMath::Sqrt(TMath::Power(cs_coh_err_stat[i],2) + TMath::Power(cs_coh_err_syst_uncr[i],2));
        double cs_err_corr = cs_coh_err_syst_corr[i];
        gr_coh_uncr->SetPoint(i,abst_coh[i],cs_coh_val[i]);
        gr_coh_uncr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_uncr,cs_err_uncr);
        gr_coh_corr->SetPoint(i,abst_coh[i],cs_coh_val[i]);
        gr_coh_corr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_corr,cs_err_corr);
    }
    // inc
    TGraphAsymmErrors* gr_inc_uncr = new TGraphAsymmErrors(5);
    TGraphAsymmErrors* gr_inc_corr = new TGraphAsymmErrors(5);
    for(int i = 0; i < 5; i++) {
        abst_inc[i] = (limits_inc[i+1] + limits_inc[i]) / 2.;
        double abst_err_low = abst_inc[i] - limits_inc[i];
        double abst_err_upp = limits_inc[i+1] - abst_inc[i];
        double cs_err_uncr = (TMath::Sqrt(TMath::Power(cs_inc_err_stat[i],2) + TMath::Power(cs_inc_err_syst_uncr[i],2))) / 1e3;
        double cs_err_corr = cs_inc_err_syst_corr[i] / 1e3;
        gr_inc_uncr->SetPoint(i,abst_inc[i],cs_inc_val[i] / 1e3);
        gr_inc_uncr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_uncr,cs_err_uncr);
        gr_inc_corr->SetPoint(i,abst_inc[i],cs_inc_val[i] / 1e3);
        gr_inc_corr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_corr,cs_err_corr);
    }
    SetStyleUncr(gr_coh_uncr, kOpenCircle);
    SetStyleUncr(gr_inc_uncr, kFullCircle);
    SetupSysErrorBox(gr_coh_corr,kGray+3);
    SetupSysErrorBox(gr_inc_corr,kGray+3);
    
    TCanvas *c = new TCanvas ("c","Cross section dependence on |t|",900,600);
    c->SetLeftMargin(0.11);
    c->SetTopMargin(0.04);
    c->SetRightMargin(0.02);
    c->SetBottomMargin(0.13);
    c->SetLogy();

    TH1F* fr = gPad->DrawFrame(0., 0.003, 1., 10.);
    SetFrame(fr,0.05);
    fr->SetTitle(";|#it{t}| (GeV^{2});d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");
    // y-axis
    fr->GetYaxis()->SetTickLength(0.025); 
    fr->GetYaxis()->SetTitleOffset(1.00);
    // x-axis
    fr->GetXaxis()->SetTickLength(0.025); 
    fr->GetXaxis()->SetTitleOffset(1.20);
    fr->GetXaxis()->SetDecimals(1);

    c->cd();
    fr->Draw("AXIS");
    gr_coh_corr->Draw("5 SAME");
    gr_inc_corr->Draw("5 SAME");
    gr_coh_uncr->Draw("PZ SAME");
    gr_inc_uncr->Draw("PZ SAME");

    //TLegend l(0.8,);
    c->Print("plot.pdf");
    c->SetLogx();
    c->Print("plot_log.pdf");
    return;
}