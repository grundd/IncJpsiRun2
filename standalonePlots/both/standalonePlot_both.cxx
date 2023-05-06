// standalonePlot_both.cxx
// David Grund, May 6, 2023

// cpp headers
#include <iostream>
#include <fstream>
// root headers
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

// data
// coh
float limits_coh[7] = { 0., 0.00072, 0.0016, 0.0026, 0.004, 0.0062, 0.0121 };
float abst_coh[6] = {};
float cs_coh_val[6] = { 8.15, 5.75, 4.23, 2.87, 1.48, 0.4 };// mb
float cs_coh_err_stat[6] = { 0.5, 0.27, 0.2, 0.15, 0.09, 0.04 };
float cs_coh_err_syst_uncr[6] = { 0.18, 0.06, 0.03, 0.04, 0.02, 0.01 };
float cs_coh_err_syst_corr[6] = { 0.46, 0.34, 0.25, 0.17, 0.09, 0.03 };
float cd_coh_err_syst_phot[6] = { 0.20, 0.16, 0.11, 0.08, 0.04, 0.03 };
// incoh
float limits_inc[6] = { 0.04, 0.08, 0.152, 0.258, 0.477, 1. }; // mub
float abst_inc[5] = {};
float cs_inc_val[5] = { 21.8, 19.1, 13.1, 8.1, 4.6 };
float cs_inc_err_stat[5] = { 2.1, 1.9, 1.6, 1.1, 0.6 };
float cs_inc_err_syst_uncr[5] = { 0.3, 0.3, 0.4, 0.1, 0.1 };
float cs_inc_err_syst_corr[5] = { 2.1, 1.5, 0.9, 0.6, 0.3 };
// models
// coh

// incoh
TGraph* gr_MS[2] = { NULL }; // (MS-hs, MS-p)
TGraph* gr_GSZ[2] = { NULL }; // (GSZ-el+diss, GSZ-el)
TGraph* gr_GSZ_err[2] = { NULL }; // error bands

// the MS (IPsat) model
void LoadGraphs_MS(bool print = false)
{
    for(Int_t i = 0; i < 2; i++) gr_MS[i] = new TGraph(183);
    ifstream ifs;
    ifs.open("incoherent_fluct");
    for(int i = 0; i < 183; i++)
    {
        // cross section values in mb
        float abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_MS[0]->SetPoint(i,abs_t_val,sigma_val); // MS-hs
    }
    ifs.close();
    ifs.open("incoherent_nofluct");
    for(int i = 0; i < 183; i++)
    {
        float abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_MS[1]->SetPoint(i,abs_t_val,sigma_val); // MS-n
    }
    ifs.close();
    Printf("MS-hs and MS-p: graphs loaded");
    if(print)
    {
        gr_MS[0]->Print();
        gr_MS[1]->Print();
    }
    return;
}

// the GSZ model
void LoadGraphs_GSZ(bool print = false)
{
    for(Int_t i = 0; i < 2; i++) gr_GSZ[i] = new TGraph(100);
    gr_GSZ_err[0] = new TGraph(2*100);
    gr_GSZ_err[1] = new TGraph(2*100);
    float abs_t_val[100] = { 0 };
    float sig_el_min[100] = { 0 };
    float sig_el_max[100] = { 0 };
    float sig_diss_min[100] = { 0 };
    float sig_diss_max[100] = { 0 };
    float sig_tot_min[100] = { 0 };
    float sig_tot_max[100] = { 0 };
    ifstream ifs;
    ifs.open("incoh_tdep_nuc_run2.dat");
    for(int i = 0; i < 100; i++)
    {
        // cross section values in nb
        ifs >> abs_t_val[i]
            >> sig_el_min[i] >> sig_el_max[i]
            >> sig_diss_min[i] >> sig_diss_max[i]
            >> sig_tot_min[i] >> sig_tot_max[i];
        // transfer to mb:
        float sig_el_mid = (sig_el_max[i] + sig_el_min[i]) / 2 / 1e6;
        float sig_tot_mid = (sig_tot_max[i] + sig_tot_min[i]) / 2 / 1e6;
        gr_GSZ[0]->SetPoint(i,abs_t_val[i],sig_tot_mid); // GSZ-el+diss
        gr_GSZ[1]->SetPoint(i,abs_t_val[i],sig_el_mid);  // GSZ-el
    }
    ifs.close();
    // fill graphs showing the error bands
    for (int i = 0; i < 100; i++)
    {
        gr_GSZ_err[0]->SetPoint(i, abs_t_val[i], sig_tot_max[i] / 1e6);
        gr_GSZ_err[0]->SetPoint(100+i, abs_t_val[100-i-1], sig_tot_min[100-i-1] / 1e6);
        gr_GSZ_err[1]->SetPoint(i, abs_t_val[i], sig_el_max[i] / 1e6);
        gr_GSZ_err[1]->SetPoint(100+i, abs_t_val[100-i-1], sig_el_min[100-i-1] / 1e6);
    }
    Printf("GSZ-el+diss and GSZ-el: graphs loaded");
    if(print)
    {
        gr_GSZ[0]->Print();
        gr_GSZ[1]->Print();
        gr_GSZ_err[0]->Print();
        gr_GSZ_err[1]->Print();
    }
    return;
}

void SetFrame(TH1* fr, float textSize)
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
    g->SetMarkerSize(0.85);
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    g->SetMarkerColor(kBlack);
    return;
}

void SetupSysErrorBox(TGraph* g, Color_t color, float transparency)
{
    g->SetLineWidth(0);
    g->SetMarkerSize(0);
    g->SetFillStyle(1001);
    g->SetFillColorAlpha(color,transparency);
    return;
}

void SetLineMarkerProperties(TGraph *gr, Color_t color, int lineStyle, int markerStyle = kCircle, Size_t markerSize = 1.)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(lineStyle);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerSize(markerSize);
    return;
}

void standalonePlot_both () 
{   
    // load graphs: coh from Roman
    TFile *f = TFile::Open("resultObjects.root", "read");
    TGraph *gr_LTA = (TGraph*)f->Get("graph_LTA");
    TGraph *gr_BKA = (TGraph*)f->Get("graph_BKA");
    // load graphs: incoh
    LoadGraphs_MS();
    LoadGraphs_GSZ();

    // coh
    TGraphAsymmErrors* gr_coh_uncr = new TGraphAsymmErrors(6);
    TGraphAsymmErrors* gr_coh_corr = new TGraphAsymmErrors(6);
    TGraphAsymmErrors* gr_coh_phot = new TGraphAsymmErrors(6);
    for(int i = 0; i < 6; i++) {
        abst_coh[i] = (limits_coh[i+1] + limits_coh[i]) / 2.;
        float abst_err_low = abst_coh[i] - limits_coh[i];
        float abst_err_upp = limits_coh[i+1] - abst_coh[i];
        float cs_err_uncr = TMath::Sqrt(TMath::Power(cs_coh_err_stat[i],2) + TMath::Power(cs_coh_err_syst_uncr[i],2));
        float cs_err_corr = cs_coh_err_syst_corr[i];
        float cs_err_phot = cd_coh_err_syst_phot[i];
        gr_coh_uncr->SetPoint(i,abst_coh[i],cs_coh_val[i]);
        gr_coh_uncr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_uncr,cs_err_uncr);
        gr_coh_corr->SetPoint(i,abst_coh[i],cs_coh_val[i]);
        gr_coh_corr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_corr,cs_err_corr);
        gr_coh_phot->SetPoint(i,abst_coh[i],cs_coh_val[i]);
        gr_coh_phot->SetPointError(i,abst_err_low,abst_err_upp,cs_err_phot,cs_err_phot);
    }
    // inc
    TGraphAsymmErrors* gr_inc_uncr = new TGraphAsymmErrors(5);
    TGraphAsymmErrors* gr_inc_corr = new TGraphAsymmErrors(5);
    for(int i = 0; i < 5; i++) {
        abst_inc[i] = (limits_inc[i+1] + limits_inc[i]) / 2.;
        float abst_err_low = abst_inc[i] - limits_inc[i];
        float abst_err_upp = limits_inc[i+1] - abst_inc[i];
        float cs_err_uncr = (TMath::Sqrt(TMath::Power(cs_inc_err_stat[i],2) + TMath::Power(cs_inc_err_syst_uncr[i],2))) / 1e3;
        float cs_err_corr = cs_inc_err_syst_corr[i] / 1e3;
        gr_inc_uncr->SetPoint(i,abst_inc[i],cs_inc_val[i] / 1e3);
        gr_inc_uncr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_uncr,cs_err_uncr);
        gr_inc_corr->SetPoint(i,abst_inc[i],cs_inc_val[i] / 1e3);
        gr_inc_corr->SetPointError(i,abst_err_low,abst_err_upp,cs_err_corr,cs_err_corr);
    }
    SetStyleUncr(gr_coh_uncr, kFullSquare);
    SetStyleUncr(gr_inc_uncr, kFullCircle);
    SetupSysErrorBox(gr_coh_corr,kGray+3,0.35);
    SetupSysErrorBox(gr_coh_phot,kYellow,0.35);
    SetupSysErrorBox(gr_inc_corr,kGray+3,0.35);
    
    TCanvas *c = new TCanvas ("c","Cross section dependence on |t|",800,700);
    c->SetLeftMargin(0.12);
    c->SetTopMargin(0.03);
    c->SetRightMargin(0.02);
    c->SetBottomMargin(0.13);
    c->SetLogy();
    c->SetLogx();

    TH1F* fr = gPad->DrawFrame(0., 0.003, 1., 10.);
    SetFrame(fr,0.048);
    fr->SetTitle(";|#it{t}| (GeV^{2});d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");
    // y-axis
    fr->GetYaxis()->SetTickLength(0.025); 
    fr->GetYaxis()->SetTitleOffset(1.08);
    // x-axis
    fr->GetXaxis()->SetTickLength(0.025); 
    fr->GetXaxis()->SetTitleOffset(1.20);
    fr->GetXaxis()->SetDecimals(1);

    // models:
    // coh
    SetLineMarkerProperties(gr_LTA,kGray+3,5);
    gStyle->SetLineStyleString(11,"40 20");
    SetLineMarkerProperties(gr_BKA,kRed+2,11);
    // incoh
    SetLineMarkerProperties(gr_MS[0],kViolet-1,7);
    SetLineMarkerProperties(gr_MS[1],kBlue,8);
    SetLineMarkerProperties(gr_GSZ[0],kGreen+2,9);
    SetLineMarkerProperties(gr_GSZ[1],kOrange+2,4);
    SetupSysErrorBox(gr_GSZ_err[0],kGreen,0.35);
    SetupSysErrorBox(gr_GSZ_err[1],kOrange,0.35);

    c->cd();
    fr->Draw("AXIS");
    gr_GSZ_err[1]->Draw("F SAME");
    gr_GSZ[1]->Draw("L SAME");
    gr_GSZ_err[0]->Draw("F SAME");
    gr_GSZ[0]->Draw("L SAME");
    gr_LTA->Draw("L SAME");
    gr_BKA->Draw("L SAME");
    gr_coh_corr->Draw("5 SAME");
    gr_coh_phot->Draw("2 SAME");
    gr_inc_corr->Draw("5 SAME");
    gr_MS[1]->Draw("L SAME");
    gr_MS[0]->Draw("L SAME");
    gr_coh_uncr->Draw("PZ SAME");
    gr_inc_uncr->Draw("PZ SAME");
    // title
    gStyle->SetTextFont(42);
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(0.048);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.925,"ALICE, Pb#minusPb UPC   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // legend
    TLegend *l = new TLegend(0.56,0.50,0.93,0.88);
    l->SetFillColor(0);
    l->SetBorderSize(0);
    l->SetTextSize(0.038);
    l->SetFillStyle(0);
    l->SetMargin(0.20);
    l->AddEntry(gr_coh_uncr,"Coherent J/#psi, |#it{y}| < 0.8", "LEP");
    l->AddEntry(gr_LTA,"LTA (nuclear shadowing)", "L");
    l->AddEntry(gr_BKA,"b-BK (gluon saturation)", "L");
    l->AddEntry(gr_inc_uncr,"Incoherent J/#psi, |#it{y}| < 0.8", "LEP");
    l->AddEntry(gr_MS[0],"MS-hs", "L");
    l->AddEntry(gr_MS[1],"MS-p", "L");
    l->AddEntry(gr_GSZ[0],"GSZ-el+diss", "L");
    l->AddEntry(gr_GSZ[1],"GSZ-el", "L");
    l->Draw();
    // save the canvas
    c->Print("plot.pdf");
    return;
}