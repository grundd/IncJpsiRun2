// CrossSec_Utilities.h
// David Grund, Sep 3, 2022

// c++ headers
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>

// root headers
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>

// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

TString str_data = "data";
TGraph* gr_data_correl = NULL;
TGraph* gr_data_uncorr = NULL;
TString str_models[7] = {"STARlight",
                         "CCK-hs",
                         "CCK-n",
                         "MS-hs",
                         "MS-p",
                         "GSZ-el+diss",
                         "GSZ-el"};
TGraph* gr_models[7] = { NULL };
TGraph* gr_GSZ_err[2] = { NULL };
TH1D* h_models[7] = { NULL };
Int_t n_models[7] = {125, 75, 75, 183, 183, 100, 100};

void InitObjects()
{
    for(Int_t i = 0; i < 7; i++) gr_models[i] = new TGraph(n_models[i]);
    gr_GSZ_err[0] = new TGraph(2*n_models[5]);
    gr_GSZ_err[1] = new TGraph(2*n_models[6]);
}

// the STARlight model
void LoadGraphs_SL(Bool_t print = kFALSE)
{
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/STARlight/IncJ_tDep_0.00-2.50.txt");
    for(Int_t i = 0; i < 125; i++)
    {
        // cross section values in mb
        Double_t abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_models[0]->SetPoint(i,abs_t_val,sigma_val);
    }
    ifs.close();
    Printf("STARlight TGraph loaded.");
    if(print) gr_models[0]->Print();

    return;
}

// the CCK (hot-spot) model
void LoadGraphs_CCK(Bool_t print = kFALSE)
{
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/HSModel/data-dtdy-y_0.6-Run1.txt");
    for(Int_t i = 0; i < 75; i++)
    {
        // cross section values in mb
        Double_t x(0.), tmp, abs_t_val(0.), 
            sig_coh_n_val(0.), sig_coh_n_err(0.), 
            sig_inc_n_val(0.), sig_inc_n_err(0.), 
            sig_coh_hs_val(0.), sig_coh_hs_err(0.), 
            sig_inc_hs_val(0.), sig_inc_hs_err(0.);
        ifs >> x >> tmp >> tmp
            >> abs_t_val
            >> sig_coh_n_val >> sig_coh_n_err
            >> sig_inc_n_val >> sig_inc_n_err
            >> sig_coh_hs_val >> sig_coh_hs_err
            >> sig_inc_hs_val >> sig_inc_hs_err;
        gr_models[1]->SetPoint(i,abs_t_val,sig_inc_hs_val); // CCK-hs
        gr_models[2]->SetPoint(i,abs_t_val,sig_inc_n_val);  // CCK-n
    }
    ifs.close();
    Printf("CCK TGraphs loaded.");
    if(print)
    {
        gr_models[1]->Print();
        gr_models[2]->Print();
    }

    return;
}

// the MS (IPsat) model
void LoadGraphs_MS(Bool_t print = kFALSE)
{
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_fluct");
    for(Int_t i = 0; i < 183; i++)
    {
        // cross section values in mb
        Double_t abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_models[3]->SetPoint(i,abs_t_val,sigma_val); // MS-hs
    }
    ifs.close();
    ifs.open("Trees/PhotoCrossSec/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_nofluct");
    for(Int_t i = 0; i < 183; i++)
    {
        Double_t abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_models[4]->SetPoint(i,abs_t_val,sigma_val); // MS-n
    }
    ifs.close();
    Printf("MS TGraphs loaded.");
    if(print)
    {
        gr_models[3]->Print();
        gr_models[4]->Print();
    }

    return;
}

// the GSZ model
void LoadGraphs_GSZ(Bool_t print = kFALSE)
{
    Double_t abs_t_val[100] = { 0 };
    Double_t sig_el_min[100] = { 0 };
    Double_t sig_el_max[100] = { 0 };
    Double_t sig_diss_min[100] = { 0 };
    Double_t sig_diss_max[100] = { 0 };
    Double_t sig_tot_min[100] = { 0 };
    Double_t sig_tot_max[100] = { 0 };
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Guzey/incoh_tdep_nuc_run2.dat");
    for(Int_t i = 0; i < 100; i++)
    {
        // cross section values in nb
        ifs >> abs_t_val[i]
            >> sig_el_min[i] >> sig_el_max[i]
            >> sig_diss_min[i] >> sig_diss_max[i]
            >> sig_tot_min[i] >> sig_tot_max[i];
        // transfer to mb:
        Double_t sig_el_mid = (sig_el_max[i] + sig_el_min[i]) / 2 / 1e6;
        Double_t sig_tot_mid = (sig_tot_max[i] + sig_tot_min[i]) / 2 / 1e6;
        gr_models[5]->SetPoint(i,abs_t_val[i],sig_tot_mid); // GSZ-el+diss
        gr_models[6]->SetPoint(i,abs_t_val[i],sig_el_mid);  // GSZ-el
    }
    ifs.close();
    // fill graphs showing the error bands
    for (Int_t i = 0; i < 100; i++)
    {
        gr_GSZ_err[0]->SetPoint(i, abs_t_val[i], sig_tot_max[i] / 1e6);
        gr_GSZ_err[0]->SetPoint(100+i, abs_t_val[100-i-1], sig_tot_min[100-i-1] / 1e6);
        gr_GSZ_err[1]->SetPoint(i, abs_t_val[i], sig_el_max[i] / 1e6);
        gr_GSZ_err[1]->SetPoint(100+i, abs_t_val[100-i-1], sig_el_min[100-i-1] / 1e6);
    }
    Printf("GSZ TGraphs loaded");
    if(print)
    {
        gr_models[5]->Print();
        gr_models[6]->Print();
        gr_GSZ_err[0]->Print();
        gr_GSZ_err[1]->Print();
    }

    return;
}

void CreateHistogramFromGraph(Int_t iM, Bool_t print = kFALSE) // iM = index of a model
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/HistogramsFromGraphs/");

    vector<Double_t> t_edges;
    // get arrays of |t| and sigma values from the graphs
    Double_t *abs_t_val = gr_models[iM]->GetX();
    Double_t *sigma_val = gr_models[iM]->GetY();
    // print the values:
    if(print)
    {
        Printf("Values loaded from the graph:");
        for(Int_t i = 0; i < n_models[iM]; i++) cout << i << "\t" << abs_t_val[i] << "\t" << sigma_val[i] << "\n";
    } 
    // first edge:
    Double_t first_edge = abs_t_val[0] - (abs_t_val[0] + abs_t_val[1]) / 2.;
    if(first_edge < 0) t_edges.push_back(0.);
    else               t_edges.push_back(first_edge);
    // up to the next-to-last one:
    for(Int_t i = 0; i < n_models[iM]-1; i++)
    {
        t_edges.push_back((abs_t_val[i] + abs_t_val[i+1]) / 2);
    } 
    // last edge:
    t_edges.push_back(abs_t_val[n_models[iM]-1] + (abs_t_val[n_models[iM]-1] - abs_t_val[n_models[iM]-2]) / 2.);
    // print the values
    if(print)
    {
        Printf("Values to be put into the histogram:");
        for(Int_t i = 0; i < n_models[iM]; i++) cout << i << "\t" << t_edges[i] << "\t" << t_edges[i+1] << "\t" << sigma_val[i] << "\n";
    } 
    // create array from vector
    Double_t *t_edges_arr = &t_edges[0]; 
    // create the histogram
    h_models[iM] = new TH1D("h_" + str_models[iM], "h_" + str_models[iM], n_models[iM], t_edges_arr);
    for(Int_t iBin = 1; iBin <= n_models[iM]; iBin++) h_models[iM]->SetBinContent(iBin,sigma_val[iBin-1]);

    // draw the graph and the histogram on the same canvas
    TCanvas *c = new TCanvas("c_"+str_models[iM],"c_"+str_models[iM],1600,600);
    c->SetLogy();
    // margins
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.02);
    c->SetLeftMargin(0.09);
    // histogram settings
    h_models[iM]->SetTitle(";|#it{t}| (GeV^{2}); d#sigma_{#gammaPb}/d|#it{t}| (mb/GeV^{2})");  
    h_models[iM]->SetLineColor(kRed);
    h_models[iM]->SetLineWidth(1);
    // vertical axis
    h_models[iM]->GetYaxis()->SetTitleSize(0.05);
    h_models[iM]->GetYaxis()->SetTitleOffset(0.9);
    h_models[iM]->GetYaxis()->SetLabelSize(0.05);
    // horizontal axis
    h_models[iM]->GetXaxis()->SetTitleSize(0.05);
    h_models[iM]->GetXaxis()->SetTitleOffset(1.2);
    h_models[iM]->GetXaxis()->SetLabelSize(0.05);
    // draw histogram
    h_models[iM]->Draw("HIST");
    // Draw graph
    gr_models[iM]->SetMarkerStyle(kFullCircle);
    gr_models[iM]->SetMarkerColor(kBlue);
    gr_models[iM]->SetMarkerSize(0.5);
    gr_models[iM]->Draw("P SAME");

    c->Print("Results/" + str_subfolder + "CrossSec/HistogramsFromGraphs/" + str_models[iM] + ".pdf");

    return;
}

Double_t IntegrateModel(Int_t iM, Double_t t_min, Double_t t_max) // iM = index of a model
// this function works for t_min, t_max in the range 0.04 to 1.0
{
    //if(t_max < t_min)  return -1e10; // undefined
    //if(TMath::Abs(t_min - t_max) < 1e-10) return 0.; // if equal, return zero
    Int_t iBinMin = h_models[iM]->FindBin(t_min);
    Int_t iBinMax = h_models[iM]->FindBin(t_max);
    Double_t integral(0.);
    //if(iBinMin == iBinMax) integral += h_models[iM]->GetBinContent(iBinMin) * (t_max - t_min);    
    integral += (h_models[iM]->GetBinLowEdge(iBinMin+1) - t_min) * h_models[iM]->GetBinContent(iBinMin);
    integral += (t_max - h_models[iM]->GetBinLowEdge(iBinMax)) * h_models[iM]->GetBinContent(iBinMax);
    // only if iBinMin and iBinMax differ by at least 1:
    //if(iBinMax > iBinMin+1) 
    integral += h_models[iM]->Integral(iBinMin+1, iBinMax-1, "width");
    Printf("***************************************");
    Printf("model %s, range (%.3f,%.3f) GeV", str_models[iM].Data(), t_min, t_max);
    Printf("iBinMin = %i, bin range = (%.3f,%.3f) GeV", iBinMin, h_models[iM]->GetBinLowEdge(iBinMin), h_models[iM]->GetBinLowEdge(iBinMin+1));
    Printf("iBinMax = %i, bin range = (%.3f,%.3f) GeV", iBinMax, h_models[iM]->GetBinLowEdge(iBinMax), h_models[iM]->GetBinLowEdge(iBinMax+1));
    Printf("integral = %.3f mub", integral * 1e3);
    Printf("***************************************");
    return integral;
}