// CrossSec_Utilities.h
// David Grund, Sep 3, 2022

// c++ headers
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>
#include <vector>

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
#include "SetPtBinning.h"

TString str_data = "data";
TGraphAsymmErrors* gr_data_uncr = NULL;
TGraphAsymmErrors* gr_data_corr = NULL;
TGraphAsymmErrors* gr_data_stat = NULL;
TGraphAsymmErrors* gr_data_syst_uncr = NULL;
TGraphAsymmErrors* gr_data_syst_corr = NULL;
// order of the models within the analysis: STARlight, CCK-hs, CCK-n, MS-hs, MS-p, GSZ-el+diss, GSZ-el
TString str_models[9] = {"STARlight",
                         "CCK-hs",
                         "CCK-n",
                         "MS-hs",
                         "MS-p",
                         "GSZ-el+diss",
                         "GSZ-el",
                         "MSS-CGC+fl",
                         "MSS-CGC"};
TGraph* gr_models[9] = { NULL };
TGraph* gr_GSZ_err[2] = { NULL };
Double_t GSZ_err_scale_upp[2] = { 0 };
Double_t GSZ_err_scale_low[2] = { 0 };
TH1D* h_models[9] = { NULL };
Int_t n_models[9] = {125, 75, 75, 183, 183, 100, 100, 696, 696};
Double_t *tBoundaries = NULL;
Double_t tBoundaries_4bins[5] = { 0 };
Double_t tBoundaries_5bins[6] = { 0 };
Int_t lineWidth = 3;
Color_t colors[9] = {
    kYellow+2, // STARlight
    kRed+1, // CCK-hs
    kCyan+2, // CCK-n
    kViolet-1, // MS-hs
    kBlue, // MS-p
    kGreen+2, // GSZ-el+diss
    kOrange+2, // GSZ-el
    kRed+1, // MSS-CGC+fl (!)
    kYellow+2, // MSS-CGC (!)
};

void InitObjects()
{
    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gr_data_uncr = new TGraphAsymmErrors(nPtBins);
    gr_data_corr = new TGraphAsymmErrors(nPtBins);
    gr_data_stat = new TGraphAsymmErrors(nPtBins);
    gr_data_syst_uncr = new TGraphAsymmErrors(nPtBins);
    gr_data_syst_corr = new TGraphAsymmErrors(nPtBins);
    gr_data_uncr->SetName("gr_data_uncr");
    gr_data_corr->SetName("gr_data_corr");
    gr_data_stat->SetName("gr_data_stat");
    gr_data_syst_uncr->SetName("gr_data_syst_uncr");
    gr_data_syst_corr->SetName("gr_data_syst_corr");
    
    for(Int_t i = 0; i < 9; i++)
    {
        gr_models[i] = new TGraph(n_models[i]);
        gr_models[i]->SetName("gr_" + str_models[i]);
    } 
    gr_GSZ_err[0] = new TGraph(2*n_models[5]);
    gr_GSZ_err[0]->SetName("gr_err_GSZ-el+diss");
    gr_GSZ_err[1] = new TGraph(2*n_models[6]);
    gr_GSZ_err[1]->SetName("gr_err_GSZ-el");

    SetPtBinning();
    if(nPtBins == 4) tBoundaries = tBoundaries_4bins;
    if(nPtBins == 5) tBoundaries = tBoundaries_5bins;
    for(Int_t i = 0; i <= nPtBins; i++) tBoundaries[i] = ptBoundaries[i] * ptBoundaries[i];
}

void LoadGraphs_data(Bool_t abstFromModel = kFALSE, Int_t abstModel = 1, Bool_t print = kFALSE)
// abstFromModel:
// = false => bin centers
// = true => use a model
// abstModel:
// = 0 => SL ... original approach
// = 1 => GSZ-el+diss
// = 2 => MS-hs
{
    Double_t abs_t_val[5] = { 0 };
    Double_t sig_val[5] = { 0 };
    Double_t sig_err_stat[5] = { 0 };
    Double_t sig_err_syst_uncr[5] = { 0 };
    Double_t sig_err_syst_corr[5] = { 0 };
    Double_t sig_err_uncr[5] = { 0 };
    Double_t sig_err_corr[5] = { 0 };
    ifstream ifs;
    ifs.open("Results/" + str_subfolder + "CrossSec/CrossSec_photo.txt");
    for(Int_t i = 0; i <= nPtBins; i++)
    {
        // cross section values in mub
        Int_t bin; Double_t t_low, t_upp, val, stat, syst_uncr, syst_corr;
        ifs >> bin >> t_low >> t_upp >> val >> stat >> syst_uncr >> syst_corr;
        if(i > 0) { // index 0 -> fiducial
            sig_val[i-1] = val;
            sig_err_stat[i-1] = stat;
            sig_err_syst_uncr[i-1] = syst_uncr;
            sig_err_syst_corr[i-1] = syst_corr;
        }
    }
    ifs.close();
    TString str_avgt;
    if(!abstFromModel) {
        for(Int_t i = 0; i < nPtBins; i++) {
            abs_t_val[i] = (tBoundaries[i+1] + tBoundaries[i]) / 2.;
            //cout << abs_t_val[i] << endl;
        }
    } else {
        if(abstModel == 0)      str_avgt = "Results/" + str_subfolder + "STARlight_tVsPt/AvgTPerBin.txt";
        else if(abstModel == 1) str_avgt = "Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/AverageT/" + str_models[5] + ".txt";
        else if(abstModel == 2) str_avgt = "Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/AverageT/" + str_models[3] + ".txt";
        else return;
        ifs.open(str_avgt.Data());
        for(Int_t i = 0; i < nPtBins; i++)
        {
            Int_t bin; Double_t tlow, tupp;
            if(abstModel == 0) ifs >> bin >> abs_t_val[i];
            else               ifs >> tlow >> tupp >> abs_t_val[i];
            //cout << abs_t_val[i] << endl;
        }
        ifs.close();
    }
    for(Int_t i = 0; i < nPtBins; i++)
    {
        // from mub to mb
        Double_t abs_t_err_low = abs_t_val[i] - tBoundaries[i];
        Double_t abs_t_err_upp = tBoundaries[i+1] - abs_t_val[i];
        sig_val[i] = sig_val[i] / 1e3;
        sig_err_stat[i] = sig_err_stat[i] / 1e3;
        sig_err_syst_uncr[i] = sig_err_syst_uncr[i] / 1e3;
        sig_err_syst_corr[i] = sig_err_syst_corr[i] / 1e3;
        sig_err_uncr[i] = TMath::Sqrt(TMath::Power(sig_err_stat[i],2) + TMath::Power(sig_err_syst_uncr[i],2));
        sig_err_corr[i] = sig_err_syst_corr[i];
        gr_data_uncr->SetPoint(i,abs_t_val[i],sig_val[i]);
        gr_data_uncr->SetPointError(i,abs_t_err_low,abs_t_err_upp,sig_err_uncr[i],sig_err_uncr[i]);
        gr_data_corr->SetPoint(i,abs_t_val[i],sig_val[i]);
        gr_data_corr->SetPointError(i,abs_t_err_low,abs_t_err_upp,sig_err_corr[i],sig_err_corr[i]);
        gr_data_stat->SetPoint(i,abs_t_val[i],sig_val[i]);
        gr_data_stat->SetPointError(i,abs_t_err_low,abs_t_err_upp,sig_err_stat[i],sig_err_stat[i]);
        gr_data_syst_uncr->SetPoint(i,abs_t_val[i],sig_val[i]);
        gr_data_syst_uncr->SetPointError(i,abs_t_err_low,abs_t_err_upp,sig_err_syst_uncr[i],sig_err_syst_uncr[i]);
        gr_data_syst_corr->SetPoint(i,abs_t_val[i],sig_val[i]);
        gr_data_syst_corr->SetPointError(i,abs_t_err_low,abs_t_err_upp,sig_err_syst_corr[i],sig_err_syst_corr[i]);

        if(print) Printf("bin %i t_low: %.4f t_upp: %.4f sig: %.3f [mub]", i+1, 
            abs_t_val[i] - abs_t_err_low, abs_t_val[i] + abs_t_err_upp, sig_val[i] * 1e3);
    }
    Printf("TGraph for data created.");
    if(print)
    {
        gr_data_uncr->Print();
        gr_data_corr->Print();
        gr_data_stat->Print();
        gr_data_syst_uncr->Print();
        gr_data_syst_corr->Print();
    } 

    return;
}

Double_t IntegrateData()
{
    Double_t integral(0.);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        Double_t t_low = gr_data_uncr->GetPointX(i) - gr_data_uncr->GetErrorXlow(i);
        Double_t t_upp = gr_data_uncr->GetPointX(i) + gr_data_uncr->GetErrorXhigh(i);
        integral += gr_data_uncr->GetPointY(i) * (t_upp - t_low);
    }
    return integral;
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
    Printf("TGraph for STARlight created.");
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
    Printf("TGraphs for CCK created.");
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
    Printf("TGraphs for MS created.");
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
    // determine the scales
    for(Int_t i = 0; i < 100; i++)
    {
        GSZ_err_scale_upp[0] += sig_tot_max[i] / ((sig_tot_max[i] + sig_tot_min[i]) / 2.);
        GSZ_err_scale_low[0] += sig_tot_min[i] / ((sig_tot_max[i] + sig_tot_min[i]) / 2.);
        GSZ_err_scale_upp[1] += sig_el_max[i] / ((sig_el_max[i] + sig_el_min[i]) / 2.);
        GSZ_err_scale_low[1] += sig_el_min[i] / ((sig_el_max[i] + sig_el_min[i]) / 2.);
    }
    for(Int_t i = 0; i < 2; i++)
    {
        GSZ_err_scale_upp[i] = GSZ_err_scale_upp[i] / 100.;
        GSZ_err_scale_low[i] = GSZ_err_scale_low[i] / 100.;
    }
    Printf("TGraphs for GSZ created");
    if(print)
    {
        gr_models[5]->Print();
        gr_models[6]->Print();
        gr_GSZ_err[0]->Print();
        gr_GSZ_err[1]->Print();
        Printf(" +++++++++++++++++++++++++++++++");
        Printf(" Scale of the upper error band:");
        Printf(" el+diss: %.5f", GSZ_err_scale_upp[0]);
        Printf(" el     : %.5f", GSZ_err_scale_upp[1]);
        Printf(" Scale of the lower error band:");
        Printf(" el+diss: %.5f", GSZ_err_scale_low[0]);
        Printf(" el     : %.5f", GSZ_err_scale_low[1]);
        Printf(" +++++++++++++++++++++++++++++++");
    }

    return;
}

// the MSS (CGC) model
void LoadGraphs_MSS(Bool_t print = kFALSE)
{
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Heikki_new/cgc_with_shapefluct.txt"); 
    for(Int_t i = 0; i < 696; i++)
    {
        // cross section values in mb
        Double_t abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_models[7]->SetPoint(i,abs_t_val,sigma_val); // MSS-CGC+fluct
    }
    ifs.close();
    ifs.open("Trees/PhotoCrossSec/Heikki_new/cgc_no_shapefluct.txt"); 
    for(Int_t i = 0; i < 696; i++)
    {
        Double_t abs_t_val(0.), sigma_val(0.);
        ifs >> abs_t_val >> sigma_val;
        gr_models[8]->SetPoint(i,abs_t_val,sigma_val); // MSS-CGC
    }
    ifs.close();
    Printf("TGraphs for MSS created.");
    if(print)
    {
        gr_models[7]->Print();
        gr_models[8]->Print();
    }

    return;
}

void CreateHistogramFromGraph(Int_t iM, Bool_t print = kTRUE) // iM = index of a model
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/HistogramsFromGraphs/");
    ofstream os;
    os.open("Results/" + str_subfolder + "CrossSec/HistogramsFromGraphs/log_" + str_models[iM] + ".txt");
    os << std::fixed << std::setprecision(4);

    vector<Double_t> t_edges;
    // get arrays of |t| and sigma values from the graphs
    Double_t *abs_t_val = gr_models[iM]->GetX();
    Double_t *sigma_val = gr_models[iM]->GetY();
    // print the values:
    if(print)
    {
        os << "values loaded from the graph:\n";
        for(Int_t i = 0; i < n_models[iM]; i++) os << i+1 << "\t" << abs_t_val[i] << "\t" << sigma_val[i] << "\n";
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
    // print the values:
    if(print)
    {
        os << "values to be put into the histogram:\n";
        for(Int_t i = 0; i < n_models[iM]; i++) os << i+1 << "\t" << t_edges[i] << "\t" << t_edges[i+1] << "\t" << sigma_val[i] << "\n";
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
    os.close();
    delete c;

    return;
}

Bool_t IntegrateModel(Int_t iM, Double_t t_min, Double_t t_max, Double_t &integral, Double_t &avgt) // iM = index of a model
// this function works for t_min, t_max in the range 0.04 to 1.0
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/Integrals/");
    integral = 0; avgt = 0;
    if(t_max < t_min) return kFALSE; // undefined, return false
    if(t_min == t_max) return kTRUE; // if the borders are equal, return zeros and true
    Int_t iBinMin = h_models[iM]->FindBin(t_min);
    Int_t iBinMax = h_models[iM]->FindBin(t_max);
    // if we integrate only within one bin:
    if(iBinMin == iBinMax){
        integral += h_models[iM]->GetBinContent(iBinMin) * (t_max - t_min);
        avgt = t_max + t_min / 2.;
    } 
    // if we integrate within at least two bins:
    else{
        // first bin
        Double_t first_integral = (h_models[iM]->GetBinLowEdge(iBinMin+1) - t_min) * h_models[iM]->GetBinContent(iBinMin);
        Double_t first_center = (t_min + h_models[iM]->GetBinLowEdge(iBinMin+1)) / 2.;
        integral += first_integral;
        avgt     += first_integral * first_center;
        // last bin
        Double_t last_integral = (t_max - h_models[iM]->GetBinLowEdge(iBinMax)) * h_models[iM]->GetBinContent(iBinMax);
        Double_t last_center = (h_models[iM]->GetBinLowEdge(iBinMax) + t_max) / 2.;
        integral += last_integral;
        avgt     += last_integral * last_center;
        // if within more than two bins => integral in between the first and last bins:
        if(iBinMax > iBinMin+1){
            for(Int_t i = iBinMin+1; i <= iBinMax-1; i++){
                integral += (h_models[iM]->GetBinLowEdge(i+1) - h_models[iM]->GetBinLowEdge(i)) * h_models[iM]->GetBinContent(i);
                avgt     += (h_models[iM]->GetBinLowEdge(i+1) - h_models[iM]->GetBinLowEdge(i)) * h_models[iM]->GetBinContent(i) * h_models[iM]->GetBinCenter(i);
            }
            //integral += h_models[iM]->Integral(iBinMin+1, iBinMax-1, "width"); // alternative
        } 
    }  
    avgt = avgt / integral; 
    // print the results
    ofstream os;
    os.open("Results/" + str_subfolder + "CrossSec/Integrals/" + Form("%s %.3f-%.3f GeV.txt", str_models[iM].Data(), t_min, t_max)); 
    os << std::fixed << std::setprecision(3);
    os << "model " << str_models[iM].Data() << "\n"
       << "range (" << t_min << "," << t_max << ") GeV\n"
       << "integral: " << integral * 1e3 << "\n"
       << "avg |t| : " << avgt << "\n";
    os.close();
    // print to the console
    Printf(" +++++++++++++++++++++++++++++++++++++++");
    Printf(" %s, range (%.3f,%.3f) GeV", str_models[iM].Data(), t_min, t_max);
    Printf(" iBinMin: %i \tbin range: (%.3f,%.3f) GeV", iBinMin, h_models[iM]->GetBinLowEdge(iBinMin), h_models[iM]->GetBinLowEdge(iBinMin+1));
    Printf(" iBinMax: %i \tbin range: (%.3f,%.3f) GeV", iBinMax, h_models[iM]->GetBinLowEdge(iBinMax), h_models[iM]->GetBinLowEdge(iBinMax+1));
    Printf(" integral: %.3f mub", integral * 1e3);
    Printf(" avg |t| : %.3f GeV", avgt);
    Printf(" +++++++++++++++++++++++++++++++++++++++");

    return kTRUE;
}

TLegend *SetLegend(Double_t x_leftdown, Double_t y_leftdown, Double_t x_rightup, Double_t y_rightup)
{
    TLegend *l = new TLegend(x_leftdown, y_leftdown, x_rightup, y_rightup);
    l->SetFillColor(0);
    //l->SetFillStyle(0);
    l->SetBorderSize(0);
    return l;
}

void SetStyle(TGraph* g, Color_t color, Style_t style, Width_t width = 3)
{
    g->SetLineColor(color);
    g->SetLineStyle(style);
    g->SetLineWidth(width);
}

void SetLineMarkerProperties(TGraph *gr, Color_t color, Int_t lineStyle, Int_t markerStyle = kCircle, Size_t markerSize = 1.)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(lineStyle);
    gr->SetLineWidth(lineWidth);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerSize(markerSize);
    return;
}

void SetupSysErrorBox(TGraph* g, Color_t color)
{
    g->SetMarkerSize(0);
    g->SetFillStyle(1001);
    g->SetFillColorAlpha(color,0.35);
    SetStyle(g,color,1,0);
}