// PhotoCrossSec_Utilities.h
// David Grund, Apr 24, 2022

// c++ headers
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>

// root headers
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

// To read the values from the files:
Double_t sig_val[5] = { 0 };
Double_t sig_err_stat[5] = { 0 };
Double_t sig_err_syst[5] = { 0 };
Double_t abs_t_val[5] = { 0 };
Double_t t_boundaries[5+1] = { 0 };

// For TGraphAsymmErrors:
Double_t sig_err_stat_upp[5] = { 0 };
Double_t sig_err_stat_low[5] = { 0 };
Double_t sig_err_syst_upp[5] = { 0 };
Double_t sig_err_syst_low[5] = { 0 };
Double_t abs_t_err_low[5] = { 0 };
Double_t abs_t_err_upp[5] = { 0 };

// HS predictions
// reserve space for data to be read
const Int_t nData_HS = 75;
Double_t abs_t_HS[nData_HS];
Double_t sig_HS_coh_n[nData_HS];
Double_t sig_HS_coh_n_err[nData_HS];
Double_t sig_HS_inc_n[nData_HS];
Double_t sig_HS_inc_n_err[nData_HS];
Double_t sig_HS_coh_hs[nData_HS];
Double_t sig_HS_coh_hs_err[nData_HS];
Double_t sig_HS_inc_hs[nData_HS];
Double_t sig_HS_inc_hs_err[nData_HS];

// Guzey predictions
const Int_t nData_GZ = 100;
Double_t abs_t_GZ[nData_GZ];
Double_t sig_GZ_el_min[nData_GZ];
Double_t sig_GZ_el_max[nData_GZ];
Double_t sig_GZ_diss_min[nData_GZ];
Double_t sig_GZ_diss_max[nData_GZ];
Double_t sig_GZ_tot_min[nData_GZ];
Double_t sig_GZ_tot_max[nData_GZ];

// Heikki predictions
const Int_t nData_HM = 183;
Double_t abs_t_HM[nData_HM];
Double_t sig_HM_fluct[nData_HM];
Double_t sig_HM_noflu[nData_HM];

// STARlight predictions
const Int_t nData_SL = 125;
Double_t abs_t_SL[nData_SL];
Double_t sig_SL[nData_SL];

Int_t lineWidth = 3;

//#####################################################################################################
// Functions to read the input

void ReadInput_Measurement()
{
    // read the input file for measured cross section
    ifstream ifs;
    t_boundaries[0] = 0.04;
    TString str = "Results/" + str_subfolder + "PhotoCrossSec/CrossSec_photo.txt";
    ifs.open(str.Data());
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            // skip first line
            Int_t bin;
            Double_t tLow;
            if(i > 0) istr >> bin >> tLow >> t_boundaries[i] >> sig_val[i-1] >> sig_err_stat[i-1] >> sig_err_syst[i-1];
            i++;   
        }
        ifs.close();
    }
    Printf("Values of the photonuclear cross section loaded.");

    str = "Results/" + str_subfolder + "STARlight_tVsPt/AvgTPerBin.txt";
    ifs.open(str.Data()); 
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            Int_t bin;
            istr >> bin >> abs_t_val[i];
            i++;   
        }
        ifs.close();
    }  
    Printf("Values of an avg |t| value per bin loaded.");

    return;
}

void ReadInput_HSModel()
{
    // read the input file for hot-spot model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/HSModel/data-dtdy-y_0.6-Run1.txt");
    for(Int_t i = 0; i < nData_HS; i++){
        Double_t x;
        ifs >> x;
        Double_t tmp;
        ifs >> tmp;
        ifs >> tmp;
        ifs >> abs_t_HS[i];
        ifs >> sig_HS_coh_n[i];
        ifs >> sig_HS_coh_n_err[i];        
        ifs >> sig_HS_inc_n[i];
        ifs >> sig_HS_inc_n_err[i];        
        ifs >> sig_HS_coh_hs[i];
        ifs >> sig_HS_coh_hs_err[i];      
        ifs >> sig_HS_inc_hs[i];
        ifs >> sig_HS_inc_hs_err[i];
        //std::cout << i << " " << abs_t[i] << " " <<sig_HS_inc_n[i]<< " " <<sig_HS_inc_hs[i] << endl;
    }
    ifs.close();
    Printf("Predictions of the HS model loaded.");

    return;
}

void ReadInput_Guzey()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Guzey/incoh_tdep_nuc_run2.dat");
    for(Int_t i = 0; i < nData_GZ; i++){
        ifs >> abs_t_GZ[i];
        ifs >> sig_GZ_el_min[i];
        ifs >> sig_GZ_el_max[i];
        ifs >> sig_GZ_diss_min[i];
        ifs >> sig_GZ_diss_max[i];
        ifs >> sig_GZ_tot_min[i];
        ifs >> sig_GZ_tot_max[i];
        //std::cout << i << " " << abs_t_GZ[i] << " " << sig_GZ_diss_min[i]<< " " << sig_GZ_diss_max[i] << endl;
    }
    ifs.close();
    Printf("Predictions of Guzey's model loaded.");

    return;
}

void ReadInput_Heikki()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_fluct");
    for(Int_t i = 0; i < nData_HM; i++){
        ifs >> abs_t_HM[i];
        ifs >> sig_HM_fluct[i];
        //std::cout << i << " " << abs_t_HM[i] << " " << sig_HM_fluct[i] << endl;
    }
    ifs.close();
    ifs.open("Trees/PhotoCrossSec/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_nofluct");
    for(Int_t i = 0; i < nData_HM; i++){
        ifs >> abs_t_HM[i];
        ifs >> sig_HM_noflu[i];
        //std::cout << i << " " << abs_t_HM[i] << " " << sig_HM_noflu[i] << endl;
    }
    ifs.close();
    Printf("Predictions of Heikki's model loaded.");

    return;
}

void ReadInput_STARlight()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/STARlight/IncJ_tDep_0.00-2.50.txt");
    for(Int_t i = 0; i < nData_SL; i++){
        ifs >> abs_t_SL[i];
        ifs >> sig_SL[i];
        //std::cout << i << " " << abs_t_SL[i] << " " << sig_SL[i] << endl;
    }
    ifs.close();
    Printf("STARlight predictions loaded.");

    return;
}

//#####################################################################################################
// Function to integrate a graph in a given range:

Double_t GraphIntegral_Calculate(TCanvas *c, TString str_name, Int_t n_data, Double_t *abs_t_val, Double_t *sig_val, Double_t t_min, Double_t t_max)
{
    vector<Double_t> t_edges; 
    vector<Double_t> sigmas;
    t_edges.push_back(t_min);
    Int_t iPoint = 0;
    Double_t new_t_edge = 0;
    // while below desired range of |t|
    while(abs_t_val[iPoint] <= t_min) iPoint++;
    // while within desired range of |t|
    while(abs_t_val[iPoint] <= t_max)
    {
        // save the value of the cross section at this point
        sigmas.push_back(sig_val[iPoint]);
        // if last point from the dataset
        if(iPoint == n_data-1)
        {
            /*
            // Previous implementation:
            new_t_edge = abs_t_val[iPoint] + (abs_t_val[iPoint] - t_edges.back());
            t_edges.push_back(new_t_edge);
            */
            t_edges.push_back(t_max);
            //Printf("%i \t%.4f \t%.5f", iPoint, new_t_edge, sig_val[iPoint]);
            break;
        } 
        if((abs_t_val[iPoint] + abs_t_val[iPoint+1]) / 2. < t_max)
        {
            new_t_edge = (abs_t_val[iPoint] + abs_t_val[iPoint+1]) / 2.;
            t_edges.push_back(new_t_edge);
            //Printf("%i \t%.4f \t%.5f", iPoint, new_t_edge, sig_val[iPoint]);
            iPoint++;
        } 
        // if |t| extends beyond t_max before the dataset ends
        else 
        {
            new_t_edge = t_max;
            t_edges.push_back(new_t_edge);
            //Printf("%i \t%.4f \t%.5f", iPoint, new_t_edge, sig_val[iPoint]);
            break;
        }
    }

    // Histograms
    Double_t *t_edges_ptr;
    t_edges_ptr = &t_edges[0];
    TH1D *hist = new TH1D("hist","hist",sigmas.size(),t_edges_ptr);
    for(unsigned i = 1; i <= sigmas.size(); i++) hist->SetBinContent(i, sigmas[i-1]);
    // Graph
    TGraph *graph = new TGraph(n_data, abs_t_val, sig_val);

    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Plots
    c->cd();
    c->SetLogy(); 
    // Margins
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.12);
    // Histogram settings
    hist->SetTitle(";|#it{t}| (GeV^{2}); d#sigma_{#gammaPb}/d|#it{t}| (mb/GeV^{2})");    
    hist->SetLineColor(kBlue);
    hist->SetLineWidth(2);
    // Vertical axis
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetLabelSize(0.05);
    //hist->GetYaxis()->SetRangeUser(1e-8,1e-1);
    // Horizontal axis
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.05);
    //hist->GetXaxis()->SetRangeUser(t_min,t_max);
    // Draw histogram
    hist->Draw("");
    // Draw graph
    graph->SetLineStyle(2);
    graph->SetLineColor(kRed);
    graph->SetLineWidth(2);  
    graph->Draw("C SAME");
    // Calculate integral
    Double_t integral_histo = 0;
    Double_t integral_graph = 0;
    for(unsigned i = 1; i <= sigmas.size(); i++){
        integral_graph += sigmas[i-1] * (t_edges[i] - t_edges[i-1]);
        integral_histo += hist->GetBinContent(i) * (hist->GetBinLowEdge(i+1) - hist->GetBinLowEdge(i));
    }
    if(TMath::Abs(integral_graph - integral_histo) > 1e-5){
        Printf("Error calculating integral.");
        return -1;
    }
    // Legend
    TLegend *l = new TLegend(0.55,0.75,0.80,0.95);
    l->SetMargin(0.);
    l->AddEntry((TObject*)0,Form("%s",str_name.Data()), "");
    l->AddEntry((TObject*)0,Form("In range |#it{t}| #in (%.2f,%.2f) GeV^{2} #it{c}^{-2}:", t_min, t_max), "");
    l->AddEntry((TObject*)0,Form("integral = %.3f #mub", integral_graph * 1000.), "");
    l->SetTextSize(0.045);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    

    return integral_graph;
}

//#####################################################################################################
// Roman's functions:

TLegend *SetLegend(Double_t x_leftdown, Double_t y_leftdown, Double_t x_rightup, Double_t y_rightup)
{
    TLegend *leg = new TLegend(x_leftdown, y_leftdown, x_rightup, y_rightup);
    //leg->SetFillColor(0);
    //leg->SetFillStyle(0);
    //leg->SetBorderSize(0);
    return leg;
}

void SetPadMargins(TVirtualPad* pad, Double_t left, Double_t top, Double_t right, Double_t bottom)
{
    pad->SetTopMargin(top);
    pad->SetBottomMargin(bottom);
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
}

void SetFrame(TH1* frame)
{
    frame->GetXaxis()->SetTitleOffset(1.);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetNdivisions(306);
    frame->GetXaxis()->SetNoExponent();
    //frame->GetYaxis()->SetNoExponent();
}

void SetStyle(TGraph* g, Color_t color, Style_t style, Width_t width = 3)
{
    g->SetLineColor(color);
    g->SetLineStyle(style);
    g->SetLineWidth(width);
}

void SetupSysErrorBox(TGraph* g, Color_t color)
{
    g->SetMarkerSize(0);
    g->SetFillStyle(1001);
    g->SetFillColorAlpha(color,0.35);
    SetStyle(g,color,1,0);
}