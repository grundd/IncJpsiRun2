// PhotoCrossSec_Total.C
// David Grund, May 11, 2022

// cpp headers
#include <vector>
// root headers
#include "TSystem.h"
#include "TAxis.h"
// my headers
#include "PhotoCrossSec_Utilities.h"

Double_t integral_data(0), 
    integral_SL(0), 
    integral_HS_hs(0), integral_HS_n(0), 
    integral_MS_fl(0), integral_MS_nf(0),
    integral_GZ_up(0), integral_GZ_lo(0);
Double_t err_stat_low(0), err_syst_low(0), 
    err_stat_upp(0), err_syst_upp(0),
    err_tot(0);

Double_t GraphIntegral(TString str_name, Int_t n_data, Double_t *abs_t_val, Double_t *sig_val, Double_t t_min, Double_t t_max);
void GraphIntegral_All(Double_t t_min, Double_t t_max);
void PlotTotal();
void ModelsRatios(Double_t t_low_1, Double_t t_upp_1, Double_t t_low_2, Double_t t_upp_2);

void PhotoCrossSec_Total(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/Total/");

    PlotTotal();

    //ModelsRatios(0.04,1.00,0.00,2.00);

    return;
}

Double_t GraphIntegral(TString str_name, Int_t n_data, Double_t *abs_t_val, Double_t *sig_val, Double_t t_min = 0.04, Double_t t_max = 1.00)
{
    TCanvas *c = new TCanvas("c","c",900,600);

    Double_t integral = GraphIntegral_Calculate(c,str_name,n_data,abs_t_val,sig_val,t_min,t_max);

    Int_t oldLevel = gErrorIgnoreLevel; 
    gErrorIgnoreLevel = kWarning; 
    gSystem->Exec("mkdir -p Results/" + str_subfolder + Form("PhotoCrossSec/Total/%.2f-%.2f/", t_min, t_max));
    TString path = "Results/" + str_subfolder + Form("PhotoCrossSec/Total/%.2f-%.2f/%s", t_min, t_max, str_name.Data());
    c->Print((path + ".pdf").Data());
    c->Print((path + ".png").Data());
    gErrorIgnoreLevel = oldLevel; 

    Printf("%s: %.3f micro barns.", str_name.Data(), integral * 1e3);

    delete c;

    return integral;
}

void GraphIntegral_All(Double_t t_min, Double_t t_max)
{
    // STARlight
    ReadInput_STARlight();
    TString str_SL = "STARlight";
    integral_SL = GraphIntegral(str_SL,nData_SL,abs_t_SL,sig_SL, t_min, t_max);    

    // HS model
    ReadInput_HSModel();
    // GG-hs
    TString str_HS_hs = "CCK GG-hs";
    integral_HS_hs = GraphIntegral(str_HS_hs,nData_HS,abs_t_HS,sig_HS_inc_hs, t_min, t_max);
    // GG-n
    TString str_HS_n = "CCK GG-n";
    integral_HS_n = GraphIntegral(str_HS_n,nData_HS,abs_t_HS,sig_HS_inc_n, t_min, t_max);

    // Heikki's model
    ReadInput_Heikki();
    // IPsat fluctuations
    TString str_MS_fl = "MS IPsat flu";
    integral_MS_fl = GraphIntegral(str_MS_fl,nData_HM,abs_t_HM,sig_HM_fluct, t_min, t_max);
    // IPsat no fluctuations
    TString str_MS_nf = "MS IPsat no flu";
    integral_MS_nf = GraphIntegral(str_MS_nf,nData_HM,abs_t_HM,sig_HM_noflu, t_min, t_max);

    // Guzey's model
    ReadInput_Guzey();
    for (Int_t i = 0; i < nData_GZ; i++){
        sig_GZ_tot_min[i] = sig_GZ_tot_min[i] / 1e6;
        sig_GZ_tot_max[i] = sig_GZ_tot_max[i] / 1e6;
    }
    // Upper error
    TString str_GZ_up = "GSZ upp";
    integral_GZ_up = GraphIntegral(str_GZ_up,nData_GZ,abs_t_GZ,sig_GZ_tot_max, t_min, t_max);
    // Lower error
    TString str_GZ_lo = "GSZ low";
    integral_GZ_lo = GraphIntegral(str_GZ_lo,nData_GZ,abs_t_GZ,sig_GZ_tot_min, t_min, t_max);

    return;
}

void PlotTotal()
{
    Double_t MarkerSize = 2.;
    // Integrate data in 0.04 < |t| < 1.0 GeV^2
    ReadInput_Measurement();
    Double_t integral_stat_low(0), integral_stat_upp(0);
    Double_t integral_syst_low(0), integral_syst_upp(0);
    for(Int_t i = 0; i < nPtBins; i++){
        sig_val[i] = sig_val[i] / 1e3;
        sig_err_stat[i] = sig_err_stat[i] / 1e3;
        sig_err_syst[i] = sig_err_syst[i] / 1e3;
        integral_data += sig_val[i] * (t_boundaries[i+1] - t_boundaries[i]);
        integral_stat_low += (sig_val[i] - sig_err_stat[i]) * (t_boundaries[i+1] - t_boundaries[i]);
        integral_stat_upp += (sig_val[i] + sig_err_stat[i]) * (t_boundaries[i+1] - t_boundaries[i]);
        integral_syst_low += (sig_val[i] - sig_err_syst[i]) * (t_boundaries[i+1] - t_boundaries[i]);
        integral_syst_upp += (sig_val[i] + sig_err_syst[i]) * (t_boundaries[i+1] - t_boundaries[i]);
    }
    err_stat_low = integral_data - integral_stat_low;
    err_syst_low = integral_data - integral_syst_low;
    err_stat_upp = integral_stat_upp - integral_data;
    err_syst_upp = integral_syst_upp - integral_data;
    err_tot = TMath::Sqrt(TMath::Power(err_stat_low, 2) + TMath::Power(err_syst_low, 2));
    Printf("Err stat low: %.3f", err_stat_low * 1e3);
    Printf("Err stat upp: %.3f", err_stat_upp * 1e3);
    Printf("Err syst low: %.3f", err_syst_low * 1e3);
    Printf("Err syst upp: %.3f", err_syst_upp * 1e3);
    Printf("Data: (%.3f pm %.3f (stat.) pm %.3f (syst.)) micro barns.", 
        integral_data * 1e3, err_stat_low * 1e3, err_syst_low * 1e3);
    Printf("Data: (%.3f pm %.3f) micro barns (stat. and syst. added in quadrature).", 
        integral_data * 1e3, err_tot * 1e3);

    GraphIntegral_All(0.04, 1.00);

    // Graph with data point and stat uncertainty
    TGraphErrors *gr_data = new TGraphErrors(); 
    gr_data->SetPoint(0, integral_data * 1e3, 8.);
    gr_data->SetPointError(0, err_stat_low * 1e3, 0.);
    gr_data->SetMarkerStyle(kFullSquare);
    gr_data->SetMarkerColor(215);
    gr_data->SetMarkerSize(MarkerSize);
    gr_data->SetLineColor(215);
    gr_data->SetLineWidth(3.);
    // Graph with the systematic uncertainty (colored box)
    Double_t arr_err_syst_x[4] = {integral_syst_low * 1e3, integral_syst_low * 1e3, 
        integral_syst_upp * 1e3, integral_syst_upp * 1e3};
    Double_t arr_err_syst_y[4] = {0., 9., 9., 0.};
    TGraph *gr_err_syst = new TGraph(4,arr_err_syst_x,arr_err_syst_y);
    gr_err_syst->SetFillStyle(3004);
    gr_err_syst->SetFillColor(215);
    // Graph with the total uncertainty (colored box)
    Double_t arr_err_tot_x[4] = {(integral_data - err_tot) * 1e3, (integral_data - err_tot) * 1e3, 
        (integral_data + err_tot) * 1e3, (integral_data + err_tot) * 1e3};
    Double_t arr_err_tot_y[4] = {0., 9., 9., 0.};
    TGraph *gr_err_tot = new TGraph(4,arr_err_tot_x,arr_err_tot_y);
    gr_err_tot->SetFillStyle(3005);
    gr_err_tot->SetFillColor(kGray+2);
    // Models
    Double_t integrals[7] = {
        integral_HS_n,
        integral_HS_hs, 
        integral_GZ_lo,
        integral_GZ_up, 
        integral_MS_nf,
        integral_MS_fl, 
        integral_SL
    };
    TGraph *gr_models = new TGraph();
    Double_t y = 7.;
    for(Int_t i = 6; i >= 0; i--){
        gr_models->SetPoint(i,integrals[i] * 1e3,y);
        gr_models->SetMarkerStyle(kFullCircle);
        gr_models->SetMarkerColor(kBlack);
        gr_models->SetMarkerSize(MarkerSize);
        y = y - 1.;
    }  
    gr_err_tot->GetYaxis()->SetTickLength(0.0);
    gr_err_tot->GetYaxis()->SetRangeUser(0.,9.);
    gr_err_tot->GetXaxis()->SetTitle("#sigma_{#gammaPb} (#mub)");
    gr_err_tot->GetXaxis()->SetTitleSize(0.06);
    gr_err_tot->GetXaxis()->SetTitleOffset(1.05);
    gr_err_tot->GetXaxis()->SetLabelSize(0.06);
    // Set range on x-axis
    // https://root-forum.cern.ch/t/setrangeuser-on-tgraphs/8213
    TAxis *axis = gr_err_tot->GetXaxis();
    Double_t x_min = 1.0;
    Double_t x_max = 20.;
    axis->SetLimits(x_min,x_max);
    // Make the plot 
    TCanvas *c2 = new TCanvas("c2","c2",900,600);
    TPad *pL = new TPad("pL","pL",0.0,0.0,0.03,1.0);
    pL->Draw();
    TPad *pR = new TPad("pR","pR",0.03,0.0,1.0,1.0);
    pR->Draw();
    pR->cd();
    // Margins
    pR->SetTopMargin(0.03);
    pR->SetBottomMargin(0.14);
    pR->SetRightMargin(0.03);
    pR->SetLeftMargin(0.0);   
    // Draw points
    gr_err_tot->Draw("AF");
    gr_err_syst->Draw("F SAME");
    gr_models->Draw("P SAME");
    for(Int_t i = 0; i < 7; i++) gr_data->Draw("P SAME");

    Double_t y_line = 0;
    TLine *line[4] = { NULL };
    for(Int_t i = 0; i < 4; i++){
        if(i == 0) y_line = 7.5;
        if(i == 1) y_line = 6.5;
        if(i == 2) y_line = 4.5;
        if(i == 3) y_line = 2.5;
        line[i] = new TLine(x_min,y_line,x_max,y_line);
        line[i]->SetLineColor(kBlack);
        line[i]->SetLineWidth(1);
        line[i]->SetLineStyle(7);
        line[i]->Draw("SAME");
    }
    TLatex *latex[8] = { NULL };
    TString names[8] = {"ALICE","STARlight","MS IPsat flu.","MS IPsat no flu.","GSZ upper","GSZ lower","CCK GG-hs","CCK GG-n"};
    Double_t y_step = 0.08;
    for(Int_t i = 0; i < 8; i++){
        latex[i] = new TLatex(); 
        latex[i]->SetTextSize(0.06);
        // https://root-forum.cern.ch/t/settextalign/7458
        latex[i]->SetTextAlign(12);
        if(i == 0) latex[i]->DrawLatex(14.5,8.0-i,Form("#bf{#color[215]{%s}}", names[i].Data()));
        else latex[i]->DrawLatex(14.5,8.0-i,Form("#bf{%s}", names[i].Data()));
    }    

    TString path = "Results/" + str_subfolder + "PhotoCrossSec/Total/TotalCrossSection";
    c2->Print((path + ".pdf").Data());
    c2->Print((path + ".png").Data());

    return;
}

void ModelsRatios(Double_t t_low_1, Double_t t_upp_1, Double_t t_low_2, Double_t t_upp_2)
{
    GraphIntegral_All(t_low_1, t_upp_1);
    Double_t int_SL_1 = integral_SL;
    Double_t int_HS_hs_1 = integral_HS_hs;
    Double_t int_HS_n_1 = integral_HS_n;
    Double_t int_MS_fl_1 = integral_MS_fl;
    Double_t int_MS_nf_1 = integral_MS_nf;
    Double_t int_GZ_up_1 = integral_GZ_up;
    Double_t int_GZ_lo_1 = integral_GZ_lo;
    GraphIntegral_All(t_low_2, t_upp_2);
    Double_t int_SL_2 = integral_SL;
    Double_t int_HS_hs_2 = integral_HS_hs;
    Double_t int_HS_n_2 = integral_HS_n;
    Double_t int_MS_fl_2 = integral_MS_fl;
    Double_t int_MS_nf_2 = integral_MS_nf;
    Double_t int_GZ_up_2 = integral_GZ_up;
    Double_t int_GZ_lo_2 = integral_GZ_lo;

    Printf("Ratio SL: %.3f", int_SL_1/int_SL_2);
    Printf("Ratio HS hs: %.3f", int_HS_hs_1/int_HS_hs_2);
    Printf("Ratio HS n: %.3f", int_HS_n_1/int_HS_n_2);
    Printf("Ratio MS fl: %.3f", int_MS_fl_1/int_MS_fl_2);
    Printf("Ratio MS nf: %.3f", int_MS_nf_1/int_MS_nf_2);
    Printf("Ratio GZ up: %.3f", int_GZ_up_1/int_GZ_up_2);
    Printf("Ratio GZ lo: %.3f", int_GZ_lo_1/int_GZ_lo_2);

    TString path = "Results/" + str_subfolder + "PhotoCrossSec/Total/ratios.txt";
    ofstream fout_ratios(path.Data());
    fout_ratios << std::fixed << std::setprecision(3);
    fout_ratios << Form("Ratio STARlight: %.3f\n", int_SL_1/int_SL_2);
    fout_ratios << Form("Ratio HS GG-hs:  %.3f\n", int_HS_hs_1/int_HS_hs_2);
    fout_ratios << Form("Ratio HS GG-n:   %.3f\n", int_HS_n_1/int_HS_n_2);
    fout_ratios << Form("Ratio MS fluct:  %.3f\n", int_MS_fl_1/int_MS_fl_2);
    fout_ratios << Form("Ratio MS noflu:  %.3f\n", int_MS_nf_1/int_MS_nf_2);
    fout_ratios << Form("Ratio GZ upp:    %.3f\n", int_GZ_up_1/int_GZ_up_2);
    fout_ratios << Form("Ratio GZ low:    %.3f\n", int_GZ_lo_1/int_GZ_lo_2);
    fout_ratios.close();

    return;
}