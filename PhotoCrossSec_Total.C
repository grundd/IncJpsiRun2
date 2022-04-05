// PhotoCrossSec_Total.C
// David Grund, Apr 05, 2022

// cpp headers
#include <vector>
// root headers
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

Double_t GraphIntegral(TString str_name, Int_t n_data, Double_t *abs_t_val, Double_t *sig_val, Double_t t_min = 0.04, Double_t t_max = 1.00)
{
    vector<Double_t> t_edges; 
    vector<Double_t> sigmas;
    t_edges.push_back(t_min);
    Int_t iPoint = 0;
    while(abs_t_val[iPoint] <= t_min) iPoint++;
    while(abs_t_val[iPoint] < t_max){
        if(iPoint == n_data-1){
            t_edges.push_back(abs_t_val[iPoint] + (abs_t_val[iPoint] - t_edges.back()));
            sigmas.push_back(sig_val[iPoint]);
            break;
        } 
        if((abs_t_val[iPoint] + abs_t_val[iPoint+1]) / 2. < t_max){
            t_edges.push_back((abs_t_val[iPoint] + abs_t_val[iPoint+1]) / 2.);
            sigmas.push_back(sig_val[iPoint]);
            //Printf("%.4f \t%.5f", (abs_t_val[iPoint] + abs_t_val[iPoint+1]) / 2., sig_val[iPoint]);
            iPoint++;
        } else {
            t_edges.push_back(t_max);
            sigmas.push_back(sig_val[iPoint]);
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
    TCanvas *c = new TCanvas("c","c",900,600);
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
    hist->GetXaxis()->SetRangeUser(t_min,t_max);
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

    Int_t oldLevel = gErrorIgnoreLevel; 
    gErrorIgnoreLevel = kWarning; 
    c->Print(Form("PhotoCrossSec/img_Total/%.2f-%.2f/%s.pdf", t_min, t_max, str_name.Data()));
    c->Print(Form("PhotoCrossSec/img_Total/%.2f-%.2f/%s.png", t_min, t_max, str_name.Data()));
    gErrorIgnoreLevel = oldLevel; 

    Printf("%s: %.3f micro barns.", str_name.Data(), integral_graph * 1e3);

    delete c;
    delete hist;

    return integral_graph;
}

void GraphIntegralAll(Double_t t_min, Double_t t_max)
{
    // STARlight
    ReadInputSTARlight();
    TString str_SL = "STARlight";
    integral_SL = GraphIntegral(str_SL,nData_SL,abs_t_SL,sig_SL, t_min, t_max);    

    // HS model
    ReadInputHSModel();
    // GG-hs
    TString str_HS_hs = "CCK GG-hs";
    integral_HS_hs = GraphIntegral(str_HS_hs,nData_HS,abs_t_HS,sig_HS_inc_hs, t_min, t_max);
    // GG-n
    TString str_HS_n = "CCK GG-n";
    integral_HS_n = GraphIntegral(str_HS_n,nData_HS,abs_t_HS,sig_HS_inc_n, t_min, t_max);

    // Heikki's model
    ReadInputHeikki();
    // IPsat fluctuations
    TString str_MS_fl = "MS IPsat flu";
    integral_MS_fl = GraphIntegral(str_MS_fl,nData_HM,abs_t_HM,sig_HM_fluct, t_min, t_max);
    // IPsat no fluctuations
    TString str_MS_nf = "MS IPsat no flu";
    integral_MS_nf = GraphIntegral(str_MS_nf,nData_HM,abs_t_HM,sig_HM_noflu, t_min, t_max);

    // Guzey's model
    ReadInputGuzey();
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
    ReadInputMeasurement();
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

    GraphIntegralAll(0.04, 1.00);

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
    TCanvas *c = new TCanvas("c","c",900,600);
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

    c->Print("PhotoCrossSec/img_Total/TotalCrossSection.pdf");
    c->Print("PhotoCrossSec/img_Total/TotalCrossSection.png");

    return;
}

void ModelsRatios(Double_t t_low_1, Double_t t_upp_1, Double_t t_low_2, Double_t t_upp_2)
{
    GraphIntegralAll(t_low_1, t_upp_1);
    Double_t int_SL_1 = integral_SL;
    Double_t int_HS_hs_1 = integral_HS_hs;
    Double_t int_HS_n_1 = integral_HS_n;
    Double_t int_MS_fl_1 = integral_MS_fl;
    Double_t int_MS_nf_1 = integral_MS_nf;
    Double_t int_GZ_up_1 = integral_GZ_up;
    Double_t int_GZ_lo_1 = integral_GZ_lo;
    GraphIntegralAll(t_low_2, t_upp_2);
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

    ofstream fout_ratios("PhotoCrossSec/img_Total/ratios.txt");
    fout_ratios << std::fixed << std::setprecision(3);
    fout_ratios << Form("Ratio STARlight: %.3f\n", int_SL_1/int_SL_2);
    fout_ratios << Form("Ratio HS GG-hs:\t %.3f\n", int_HS_hs_1/int_HS_hs_2);
    fout_ratios << Form("Ratio HS GG-n:\t %.3f\n", int_HS_n_1/int_HS_n_2);
    fout_ratios << Form("Ratio MS fluct:\t %.3f\n", int_MS_fl_1/int_MS_fl_2);
    fout_ratios << Form("Ratio MS noflu:\t %.3f\n", int_MS_nf_1/int_MS_nf_2);
    fout_ratios << Form("Ratio GZ upp: \t %.3f\n", int_GZ_up_1/int_GZ_up_2);
    fout_ratios << Form("Ratio GZ low: \t %.3f\n", int_GZ_lo_1/int_GZ_lo_2);
    fout_ratios.close();

    return;
}

void PhotoCrossSec_Total()
{
    PlotTotal();

    //ModelsRatios(0.04,1.00,0.00,2.00);

    return;
}