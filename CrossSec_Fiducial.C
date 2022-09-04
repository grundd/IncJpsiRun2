// CrossSec_Fiducial.C
// David Grund, Sep 04, 2022

// root headers
#include "TSystem.h"
#include "TAxis.h"
// my headers
#include "CrossSec_Utilities.h"

void PlotFiducialCrossSection();

void CrossSec_Fiducial(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    InitObjects();

    // load graphs
    LoadGraphs_data();
    LoadGraphs_SL();
    LoadGraphs_CCK();
    LoadGraphs_MS();
    LoadGraphs_GSZ();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/Fiducial/");

    PlotFiducialCrossSection();

    return;
}

void PlotFiducialCrossSection()
{
    Double_t MarkerSize = 2.;
    // integrate data in 0.04 < |t| < 1.0 GeV^2
    Double_t int_data(0.);
    Double_t int_uncr_low(0), int_uncr_upp(0);
    Double_t int_corr_low(0), int_corr_upp(0);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        Double_t t_low = gr_data_uncr->GetPointX(i) - gr_data_uncr->GetErrorXlow(i);
        Double_t t_upp = gr_data_uncr->GetPointX(i) + gr_data_uncr->GetErrorXhigh(i);
        int_data     += gr_data_uncr->GetPointY(i) * (t_upp - t_low);
        int_uncr_low += (gr_data_uncr->GetPointY(i) - gr_data_uncr->GetErrorYlow(i)) * (t_upp - t_low);
        int_uncr_upp += (gr_data_uncr->GetPointY(i) + gr_data_uncr->GetErrorYhigh(i)) * (t_upp - t_low);
        int_corr_low += (gr_data_corr->GetPointY(i) - gr_data_corr->GetErrorYlow(i)) * (t_upp - t_low);
        int_corr_upp += (gr_data_corr->GetPointY(i) + gr_data_corr->GetErrorYhigh(i)) * (t_upp - t_low);
    }
    Double_t err_uncr_low = int_data - int_uncr_low;
    Double_t err_corr_low = int_data - int_corr_low;
    Double_t err_uncr_upp = int_uncr_upp - int_data;
    Double_t err_corr_upp = int_corr_upp - int_data;
    Double_t err_tot = TMath::Sqrt(TMath::Power(err_uncr_low, 2) + TMath::Power(err_corr_low, 2));
    Printf(" +++++++++++++++++++++++++++++++++++++++");
    Printf(" err uncr low: %.3f", err_uncr_low * 1e3);
    Printf(" err uncr upp: %.3f", err_uncr_upp * 1e3);
    Printf(" err corr low: %.3f", err_corr_low * 1e3);
    Printf(" err corr upp: %.3f", err_corr_upp * 1e3);
    Printf(" data: (%.3f pm %.3f (uncr.) pm %.3f (corr.)) mub", int_data * 1e3, err_uncr_low * 1e3, err_corr_low * 1e3);
    Printf(" data: (%.3f pm %.3f) mub (uncertainties added in quadrature)", int_data * 1e3, err_tot * 1e3);
    Printf(" +++++++++++++++++++++++++++++++++++++++");

    // graph with the data point and uncorrelated uncertainty (blue bar)
    TGraphErrors *gr_data = new TGraphErrors(); 
    gr_data->SetPoint(0, int_data * 1e3, 8.);
    gr_data->SetPointError(0, err_uncr_low * 1e3, 0.);
    gr_data->SetMarkerStyle(kFullSquare);
    gr_data->SetMarkerColor(215);
    gr_data->SetMarkerSize(MarkerSize);
    gr_data->SetLineColor(215);
    gr_data->SetLineWidth(3.);
    // graph with the correlated uncertainty (blue colored box)
    Double_t arr_err_corr_x[4] = {int_corr_low * 1e3, int_corr_low * 1e3, int_corr_upp * 1e3, int_corr_upp * 1e3};
    Double_t arr_err_corr_y[4] = {0., 9., 9., 0.};
    TGraph *gr_err_corr = new TGraph(4,arr_err_corr_x,arr_err_corr_y);
    gr_err_corr->SetFillStyle(3345);
    gr_err_corr->SetFillColor(215);
    // graph with the total uncertainty (gray colored box)
    Double_t arr_err_tot_x[4] = {(int_data - err_tot) * 1e3, (int_data - err_tot) * 1e3, (int_data + err_tot) * 1e3, (int_data + err_tot) * 1e3};
    Double_t arr_err_tot_y[4] = {0., 9., 9., 0.};
    TGraph *gr_err_tot = new TGraph(4,arr_err_tot_x,arr_err_tot_y);
    gr_err_tot->SetFillStyle(3354);
    gr_err_tot->SetFillColor(kGray+2);
    // models
    Double_t integrals[7] = { 0 };
    for(Int_t i = 0; i < 7; i++)
    {
        CreateHistogramFromGraph(i);
        Double_t integral(0.), avgt(0.);
        IntegrateModel(i,0.04,1.00,integral,avgt);
        integrals[i] = integral * 1e3;
    }
    TGraph *gr_models = new TGraph();
    Double_t y = 7.;
    for(Int_t i = 0; i < 7; i++)
    {
        gr_models->SetPoint(i,integrals[i],y);
        gr_models->SetMarkerStyle(kFullCross);
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
    // set range on x-axis
    // https://root-forum.cern.ch/t/setrangeuser-on-tgraphs/8213
    TAxis *axis = gr_err_tot->GetXaxis();
    Double_t x_min = 1.0;
    Double_t x_max = 19.;
    axis->SetLimits(x_min,x_max);
    // make the plot 
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
    gr_err_corr->Draw("F SAME");
    gr_models->Draw("P SAME");
    for(Int_t i = 0; i < 7; i++) gr_data->Draw("P SAME");

    Double_t y_line[4] = {7.5, 6.5, 4.5, 2.5};
    TLine *line[4] = { NULL };
    for(Int_t i = 0; i < 4; i++)
    {
        line[i] = new TLine(x_min,y_line[i],x_max,y_line[i]);
        line[i]->SetLineColor(kBlack);
        line[i]->SetLineWidth(1);
        line[i]->SetLineStyle(7);
        line[i]->Draw("SAME");
    }

    TLatex *latex[8] = { NULL };
    Double_t y_step = 0.08;
    for(Int_t i = 0; i < 8; i++){
        latex[i] = new TLatex(); 
        latex[i]->SetTextSize(0.055);
        // https://root-forum.cern.ch/t/settextalign/7458
        latex[i]->SetTextAlign(12);
        if(i == 0) latex[i]->DrawLatex(15.0,8.0-i,Form("#bf{#color[215]{%s}}", "ALICE"));
        else latex[i]->DrawLatex(15.0,8.0-i,Form("#bf{%s}", str_models[i-1].Data()));
    }  

    TString path = "Results/" + str_subfolder + "CrossSec/Fiducial/fiducial";
    c->Print((path + ".pdf").Data());

    gr_err_corr->Print();

    return;
}