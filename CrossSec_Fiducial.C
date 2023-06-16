// CrossSec_Fiducial.C
// David Grund, Sep 04, 2022

// root headers
#include "TSystem.h"
#include "TAxis.h"
// my headers
#include "CrossSec_Utilities.h"

void PlotFiducialCrossSection(bool onlyPaperModels);

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

    PlotFiducialCrossSection(kFALSE);
    PlotFiducialCrossSection(kTRUE);

    return;
}

void PlotFiducialCrossSection(bool onlyPaperModels)
{
    // load the value of the fiducial cross section calculated directly:
    Double_t val, err_stat, err_syst_uncr, err_syst_corr;
    ifstream ifs("Results/" + str_subfolder + "CrossSec/CrossSec_fiducial_dir.txt");
    ifs >> val >> err_stat >> err_syst_uncr >> err_syst_corr;
    ifs.close();
    cout << "Loaded values from the direct calculation:\n"
         << Form("%.3f pm %.3f (stat.) pm %.3f (syst. uncr.) pm %.3f (syst. corr.)\n",
            val, err_stat, err_syst_uncr, err_syst_corr);

    Double_t MarkerSize = 1.5;
    // integrate the measurement in 0.04 < |t| < 1.0 GeV^2:
    Double_t int_data(0.);
    Double_t int_uncr_low(0), int_uncr_upp(0);
    Double_t int_corr_low(0), int_corr_upp(0);
    Double_t int_stat_low(0), int_stat_upp(0);
    Double_t int_syst_uncr_low(0), int_syst_uncr_upp(0);
    Double_t int_syst_corr_low(0), int_syst_corr_upp(0);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        Double_t t_low = gr_data_uncr->GetPointX(i) - gr_data_uncr->GetErrorXlow(i);
        Double_t t_upp = gr_data_uncr->GetPointX(i) + gr_data_uncr->GetErrorXhigh(i);
        int_data     +=  gr_data_uncr->GetPointY(i) * (t_upp - t_low);
        int_uncr_low += (gr_data_uncr->GetPointY(i) - gr_data_uncr->GetErrorYlow(i)) * (t_upp - t_low);
        int_uncr_upp += (gr_data_uncr->GetPointY(i) + gr_data_uncr->GetErrorYhigh(i)) * (t_upp - t_low);
        int_corr_low += (gr_data_corr->GetPointY(i) - gr_data_corr->GetErrorYlow(i)) * (t_upp - t_low);
        int_corr_upp += (gr_data_corr->GetPointY(i) + gr_data_corr->GetErrorYhigh(i)) * (t_upp - t_low);
        int_stat_low += (gr_data_stat->GetPointY(i) - gr_data_stat->GetErrorYlow(i)) * (t_upp - t_low);
        int_stat_upp += (gr_data_stat->GetPointY(i) + gr_data_stat->GetErrorYhigh(i)) * (t_upp - t_low);
        int_syst_uncr_low += (gr_data_syst_uncr->GetPointY(i) - gr_data_syst_uncr->GetErrorYlow(i)) * (t_upp - t_low);
        int_syst_uncr_upp += (gr_data_syst_uncr->GetPointY(i) + gr_data_syst_uncr->GetErrorYhigh(i)) * (t_upp - t_low);
        int_syst_corr_low += (gr_data_syst_corr->GetPointY(i) - gr_data_syst_corr->GetErrorYlow(i)) * (t_upp - t_low);
        int_syst_corr_upp += (gr_data_syst_corr->GetPointY(i) + gr_data_syst_corr->GetErrorYhigh(i)) * (t_upp - t_low);
    }
    // errors from the integration:
    Double_t err_uncr_low = int_data - int_uncr_low;
    Double_t err_uncr_upp = int_uncr_upp - int_data;
    Double_t err_corr_low = int_data - int_corr_low;
    Double_t err_corr_upp = int_corr_upp - int_data;
    Double_t err_stat_low = int_data - int_stat_low;
    Double_t err_stat_upp = int_stat_upp - int_data;
    Double_t err_syst_uncr_low = int_data - int_syst_uncr_low;
    Double_t err_syst_uncr_upp = int_syst_uncr_upp - int_data;
    Double_t err_syst_corr_low = int_data - int_syst_corr_low;
    Double_t err_syst_corr_upp = int_syst_corr_upp - int_data;
    Double_t err_tot_int = TMath::Sqrt(TMath::Power(err_uncr_low, 2) + TMath::Power(err_corr_low, 2));
    // errors from the direct calculation -- scale them to the integrated value:
    err_stat = err_stat * int_data * 1e3 / val;
    err_syst_uncr = err_syst_uncr * int_data * 1e3 / val;
    err_syst_corr = err_syst_corr * int_data * 1e3 / val;
    Double_t err_uncr = TMath::Sqrt(TMath::Power(err_stat, 2) + TMath::Power(err_syst_uncr, 2));
    Double_t err_corr = err_syst_corr;
    Double_t err_tot_dir = TMath::Sqrt(TMath::Power(err_uncr, 2) + TMath::Power(err_corr, 2));
    Double_t err_syst = TMath::Sqrt(TMath::Power(err_syst_uncr, 2) + TMath::Power(err_syst_corr, 2));
    // print the results
    cout << " +++++++++++++++++++++++++++++++++++++++ \n"
         << " errors from the integration: [mub] \n"
         << Form(" uncr low: %.3f, upp: %.3f \n", err_uncr_low * 1e3, err_uncr_upp * 1e3)
         << Form(" corr low: %.3f, upp: %.3f \n", err_corr_low * 1e3, err_corr_upp * 1e3)
         << Form(" stat low: %.3f, upp: %.3f \n", err_stat_low * 1e3, err_stat_upp * 1e3)
         << Form(" syst uncr low: %.3f, upp: %.3f \n", err_syst_uncr_low * 1e3, err_syst_uncr_upp * 1e3) 
         << Form(" syst corr low: %.3f, upp: %.3f \n", err_syst_corr_low * 1e3, err_syst_corr_upp * 1e3)
         << " results - errors from the integration: \n"
         << Form(" sigma_gPb = (%.3f pm %.3f (stat.) pm %.3f (uncr. syst.) pm %.3f (corr. syst.)) mub \n", 
                int_data * 1e3, err_stat_low * 1e3, err_syst_uncr_low * 1e3, err_syst_corr_low * 1e3)
         << Form(" sigma_gPb = (%.3f pm %.3f (uncr.) pm %.3f (corr.)) mub \n", int_data * 1e3, err_uncr_low * 1e3, err_corr_low * 1e3)
         << Form(" sigma_gPb = (%.3f pm %.3f) mub (uncertainties added in quadrature) \n", int_data * 1e3, err_tot_int * 1e3)
         << " results - errors from the direct calculation: \n"
         << Form(" sigma_gPb = (%.3f pm %.3f (stat.) pm %.3f (uncr. syst.) pm %.3f (corr. syst.)) mub \n", 
                int_data * 1e3, err_stat, err_syst_uncr, err_syst_corr)
         << Form(" sigma_gPb = (%.3f pm %.3f (stat.) pm %.3f (syst.)) mub \n", 
                int_data * 1e3, err_stat, err_syst)
         << Form(" sigma_gPb = (%.3f pm %.3f (uncr.) pm %.3f (corr.)) mub \n", int_data * 1e3, err_uncr, err_corr)
         << Form(" sigma_gPb = (%.3f pm %.3f) mub (uncertainties added in quadrature) \n", int_data * 1e3, err_tot_dir)
         << " +++++++++++++++++++++++++++++++++++++++ \n";

    // graph with the data point and uncorrelated uncertainty (blue bar)
    TGraphErrors *gr_data = new TGraphErrors(); 
    gr_data->SetPoint(0, int_data * 1e3, 8.);
    gr_data->SetPointError(0, err_stat, 0.);
    gr_data->SetMarkerStyle(kFullSquare);
    gr_data->SetMarkerColor(215);
    gr_data->SetMarkerSize(MarkerSize);
    gr_data->SetLineColor(215);
    gr_data->SetLineWidth(3);
    // graph with the systematic uncertainty (blue colored box)
    Double_t arr_err_corr_x[4] = {int_data * 1e3 - err_syst, int_data * 1e3 - err_syst, int_data * 1e3 + err_syst, int_data * 1e3 + err_syst};
    Double_t arr_err_corr_y[4] = {0., 9., 9., 0.};
    TGraph *gr_err_syst = new TGraph(4,arr_err_corr_x,arr_err_corr_y);
    gr_err_syst->SetFillStyle(3345);
    gr_err_syst->SetFillColor(215);
    // graph with the total uncertainty (gray colored box)
    Double_t arr_err_tot_x[4] = {int_data * 1e3 - err_tot_dir, int_data * 1e3 - err_tot_dir, int_data * 1e3 + err_tot_dir, int_data * 1e3 + err_tot_dir};
    Double_t arr_err_tot_y[4] = {0., 9., 9., 0.};
    TGraph *gr_err_tot = new TGraph(4,arr_err_tot_x,arr_err_tot_y);
    gr_err_tot->SetFillStyle(1001);
    gr_err_tot->SetFillColorAlpha(kGray+2,0.35);
    // models
    Double_t integrals[7] = { 0 };
    for(Int_t i = 0; i < 7; i++)
    {
        CreateHistogramFromGraph(i);
        Double_t integral(0.), avgt(0.);
        IntegrateModel(i,0.04,1.00,integral,avgt);
        integrals[i] = integral * 1e3;
    }
    gStyle->SetEndErrorSize(4); 
    TGraphErrors *gr_models = new TGraphErrors();
    Double_t y = 7.;
    for(Int_t i = 0; i < 7; i++)
    {
        if(onlyPaperModels) if(i < 3) continue;
        gr_models->SetPoint(i,integrals[i],y);
        if(i == 5 || i == 6) gr_models->SetPointError(i,integrals[i]*(1 - GSZ_err_scale_low[i-5]),0.);
        else                 gr_models->SetPointError(i,0.,0.);
        gr_models->SetMarkerStyle(kFullCircle);
        gr_models->SetMarkerColor(kBlack);
        gr_models->SetMarkerSize(MarkerSize);
        gr_models->SetLineWidth(3);
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
    int ySize = 600;
    double textSize = 0.05;
    if(onlyPaperModels) {
        ySize = 400;
        textSize = 0.07;
    }
    TCanvas *c = new TCanvas("c","c",700,ySize);
    TPad *pL = new TPad("pL","pL",0.0,0.0,0.03,1.0);
    pL->Draw();
    TPad *pR = new TPad("pR","pR",0.03,0.0,1.0,1.0);
    pR->Draw();
    pR->cd();
    // Margins
    pR->SetTopMargin(0.03);
    pR->SetBottomMargin(0.14);
    if(onlyPaperModels) {
        pR->SetBottomMargin(0.18);
        gr_err_tot->GetXaxis()->SetTitleSize(textSize+0.01);
        gr_err_tot->GetXaxis()->SetLabelSize(textSize+0.01);
    }
    pR->SetRightMargin(0.03);
    pR->SetLeftMargin(0.0); 
    // Draw points
    if(onlyPaperModels) gr_err_tot->GetYaxis()->SetRangeUser(3.,9.);
    gr_err_tot->Draw("AF");
    gr_err_syst->Draw("F SAME");
    gr_models->Draw("PZ SAME");
    gr_data->Draw("PZ SAME");

    Double_t y_line_1[4] = {7.5, 6.5, 4.5, 2.5};
    Double_t y_line_2[4] = {7.5, 5.5, 0., 0.};
    TLine *line[4] = { NULL };
    for(Int_t i = 0; i < 4; i++)
    {
        if(onlyPaperModels) line[i] = new TLine(x_min,y_line_2[i],x_max,y_line_2[i]);
        else line[i] = new TLine(x_min,y_line_1[i],x_max,y_line_1[i]);
        line[i]->SetLineColor(kBlack);
        line[i]->SetLineWidth(1);
        line[i]->SetLineStyle(7);
        line[i]->Draw("SAME");
    }

    TLatex *latex[8] = { NULL };
    y = 8.0;
    for(Int_t i = 0; i < 8; i++){
        if(onlyPaperModels) if(i < 4 && i != 0) continue;
        latex[i] = new TLatex(); 
        latex[i]->SetTextSize(textSize+0.01);
        // https://root-forum.cern.ch/t/settextalign/7458
        latex[i]->SetTextAlign(12);
        if(i == 0) latex[i]->DrawLatex(14.2,y,Form("#bf{#color[215]{%s}}", "ALICE"));
        else latex[i]->DrawLatex(14.2,y,Form("#bf{%s}", str_models[i-1].Data()));
        y = y - 1;
    }  

    TString path = "Results/" + str_subfolder + "CrossSec/Fiducial/";
    if(onlyPaperModels) path += "fiducial_paperModels.pdf";
    else path += "fiducial_all.pdf";
    c->Print(path.Data());

    if(onlyPaperModels) {
        TLegend *ltw = new TLegend(0.02,0.87,0.20,0.93);
        ltw->AddEntry((TObject*)0,"#bf{This work}","");
        ltw->SetMargin(0.);
        ltw->SetTextSize(textSize);
        ltw->SetBorderSize(0);
        ltw->SetFillStyle(0);
        ltw->Draw();

        path = "Results/" + str_subfolder + "_rozprava/fiducial.pdf";
        c->Print(path.Data());
    }

    gr_err_syst->Print();

    path = "Results/" + str_subfolder + "CrossSec/CrossSec_fiducial_int.txt";
    ofstream outfile(path.Data());
    outfile << std::fixed << std::setprecision(3)
            << int_data * 1e3 << "\t"
            << err_stat_low * 1e3 << "\t"
            << err_syst_uncr_low * 1e3 << "\t"
            << err_syst_corr_low * 1e3 << "\n";
    outfile.close();
    Printf("Results printed to %s.", path.Data());

    return;
}