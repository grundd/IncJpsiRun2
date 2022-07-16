// PhotoCrossSec_Plot.C
// David Grund, Apr 24, 2022

// root headers
#include "TSystem.h"
// my headers
#include "PhotoCrossSec_Utilities.h"

void Plot();

void PhotoCrossSec_Plot(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/Plot/");

    Plot();

    return;
}

void Plot()
{
    ReadInput_Measurement();

    ReadInput_HSModel();

    ReadInput_Guzey();

    ReadInput_Heikki();

    ReadInput_STARlight();

    // Fill the histogram
    for(Int_t i = 0; i < nPtBins; i++){
        // from microbarns to milibarns
        sig_val[i] = sig_val[i] / 1000.;
        sig_err_stat[i] = sig_err_stat[i] / 1000.;
        sig_err_syst[i] = sig_err_syst[i] / 1000.;
        sig_err_stat_low[i] = sig_err_stat[i];
        sig_err_stat_upp[i] = sig_err_stat[i];
        sig_err_syst_low[i] = sig_err_syst[i];
        sig_err_syst_upp[i] = sig_err_syst[i];
        Printf("Cross section in bin %i: %.5f pm %.5f(stat.) pm %.5f(syst.).", i+1, sig_val[i], sig_err_stat_low[i], sig_err_syst_low[i]);
        abs_t_err_low[i] = abs_t_val[i] - t_boundaries[i];
        abs_t_err_upp[i] = t_boundaries[i+1] - abs_t_val[i];
    }
    TGraphAsymmErrors *grData_stat = new TGraphAsymmErrors(nPtBins,abs_t_val,sig_val,abs_t_err_low,abs_t_err_upp,sig_err_stat_low,sig_err_stat_upp);
    TGraphAsymmErrors *grData_syst = new TGraphAsymmErrors(nPtBins,abs_t_val,sig_val,abs_t_err_low,abs_t_err_upp,sig_err_syst_low,sig_err_syst_upp);
    // with stat errors
    gStyle->SetEndErrorSize(4);         
    grData_stat->SetMarkerStyle(kFullCircle);
    grData_stat->SetMarkerSize(0.7);
    grData_stat->SetLineColor(kBlack);
    grData_stat->SetLineWidth(2);
    grData_stat->SetMarkerColor(kBlack);
    // with syst errors 
    grData_syst->SetFillColor(17);
    grData_syst->SetMarkerSize(0);
    grData_syst->SetMarkerStyle(1);
    grData_syst->SetMarkerColor(kBlack);
    grData_syst->SetLineWidth(0); // to have no line around the box

    // STARlight
    TGraph *gr_SL = new TGraph(nData_SL, abs_t_SL, sig_SL);
    gr_SL->SetLineColor(kBlue);
    gr_SL->SetLineStyle(1);
    gr_SL->SetLineWidth(lineWidth); 

    // Hot-spot model
    // With subnucleonic degrees of freedom (hot spots):
    TGraphErrors *gr_HS_hs = new TGraphErrors(nData_HS,abs_t_HS,sig_HS_inc_hs,NULL,sig_HS_inc_hs_err);
    // https://root.cern.ch/doc/master/classTGraphErrors.html#a0f51786d0f0e210869a53ab58c0a3ffb 
    // number of points; x-values; y-values; x-errors; y-errors
    gr_HS_hs->SetMarkerStyle(20);
    gr_HS_hs->SetMarkerColor(kRed+1);  
    gr_HS_hs->SetLineStyle(9);
    gr_HS_hs->SetLineColor(kRed+1);
    gr_HS_hs->SetLineWidth(lineWidth);  
    // Without subnucleonic degrees of freedom:
    TGraphErrors *gr_HS_n = new TGraphErrors(nData_HS,abs_t_HS,sig_HS_inc_n,NULL,sig_HS_inc_n_err);
    gr_HS_n->SetMarkerStyle(20);
    gr_HS_n->SetMarkerColor(kRed+1);
    gr_HS_n->SetLineStyle(8);
    gr_HS_n->SetLineColor(kRed+1);
    gr_HS_n->SetLineWidth(lineWidth);  

    // Guzey's model
    // First scale the values (Guzey uses nb instead of mb)
    for (Int_t i = 0; i < nData_GZ; i++){
        sig_GZ_tot_min[i] = sig_GZ_tot_min[i] / 1e6;
        sig_GZ_tot_max[i] = sig_GZ_tot_max[i] / 1e6;
    }
    // Then fill the graph
    // https://root.cern/doc/master/graphShade_8C.html
    TGraph *gr_GZ_min = new TGraph(nData_GZ, abs_t_GZ, sig_GZ_tot_min);
    TGraph *gr_GZ_max = new TGraph(nData_GZ, abs_t_GZ, sig_GZ_tot_max);
    TGraph *gr_GZ_area = new TGraph(2*nData_GZ);
    for (Int_t i = 0; i < nData_GZ; i++){
        gr_GZ_area->SetPoint(i, abs_t_GZ[i], sig_GZ_tot_max[i]);
        gr_GZ_area->SetPoint(nData_GZ+i, abs_t_GZ[nData_GZ-i-1], sig_GZ_tot_min[nData_GZ-i-1]);
    }
    gr_GZ_min->SetLineStyle(10);
    gr_GZ_min->SetLineColor(kGreen);
    gr_GZ_min->SetLineWidth(lineWidth);
    gr_GZ_max->SetLineStyle(10);
    gr_GZ_max->SetLineColor(kGreen);
    gr_GZ_max->SetLineWidth(lineWidth);
    SetupSysErrorBox(gr_GZ_area,kGreen);
    gr_GZ_area->SetFillStyle(3013);

    // Heikki's model
    // With fluctuations
    TGraph *gr_HM_fluct = new TGraph(nData_HM, abs_t_HM, sig_HM_fluct);
    gr_HM_fluct->SetLineStyle(9);
    gr_HM_fluct->SetLineColor(kGray+3);
    gr_HM_fluct->SetLineWidth(lineWidth);
    // Without fluctuations
    TGraph *gr_HM_noflu = new TGraph(nData_HM, abs_t_HM, sig_HM_noflu);
    gr_HM_noflu->SetLineStyle(8);
    gr_HM_noflu->SetLineColor(kGray+3);
    gr_HM_noflu->SetLineWidth(lineWidth);

    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Canvas
    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetLogy();  
    //c->SetLogx();
    // Margins
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.12);
    //Plot the graphs
    TH1 *h = (TH1*) gr_GZ_area->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    h->SetMinimum(1e-6);
    h->SetMaximum(0.1);
    // Vertical axis
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetRangeUser(0.04,1.0);
    // Draw everything
    // https://root.cern.ch/doc/master/classTGraphPainter.html
    // 1) Guzey's model (green shaded area)
    gr_GZ_area->Draw("AF");
    // 2) Systematic errors (gray shaded areas)
    grData_syst->Draw("5 SAME");
    // 3) Draw Guzey's model (lines)
    gr_GZ_min->Draw("L SAME");
    gr_GZ_max->Draw("L SAME"); 
    // 4) Draw STARlight curve
    gr_SL->Draw("L SAME");
    // 5) Draw hot-spot curves
    gr_HS_hs->Draw("CX SAME");
    gr_HS_n->Draw("CX SAME");
    // 6) Draw Heikki's curves
    gr_HM_fluct->Draw("CX SAME");
    gr_HM_noflu->Draw("CX SAME");    
    // 7) Draw data with statistic uncertainties
    grData_stat->Draw("P SAME");

    // Legend
    TLegend *l = new TLegend(0.15,0.18,0.48,0.56);
    l->AddEntry(gr_SL,"STARlight","L");
    l->AddEntry(gr_HM_fluct,"MS: IPsat flu.","L");
    l->AddEntry(gr_HM_noflu,"MS: IPsat no flu.","L");
    l->AddEntry(gr_GZ_area,"GSZ: el. + diss.","F");
    l->AddEntry(gr_HS_hs,"CCK: GG-hs","L");
    l->AddEntry(gr_HS_n,"CCK: GG-n","L");
    l->AddEntry(grData_stat,"ALICE measurement","EP");
    l->SetTextSize(0.048);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    TString path = "Results/" + str_subfolder + Form("PhotoCrossSec/Plot/plot_%ibins", nPtBins);
    c->Print((path + ".pdf").Data());
    c->Print((path + ".png").Data());

    // Plot the measurement only
    TCanvas *cMeas = new TCanvas("cMeas","cMeas",800,600);
    cMeas->SetLogy();  
    // Margins
    cMeas->SetTopMargin(0.03);
    cMeas->SetBottomMargin(0.14);
    cMeas->SetRightMargin(0.03);
    cMeas->SetLeftMargin(0.12);
    //Plot the graphs
    TH1 *h2 = (TH1*) grData_syst->GetHistogram();
    h2->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    //h->SetMinimum(1e-6);
    h->SetMaximum(0.2);
    // Vertical axis
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleOffset(1.1);
    h2->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetRangeUser(0.04,1.0);
    // Draw everything
    grData_syst->Draw("A5");
    grData_stat->Draw("P SAME");
    // legends
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.05);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.92,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // Draw legend with data+unc. description
    gStyle->SetLegendBorderSize(0);
    TLegend *leg2 = SetLegend(0.45,0.70,0.94,0.88);
    leg2->SetTextSize(0.05);
    leg2->SetMargin(0.15);
    leg2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |y|<0.8", "");
    leg2->AddEntry(grData_stat,"Experimental stat.", "EPL");
    leg2->AddEntry(grData_syst,"Experimental syst.", "F");
    leg2->Draw();

    TString pathMeas = "Results/" + str_subfolder + Form("PhotoCrossSec/Plot/measurement_%ibins", nPtBins);
    cMeas->Print((pathMeas + ".pdf").Data());
    cMeas->Print((pathMeas + ".png").Data());

    return;
}