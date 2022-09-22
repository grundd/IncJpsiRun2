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
    ReadInput_data();

    ReadInput_CCK();

    ReadInput_GSZ();

    ReadInput_MS();

    ReadInput_SL();

    // fill the histogram
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
    TGraphAsymmErrors *gr_data_stat = new TGraphAsymmErrors(nPtBins,abs_t_val,sig_val,abs_t_err_low,abs_t_err_upp,sig_err_stat_low,sig_err_stat_upp);
    TGraphAsymmErrors *gr_data_syst = new TGraphAsymmErrors(nPtBins,abs_t_val,sig_val,abs_t_err_low,abs_t_err_upp,sig_err_syst_low,sig_err_syst_upp);
    // with stat errors
    gStyle->SetEndErrorSize(4);         
    gr_data_stat->SetMarkerStyle(kFullCircle);
    gr_data_stat->SetMarkerSize(0.7);
    gr_data_stat->SetLineColor(kBlack);
    gr_data_stat->SetLineWidth(2);
    gr_data_stat->SetMarkerColor(kBlack);
    // with syst errors 
    gr_data_syst->SetFillColor(17);
    gr_data_syst->SetMarkerSize(0);
    gr_data_syst->SetMarkerStyle(1);
    gr_data_syst->SetMarkerColor(kBlack);
    gr_data_syst->SetLineWidth(0); // to have no line around the box

    // STARlight
    TGraph *gr_SL = new TGraph(n_SL, abs_t_SL, sig_SL);
    gr_SL->SetLineColor(kBlue);
    gr_SL->SetLineStyle(1);
    gr_SL->SetLineWidth(lineWidth); 

    // CCK (hot-spot) model
    // with subnucleonic degrees of freedom (hot spots):
    TGraphErrors *gr_CCK_hs = new TGraphErrors(n_CCK,abs_t_CCK,sig_CCK_inc_hs,NULL,sig_CCK_inc_hs_err);
    // https://root.cern.ch/doc/master/classTGraphErrors.html#a0f51786d0f0e210869a53ab58c0a3ffb 
    // number of points; x-values; y-values; x-errors; y-errors
    gr_CCK_hs->SetMarkerStyle(20);
    gr_CCK_hs->SetMarkerColor(kRed+1);  
    gr_CCK_hs->SetLineStyle(9);
    gr_CCK_hs->SetLineColor(kRed+1);
    gr_CCK_hs->SetLineWidth(lineWidth);  
    // without subnucleonic degrees of freedom:
    TGraphErrors *gr_CCK_n = new TGraphErrors(n_CCK,abs_t_CCK,sig_CCK_inc_n,NULL,sig_CCK_inc_n_err);
    gr_CCK_n->SetMarkerStyle(20);
    gr_CCK_n->SetMarkerColor(kRed+1);
    gr_CCK_n->SetLineStyle(8);
    gr_CCK_n->SetLineColor(kRed+1);
    gr_CCK_n->SetLineWidth(lineWidth);  

    // GSZ model
    // first scale the values (Guzey uses nb instead of mb)
    for (Int_t i = 0; i < n_GSZ; i++){
        sig_GSZ_el_min[i] = sig_GSZ_el_min[i] / 1e6;
        sig_GSZ_el_max[i] = sig_GSZ_el_max[i] / 1e6;
        sig_GSZ_tot_min[i] = sig_GSZ_tot_min[i] / 1e6;
        sig_GSZ_tot_max[i] = sig_GSZ_tot_max[i] / 1e6;
    }
    // then fill the graph
    // https://root.cern/doc/master/graphShade_8C.html
    // total cross section
    TGraph *gr_GSZ_tot_min = new TGraph(n_GSZ, abs_t_GSZ, sig_GSZ_tot_min);
    TGraph *gr_GSZ_tot_max = new TGraph(n_GSZ, abs_t_GSZ, sig_GSZ_tot_max);
    TGraph *gr_GSZ_tot_area = new TGraph(2*n_GSZ);
    for (Int_t i = 0; i < n_GSZ; i++){
        gr_GSZ_tot_area->SetPoint(i, abs_t_GSZ[i], sig_GSZ_tot_max[i]);
        gr_GSZ_tot_area->SetPoint(n_GSZ+i, abs_t_GSZ[n_GSZ-i-1], sig_GSZ_tot_min[n_GSZ-i-1]);
    }
    gr_GSZ_tot_min->SetLineStyle(10);
    gr_GSZ_tot_min->SetLineColor(kGreen);
    gr_GSZ_tot_min->SetLineWidth(lineWidth);
    gr_GSZ_tot_max->SetLineStyle(10);
    gr_GSZ_tot_max->SetLineColor(kGreen);
    gr_GSZ_tot_max->SetLineWidth(lineWidth);
    SetupSysErrorBox(gr_GSZ_tot_area,kGreen);
    // elastic only
    TGraph *gr_GSZ_el_min = new TGraph(n_GSZ, abs_t_GSZ, sig_GSZ_el_min);
    TGraph *gr_GSZ_el_max = new TGraph(n_GSZ, abs_t_GSZ, sig_GSZ_el_max);
    TGraph *gr_GSZ_el_area = new TGraph(2*n_GSZ);
    for (Int_t i = 0; i < n_GSZ; i++){
        gr_GSZ_el_area->SetPoint(i, abs_t_GSZ[i], sig_GSZ_el_max[i]);
        gr_GSZ_el_area->SetPoint(n_GSZ+i, abs_t_GSZ[n_GSZ-i-1], sig_GSZ_el_min[n_GSZ-i-1]);
    }
    gr_GSZ_el_min->SetLineStyle(10);
    gr_GSZ_el_min->SetLineColor(kCyan);
    gr_GSZ_el_min->SetLineWidth(lineWidth);
    gr_GSZ_el_max->SetLineStyle(10);
    gr_GSZ_el_max->SetLineColor(kCyan);
    gr_GSZ_el_max->SetLineWidth(lineWidth);
    SetupSysErrorBox(gr_GSZ_el_area,kCyan);

    // MS (IPsat) model
    // fluctuations
    TGraph *gr_MS_fluct = new TGraph(n_MS, abs_t_MS, sig_MS_fluct);
    gr_MS_fluct->SetLineStyle(9);
    gr_MS_fluct->SetLineColor(kGray+3);
    gr_MS_fluct->SetLineWidth(lineWidth);
    // no fluctuations
    TGraph *gr_MS_noflu = new TGraph(n_MS, abs_t_MS, sig_MS_noflu);
    gr_MS_noflu->SetLineStyle(8);
    gr_MS_noflu->SetLineColor(kGray+3);
    gr_MS_noflu->SetLineWidth(lineWidth);

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
    TH1 *h = (TH1*) gr_GSZ_tot_area->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    h->SetMinimum(0.0008);
    h->SetMaximum(0.06);
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
    // GSZ model (green shaded area)
    gr_GSZ_tot_area->Draw("AF");
    gr_GSZ_el_area->Draw("F SAME");
    // systematic errors (gray shaded areas)
    gr_data_syst->Draw("5 SAME");
    // GSZ model (lines)
    gr_GSZ_tot_min->Draw("L SAME");
    gr_GSZ_tot_max->Draw("L SAME"); 
    gr_GSZ_el_min->Draw("L SAME");
    gr_GSZ_el_max->Draw("L SAME"); 
    // STARlight
    gr_SL->Draw("L SAME");
    // CCK
    gr_CCK_hs->Draw("CX SAME");
    gr_CCK_n->Draw("CX SAME");
    // MS (IPsat)
    gr_MS_fluct->Draw("CX SAME");
    gr_MS_noflu->Draw("CX SAME");    
    // data with statistic uncertainties
    gr_data_stat->Draw("P SAME");

    // Legend
    TLegend *l = new TLegend(0.70,0.65,0.90,0.95);
    l->AddEntry(gr_SL,"STARlight","L");
    l->AddEntry(gr_MS_fluct,"MS: IPsat flu.","L");
    l->AddEntry(gr_MS_noflu,"MS: IPsat no flu.","L");
    l->AddEntry(gr_GSZ_tot_area,"GSZ: el. + diss.","F");
    l->AddEntry(gr_GSZ_el_area,"GSZ: el.","F");
    l->AddEntry(gr_CCK_hs,"CCK: GG-hs","L");
    l->AddEntry(gr_CCK_n,"CCK: GG-n","L");
    l->AddEntry(gr_data_stat,"ALICE measurement","EP");
    l->SetTextSize(0.038);
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
    cMeas->SetLeftMargin(0.13);
    //Plot the graphs
    TH1 *h2 = (TH1*) gr_data_syst->GetHistogram();
    h2->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    h2->SetMaximum(0.04);
    // Vertical axis
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleOffset(1.25);
    h2->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetRangeUser(0.04,1.0);
    // Draw everything
    gr_data_syst->Draw("A5");
    gr_data_stat->Draw("P SAME");
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
    leg2->AddEntry(gr_data_stat,"Experimental stat.", "EPL");
    leg2->AddEntry(gr_data_syst,"Experimental syst.", "F");
    leg2->Draw();

    TString pathMeas = "Results/" + str_subfolder + Form("PhotoCrossSec/Plot/measurement_%ibins", nPtBins);
    cMeas->Print((pathMeas + ".pdf").Data());
    cMeas->Print((pathMeas + ".png").Data());

    return;
}