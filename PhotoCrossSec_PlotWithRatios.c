// PhotoCrossSec_PlotWithRatios.c
// David Grund, Dec 1, 2021

// my headers
#include "PhotoCrossSec_Utilities.h"

Int_t lineWidth = 2;

void PlotWithRatios();

void PhotoCrossSec_PlotWithRatios()
{
    iFeedDown = 0; 
    PlotWithRatios();

    iFeedDown = 1; 
    PlotWithRatios();

    return;
}

void PlotWithRatios()
{
    ReadInputMeasurement();

    ReadInputHSModel();

    ReadInputGuzey();

    ReadInputHeikki();

    ReadInputSTARlight();

    // Fill the data graphs
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

    // Using Roman's settings:
    gStyle->SetErrorX(0.02);
    gStyle->SetLineScalePS(2.0);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleOffset(1.25,"XYZ");
    gStyle->SetTitleSize(0.04,"XYZ");
    gStyle->SetLabelSize(0.04,"XYZ");
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetTextSize(0.04);
    gStyle->SetTextFont(42);

    TCanvas *cCSont = new TCanvas ("cCSont","Cross section dependence on p_{t}^{2}",1050,800);
    SetPadMargins(gPad,0.13,0.03,0.03,0.11);
    gStyle->SetOptStat("0");
    gStyle->SetOptFit(0);

    // Draw frame
    TH1F* fCSont = gPad->DrawFrame(0.04,
        //0.2*TMath::MinElement(grData_stat->GetN(),grData_stat->GetY()), // GetN() = get number of points stored in TGraph
        0.00002,
        1.0,
        //1.4*TMath::MaxElement(grData_stat->GetN(),grData_stat->GetY()));
        0.08);
    SetFrame(fCSont);
    fCSont->GetYaxis()->SetTickLength(0.025); 
    fCSont->GetXaxis()->SetTickLength(0.025); 
    fCSont->GetYaxis()->SetTitleOffset(1.2);
    fCSont->SetTitle("Cross section dependence on |#it{t}|;|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    // Format: Title of the canvas, title of the X axis, title of the Y axis

    cCSont->cd();
    cCSont->SetLogy();
    cCSont->Modified();
    fCSont->Draw("AXIS");

    // Create boxes with systematic uncertainties
    SetupSysErrorBox(grData_syst,kGray);

    // Set data properties 
    gStyle->SetEndErrorSize(4);         
    grData_stat->SetMarkerStyle(kFullCircle);
    grData_stat->SetMarkerSize(0.7);
    grData_stat->SetLineColor(kBlack);
    grData_stat->SetLineWidth(2);
    grData_stat->SetMarkerColor(kBlack);

    // STARlight
    TGraph *gr_SL = new TGraph(nData_SL, abs_t_SL, sig_SL);
    gr_SL->SetLineColor(kBlue);
    gr_SL->SetLineStyle(1);
    gr_SL->SetLineWidth(lineWidth);  

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

    // Hot-spot model
    // With subnucleonic degrees of freedom (hot spots):
    TGraphErrors *gr_HS_hs = new TGraphErrors(nData_HS,abs_t_HS,sig_HS_inc_hs,NULL,sig_HS_inc_hs_err);
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

    // **********************************************************************
    // Draw everything
    // 1) Guzey's model (green shaded area)
    gr_GZ_area->Draw("F SAME");
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
    cCSont->Modified();
    cCSont->Update(); 

    // ALICE PbPb label
    //TLegend *leg0 = SetLegend(0.22,0.93,0.88,0.97);
    //leg0->Draw();
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.035);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // Draw legend with models
    TLegend *leg1 = SetLegend(0.17,0.14,0.42,0.40);
    leg1->SetTextSize(0.035);
    leg1->SetMargin(0.30);
    leg1->AddEntry(gr_SL,"STARlight", "L");
    leg1->AddEntry(gr_HM_fluct,"MS: IPsat flu.", "L");
    leg1->AddEntry(gr_HM_noflu,"MS: IPsat no flu.", "L");
    leg1->AddEntry(gr_GZ_area,"GSZ: el. + diss.", "F");
    leg1->AddEntry(gr_HS_hs,"CCK: GG-hs", "L");
    leg1->AddEntry(gr_HS_n, "CCK: GG-n", "L");
    leg1->Draw();
    cCSont->Modified();
    cCSont->Update();

    // Draw legend with data+unc. description
    TLegend *leg2 = SetLegend(0.59,0.76,0.94,0.90);
    leg2->SetTextSize(0.035);
    leg2->SetMargin(0.11);
    leg2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |y|<0.8", "");
    leg2->AddEntry(grData_stat,"Experimental stat.", "EPL");
    leg2->AddEntry(grData_syst,"Experimental syst.", "F");
    leg2->Draw();
    cCSont->Modified();
    cCSont->Update();

    cCSont->Print(Form("PhotoCrossSec/img_PlotWithRatios/Plot_FeedDown%i_%ibins.pdf", iFeedDown, nPtBins));
    cCSont->Print(Form("PhotoCrossSec/img_PlotWithRatios/Plot_FeedDown%i_%ibins.png", iFeedDown, nPtBins));

    // *****************************************************************************
    // Calculate and plot ratios
    TGraphErrors *grRatio_SL = new TGraphErrors(nPtBins);
    TGraphErrors *grRatio_GZ = new TGraphErrors(nPtBins);
    TGraphErrors *grRatio_HS_hs = new TGraphErrors(nPtBins);
    TGraphErrors *grRatio_HS_n = new TGraphErrors(nPtBins);
    TGraphErrors *grRatio_MS_fl = new TGraphErrors(nPtBins);
    TGraphErrors *grRatio_MS_nf = new TGraphErrors(nPtBins);
    TGraphAsymmErrors *gr_stat_data = new TGraphAsymmErrors(nPtBins);
    TGraphAsymmErrors *gr_syst_data = new TGraphAsymmErrors(nPtBins);

    SetupSysErrorBox(gr_syst_data,kGray);

    gr_stat_data->SetMarkerSize(0.);
    gr_stat_data->SetLineColor(kBlack);
    gr_stat_data->SetLineWidth(2);
    gr_stat_data->SetMarkerColor(kBlack);

    for(Int_t i = 0; i < nPtBins; ++i) {
        Double_t x,y,err_y_stat,err_y_syst,ey_gama,err_x_low,err_x_upp;
        grData_stat->GetPoint(i, x, y);
        err_y_stat = grData_stat->GetErrorYhigh(i);
        err_y_syst = grData_syst->GetErrorYhigh(i);
        err_x_upp = grData_syst->GetErrorXhigh(i);
        err_x_low = grData_syst->GetErrorXlow(i);
        
        grRatio_SL->SetPoint(i, x, gr_SL->Eval(x)/y);
        //grRatio_GZ -> how?
        grRatio_HS_hs->SetPoint(i, x, gr_HS_hs->Eval(x)/y);
        grRatio_HS_n->SetPoint(i, x, gr_HS_n->Eval(x)/y);
        grRatio_MS_fl->SetPoint(i, x, gr_HM_fluct->Eval(x)/y);
        grRatio_MS_nf->SetPoint(i, x, gr_HM_noflu->Eval(x)/y);

        Double_t ratio_correction = 1/y;
        err_y_stat *= ratio_correction;
        err_y_syst *= ratio_correction;

        gr_stat_data->SetPoint(i, x, 1.);
        gr_stat_data->SetPointError(i,0.,0.,err_y_stat,err_y_stat);
        gr_syst_data->SetPoint(i, x, 1.);
        gr_syst_data->SetPointError(i,err_x_low,err_x_upp,err_y_syst,err_y_syst);
    }

    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"y");
    TCanvas *cDataModel = new TCanvas ("cDataModel","",1600,300);
    SetPadMargins(gPad,0.05,0.03,0.03,0.12);
    TH1F* fCSratio = gPad->DrawFrame(0.04,0.0,1.0,4.5);
    //TH1F* fCSratio = gPad->DrawFrame(0.04,0.01,1.0,9.5);
    SetFrame(fCSratio);
    fCSratio->SetTitle("Ratios model/data;|#it{t}| (GeV^{2} #it{c}^{-2});Model / Data");
    fCSratio->GetYaxis()->SetTitleOffset(0.5);

    cDataModel->cd();
    cDataModel->Modified();
    fCSratio->Draw("AXIS");
    // Errors
    gr_syst_data->Draw("5 SAME");
    cDataModel->Modified();
    cDataModel->Update();
    // STARlight ratios
    grRatio_SL->SetLineColor(kBlue);
    grRatio_SL->SetMarkerColor(kBlue);
    grRatio_SL->SetMarkerStyle(kOpenSquare);
    grRatio_SL->Draw("SAME P");
    cDataModel->Modified();
    cDataModel->Update();
    // Hot-spot model ratios
    grRatio_HS_hs->SetLineColor(kRed+1);
    grRatio_HS_hs->SetMarkerColor(kRed+1);
    grRatio_HS_hs->SetMarkerStyle(kOpenCircle);
    grRatio_HS_hs->Draw("SAME P");
    grRatio_HS_n->SetLineColor(kRed+1);
    grRatio_HS_n->SetMarkerColor(kRed+1);
    grRatio_HS_n->SetMarkerStyle(kFullCircle);
    grRatio_HS_n->Draw("SAME P");
    cDataModel->Modified();
    cDataModel->Update();
    // IPsat model ratios
    grRatio_MS_fl->SetLineColor(kGray+3);
    grRatio_MS_fl->SetMarkerColor(kGray+3);
    grRatio_MS_fl->SetMarkerStyle(kOpenTriangleDown);
    grRatio_MS_fl->Draw("SAME P");
    grRatio_MS_nf->SetLineColor(kGray+3);
    grRatio_MS_nf->SetMarkerColor(kGray+3);
    grRatio_MS_nf->SetMarkerStyle(kFullTriangleDown);
    grRatio_MS_nf->Draw("SAME P");
    // Draw data with stat. errors
    gr_stat_data->Draw("SAME P");

    TLegend *leg3 = SetLegend(0.75,0.35,0.95,0.95);
    leg3->SetTextSize(0.035);
    leg3->AddEntry(grRatio_SL,"STARlight / Data", "P");
    leg3->AddEntry(grRatio_HS_hs,"CCK: GG-hs / Data", "P");
    leg3->AddEntry(grRatio_HS_n,"CCK: GG-n / Data", "P");
    leg3->AddEntry(grRatio_MS_fl,"MS: IPsat flu. / Data", "P");
    leg3->AddEntry(grRatio_MS_nf,"MS: IPsat no flu. / Data", "P");
    leg3->Draw();
    cDataModel->Modified();
    cDataModel->Update();

    // Prepare dashed line at y = 1
    TLine *line = new TLine(0.04,1.,1.0,1.0);
    line->SetLineColor(kOrange+2);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->Draw("SAME");

    cDataModel->Print(Form("PhotoCrossSec/img_PlotWithRatios/Ratios_FeedDown%i_%ibins.pdf", iFeedDown, nPtBins));
    cDataModel->Print(Form("PhotoCrossSec/img_PlotWithRatios/Ratios_FeedDown%i_%ibins.png", iFeedDown, nPtBins));

    // *****************************************************************************
    // Draw both
    TCanvas *cBoth = new TCanvas("cBoth","Cross section dependence on p_{t}^{2}",1050,1000);

    // Draw legend with models
    TLegend *leg1mod = SetLegend(0.17,0.04,0.42,0.30);
    leg1mod->SetTextSize(0.035);
    leg1mod->SetMargin(0.30);
    leg1mod->AddEntry(gr_SL,"STARlight", "L");
    leg1mod->AddEntry(gr_HM_fluct,"MS: IPsat flu.", "L");
    leg1mod->AddEntry(gr_HM_noflu,"MS: IPsat no flu.", "L");
    leg1mod->AddEntry(gr_GZ_area,"GSZ: el. + diss.", "F");
    leg1mod->AddEntry(gr_HS_hs,"CCK: GG-hs", "L");
    leg1mod->AddEntry(gr_HS_n, "CCK: GG-n", "L");

    TPad *pMain = new TPad("pMain","pMain",0.,0.25,1.,1.);
    SetPadMargins(pMain,0.13,0.03,0.03,0.0);
    pMain->SetLogy();
    pMain->Draw();
    pMain->cd();
    // Draw everything needed
    fCSont->GetXaxis()->SetLabelSize(0);
    fCSont->GetXaxis()->SetTitle("");
    fCSont->Draw("AXIS");
    gr_GZ_area->Draw("F SAME");
    grData_syst->Draw("5 SAME");
    gr_GZ_min->Draw("L SAME");
    gr_GZ_max->Draw("L SAME"); 
    gr_SL->Draw("L SAME");
    gr_HS_hs->Draw("CX SAME");
    gr_HS_n->Draw("CX SAME");
    gr_HM_fluct->Draw("CX SAME");
    gr_HM_noflu->Draw("CX SAME");    
    grData_stat->Draw("P SAME");
    //leg0->Draw();
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    leg1mod->Draw();
    leg2->Draw();

    cBoth->cd();
    TPad *pRatio = new TPad("pRatio","pRatio",0.,0.,1.,0.25);
    SetPadMargins(pRatio,0.13,0.0,0.03,0.33);
    pRatio->Draw();
    //pRatio->SetLogy();
    pRatio->cd();
    fCSratio->GetYaxis()->SetTickLength(0.025);
    fCSratio->GetXaxis()->SetTickLength(0.025);
    fCSratio->GetYaxis()->SetNdivisions(205);
    fCSratio->GetXaxis()->SetTitleOffset(1.);
    //fCSratio->GetYaxis()->SetTitleOffset(0.33);
    fCSratio->GetYaxis()->SetTitleOffset(0.4);
    fCSratio->GetXaxis()->SetTitleSize(0.15);
    fCSratio->GetYaxis()->SetTitleSize(0.15);
    fCSratio->GetXaxis()->SetLabelSize(0.15);
    fCSratio->GetYaxis()->SetLabelSize(0.15);
    fCSratio->Draw("AXIS");
    gr_syst_data->Draw("5 SAME");
    grRatio_SL->Draw("SAME P");
    grRatio_HS_hs->Draw("SAME P");
    grRatio_HS_n->Draw("SAME P");
    grRatio_MS_fl->Draw("SAME P");
    grRatio_MS_nf->Draw("SAME P");
    gr_stat_data->Draw("SAME P");
    line->Draw("SAME");
    TLegend *leg4 = SetLegend(0.72,0.40,0.94,0.90);
    leg4->SetMargin(0.12);
    leg4->AddEntry(grRatio_SL,"STARlight / Data", "P");
    leg4->AddEntry(grRatio_HS_hs,"GG-hs / Data", "P");
    leg4->AddEntry(grRatio_HS_n, "GG-n / Data", "P");
    leg4->AddEntry(grRatio_MS_fl,"IPsat flu. / Data", "P");
    leg4->AddEntry(grRatio_MS_nf,"IPsat no flu. / Data", "P");
    leg4->SetTextSize(0.105);
    leg4->Draw();

    cBoth->Print(Form("PhotoCrossSec/img_PlotWithRatios/RatiosPlot_FeedDown%i_%ibins.pdf", iFeedDown, nPtBins));
    cBoth->Print(Form("PhotoCrossSec/img_PlotWithRatios/RatiosPlot_FeedDown%i_%ibins.png", iFeedDown, nPtBins));

    return;
}