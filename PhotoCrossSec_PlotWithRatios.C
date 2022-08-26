// PhotoCrossSec_PlotWithRatios.C
// David Grund, Apr 24, 2022

// root headers
#include "TSystem.h"
// my headers
#include "PhotoCrossSec_Utilities.h"

Double_t scale(0.);
TGraphAsymmErrors *gr_data_stat = NULL;
TGraphAsymmErrors *gr_data_syst = NULL;
TGraph *gr_models[9] = { NULL };
// order: SL, CCK_hs, CCK_n, MS_fl, MS_nf, GSZ_tot_max, GSZ_tot_min, GSZ_el_max, GSZ_el_min
TGraph *gr_GSZ_tot_area = NULL;
TGraph *gr_GSZ_el_area = NULL;
TGraphErrors *grRat_SL = NULL;
TGraphErrors *grRat_CCK_hs = NULL;
TGraphErrors *grRat_CCK_n = NULL;
TGraphErrors *grRat_MS_fl = NULL;
TGraphErrors *grRat_MS_nf = NULL;
TGraphErrors *grRat_GSZ = NULL;

void PlotWithRatios(Int_t iNorm, Int_t iBinn);
// iNorm == 0 => original normalizations of the models
//       == 1 => everything scaled to the measured cross section
// iBinn == 0 => original binning of the models
//       == 1 => everything in 5 bins as data
void NormalizeHisto(TH1D *h){ h->Scale(scale / h->Integral("width")); return; }
void SetLineColorStyleWidth(TGraph *gr, Color_t col, Int_t stl);
void DrawLegend1(Int_t iNorm, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize);
void DrawLegend2(Int_t iNorm, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize);
void DrawLegend3(Int_t iNorm, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize);

void PhotoCrossSec_PlotWithRatios(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PlotWithRatios/");

    // original normalization of the models, original binning ("continuous")
    PlotWithRatios(0,0);

    return;
}

void PlotWithRatios(Int_t iNorm, Int_t iBinn)
{
    // load histograms with pheno predictions
    // these are created in PhotoCrossSec_PrepareHistograms.C
    TString path;
    TString binning;
    if(iBinn == 0)
    {
        path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/continuous/all.root";
        binning = "continuous/";
    }      
    else if(iBinn == 1)
    {
        path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/all.root";
        binning = "binned/";
    } 
    else return;  
    // open the file with histograms
    TFile *f = TFile::Open(path.Data(),"read");
    if(f) Printf("Input file %s loaded.", f->GetName()); 
    TList *l = (TList*) f->Get("HistList");
    if(l) Printf("List %s loaded.", l->GetName()); 
    // load histograms
    TH1D *h_models[9] = { NULL };
    for(Int_t i = 0; i < 9; i++) h_models[i] = (TH1D*)l->FindObject(binning + str_models[i]);
    // load the measured cross section
    ifstream ifs;
    path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/integrals/data.txt";
    ifs.open(path.Data());
    if(!ifs.fail())
    {
        ifs >> scale;
        ifs.close();
        Printf("Loaded measured cross section: %.3f mub.", scale);
    }
    else 
    {
        Printf("File %s not found. Terminating...", path.Data());
        return;
    }
    scale = scale / 1e3; // from mub to mb
    // change normalizations
    if(iNorm == 0) ; // nothing
    else if(iNorm == 1)
    {
        // we don't care about the errors of the histograms because we don't use them
        for(Int_t i = 0; i < 9; i++) NormalizeHisto(h_models[i]);
    }
    else return;

    // fill data graphs
    ReadInput_data();
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
    gr_data_stat = new TGraphAsymmErrors(nPtBins,abs_t_val,sig_val,abs_t_err_low,abs_t_err_upp,sig_err_stat_low,sig_err_stat_upp);
    gr_data_syst = new TGraphAsymmErrors(nPtBins,abs_t_val,sig_val,abs_t_err_low,abs_t_err_upp,sig_err_syst_low,sig_err_syst_upp);

    // using Roman's settings:
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

    TCanvas *cPlot = new TCanvas ("cPlot","Cross section dependence on p_{t}^{2}",1050,800);
    SetPadMargins(gPad,0.13,0.03,0.03,0.11);
    gStyle->SetOptStat("0");
    gStyle->SetOptFit(0);

    // draw frame
    TH1F* fCSont = NULL;
    fCSont = gPad->DrawFrame(0.04, 0.0002, 1.0, 0.04);
    SetFrame(fCSont);
    fCSont->GetYaxis()->SetTickLength(0.025); 
    fCSont->GetXaxis()->SetTickLength(0.025); 
    fCSont->GetYaxis()->SetTitleOffset(1.2);
    fCSont->SetTitle("Cross section dependence on |#it{t}|;|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    // format: title of the canvas, title of the X axis, title of the Y axis

    cPlot->cd();
    cPlot->SetLogy();
    cPlot->Modified();
    fCSont->Draw("AXIS");

    // create boxes with systematic uncertainties
    SetupSysErrorBox(gr_data_syst,kGray+3);

    // Set data properties 
    gStyle->SetEndErrorSize(4);         
    gr_data_stat->SetMarkerStyle(kFullCircle);
    gr_data_stat->SetMarkerSize(0.7);
    gr_data_stat->SetLineColor(kBlack);
    gr_data_stat->SetLineWidth(2);
    gr_data_stat->SetMarkerColor(kBlack);

    // STARlight
    gr_models[0] = new TGraph(h_models[0]);
    SetLineColorStyleWidth(gr_models[0],kBlue,1);

    // CCK model
    // with hot spots:
    gr_models[1] = new TGraph(h_models[1]);
    SetLineColorStyleWidth(gr_models[1],kRed+1,9);
    // without hot spots:
    gr_models[2] = new TGraph(h_models[2]);
    SetLineColorStyleWidth(gr_models[2],kRed+1,8);

    // MS model
    // with subnucleonic fluctuations
    gr_models[3] = new TGraph(h_models[3]);
    SetLineColorStyleWidth(gr_models[3],kGray+3,9);
    // without subnucleonic fluctuations
    gr_models[4] = new TGraph(h_models[4]);
    SetLineColorStyleWidth(gr_models[4],kGray+3,8);

    // GSZ model
    // total cross section
    gr_models[5] = new TGraph(h_models[5]);
    gr_models[6] = new TGraph(h_models[6]);
    SetLineColorStyleWidth(gr_models[6],kGreen,10);
    SetLineColorStyleWidth(gr_models[5],kGreen,10);
    // elastic only
    gr_models[7] = new TGraph(h_models[7]);
    gr_models[8] = new TGraph(h_models[8]);
    SetLineColorStyleWidth(gr_models[8],kCyan,10);
    SetLineColorStyleWidth(gr_models[7],kCyan,10);
    // areas
    ReadInput_GSZ();
    // scale the values (GSZ uses nb instead of mb)
    for (Int_t i = 0; i < n_GSZ; i++){
        sig_GSZ_tot_min[i] = sig_GSZ_tot_min[i] / 1e6;
        sig_GSZ_tot_max[i] = sig_GSZ_tot_max[i] / 1e6;
        sig_GSZ_el_min[i] = sig_GSZ_el_min[i] / 1e6;
        sig_GSZ_el_max[i] = sig_GSZ_el_max[i] / 1e6;
    }
    if(iNorm == 0)
    {
        gr_GSZ_tot_area = new TGraph(2*n_GSZ);
        for (Int_t i = 0; i < n_GSZ; i++){
            gr_GSZ_tot_area->SetPoint(i, abs_t_GSZ[i], sig_GSZ_tot_max[i]);
            gr_GSZ_tot_area->SetPoint(n_GSZ+i, abs_t_GSZ[n_GSZ-i-1], sig_GSZ_tot_min[n_GSZ-i-1]);
        }
        gr_GSZ_el_area = new TGraph(2*n_GSZ);
        for (Int_t i = 0; i < n_GSZ; i++){
            gr_GSZ_el_area->SetPoint(i, abs_t_GSZ[i], sig_GSZ_el_max[i]);
            gr_GSZ_el_area->SetPoint(n_GSZ+i, abs_t_GSZ[n_GSZ-i-1], sig_GSZ_el_min[n_GSZ-i-1]);
        }
        SetupSysErrorBox(gr_GSZ_tot_area,kGreen);
        SetupSysErrorBox(gr_GSZ_el_area,kCyan);
    }
    else if(iNorm == 1) // if everything scaled to the measured cross section => no area displayed
    {
        gr_GSZ_tot_area = new TGraph(h_models[5]);
        SetLineColorStyleWidth(gr_GSZ_tot_area,kGreen,10);
        gr_GSZ_el_area = new TGraph(h_models[7]);
        SetLineColorStyleWidth(gr_GSZ_el_area,kGreen,10);
    }
    else return;

    // **********************************************************************
    // draw everything
    if(iNorm == 0)
    {
        gr_GSZ_tot_area->Draw("F SAME");
        gr_GSZ_el_area->Draw("F SAME");
        gr_data_syst->Draw("5 SAME");
        gr_models[6]->Draw("L SAME");
        gr_models[5]->Draw("L SAME");
        gr_models[8]->Draw("L SAME");
        gr_models[7]->Draw("L SAME"); 
    }
    else if(iNorm == 1)
    {
        gr_GSZ_tot_area->Draw("L SAME"); // here it is not an area but a line
        gr_GSZ_el_area->Draw("L SAME");  // here it is not an area but a line
        gr_data_syst->Draw("5 SAME");
    }
    else return;
    gr_models[0]->Draw("L SAME");
    gr_models[1]->Draw("CX SAME");
    gr_models[2]->Draw("CX SAME");
    gr_models[3]->Draw("CX SAME");
    gr_models[4]->Draw("CX SAME");
    gr_data_stat->Draw("P SAME");
    cPlot->Modified();
    cPlot->Update(); 

    // ALICE PbPb label
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.035);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // draw the legend with the models
    DrawLegend1(iNorm,0.17,0.15,0.40,0.40,0.032);
    cPlot->Modified();
    cPlot->Update();

    // draw the legend with the data + uncertainties
    DrawLegend2(iNorm,0.62,0.76,0.96,0.90,0.032);
    cPlot->Modified();
    cPlot->Update();

    path = "Results/" + str_subfolder + Form("PhotoCrossSec/PlotWithRatios/plot");
    if(iNorm == 1) path += "_scaled";
    cPlot->Print((path + ".pdf").Data());

    // *****************************************************************************
    // calculate and plot ratios
    grRat_SL = new TGraphErrors(nPtBins);
    grRat_GSZ = new TGraphErrors(nPtBins);
    grRat_CCK_hs = new TGraphErrors(nPtBins);
    grRat_CCK_n = new TGraphErrors(nPtBins);
    grRat_MS_fl = new TGraphErrors(nPtBins);
    grRat_MS_nf = new TGraphErrors(nPtBins);
    TGraphAsymmErrors *gr_stat_data = new TGraphAsymmErrors(nPtBins);
    TGraphAsymmErrors *gr_syst_data = new TGraphAsymmErrors(nPtBins);

    SetupSysErrorBox(gr_syst_data,kGray+3);

    gr_stat_data->SetMarkerSize(0.);
    gr_stat_data->SetLineColor(kBlack);
    gr_stat_data->SetLineWidth(2);
    gr_stat_data->SetMarkerColor(kBlack);

    for(Int_t i = 0; i < nPtBins; ++i) {
        Double_t x,y,err_y_stat,err_y_syst,ey_gama,err_x_low,err_x_upp;
        gr_data_stat->GetPoint(i, x, y);
        err_y_stat = gr_data_stat->GetErrorYhigh(i);
        err_y_syst = gr_data_syst->GetErrorYhigh(i);
        err_x_upp = gr_data_syst->GetErrorXhigh(i);
        err_x_low = gr_data_syst->GetErrorXlow(i);
        
        grRat_SL->SetPoint(i, x, gr_models[0]->Eval(x)/y);
        if(iNorm == 1) grRat_GSZ->SetPoint(i, x, gr_GSZ_tot_area->Eval(x)/y);
        grRat_CCK_hs->SetPoint(i, x, gr_models[1]->Eval(x)/y);
        grRat_CCK_n->SetPoint(i, x, gr_models[2]->Eval(x)/y);
        grRat_MS_fl->SetPoint(i, x, gr_models[3]->Eval(x)/y);
        grRat_MS_nf->SetPoint(i, x, gr_models[4]->Eval(x)/y);

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
    TCanvas *cRatios = new TCanvas ("cRatios","",1050,300);
    SetPadMargins(gPad,0.05,0.03,0.03,0.12);
    TH1F* fCSratio = NULL;
    if(iNorm == 0)      fCSratio = gPad->DrawFrame(0.04,0.0,1.0,2.3);
    else if(iNorm == 1) fCSratio = gPad->DrawFrame(0.04,0.0,1.0,3.3);
    else return;
    SetFrame(fCSratio);
    fCSratio->SetTitle("Ratios model/data;|#it{t}| (GeV^{2} #it{c}^{-2});Model / Data");
    fCSratio->GetYaxis()->SetTitleOffset(0.5);

    cRatios->cd();
    cRatios->Modified();
    fCSratio->Draw("AXIS");
    // Errors
    gr_syst_data->Draw("5 SAME");
    cRatios->Modified();
    cRatios->Update();
    // STARlight ratios
    grRat_SL->SetLineColor(kBlue);
    grRat_SL->SetMarkerColor(kBlue);
    grRat_SL->SetMarkerStyle(kFullSquare);
    grRat_SL->Draw("SAME P");
    cRatios->Modified();
    cRatios->Update();
    // GSZ model ratios
    grRat_GSZ->SetLineColor(kGreen);
    grRat_GSZ->SetMarkerColor(kGreen);
    grRat_GSZ->SetMarkerStyle(kFullCross);
    grRat_GSZ->Draw("SAME P");
    cRatios->Modified();
    cRatios->Update();    
    // Hot-spot model ratios
    grRat_CCK_hs->SetLineColor(kRed+1);
    grRat_CCK_hs->SetMarkerColor(kRed+1);
    grRat_CCK_hs->SetMarkerStyle(kOpenCircle);
    grRat_CCK_hs->Draw("SAME P");
    grRat_CCK_n->SetLineColor(kRed+1);
    grRat_CCK_n->SetMarkerColor(kRed+1);
    grRat_CCK_n->SetMarkerStyle(kFullCircle);
    grRat_CCK_n->Draw("SAME P");
    cRatios->Modified();
    cRatios->Update();
    // IPsat model ratios
    grRat_MS_fl->SetLineColor(kGray+3);
    grRat_MS_fl->SetMarkerColor(kGray+3);
    grRat_MS_fl->SetMarkerStyle(kOpenTriangleDown);
    grRat_MS_fl->Draw("SAME P");
    grRat_MS_nf->SetLineColor(kGray+3);
    grRat_MS_nf->SetMarkerColor(kGray+3);
    grRat_MS_nf->SetMarkerStyle(kFullTriangleDown);
    grRat_MS_nf->Draw("SAME P");
    // Draw data with stat. errors
    gr_stat_data->Draw("SAME P");

    // prepare dashed line at y = 1
    TLine *line = new TLine(0.04,1.,1.0,1.0);
    line->SetLineColor(kOrange+2);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // draw a legend
    DrawLegend3(iNorm,0.69,0.28,0.95,0.88,0.08);
    cRatios->Modified();
    cRatios->Update();

    path = "Results/" + str_subfolder + Form("PhotoCrossSec/PlotWithRatios/ratios");
    if(iNorm == 1) path += "_scaled";
    cRatios->Print((path + ".pdf").Data());

    // *****************************************************************************
    // Draw both
    TCanvas *cBoth = new TCanvas("cBoth","Cross section dependence on p_{t}^{2}",1050,1000);
    TPad *pMain = new TPad("pMain","pMain",0.,0.25,1.,1.);
    SetPadMargins(pMain,0.13,0.03,0.03,0.0);
    pMain->SetLogy();
    pMain->Draw();
    pMain->cd();
    // Draw everything needed
    fCSont->GetYaxis()->SetTitleOffset(1.02);
    fCSont->GetXaxis()->SetTitleSize(0.06);
    fCSont->GetYaxis()->SetTitleSize(0.06);
    fCSont->GetXaxis()->SetLabelSize(0.06);
    fCSont->GetYaxis()->SetLabelSize(0.06);
    fCSont->GetXaxis()->SetLabelSize(0);
    fCSont->GetXaxis()->SetTitle("");
    fCSont->Draw("AXIS");
    if(iNorm == 0)
    {
        gr_GSZ_tot_area->Draw("F SAME");
        gr_GSZ_el_area->Draw("F SAME");
        gr_data_syst->Draw("5 SAME");
        gr_models[6]->Draw("L SAME");
        gr_models[5]->Draw("L SAME");
        gr_models[8]->Draw("L SAME");
        gr_models[7]->Draw("L SAME"); 
    }
    else if(iNorm == 1)
    {
        gr_GSZ_tot_area->Draw("L SAME"); // here it is not an area but a line
        gr_GSZ_el_area->Draw("L SAME");  // here it is not an area but a line
        gr_data_syst->Draw("5 SAME");
    }
    else return;
    gr_models[0]->Draw("L SAME");
    gr_models[1]->Draw("CX SAME");
    gr_models[2]->Draw("CX SAME");
    gr_models[3]->Draw("CX SAME");
    gr_models[4]->Draw("CX SAME");
    gr_data_stat->Draw("P SAME");
    latex->SetTextSize(0.044);
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // draw legends
    DrawLegend1(iNorm,0.17,0.04,0.42,0.34,0.038);
    DrawLegend2(iNorm,0.62,0.76,0.96,0.90,0.038);
    // draw pad with ratios
    cBoth->cd();
    TPad *pRatio = new TPad("pRatio","pRatio",0.,0.,1.,0.25);
    SetPadMargins(pRatio,0.13,0.0,0.03,0.33);
    pRatio->Draw();
    pRatio->cd();
    fCSratio->GetYaxis()->SetTickLength(0.025);
    fCSratio->GetXaxis()->SetTickLength(0.025);
    fCSratio->GetYaxis()->SetNdivisions(205);
    fCSratio->GetXaxis()->SetTitleOffset(1.);
    fCSratio->GetYaxis()->SetTitleOffset(0.4);
    fCSratio->GetXaxis()->SetTitleSize(0.15);
    fCSratio->GetYaxis()->SetTitleSize(0.15);
    fCSratio->GetXaxis()->SetLabelSize(0.15);
    fCSratio->GetYaxis()->SetLabelSize(0.15);
    fCSratio->Draw("AXIS");
    gr_syst_data->Draw("5 SAME");
    grRat_SL->Draw("SAME P");
    grRat_CCK_hs->Draw("SAME P");
    grRat_CCK_n->Draw("SAME P");
    if(iNorm == 1) grRat_GSZ->Draw("SAME P");
    grRat_MS_fl->Draw("SAME P");
    grRat_MS_nf->Draw("SAME P");
    gr_stat_data->Draw("SAME P");
    line->Draw("SAME");
    DrawLegend3(iNorm,0.72,0.38,0.945,0.94,0.100);

    path = "Results/" + str_subfolder + "PhotoCrossSec/PlotWithRatios/plotRatios";
    if(iNorm == 1) path += "_scaled";
    cBoth->Print((path + ".pdf").Data());

    return;
}

void DrawLegend1(Int_t iNorm, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize)
{
    TLegend *l1 = SetLegend(x1,y1,x2,y2);
    l1->SetFillStyle(0);
    l1->SetTextSize(textsize);
    l1->SetMargin(0.30);
    l1->AddEntry(gr_models[0],"STARlight", "L");
    l1->AddEntry(gr_models[3],"MS: IPsat flu.", "L");
    l1->AddEntry(gr_models[4],"MS: IPsat no flu.", "L");
    if(iNorm == 0)
    {
        l1->AddEntry(gr_GSZ_tot_area,"GSZ: el. + diss.", "F");
        l1->AddEntry(gr_GSZ_el_area,"GSZ: el.", "F");
    }
    else if(iNorm == 1)
    {
        l1->AddEntry(gr_GSZ_tot_area,"GSZ: el. + diss.", "L");
        l1->AddEntry(gr_GSZ_el_area,"GSZ: el.", "L");
    }
    else return;
    l1->AddEntry(gr_models[1],"CCK: GG-hs", "L");
    l1->AddEntry(gr_models[2], "CCK: GG-n", "L");
    l1->Draw();
    return;
}

void DrawLegend2(Int_t iNorm, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize)
{
    TLegend *l2 = SetLegend(x1,y1,x2,y2);
    l2->SetTextSize(textsize);
    l2->SetMargin(0.11);
    l2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |y| < 0.8", "");
    l2->AddEntry(gr_data_stat,"Experimental stat.", "EPL");
    l2->AddEntry(gr_data_syst,"Experimental syst.", "F");
    l2->Draw();
    return;
}

void DrawLegend3(Int_t iNorm, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize)
{
    TLegend *l3 = SetLegend(x1,y1,x2,y2);
    l3->SetTextSize(textsize);
    l3->SetMargin(0.13);
    l3->AddEntry(grRat_SL,"STARlight / Data", "P");
    l3->AddEntry(grRat_MS_fl,"MS: IPsat flu. / Data", "P");
    l3->AddEntry(grRat_MS_nf,"MS: IPsat no flu. / Data", "P");
    if(iNorm == 1) l3->AddEntry(grRat_GSZ,"GSZ: el. + diss. / Data", "P");
    l3->AddEntry(grRat_CCK_hs,"CCK: GG-hs / Data", "P");
    l3->AddEntry(grRat_CCK_n,"CCK: GG-n / Data", "P");
    l3->Draw();
    return;
}