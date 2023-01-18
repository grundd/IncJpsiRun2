// CrossSec_PlotWithRatios.C
// David Grund, Sep 04, 2022

// my headers
#include "CrossSec_Utilities.h"

TGraphErrors *gr_ratios[7] = { NULL };
TGraph *gr_binned[7] = { NULL };
Double_t textSize1 = 0.044; // for plot only
Double_t textSize2 = 0.100; // for ratios only
Double_t textSize3 = 0.050; // for plot in plot with ratios
Double_t textSize4 = 0.140; // for ratios in plot with ratios

void PlotWithRatios(Int_t iBinn, Int_t iModels = 0);
// iBinn == 0 => original binning of the models
//       == 1 => everything in 5 bins as the data
// iModels == 0 => all the models plotted in c3
//         == 1 => only the models with fluctuations (CCK-hs, MS-hs, GSZ-el+diss)
//         == 2 => only those without (SL, CCK-n, MS-p, GSZ-el)
//         == 3 => only MS and GSZ (both versions of them)

void SetPadMarginsLTRB(TVirtualPad* pad, Double_t l, Double_t t, Double_t r, Double_t b)
{
    pad->SetLeftMargin(l);
    pad->SetTopMargin(t);
    pad->SetRightMargin(r);
    pad->SetBottomMargin(b);
}

void SetFrame(TH1* fr, Double_t textSize)
{
    // title and label sizes
    fr->GetXaxis()->SetTitleSize(textSize);
    fr->GetYaxis()->SetTitleSize(textSize);
    fr->GetXaxis()->SetLabelSize(textSize);
    fr->GetYaxis()->SetLabelSize(textSize);
    // font types
    fr->GetXaxis()->SetTitleFont(42);
    fr->GetYaxis()->SetTitleFont(42);
    fr->GetXaxis()->SetLabelFont(42);
    fr->GetYaxis()->SetLabelFont(42);
    // divisions on the x-axis
    fr->GetXaxis()->SetNdivisions(505);
}

void DrawLegend1(Int_t iModels, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textSize, TString draw_leg)
{
    TLegend *l1 = SetLegend(x1,y1,x2,y2);
    l1->SetTextSize(textSize);
    l1->SetFillStyle(0);
    l1->SetMargin(0.30);
    //for(Int_t i = 0; i < 7; i++) l1->AddEntry(gr_models[i],str_models[i],"L");
    if(iModels != 1 && iModels != 3) l1->AddEntry(gr_models[0],str_models[0],draw_leg.Data());
    if(iModels != 2 && iModels != 3) l1->AddEntry(gr_models[1],str_models[1],draw_leg.Data());
    if(iModels != 1 && iModels != 3) l1->AddEntry(gr_models[2],str_models[2],draw_leg.Data());
    if(iModels != 2) l1->AddEntry(gr_models[3],str_models[3],draw_leg.Data());
    if(iModels != 1) l1->AddEntry(gr_models[4],str_models[4],draw_leg.Data());
    if(iModels != 2) l1->AddEntry(gr_models[5],str_models[5],draw_leg.Data());
    if(iModels != 1) l1->AddEntry(gr_models[6],str_models[6],draw_leg.Data());
    l1->Draw();

    return;
}

void DrawLegend2(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textSize)
{
    TLegend *l2 = SetLegend(x1,y1,x2,y2);
    l2->SetTextSize(textSize);
    l2->SetFillStyle(0);
    l2->SetMargin(0.12);
    l2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |y| < 0.8", "");
    l2->AddEntry(gr_data_uncr,"Uncorrelated stat. + syst.", "EPL");
    l2->AddEntry(gr_data_corr,"Correlated syst.", "F");
    l2->Draw();

    return;
}

void DrawLegend3(Int_t iModels, Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textSize)
{
    TLegend *l3 = SetLegend(x1,y1,x2,y2);
    l3->SetTextSize(textSize);
    // here we do not want the legend to be transparent
    l3->SetMargin(0.20);
    //for(Int_t i = 0; i < 7; i++) l3->AddEntry(gr_ratios[i],str_models[i],"P");
    if(iModels != 1 && iModels != 3) l3->AddEntry(gr_ratios[0],str_models[0],"P");
    if(iModels != 2 && iModels != 3) l3->AddEntry(gr_ratios[1],str_models[1],"P");
    if(iModels != 1 && iModels != 3) l3->AddEntry(gr_ratios[2],str_models[2],"P");
    if(iModels != 2) l3->AddEntry(gr_ratios[3],str_models[3],"P");
    if(iModels != 1) l3->AddEntry(gr_ratios[4],str_models[4],"P");
    if(iModels != 2) l3->AddEntry(gr_ratios[5],str_models[5],"P");
    if(iModels != 1) l3->AddEntry(gr_ratios[6],str_models[6],"P");
    l3->Draw();

    return;
}

void CrossSec_PlotWithRatios(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/PlotWithRatios/");

    // original binning of the models ("continuous")
    PlotWithRatios(0,0);
    PlotWithRatios(0,1);
    PlotWithRatios(0,2);
    PlotWithRatios(0,3);
    // binning from the measurement (5 bins)
    PlotWithRatios(1,0);
    PlotWithRatios(1,1);
    PlotWithRatios(1,2);
    PlotWithRatios(1,3);

    return;
}

void PlotWithRatios(Int_t iBinn, Int_t iModels)
{
    // open the file with histograms and graphs
    TFile *f = TFile::Open("Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/histograms_and_graphs.root","read");
    if(f) Printf("Input file %s loaded.", f->GetName()); 
    TList *lh = (TList*) f->Get("histograms");
    if(lh) Printf("List %s loaded.", lh->GetName());
    TList *lg = (TList*) f->Get("graphs");
    if(lg) Printf("List %s loaded.", lg->GetName());
    // load the graphs
    gr_data_uncr = (TGraphAsymmErrors*)lg->FindObject("gr_data_uncr");
    gr_data_corr = (TGraphAsymmErrors*)lg->FindObject("gr_data_corr");
    TString bin_load = "";
    TString bin_save = "";
    TString draw_opt = "";
    TString draw_leg = "";
    if(iBinn == 0)
    {
        bin_load = "gr_";
        bin_save = "";
        draw_opt = "L SAME";
        draw_leg = "LP";
    }    
    else if(iBinn == 1)
    {
        bin_load = "grBinned_";
        bin_save = "binned_";
        draw_opt = "P SAME";
        draw_leg = "P";
    } 
    else return;
    for(Int_t i = 0; i < 9; i++) gr_models[i] = (TGraph*)lg->FindObject(bin_load + str_models[i]);
    gr_GSZ_err[0] = (TGraph*)lg->FindObject(bin_load + "err_GSZ-el+diss");
    gr_GSZ_err[1] = (TGraph*)lg->FindObject(bin_load + "err_GSZ-el");

    // create boxes with correlated syst. uncertainties
    SetupSysErrorBox(gr_data_corr,kGray+3);

    // set properties of the data graph with uncorrelated uncertainties
    gr_data_uncr->SetMarkerStyle(kFullCircle);
    gr_data_uncr->SetMarkerSize(0.7);
    gr_data_uncr->SetLineColor(kBlack);
    gr_data_uncr->SetLineWidth(2);
    gr_data_uncr->SetMarkerColor(kBlack);

    // set up properties of the graphs
    // STARlight
    SetLineMarkerProperties(gr_models[0],colors[0],1,kFullDiamond,1.8);
    // CCK-hs
    SetLineMarkerProperties(gr_models[1],colors[1],2,kOpenCross,1.2);
    // CCK-n
    SetLineMarkerProperties(gr_models[2],colors[2],6,kFullCross,1.2);
    // MS-hs
    SetLineMarkerProperties(gr_models[3],colors[3],7,kOpenSquare,1.2);
    // MS-p
    SetLineMarkerProperties(gr_models[4],colors[4],8,kFullSquare,1.2);
    // GSZ-el+diss
    SetLineMarkerProperties(gr_models[5],colors[5],9,kOpenCircle,1.2);
    // GSZ-el
    SetLineMarkerProperties(gr_models[6],colors[6],4,kFullCircle,1.2);
    // GSZ error bands:
    SetupSysErrorBox(gr_GSZ_err[0],kGreen);
    SetupSysErrorBox(gr_GSZ_err[1],kOrange);

    // global style settings:
    gStyle->SetTextFont(42); // so that latex text is not bold
    /*
    gStyle->SetEndErrorSize(4); 
    gStyle->SetErrorX(0.02);
    gStyle->SetLineScalePS(2.0);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleOffset(1.25,"XYZ");
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetOptStat("0");
    gStyle->SetOptFit(0);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"y");
    */

    // ********************************************************************************
    // Draw the canvas with "plot"

    TCanvas *c1 = new TCanvas ("c1","Cross section dependence on |t|",880,1000);
    SetPadMarginsLTRB(gPad,0.14,0.02,0.03,0.12);
    c1->SetLogy();

    // draw frame
    TH1F* frPlot;
    if(iModels == 3) frPlot = gPad->DrawFrame(0.04, 0.0004, 1.0, 0.04);
    else             frPlot = gPad->DrawFrame(0.04, 0.0002, 1.0, 0.04);
    SetFrame(frPlot,textSize1);
    frPlot->SetTitle("Cross section dependence on |#it{t}|;|#it{t}| (GeV^{2});d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");
    // y-axis
    frPlot->GetYaxis()->SetTickLength(0.025); 
    frPlot->GetYaxis()->SetTitleOffset(1.50);
    // x-axis
    frPlot->GetXaxis()->SetTickLength(0.025); 
    frPlot->GetXaxis()->SetTitleOffset(1.10);
    frPlot->GetXaxis()->SetDecimals(1);

    c1->cd();
    frPlot->Draw("AXIS");

    gr_GSZ_err[1]->Draw("F SAME");
    gr_models[6]->Draw(draw_opt.Data());
    gr_GSZ_err[0]->Draw("F SAME");
    gr_models[5]->Draw(draw_opt.Data());
    gr_data_corr->Draw("5 SAME");
    for(Int_t i = 0; i < 5; i++) gr_models[i]->Draw(draw_opt.Data());
    gr_data_uncr->Draw("PZ SAME");

    // ALICE PbPb label
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(textSize1*0.8);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.94,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // draw the legend with the models
    DrawLegend1(0,0.17,0.15,0.40,0.40,textSize1*0.7,draw_leg);

    // draw the legend with the data + uncertainties
    DrawLegend2(0.54,0.80,0.92,0.92,textSize1*0.7);

    c1->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/" + bin_save + "plot.pdf");

    // ********************************************************************************
    // calculate and plot the ratios
    for(Int_t i = 0; i < 7; i++) {
        gr_ratios[i] = new TGraphErrors(nPtBins);
        gr_binned[i] = (TGraph*)lg->FindObject("grBinned_" + str_models[i]);
    } 
    
    TGraphAsymmErrors *gr_err_uncr = new TGraphAsymmErrors(nPtBins);
    TGraphAsymmErrors *gr_err_corr = new TGraphAsymmErrors(nPtBins);

    SetupSysErrorBox(gr_err_corr,kGray+3);

    gr_err_uncr->SetMarkerSize(0.);
    gr_err_uncr->SetLineColor(kBlack);
    gr_err_uncr->SetLineWidth(2);
    gr_err_uncr->SetMarkerColor(kBlack);

    for(Int_t iBin = 0; iBin < nPtBins; iBin++) 
    {
        Double_t x, y, err_y_uncr, err_y_corr, err_x_low, err_x_upp;
        gr_data_uncr->GetPoint(iBin, x, y);
        err_y_uncr = gr_data_uncr->GetErrorYhigh(iBin);
        err_y_corr = gr_data_corr->GetErrorYhigh(iBin);
        err_x_upp = gr_data_uncr->GetErrorXhigh(iBin);
        err_x_low = gr_data_uncr->GetErrorXlow(iBin);
        
        for(Int_t i = 0; i < 7; i++) gr_ratios[i]->SetPoint(iBin, x, gr_binned[i]->GetPointY(iBin) / gr_data_uncr->GetPointY(iBin));

        Double_t ratio_correction = 1/y;
        err_y_uncr *= ratio_correction;
        err_y_corr *= ratio_correction;

        gr_err_uncr->SetPoint(iBin, x, 1.);
        gr_err_uncr->SetPointError(iBin,0.,0.,err_y_uncr,err_y_uncr);
        gr_err_corr->SetPoint(iBin, x, 1.);
        gr_err_corr->SetPointError(iBin,err_x_low,err_x_upp,err_y_corr,err_y_corr);
    }

    // STARlight
    SetLineMarkerProperties(gr_ratios[0],colors[0],1,kFullDiamond,1.8);
    // CCK-hs
    SetLineMarkerProperties(gr_ratios[1],colors[1],1,kOpenCross,1.2);
    // CCK-n
    SetLineMarkerProperties(gr_ratios[2],colors[2],1,kFullCross,1.2);
    // MS-hs
    SetLineMarkerProperties(gr_ratios[3],colors[3],1,kOpenSquare,1.2);
    // MS-p
    SetLineMarkerProperties(gr_ratios[4],colors[4],1,kFullSquare,1.2);
    // GSZ-el+diss
    SetLineMarkerProperties(gr_ratios[5],colors[5],1,kOpenCircle,1.2);
    // GSZ-el
    SetLineMarkerProperties(gr_ratios[6],colors[6],1,kFullCircle,1.2);

    // ********************************************************************************
    // draw the canvas with "ratios"

    TCanvas *c2 = new TCanvas ("c2","Ratios models / data",1050,400);
    SetPadMarginsLTRB(gPad,0.10,0.03,0.03,0.23);

    // draw frame
    Float_t upperY = 2.4;
    if(iModels == 3) upperY = 2.2;
    TH1F* frRatio = gPad->DrawFrame(0.04,0.0,1.0,upperY);
    SetFrame(frRatio,textSize2);
    frRatio->SetTitle("Ratios model/data;|#it{t}| (GeV^{2});Model / Data");
    // y-axis
    frRatio->GetYaxis()->SetTitleOffset(0.5);
    frRatio->GetYaxis()->SetNdivisions(205);
    // x-axis
    frRatio->GetXaxis()->SetDecimals(1);

    c2->cd();
    frRatio->Draw("AXIS");
    gr_err_corr->Draw("5 SAME");
    for(Int_t i = 0; i < 7; i++) gr_ratios[i]->Draw("SAME P");
    gr_err_uncr->Draw("SAME PZ");

    // dashed line at y = 1
    TLine *line = new TLine(0.04,1.,1.0,1.0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // draw a legend
    DrawLegend3(0,0.80,0.30,0.96,0.90,textSize2*0.65);
    c2->Modified();
    c2->Update();
    
    c2->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/ratios.pdf");

    // ********************************************************************************
    // draw both

    TCanvas *c3 = new TCanvas("c3","Cross section dependence on |t|",880,1000);
    // pad with the plot
    TPad *pUpp = new TPad("pUpp","pUpp",0.,0.25,1.,1.);
    SetPadMarginsLTRB(pUpp,0.13,0.03,0.03,0.0);
    pUpp->SetLogy();
    pUpp->Draw();
    pUpp->cd();
    // x-axis
    frPlot->GetYaxis()->SetTitleOffset(1.25);
    frPlot->GetYaxis()->SetTitleSize(textSize3);
    frPlot->GetYaxis()->SetLabelSize(textSize3);
    // y-axis
    frPlot->GetXaxis()->SetLabelSize(0);
    frPlot->GetXaxis()->SetTitleSize(0);
    frPlot->GetXaxis()->SetTitle("");

    frPlot->Draw("AXIS");
    if(iModels != 1) gr_GSZ_err[1]->Draw("F SAME");
    if(iModels != 1) gr_models[6]->Draw(draw_opt.Data());
    if(iModels != 2) gr_GSZ_err[0]->Draw("F SAME");
    if(iModels != 2) gr_models[5]->Draw(draw_opt.Data());
    gr_data_corr->Draw("5 SAME");
    //for(Int_t i = 0; i < 5; i++) gr_models[i]->Draw("L SAME");
    if(iModels != 1 && iModels != 3) gr_models[0]->Draw(draw_opt.Data());
    if(iModels != 2 && iModels != 3) gr_models[1]->Draw(draw_opt.Data());
    if(iModels != 1 && iModels != 3) gr_models[2]->Draw(draw_opt.Data());
    if(iModels != 2) gr_models[3]->Draw(draw_opt.Data());
    if(iModels != 1) gr_models[4]->Draw(draw_opt.Data());
    gr_data_uncr->Draw("PZ SAME");

    // draw latex label
    ltx->SetTextSize(textSize3*0.88);
    ltx->DrawLatex(0.55,0.92,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // draw legends
    if(iModels == 0) DrawLegend1(iModels,0.17,0.04,0.42,0.34,textSize3*0.8,draw_leg);
    else             DrawLegend1(iModels,0.17,0.04,0.42,0.28,textSize3*0.8,draw_leg);
    DrawLegend2(0.54,0.74,0.98,0.88,textSize3*0.8);

    // pad with the ratios
    c3->cd();
    TPad *pLow = new TPad("pLow","pLow",0.,0.,1.,0.25);
    SetPadMarginsLTRB(pLow,0.13,0.0,0.03,0.33);
    pLow->Draw();
    pLow->cd();
    // y-axis
    frRatio->GetYaxis()->SetTitleOffset(0.45);
    frRatio->GetYaxis()->SetTickLength(0.025);
    frRatio->GetYaxis()->SetNdivisions(205);
    frRatio->GetYaxis()->SetTitleSize(textSize4);
    frRatio->GetYaxis()->SetLabelSize(textSize4);
    // x-axis
    frRatio->GetXaxis()->SetTickLength(0.07);
    frRatio->GetXaxis()->SetTitleOffset(1.0);
    frRatio->GetXaxis()->SetTitleSize(textSize4);
    frRatio->GetXaxis()->SetLabelSize(textSize4);
    
    frRatio->Draw("AXIS");
    gr_err_corr->Draw("5 SAME");
    //for(Int_t i = 0; i < 7; i++) gr_ratios[i]->Draw("SAME P");
    if(iModels != 1 && iModels != 3) gr_ratios[0]->Draw("P SAME");
    if(iModels != 2 && iModels != 3) gr_ratios[1]->Draw("P SAME");
    if(iModels != 1 && iModels != 3) gr_ratios[2]->Draw("P SAME");
    if(iModels != 2) gr_ratios[3]->Draw("P SAME");
    if(iModels != 1) gr_ratios[4]->Draw("P SAME");
    if(iModels != 2) gr_ratios[5]->Draw("P SAME");
    if(iModels != 1) gr_ratios[6]->Draw("P SAME");
    gr_err_uncr->Draw("SAME PZ");
    line->Draw("SAME");
    /*
    if(iModels == 0)      DrawLegend3(iModels,0.72,0.385,0.93,0.945,textSize4*0.60);
    else if(iModels == 3) DrawLegend3(iModels,0.71,0.385,0.96,0.945,textSize4*0.80);
    else                  DrawLegend3(iModels,0.72,0.385,0.93,0.945,textSize4*0.65); 
    */

    TString append = "";
    if(iModels == 0) append += "_allModels";
    if(iModels == 1) append += "_fluctOnly";
    if(iModels == 2) append += "_nofluOnly";
    if(iModels == 3) append += "_onlyMSGSZ";

    c3->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/" + bin_save + "plotWithRatios" + append + ".pdf");
    if(iBinn == 0 && iModels == 3) c3->Print("Results/" + str_subfolder + "_PaperFigures/crossSection.pdf");

    return;
}