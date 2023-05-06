// standalonePlot_incoh.cxx
// David Grund, May 6, 2023

// root headers
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TList.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"

TGraphAsymmErrors* gr_data_uncr = NULL;
TGraphAsymmErrors* gr_data_corr = NULL;
TGraph* gr_models[4] = { NULL };
TGraph *gr_binned[4] = { NULL };
TGraph* gr_GSZ_err[2] = { NULL };
TGraphErrors *gr_ratios[4] = { NULL };
TString str_models[4] = {
    "MS-hs",
    "MS-p",
    "GSZ-el+diss",
    "GSZ-el"};
Color_t colors[7] = {
    kViolet-1, // MS-hs
    kBlue, // MS-p
    kGreen+2, // GSZ-el+diss
    kOrange+2 // GSZ-el
};
Float_t sizeUpp = 0.050; // text size for the plot
Float_t sizeLow = 0.140; // text size for the bottom panel (ratios)

void SetupSysErrorBox(TGraph* g, Color_t color)
{
    g->SetMarkerSize(0);
    g->SetFillStyle(1001);
    g->SetFillColorAlpha(color,0.35);
    g->SetLineWidth(0);
    return;
}

void SetLineMarkerProperties(TGraph *g, Color_t color, Int_t lStyle, Int_t mStyle = kCircle, Size_t mSize = 1.)
{
    g->SetLineColor(color);
    g->SetLineStyle(lStyle);
    g->SetLineWidth(3);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(mStyle);
    g->SetMarkerSize(mSize);
    return;
}

void SetMarginsLeftTopRightBottom(TVirtualPad* p, Float_t l, Float_t t, Float_t r, Float_t b)
{
    p->SetLeftMargin(l);
    p->SetTopMargin(t);
    p->SetRightMargin(r);
    p->SetBottomMargin(b);
    return;
}

void SetFrame(TH1* f, Float_t textSize)
{
    // title and label sizes
    f->GetXaxis()->SetTitleSize(textSize);
    f->GetYaxis()->SetTitleSize(textSize);
    f->GetXaxis()->SetLabelSize(textSize);
    f->GetYaxis()->SetLabelSize(textSize);
    // font types
    f->GetXaxis()->SetTitleFont(42);
    f->GetYaxis()->SetTitleFont(42);
    f->GetXaxis()->SetLabelFont(42);
    f->GetYaxis()->SetLabelFont(42);
    // divisions on the x-axis
    f->GetXaxis()->SetNdivisions(505);
    return;
}

TLegend *SetLegend(Double_t x_leftdown, Double_t y_leftdown, Double_t x_rightup, Double_t y_rightup)
{
    TLegend *l = new TLegend(x_leftdown, y_leftdown, x_rightup, y_rightup);
    l->SetFillColor(0);
    l->SetBorderSize(0);
    return l;
}

void PlotWithRatios(string config = "")
// "a" -> no data, only models without fluctuations
// "b" -> no data, all four models
// "c" -> only data
// "d" -> everything (paper plot)
{
    bool plotData = false;
    bool plotFlu = false;
    bool plotNoflu = false;
    if(config == "c" || config == "d") plotData = true;
    if(config == "b" || config == "d") plotFlu = true;
    if(config == "a" || config == "d") plotNoflu = true;
    // open the file with histograms and graphs
    TFile *f = TFile::Open("histograms_and_graphs.root","read");
    if(f) Printf("File %s loaded.", f->GetName());
    TList *l = (TList*) f->Get("graphs");
    if(l) Printf("List %s loaded.", l->GetName());
    // load the graphs
    gr_data_uncr = (TGraphAsymmErrors*)l->FindObject("gr_data_uncr");
    gr_data_corr = (TGraphAsymmErrors*)l->FindObject("gr_data_corr");
    for(Int_t i = 0; i < 4; i++) gr_models[i] = (TGraph*)l->FindObject("gr_" + str_models[i]);
    gr_GSZ_err[0] = (TGraph*)l->FindObject("gr_err_GSZ-el+diss");
    gr_GSZ_err[1] = (TGraph*)l->FindObject("gr_err_GSZ-el");

    // create boxes with correlated syst. uncertainties
    SetupSysErrorBox(gr_data_corr,kGray+3);

    // set properties of the data graph with uncorrelated uncertainties
    gr_data_uncr->SetMarkerStyle(kFullCircle);
    gr_data_uncr->SetMarkerSize(1.0);
    gr_data_uncr->SetLineColor(kBlack);
    gr_data_uncr->SetLineWidth(2);
    gr_data_uncr->SetMarkerColor(kBlack);

    // set up properties of the graphs
    // MS-hs
    SetLineMarkerProperties(gr_models[0],colors[0],7,kOpenSquare,1.2);
    // MS-p
    SetLineMarkerProperties(gr_models[1],colors[1],8,kFullSquare,1.2);
    // GSZ-el+diss
    SetLineMarkerProperties(gr_models[2],colors[2],9,kOpenCircle,1.2);
    // GSZ-el
    SetLineMarkerProperties(gr_models[3],colors[3],4,kFullCircle,1.2);
    // GSZ error bands:
    SetupSysErrorBox(gr_GSZ_err[0],kGreen);
    SetupSysErrorBox(gr_GSZ_err[1],kOrange);

    // global style settings:
    gStyle->SetTextFont(42); // latex text not bold

    // ********************************************************************************
    // calculate the ratios
    for(Int_t i = 0; i < 4; i++) {
        gr_ratios[i] = new TGraphErrors(5);
        gr_binned[i] = (TGraph*)l->FindObject("grBinned_" + str_models[i]);
    } 
    
    TGraphAsymmErrors *gr_err_uncr = new TGraphAsymmErrors(5);
    TGraphAsymmErrors *gr_err_corr = new TGraphAsymmErrors(5);

    SetupSysErrorBox(gr_err_corr,kGray+3);

    gr_err_uncr->SetMarkerSize(0.);
    gr_err_uncr->SetLineColor(kBlack);
    gr_err_uncr->SetLineWidth(2);
    gr_err_uncr->SetMarkerColor(kBlack);

    for(Int_t iBin = 0; iBin < 5; iBin++) 
    {
        Double_t x, y, err_y_uncr, err_y_corr, err_x_low, err_x_upp;
        gr_data_uncr->GetPoint(iBin, x, y);
        err_y_uncr = gr_data_uncr->GetErrorYhigh(iBin);
        err_y_corr = gr_data_corr->GetErrorYhigh(iBin);
        err_x_upp = gr_data_uncr->GetErrorXhigh(iBin);
        err_x_low = gr_data_uncr->GetErrorXlow(iBin);
        
        for(Int_t i = 0; i < 4; i++) gr_ratios[i]->SetPoint(iBin, x, gr_binned[i]->GetPointY(iBin) / gr_data_uncr->GetPointY(iBin));

        Double_t ratio_correction = 1/y;
        err_y_uncr *= ratio_correction;
        err_y_corr *= ratio_correction;

        gr_err_uncr->SetPoint(iBin, x, 1.);
        gr_err_uncr->SetPointError(iBin,0.,0.,err_y_uncr,err_y_uncr);
        gr_err_corr->SetPoint(iBin, x, 1.);
        gr_err_corr->SetPointError(iBin,err_x_low,err_x_upp,err_y_corr,err_y_corr);
    }

    // MS-hs
    SetLineMarkerProperties(gr_ratios[0],colors[0],1,kOpenSquare,1.2);
    // MS-p
    SetLineMarkerProperties(gr_ratios[1],colors[1],1,kFullSquare,1.2);
    // GSZ-el+diss
    SetLineMarkerProperties(gr_ratios[2],colors[2],1,kOpenCircle,1.2);
    // GSZ-el
    SetLineMarkerProperties(gr_ratios[3],colors[3],1,kFullCircle,1.2);

    // ********************************************************************************
    // draw the cross section plot with the ratio panel

    TCanvas *c = new TCanvas("c","Cross section dependence on |t|",880,1000);
    // pad with the plot
    TPad *pUpp = new TPad("pUpp","pUpp",0.,0.25,1.,1.);
    SetMarginsLeftTopRightBottom(pUpp,0.13,0.03,0.03,0.0);
    pUpp->SetLogy();
    pUpp->Draw();
    pUpp->cd();

    TH1F* frPlot = gPad->DrawFrame(0.04, 0.0004, 1.0, 0.04);
    SetFrame(frPlot,sizeUpp);
    frPlot->SetTitle("Cross section dependence on |#it{t}|;|#it{t}| (GeV^{2});d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");    
    // x-axis
    frPlot->GetXaxis()->SetTickLength(0.025); 
    frPlot->GetXaxis()->SetTitleOffset(1.25);
    frPlot->GetXaxis()->SetDecimals(1);
    // y-axis
    frPlot->GetYaxis()->SetTickLength(0.025); 
    frPlot->GetYaxis()->SetTitleOffset(1.25);
    frPlot->GetXaxis()->SetLabelSize(0);
    frPlot->GetXaxis()->SetTitleSize(0);
    frPlot->GetXaxis()->SetTitle("");

    frPlot->Draw("AXIS");
    if(plotNoflu) gr_GSZ_err[1]->Draw("F SAME");
    if(plotNoflu) gr_models[3]->Draw("L SAME");
    if(plotFlu) gr_GSZ_err[0]->Draw("F SAME");
    if(plotFlu) gr_models[2]->Draw("L SAME");
    if(plotData) gr_data_corr->Draw("5 SAME");
    if(plotFlu) gr_models[0]->Draw("L SAME");
    if(plotNoflu) gr_models[1]->Draw("L SAME");
    if(plotData) gr_data_uncr->Draw("PZ SAME");

    // draw latex label
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(sizeUpp*0.88);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.92,"ALICE, Pb#minusPb UPC   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // draw legends
    // legend for the measurement
    TLegend *l1 = SetLegend(0.54,0.74,0.98,0.88);
    l1->SetTextSize(sizeUpp*0.8);
    l1->SetFillStyle(0);
    l1->SetMargin(0.12);
    l1->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |#it{y}| < 0.8", "");
    l1->AddEntry(gr_data_uncr,"Uncorrelated stat. + syst.", "EPL");
    l1->AddEntry(gr_data_corr,"Correlated syst.", "F");
    l1->Draw();
    // legend for the models
    TLegend *l2 = SetLegend(0.17,0.04,0.42,0.28);
    l2->SetTextSize(sizeUpp*0.8);
    l2->SetFillStyle(0);
    l2->SetMargin(0.30);
    if(plotFlu) l2->AddEntry(gr_models[0],str_models[0],"LP");
    if(plotNoflu) l2->AddEntry(gr_models[1],str_models[1],"LP");
    if(plotFlu) l2->AddEntry(gr_models[2],str_models[2],"LP");
    if(plotNoflu) l2->AddEntry(gr_models[3],str_models[3],"LP");
    l2->Draw();

    // pad with the ratios
    c->cd();
    TPad *pLow = new TPad("pLow","pLow",0.,0.,1.,0.25);
    SetMarginsLeftTopRightBottom(pLow,0.13,0.0,0.03,0.33);
    pLow->Draw();
    pLow->cd();

    TH1F* frRatio = gPad->DrawFrame(0.04,0.0,1.0,2.2);
    SetFrame(frRatio,sizeLow);
    frRatio->SetTitle("Ratios model/data;|#it{t}| (GeV^{2});Model / Data");
    // y-axis
    frRatio->GetYaxis()->SetTitleOffset(0.5);
    frRatio->GetYaxis()->SetNdivisions(205);
    // x-axis
    frRatio->GetXaxis()->SetDecimals(1);
    // y-axis
    frRatio->GetYaxis()->SetTitleOffset(0.45);
    frRatio->GetYaxis()->SetTickLength(0.025);
    frRatio->GetYaxis()->SetNdivisions(205);
    // x-axis
    frRatio->GetXaxis()->SetTickLength(0.07);
    frRatio->GetXaxis()->SetTitleOffset(1.0);
    
    frRatio->Draw("AXIS");
    if(plotData) {
        gr_err_corr->Draw("5 SAME");
        if(plotFlu) gr_ratios[0]->Draw("P SAME");
        if(plotNoflu) gr_ratios[1]->Draw("P SAME");
        if(plotFlu) gr_ratios[2]->Draw("P SAME");
        if(plotNoflu) gr_ratios[3]->Draw("P SAME");
        gr_err_uncr->Draw("SAME PZ");
    }
    
    // dashed line at y = 1
    TLine *line = new TLine(0.04,1.,1.0,1.0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->Draw("SAME");

    TString append = ""; 
    if(config == "a") append += "_modelsNoflu";
    if(config == "b") append += "_modelsAll";
    if(config == "c") append += "_dataOnly";
    if(config == "d") append += "_all";
    c->Print("crossSection" + append + ".pdf");

    delete c;
    return;
}

void standalonePlot_incoh()
{
    PlotWithRatios("a");
    PlotWithRatios("b");
    PlotWithRatios("c");
    PlotWithRatios("d");
    return;
}