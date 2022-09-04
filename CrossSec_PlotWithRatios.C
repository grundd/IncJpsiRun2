// CrossSec_PlotWithRatios.C
// David Grund, Sep 04, 2022

// my headers
#include "CrossSec_Utilities.h"

Int_t lineWidth = 3;

TGraphErrors *gr_ratios[7] = { NULL };
TGraph *gr_binned[7] = { NULL };

void PlotWithRatios(Int_t iBinn);
// iBinn == 0 => original binning of the models
//       == 1 => everything in 5 bins as the data

void SetLineColorStyleWidth(TGraph *gr, Color_t col, Int_t stl)
{
    gr->SetLineColor(col);
    gr->SetLineStyle(stl);
    gr->SetLineWidth(lineWidth);
    return;
}

void SetMarkerColorStyle(TGraph *gr, Color_t col, Int_t stl)
{
    gr->SetMarkerColor(col);
    gr->SetMarkerStyle(stl);
    return;
}

TLegend *SetLegend(Double_t x_leftdown, Double_t y_leftdown, Double_t x_rightup, Double_t y_rightup)
{
    TLegend *leg = new TLegend(x_leftdown, y_leftdown, x_rightup, y_rightup);
    //leg->SetFillColor(0);
    //leg->SetFillStyle(0);
    //leg->SetBorderSize(0);
    return leg;
}

void SetPadMargins(TVirtualPad* pad, Double_t left, Double_t top, Double_t right, Double_t bottom)
{
    pad->SetTopMargin(top);
    pad->SetBottomMargin(bottom);
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
}

void SetFrame(TH1* frame)
{
    frame->GetXaxis()->SetTitleOffset(1.);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetNdivisions(306);
    frame->GetXaxis()->SetNoExponent();
    //frame->GetYaxis()->SetNoExponent();
}

void SetStyle(TGraph* g, Color_t color, Style_t style, Width_t width = 3)
{
    g->SetLineColor(color);
    g->SetLineStyle(style);
    g->SetLineWidth(width);
}

void SetupSysErrorBox(TGraph* g, Color_t color)
{
    g->SetMarkerSize(0);
    g->SetFillStyle(1001);
    g->SetFillColorAlpha(color,0.35);
    SetStyle(g,color,1,0);
}

void DrawLegend1(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize)
{
    TLegend *l1 = SetLegend(x1,y1,x2,y2);
    l1->SetFillStyle(0);
    l1->SetTextSize(textsize);
    l1->SetMargin(0.30);
    for(Int_t i = 0; i < 7; i++) l1->AddEntry(gr_models[i],str_models[i],"L");
    l1->Draw();
    return;
}

void DrawLegend2(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize)
{
    TLegend *l2 = SetLegend(x1,y1,x2,y2);
    l2->SetTextSize(textsize);
    l2->SetMargin(0.11);
    l2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |y| < 0.8", "");
    l2->AddEntry(gr_data_uncr,"Exp. uncorrelated stat. + syst.", "EPL");
    l2->AddEntry(gr_data_corr,"Exp. correlated syst.", "F");
    l2->Draw();
    return;
}

void DrawLegend3(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textsize)
{
    TLegend *l3 = SetLegend(x1,y1,x2,y2);
    l3->SetTextSize(textsize);
    l3->SetMargin(0.13);
    for(Int_t i = 0; i < 7; i++) l3->AddEntry(gr_ratios[i],str_models[i] + " / Data","P");
    l3->Draw();
    return;
}

void CrossSec_PlotWithRatios(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/PlotWithRatios/");

    // original normalization of the models, original binning ("continuous")
    PlotWithRatios(0);

    return;
}

void PlotWithRatios(Int_t iBinn)
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
    if(iBinn == 0)
    {
        for(Int_t i = 0; i < 9; i++) gr_models[i] = (TGraph*)lg->FindObject("gr_" + str_models[i]);
        gr_GSZ_err[0] = (TGraph*)lg->FindObject("gr_err_GSZ-el+diss");
        gr_GSZ_err[1] = (TGraph*)lg->FindObject("gr_err_GSZ-el");
    }

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

    TCanvas *cPlot = new TCanvas ("cPlot","Cross section dependence on p_{t}^{2}",1000,900);
    SetPadMargins(gPad,0.14,0.02,0.02,0.12);
    gStyle->SetOptStat("0");
    gStyle->SetOptFit(0);

    // draw frame
    TH1F* fCSont = gPad->DrawFrame(0.04, 0.0002, 1.0, 0.04);
    SetFrame(fCSont);
    fCSont->GetYaxis()->SetTickLength(0.025); 
    fCSont->GetXaxis()->SetTickLength(0.025); 
    fCSont->GetYaxis()->SetTitleOffset(1.25);
    fCSont->SetTitle("Cross section dependence on |#it{t}|;|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    // format: title of the canvas, title of the X axis, title of the Y axis

    cPlot->cd();
    cPlot->SetLogy();
    cPlot->Modified();
    fCSont->Draw("AXIS");

    // create boxes with correlated syst. uncertainties
    SetupSysErrorBox(gr_data_corr,kGray+3);

    // set properties of the data graph with uncorrelated uncertainties
    gStyle->SetEndErrorSize(4);         
    gr_data_uncr->SetMarkerStyle(kFullCircle);
    gr_data_uncr->SetMarkerSize(0.7);
    gr_data_uncr->SetLineColor(kBlack);
    gr_data_uncr->SetLineWidth(2);
    gr_data_uncr->SetMarkerColor(kBlack);

    // set up properties of the graphs
    // STARlight
    SetLineColorStyleWidth(gr_models[0],kBlue,1);
    // CCK-hs
    SetLineColorStyleWidth(gr_models[1],kRed+1,9);
    // CCK-n
    SetLineColorStyleWidth(gr_models[2],kRed+1,5);
    // MS-hs
    SetLineColorStyleWidth(gr_models[3],kGray+3,9);
    // MS-p
    SetLineColorStyleWidth(gr_models[4],kGray+3,5);
    // GSZ-el+diss
    SetLineColorStyleWidth(gr_models[5],kGreen+2,9);
    // GSZ-el
    SetLineColorStyleWidth(gr_models[6],kCyan+2,5);
    // GSZ error bands:
    SetupSysErrorBox(gr_GSZ_err[0],kGreen);
    SetupSysErrorBox(gr_GSZ_err[1],kCyan);

    // **********************************************************************
    // draw everything
    gr_GSZ_err[0]->Draw("F SAME");
    gr_GSZ_err[1]->Draw("F SAME");
    gr_data_corr->Draw("5 SAME");
    gr_models[5]->Draw("L SAME");
    gr_models[6]->Draw("L SAME");
    for(Int_t i = 0; i < 5; i++) gr_models[i]->Draw("L SAME");
    gr_data_uncr->Draw("P SAME");
    cPlot->Modified();
    cPlot->Update(); 

    // ALICE PbPb label
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.035);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // draw the legend with the models
    DrawLegend1(0.17,0.15,0.40,0.40,0.032);
    cPlot->Modified();
    cPlot->Update();

    // draw the legend with the data + uncertainties
    DrawLegend2(0.56,0.76,0.90,0.90,0.032);
    cPlot->Modified();
    cPlot->Update();

    cPlot->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/plot.pdf");

    // *****************************************************************************
    // calculate and plot ratios
    for(Int_t i = 0; i < 7; i++)
    {
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

    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"y");
    TCanvas *cRat = new TCanvas ("cRat","",1050,300);
    SetPadMargins(gPad,0.05,0.03,0.03,0.12);
    TH1F* fCSratio = gPad->DrawFrame(0.04,0.0,1.0,2.4);
    SetFrame(fCSratio);
    fCSratio->SetTitle("Ratios model/data;|#it{t}| (GeV^{2} #it{c}^{-2});Model / Data");
    fCSratio->GetYaxis()->SetTitleOffset(0.5);

    // STARlight
    SetMarkerColorStyle(gr_ratios[0],kBlue,kFullSquare);
    // CCK-hs
    SetMarkerColorStyle(gr_ratios[1],kRed+1,kOpenCircle);
    // CCK-n
    SetMarkerColorStyle(gr_ratios[2],kRed+1,kFullCircle);
    // MS-hs
    SetMarkerColorStyle(gr_ratios[3],kGray+3,kOpenTriangleDown);
    // MS-p
    SetMarkerColorStyle(gr_ratios[4],kGray+3,kFullTriangleDown);
    // GSZ-el+diss
    SetMarkerColorStyle(gr_ratios[5],kGreen+2,kOpenCross);
    // GSZ-el
    SetMarkerColorStyle(gr_ratios[6],kCyan+2,kFullCross);

    // draw everything
    cRat->cd();
    cRat->Modified();
    fCSratio->Draw("AXIS");
    // errors
    gr_err_corr->Draw("5 SAME");
    for(Int_t i = 0; i < 7; i++) gr_ratios[i]->Draw("SAME P");
    cRat->Modified();
    cRat->Update();
    // Draw data with stat. errors
    gr_err_uncr->Draw("SAME P");

    // prepare dashed line at y = 1
    TLine *line = new TLine(0.04,1.,1.0,1.0);
    line->SetLineColor(kOrange+2);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // draw a legend
    DrawLegend3(0.69,0.28,0.95,0.88,0.08);
    cRat->Modified();
    cRat->Update();
    
    cRat->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/ratios.pdf");

    // *****************************************************************************
    // Draw both
    TCanvas *cBoth = new TCanvas("cBoth","Cross section dependence on p_{t}^{2}",900,1000);
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
    /*
    gr_GSZ_err[0]->Draw("F SAME");
    gr_GSZ_err[1]->Draw("F SAME");
    gr_data_corr->Draw("5 SAME");
    gr_models[5]->Draw("L SAME");
    gr_models[6]->Draw("L SAME");
    for(Int_t i = 0; i < 5; i++) gr_models[i]->Draw("L SAME");
    gr_data_uncr->Draw("P SAME");
    */
    gr_GSZ_err[0]->Draw("F SAME");
    //gr_GSZ_err[1]->Draw("F SAME");
    gr_data_corr->Draw("5 SAME");
    gr_models[5]->Draw("L SAME");
    //gr_models[6]->Draw("L SAME");
    //for(Int_t i = 0; i < 5; i++) gr_models[i]->Draw("L SAME");
    gr_models[1]->Draw("L SAME");
    gr_models[3]->Draw("L SAME");
    gr_data_uncr->Draw("P SAME");
    latex->SetTextSize(0.044);
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // draw legends
    DrawLegend1(0.17,0.04,0.42,0.34,0.038);
    DrawLegend2(0.58,0.76,0.84,0.90,0.038);
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
    gr_err_corr->Draw("5 SAME");
    for(Int_t i = 0; i < 7; i++) gr_ratios[i]->Draw("SAME P");
    gr_err_uncr->Draw("SAME P");
    line->Draw("SAME");
    DrawLegend3(0.72,0.38,0.945,0.94,0.09);

    cBoth->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/plotRatios.pdf");

    return;
}