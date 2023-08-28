// CrossSec_PlotWithRatios.C
// David Grund, Sep 04, 2022

// my headers
#include "CrossSec_Utilities.h"

TGraphErrors *gr_ratios[9] = { NULL };
TGraph *gr_binned[9] = { NULL };
float textSize1 = 0.044; // for plot only
float textSize2 = 0.100; // for ratios only
float textSize3 = 0.050; // for plot in plot with ratios
float textSize4 = 0.140; // for ratios in plot with ratios
// horizontal range of the plots
float t_low = 0.04; // GeV^2
float t_upp = 1.2; // Gev^2

string SelectModels(int sel) {
    string mods = ""; // which models to plot
    if(sel == 0) mods = "012345678"; // all models
    else if(sel == 1) mods = "1357"; // models with fluctuations (CCK-hs,MS-hs,GSZ-el+diss,MSS-CGC+fl)
    else if(sel == 2) mods = "02468"; // models without flu. (CCK-n,MS-p,GSZ-el,MSS-CGC)
    else if(sel == 3) mods = "3456"; // MS and GSZ
    else if(sel == 4) mods = "345678"; // MS, GSZ and MSS
    else if(sel == 5) mods = "5678"; // GSZ and MSS
    else if(sel == 6) mods = "12"; // only CCK
    return mods;
}

void SetPadMarginsLTRB(TVirtualPad* pad, float l, float t, float r, float b)
{
    pad->SetLeftMargin(l);
    pad->SetTopMargin(t);
    pad->SetRightMargin(r);
    pad->SetBottomMargin(b);
}

void SetFrame(TH1* fr, float textSize)
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
    //fr->GetXaxis()->SetNdivisions(505);
}

void DrawLegend1(int modSel, float x1, float y1, float x2, float y2, float textSize, TString draw_leg)
{
    string mods = SelectModels(modSel);
    TLegend *l1 = SetLegend(x1,y1,x2,y2);
    l1->SetTextSize(textSize);
    l1->SetFillStyle(0);
    l1->SetMargin(0.30);
    for(int i = 0; i < 5; i++) if(mods.find(std::to_string(i)) != string::npos) l1->AddEntry(gr_models[i],str_models[i],draw_leg.Data());
    for(int i = 7; i < 9; i++) if(mods.find(std::to_string(i)) != string::npos) l1->AddEntry(gr_models[i],str_models[i],draw_leg.Data());
    for(int i = 5; i < 7; i++) if(mods.find(std::to_string(i)) != string::npos) l1->AddEntry(gr_models[i],str_models[i],draw_leg.Data());
    l1->Draw();

    return;
}

void DrawLegend2(float x1, float y1, float x2, float y2, float textSize)
{
    TLegend *l2 = SetLegend(x1,y1,x2,y2);
    l2->SetTextSize(textSize);
    l2->SetFillStyle(0);
    l2->SetMargin(0.12);
    l2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |#it{y}| < 0.8", "");
    l2->AddEntry(gr_data_uncr,"Uncorrelated stat. + syst.", "EPL");
    l2->AddEntry(gr_data_corr,"Correlated syst.", "F");
    l2->Draw();

    return;
}

void DrawLegend3(int modSel, float x1, float y1, float x2, float y2, float textSize)
{
    string mods = SelectModels(modSel);
    TLegend *l3 = SetLegend(x1,y1,x2,y2);
    l3->SetTextSize(textSize);
    // here we do not want the legend to be transparent
    l3->SetMargin(0.20);
    for(int i = 0; i < 9; i++) if(mods.find(std::to_string(i)) != string::npos) l3->AddEntry(gr_ratios[i],str_models[i],"P");
    l3->Draw();

    return;
}

void PlotWithRatios(int modSel, bool binned = false)
{
    string mods = SelectModels(modSel);
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
    if(!binned) {
        bin_load = "gr_";
        bin_save = "";
        draw_opt = "L SAME";
        draw_leg = "LP";
    } else {
        bin_load = "grBinned_";
        bin_save = "binned_";
        draw_opt = "P SAME";
        draw_leg = "P";
    }
    for(int i = 0; i < 9; i++) gr_models[i] = (TGraph*)lg->FindObject(bin_load + str_models[i]);
    gr_GSZ_err[0] = (TGraph*)lg->FindObject(bin_load + "err_GSZ-el+diss");
    gr_GSZ_err[1] = (TGraph*)lg->FindObject(bin_load + "err_GSZ-el");

    // create boxes with correlated syst. uncertainties
    SetupSysErrorBox(gr_data_corr,kGray+3);

    // set properties of the data graph with uncorrelated uncertainties
    gr_data_uncr->SetMarkerStyle(kFullCircle);
    gr_data_uncr->SetMarkerSize(1.0);
    gr_data_uncr->SetLineColor(kBlack);
    gr_data_uncr->SetLineWidth(2);
    gr_data_uncr->SetMarkerColor(kBlack);
    gr_data_uncr->Print();

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
    // MSS-CGC-fl
    SetLineMarkerProperties(gr_models[7],colors[7],1,kOpenDiamond,1.8); //(!)
    // MSS-CGC
    SetLineMarkerProperties(gr_models[8],colors[8],2,kFullDiamond,1.8); //(!)

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
    if(modSel >= 3) frPlot = gPad->DrawFrame(t_low, 0.0004, t_upp, 0.04);
    else            frPlot = gPad->DrawFrame(t_low, 0.0002, t_upp, 0.04);
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
    for(int i = 3; i < 9; i++) if(i != 5 && i != 6) gr_models[i]->Draw(draw_opt.Data());
    gr_data_uncr->Draw("PZ SAME");

    // ALICE PbPb label
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(textSize1*0.8);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.94,"ALICE, Pb#minusPb UPC   #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // draw the legend with the models
    DrawLegend1(0,0.17,0.15,0.40,0.44,textSize1*0.7,draw_leg);
    //DrawLegend1(modSel,0.70,0.65,0.90,0.90,textSize1*0.7,draw_leg);

    // draw the legend with the data + uncertainties
    DrawLegend2(0.54,0.80,0.92,0.92,textSize1*0.7);

    c1->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/" + bin_save + "plot.pdf");

    // ********************************************************************************
    // calculate and plot the ratios
    for(int i = 0; i < 9; i++) {
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

    for(int iBin = 0; iBin < nPtBins; iBin++) 
    {
        Double_t x, y, err_y_uncr, err_y_corr, err_x_low, err_x_upp;
        gr_data_uncr->GetPoint(iBin, x, y);
        err_y_uncr = gr_data_uncr->GetErrorYhigh(iBin);
        err_y_corr = gr_data_corr->GetErrorYhigh(iBin);
        err_x_upp = gr_data_uncr->GetErrorXhigh(iBin);
        err_x_low = gr_data_uncr->GetErrorXlow(iBin);
        
        for(int i = 0; i < 9; i++) gr_ratios[i]->SetPoint(iBin, x, gr_binned[i]->GetPointY(iBin) / gr_data_uncr->GetPointY(iBin));

        Double_t ratio_correction = 1/y;
        err_y_uncr *= ratio_correction;
        err_y_corr *= ratio_correction;

        gr_err_uncr->SetPoint(iBin, x, 1.);
        gr_err_uncr->SetPointError(iBin,0.,0.,err_y_uncr,err_y_uncr);
        gr_err_corr->SetPoint(iBin, x, 1.);
        gr_err_corr->SetPointError(iBin,err_x_low,err_x_upp,err_y_corr,err_y_corr);
    }
    // error bars for GSZ in the ratio panel:
    /*
    LoadGraphs_GSZ();
    for(int i = 5; i < 7; i++) {
        cout << str_models[i] << "\n"; 
        for(int iBin = 0; iBin < nPtBins; iBin++) {
            cout << Form("bin %i: rat: %.2f, upp: %.2f, low: %.2f, diff low-upp: %.2f\n",
                iBin+1,
                gr_ratios[i]->GetPointY(iBin),
                gr_ratios[i]->GetPointY(iBin)*GSZ_err_scale_upp[i-5],
                gr_ratios[i]->GetPointY(iBin)*GSZ_err_scale_low[i-5],
                gr_ratios[i]->GetPointY(iBin)*GSZ_err_scale_upp[i-5]-gr_ratios[i]->GetPointY(iBin)*GSZ_err_scale_low[i-5]);
        }
    }
    */

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
    // MSS-CGC-fl
    SetLineMarkerProperties(gr_ratios[7],colors[7],1,kOpenDiamond,1.8); //(!)
    // MSS-CGC
    SetLineMarkerProperties(gr_ratios[8],colors[8],1,kFullDiamond,1.8); //(!)

    // ********************************************************************************
    // draw the canvas with "ratios"

    TCanvas *c2 = new TCanvas ("c2","Ratios models / data",1050,400);
    SetPadMarginsLTRB(gPad,0.10,0.03,0.03,0.23);

    // draw frame
    Float_t rat_upp = 2.4;
    if(modSel >= 3) rat_upp = 2.2;
    if(modSel == 5) rat_upp = 1.65;
    TH1F* frRatio = gPad->DrawFrame(t_low,0.0,t_upp,rat_upp);
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
    gr_ratios[0]->Draw("SAME P");
    for(int i = 2; i < 9; i=i+2) gr_ratios[i]->Draw("SAME P");
    for(int i = 1; i < 9; i=i+2) gr_ratios[i]->Draw("SAME P");
    gr_err_uncr->Draw("SAME PZ");

    // dashed line at y = 1
    TLine *line = new TLine(t_low,1.,t_upp,1.);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // draw a legend
    DrawLegend3(0,0.78,0.30,0.96,0.90,textSize2*0.65);
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
    // y-axis    
    frPlot->GetYaxis()->SetTitleOffset(1.25);
    frPlot->GetYaxis()->SetTitleSize(textSize3);
    frPlot->GetYaxis()->SetLabelSize(textSize3);
    // x-axis
    frPlot->GetXaxis()->SetLabelSize(0);
    frPlot->GetXaxis()->SetTitleSize(0);
    frPlot->GetXaxis()->SetTitle("");

    frPlot->Draw("AXIS");
    if(mods.find(std::to_string(6)) != string::npos) {
        gr_GSZ_err[1]->Draw("F SAME");
        gr_models[6]->Draw(draw_opt.Data());
    }
    if(mods.find(std::to_string(5)) != string::npos) {
        gr_GSZ_err[0]->Draw("F SAME");
        gr_models[5]->Draw(draw_opt.Data());
    }
    gr_data_corr->Draw("5 SAME");
    if(mods.find(std::to_string(7)) != string::npos) gr_models[7]->Draw(draw_opt.Data());
    for(int i = 0; i < 9; i++) if(mods.find(std::to_string(i)) != string::npos && i != 5 && i != 6 && i != 7) gr_models[i]->Draw(draw_opt.Data());
    gr_data_uncr->Draw("PZ SAME");

    // draw latex label
    ltx->SetTextSize(textSize3*0.88);
    Bool_t preliminary = kFALSE;
    if(preliminary) ltx->DrawLatex(0.55,0.92,"ALICE Preliminary, Pb#minusPb UPC   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    else            ltx->DrawLatex(0.55,0.92,"ALICE, Pb#minusPb UPC   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // draw legends
    if(modSel == 0) DrawLegend1(modSel,0.155,0.03,0.405,0.37,textSize3*0.8,draw_leg);
    else if(modSel == 4) DrawLegend1(modSel,0.155,0.03,0.405,0.32,textSize3*0.8,draw_leg);
    else DrawLegend1(modSel,0.155,0.03,0.405,0.28,textSize3*0.8,draw_leg);
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
    if(mods.find(std::to_string(0)) != string::npos) gr_ratios[0]->Draw("SAME P");
    for(int i = 2; i < 9; i=i+2) if(mods.find(std::to_string(i)) != string::npos) gr_ratios[i]->Draw("SAME P");
    for(int i = 1; i < 9; i=i+2) if(mods.find(std::to_string(i)) != string::npos) gr_ratios[i]->Draw("SAME P");
    gr_err_uncr->Draw("SAME PZ");
    line->Draw("SAME");
    /*
    if(iModels == 0)      DrawLegend3(iModels,0.72,0.385,0.93,0.945,textSize4*0.60);
    else if(iModels == 3) DrawLegend3(iModels,0.71,0.385,0.96,0.945,textSize4*0.80);
    else                  DrawLegend3(iModels,0.72,0.385,0.93,0.945,textSize4*0.65); 
    */

    TString append = "";
    if(modSel == 0) append += "_all";
    if(modSel == 1) append += "_flu";
    if(modSel == 2) append += "_noflu";
    if(modSel == 3) append += "_MS-GSZ";
    if(modSel == 4) append += "_MS-GSZ-MSS";
    if(modSel == 5) append += "_MSS-GSZ";
    if(modSel == 6) append += "_CCK";

    append += Form("_upto%.1f",t_upp);

    c3->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/" + bin_save + "plotWithRatios" + append + ".pdf");
    c3->Print("Results/" + str_subfolder + "CrossSec/PlotWithRatios/" + bin_save + "plotWithRatios" + append + ".C");
    if(!binned && modSel == 3) {
        if(preliminary) {
            c3->Print("Results/" + str_subfolder + "_PreliminaryFigures/crossSection.pdf");
            c3->Print("Results/" + str_subfolder + "_PreliminaryFigures/crossSection.eps");
        } else c3->Print("Results/" + str_subfolder + "_PaperFigures/crossSection.pdf");
    }

    return;
}

void CrossSec_PlotWithRatios(int iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/PlotWithRatios/");

    // original binning of the models ("continuous")
    t_upp = 1.0;
    for(int i = 0; i < 7; i++) PlotWithRatios(i);
    t_upp = 1.2;
    for(int i = 0; i < 7; i++) PlotWithRatios(i);
    //PlotWithRatios(4);
    // binning of the measurement (5 bins)
    //for(int i = 0; i < 6; i++) PlotWithRatios(i,true);

    return;
}