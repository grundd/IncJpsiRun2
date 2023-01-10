// CrossSec_Plot.C
// David Grund, Apr 24, 2022

// my headers
#include "CrossSec_Utilities.h"

void Plot();

void CrossSec_Plot(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/Plot/");

    Plot();

    return;
}

void Plot()
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
    for(Int_t i = 0; i < 9; i++) gr_models[i] = (TGraph*)lg->FindObject("gr_" + str_models[i]);
    gr_GSZ_err[0] = (TGraph*)lg->FindObject("gr_err_GSZ-el+diss");
    gr_GSZ_err[1] = (TGraph*)lg->FindObject("gr_err_GSZ-el");

    // with stat errors
    gStyle->SetEndErrorSize(4);         
    gr_data_uncr->SetMarkerStyle(kFullCircle);
    gr_data_uncr->SetMarkerSize(0.7);
    gr_data_uncr->SetLineColor(kBlack);
    gr_data_uncr->SetLineWidth(2);
    gr_data_uncr->SetMarkerColor(kBlack);
    // with syst errors 
    gr_data_corr->SetFillColor(17);
    gr_data_corr->SetMarkerSize(0);
    gr_data_corr->SetMarkerStyle(1);
    gr_data_corr->SetMarkerColor(kBlack);
    gr_data_corr->SetLineWidth(0); // to have no line around the box

    // set up properties of the graphs
    // STARlight
    SetLineMarkerProperties(gr_models[0],colors[0],1);
    // CCK-hs
    SetLineMarkerProperties(gr_models[1],colors[1],2);
    // CCK-n
    SetLineMarkerProperties(gr_models[2],colors[2],6);
    // MS-hs
    SetLineMarkerProperties(gr_models[3],colors[3],7);
    // MS-p
    SetLineMarkerProperties(gr_models[4],colors[4],8);
    // GSZ-el+diss
    SetLineMarkerProperties(gr_models[5],colors[5],9);
    // GSZ-el
    SetLineMarkerProperties(gr_models[6],colors[6],4);
    // GSZ error bands:
    SetupSysErrorBox(gr_GSZ_err[0],kGreen);
    SetupSysErrorBox(gr_GSZ_err[1],kOrange); 

    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // canvas
    TCanvas *c1 = new TCanvas("c1","c1",900,800);
    c1->SetLogy();  
    // margins
    c1->SetTopMargin(0.02);
    c1->SetBottomMargin(0.14);
    c1->SetRightMargin(0.02);
    c1->SetLeftMargin(0.14);
    // plot the graphs
    TH1 *h1 = (TH1*) gr_GSZ_err[1]->GetHistogram();
    h1->SetTitle(";|#it{t}| (GeV^{2}); d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");
    h1->SetMinimum(0.0005);
    h1->SetMaximum(0.05);
    // vertical axis
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetTitleOffset(1.52);
    h1->GetYaxis()->SetLabelSize(0.045);
    // horizontal axis
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetXaxis()->SetTitleOffset(1.2);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetRangeUser(0.04,1.0);
    // draw everything
    gr_GSZ_err[1]->Draw("AF");
    gr_models[6]->Draw("L SAME");
    gr_GSZ_err[0]->Draw("F SAME");
    gr_models[5]->Draw("L SAME");
    //gr_data_corr->Draw("5 SAME");
    for(Int_t i = 0; i < 5; i++) gr_models[i]->Draw("L SAME");
    //gr_data_uncr->Draw("P SAME");

    // legend
    TLegend *l1 = SetLegend(0.72,0.65,0.99,0.96);
    for(Int_t i = 0; i < 7; i++) l1->AddEntry(gr_models[i],str_models[i],"L");
    l1->SetTextSize(0.040);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    c1->Print("Results/" + str_subfolder + "CrossSec/Plot/models_only.pdf");

    // plot the measurement only
    TCanvas *c2 = new TCanvas("c2","c2",900,800);
    c2->SetLogy();  
    // margins
    c2->SetTopMargin(0.02);
    c2->SetBottomMargin(0.14);
    c2->SetRightMargin(0.02);
    c2->SetLeftMargin(0.14);
    // plot the graphs
    TH1 *h2 = (TH1*) gr_data_corr->GetHistogram();
    h2->SetTitle(";|#it{t}| (GeV^{2}); d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");
    h2->SetMaximum(0.04);
    // vertical axis
    h2->GetYaxis()->SetTitleSize(0.045);
    h2->GetYaxis()->SetTitleOffset(1.52);
    h2->GetYaxis()->SetLabelSize(0.045);
    // horizontal axis
    h2->GetXaxis()->SetTitleSize(0.045);
    h2->GetXaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->SetLabelSize(0.045);
    h2->GetXaxis()->SetRangeUser(0.04,1.0);
    // draw everything
    gr_data_corr->Draw("A5");
    gr_data_uncr->Draw("PZ SAME");
    // legends
    gStyle->SetTextFont(42);
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.045);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.56,0.92,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // draw legend with data+unc. description
    gStyle->SetLegendBorderSize(0);
    TLegend *l2 = SetLegend(0.44,0.70,0.76,0.88);
    l2->SetTextSize(0.042);
    l2->SetMargin(0.15);
    l2->AddEntry((TObject*)0,"ALICE incoherent J/#psi, |y| < 0.8", "");
    l2->AddEntry(gr_data_uncr,"Uncorrelated stat. + syst.", "EPL");
    l2->AddEntry(gr_data_corr,"Correlated syst.", "F");
    l2->Draw();

    c2->Print("Results/" + str_subfolder + "CrossSec/Plot/measurement.pdf");

    return;
}