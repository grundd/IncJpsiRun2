// CrossSec_KSTest.C
// David Grund, Aug 28, 2023

// my headers
#include "CrossSec_Utilities.h"

TH1F* hBinned_data = NULL;
TH1F* hBinned_models[9] = { NULL };

void DoKSTest (int iModel)
{
    TGraph *gr_data = new TGraph();
    TGraph *gr_model = new TGraph();
    TGraph *gr_KS = new TGraph(); // K-S statistics
    float KS_max = 0.;

    // total integrals
    float int_tot_data = hBinned_data->Integral(1,nPtBins);
    float int_tot_model = hBinned_models[iModel]->Integral(1,nPtBins);

    // integrals at bin edges
    for(int i = 0; i <= nPtBins; i++)
    {
        float t_low = tBoundaries[0]; // 0.04 GeV^2
        float t_upp = tBoundaries[i]; // up to 1 GeV^2
        float int_data(0.), int_model(0.);
        if(i > 0) {
            int_data = hBinned_data->Integral(1,i);
            int_model = hBinned_models[iModel]->Integral(1,i);
        }
        cout << "|t|: " << t_upp << ", data: " << int_data << ", model: " << int_model << "\n";
        gr_data->AddPoint(t_upp,int_data / int_tot_data);
        gr_model->AddPoint(t_upp,int_model / int_tot_model);
        float ks = TMath::Abs(int_data / int_tot_data - int_model / int_tot_model); // value of the KS statistics
        gr_KS->AddPoint(t_upp,ks);
        if(i == 0) KS_max = ks;
        else if(KS_max < ks) KS_max = ks;
    }

    TCanvas *c = new TCanvas("c","c",900,800);
    // set the canvas
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.12);
    c->cd();
    // get the histogram
    TH1 *h = (TH1*) gr_data->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2}); Cumulative Probability");
    // Vertical axis
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetRangeUser(0.04,1.0);
    // print the graph points
    /*
    gr_data->Print();
    gr_model->Print();
    */
    // set the graph properties
    // data
    gr_data->SetLineStyle(9);
    gr_data->SetLineColor(kRed);
    gr_data->SetLineWidth(3);  
    gr_data->Draw("AL");
    // model
    gr_model->SetLineStyle(1);
    gr_model->SetLineColor(kBlue);
    gr_model->SetLineWidth(3);  
    gr_model->Draw("L SAME");
    // KS test
    gr_KS->SetMarkerStyle(70);
    gr_KS->SetMarkerColor(1);
    gr_KS->SetMarkerSize(2.5);
    gr_KS->SetLineWidth(2);
    gr_KS->Draw("P SAME");
    // legend
    int nRows = 5;
    TLegend *l = new TLegend(0.50,0.70-nRows*0.06,0.97,0.70);
    l->AddEntry((TObject*)0,Form("#bf{%s}",str_models[iModel].Data()),"");
    l->AddEntry(gr_data,"CDF_{data}","L");
    l->AddEntry(gr_model,"CDF_{model}","L");
    l->AddEntry(gr_KS,"abs(CDF_{data}#minus CDF_{model})","P");
    l->AddEntry((TObject*)0,Form("max difference: %.3f",KS_max),"");
    l->SetTextSize(0.045);
    l->SetMargin(0.18);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();
    // print the canvas
    c->Print("Results/" + str_subfolder + "CrossSec/KSTest/KS_" + str_models[iModel] + ".pdf");
    delete c;

    return;
}

void CrossSec_KSTest(int iAnalysis)
{
    InitAnalysis(iAnalysis);
    InitObjects();
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/KSTest/");

    // open the file with histograms and graphs
    TFile *f = TFile::Open("Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/histograms_and_graphs.root","read");
    if(f) Printf("Input file %s loaded.", f->GetName());
    TList *lh = (TList*) f->Get("histograms");
    if(lh) Printf("List %s loaded.", lh->GetName());
    TList *lg = (TList*) f->Get("graphs");
    if(lg) Printf("List %s loaded.", lg->GetName());
    // model histograms
    for(int i = 0; i < 9; i++) hBinned_models[i] = (TH1F*)lh->FindObject("hBinned_" + str_models[i]);
    // data graph & histo
    gr_data_uncr = (TGraphAsymmErrors*)lg->FindObject("gr_data_uncr");
    hBinned_data = new TH1F("hBinned_data","",nPtBins,tBoundaries);
    for(int i = 1; i <= nPtBins; i++) {
        hBinned_data->SetBinContent(i, gr_data_uncr->GetPointY(i-1));
        hBinned_data->SetBinError(i, gr_data_uncr->GetErrorYhigh(i-1));
    }

    for(int i = 0; i < 9; i++) DoKSTest(i);

    return;
}