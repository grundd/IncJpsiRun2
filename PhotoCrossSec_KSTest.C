// PhotoCrossSec_KSTest.C
// David Grund, July 14, 2022

// root headers
#include "TGraph.h"
// my headers
#include "PhotoCrossSec_Utilities.h"

// function to perform Kilmogorov-Smirnov test
Double_t IntegrateData(Double_t t_min, Double_t t_max);
void DoKSTest(TString str, Int_t n_data, Double_t *t_val, Double_t *sigma_val);

void PhotoCrossSec_KSTest(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/KSTest/");

    ReadInput_HSModel();
    TString str_CCK_hs = "CCK_GG-hs"; 
    DoKSTest(str_CCK_hs,nData_HS,abs_t_HS,sig_HS_inc_hs);
    TString str_CCK_n = "CCK_GG-n"; 
    DoKSTest(str_CCK_n,nData_HS,abs_t_HS,sig_HS_inc_n);

    ReadInput_Heikki();
    TString str_MS_fl = "MS_IPsat_flu"; 
    DoKSTest(str_MS_fl,nData_HM,abs_t_HM,sig_HM_fluct);
    TString str_MS_nf = "MS_IPsat_no_flu"; 
    DoKSTest(str_MS_nf,nData_HM,abs_t_HM,sig_HM_noflu);

    ReadInput_Guzey();
    for(Int_t i = 0; i < nData_GZ; i++){
        sig_GZ_tot_min[i] = sig_GZ_tot_min[i] / 1e6;
        sig_GZ_tot_max[i] = sig_GZ_tot_max[i] / 1e6;
    }
    TString str_GZ_up = "GSZ_upp";
    DoKSTest(str_GZ_up,nData_GZ,abs_t_GZ,sig_GZ_tot_max);
    TString str_GZ_lo = "GSZ_low";
    DoKSTest(str_GZ_lo,nData_GZ,abs_t_GZ,sig_GZ_tot_min);

    return;
}

Double_t IntegrateData(Double_t t_max)
{
    // t_min always = 0.04 GeV
    Double_t integral = 0.;
    Int_t tMaxBin = 1;
    while(t_max > t_boundaries[tMaxBin]) tMaxBin++;
    //Printf("t_max = %.3f within bin %i: %.3f to %.3f", t_max, tMaxBin, t_boundaries[tMaxBin-1], t_boundaries[tMaxBin]);
    for(Int_t i = 0; i < tMaxBin; i++) integral += sig_val[i] * (t_boundaries[i+1] - t_boundaries[i]);
    integral += sig_val[tMaxBin-1] * (t_max - t_boundaries[tMaxBin]);

    return integral;
}

void DoKSTest(TString str, Int_t n_data, Double_t *t_val, Double_t *sigma_val)
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/KSTest/" + str + "/");

    Double_t t_low, t_upp;
    Double_t t_step = 0.04;
    Int_t nPoints = 24;
    ReadInput_Measurement();
    
    TGraph *grData = new TGraph();
    TGraph *grModl = new TGraph();
    TGraph *grKS_statistics = new TGraph();
    Double_t KS_statistics_max = 0.;

    Double_t integral_data_tot = IntegrateData(1.0);
    TCanvas *C = new TCanvas("C","C",900,600);
    Double_t integral_modl_tot = GraphIntegral_Calculate(C,str,n_data,t_val,sigma_val,0.04,1.0) * 1e3; // in micro barns
    C->Print("Results/" + str_subfolder + "PhotoCrossSec/KSTest/" + str + "/integral_total.pdf");

    Double_t integral_data, integral_modl;
    for(Int_t i = 0; i <= nPoints; i++)
    {
        t_low = 0.04;
        t_upp = 0.04 + i * t_step;
        TCanvas *c = new TCanvas(Form("cBin%i",i),Form("cBin%i",i),900,600);
        integral_data = IntegrateData(t_upp);
        integral_modl = GraphIntegral_Calculate(c,str,n_data,t_val,sigma_val,t_low,t_upp) * 1e3; // in micro barns
        Printf("|t|: %.2f data: %.2f model: %.2f", t_upp, integral_data, integral_modl);
        grData->AddPoint(t_upp,integral_data / integral_data_tot);
        grModl->AddPoint(t_upp,integral_modl / integral_modl_tot);
        Double_t statistics = TMath::Abs(integral_data / integral_data_tot - integral_modl / integral_modl_tot);
        grKS_statistics->AddPoint(t_upp,statistics);
        if(i == 0) KS_statistics_max = statistics;
        else       if(KS_statistics_max < statistics) KS_statistics_max = statistics;
        c->Print("Results/" + str_subfolder + "PhotoCrossSec/KSTest/" + str + "/" + Form("%.2f-%.2f.pdf",t_low,t_upp));
        delete c;
    }

    TCanvas *c = new TCanvas("c","c",900,600);
    // set the canvas
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.1);
    c->cd();
    // get the histogram
    TH1 *h = (TH1*) grData->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); Cumulative Probability");
    // Vertical axis
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetDecimals(1);
    // Horizontal axis
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetRangeUser(0.04,1.0);
    // print the graph points
    grData->Print();
    grModl->Print();
    // set the graph properties
    // data
    grData->SetLineStyle(9);
    grData->SetLineColor(kRed);
    grData->SetLineWidth(3);  
    grData->Draw("AL");
    // model
    grModl->SetLineStyle(1);
    grModl->SetLineColor(kBlue);
    grModl->SetLineWidth(3);  
    grModl->Draw("L SAME");
    // KS test
    grKS_statistics->SetMarkerStyle(2);
    grKS_statistics->SetMarkerColor(1);
    grKS_statistics->SetMarkerSize(2);
    grKS_statistics->Draw("P SAME");
    // legend
    TLegend *l = new TLegend(0.65,0.35,0.95,0.70);
    l->SetMargin(0.15);
    l->AddEntry((TObject*)0,Form("%s",str.Data()),"");
    l->AddEntry(grData,"CDF_{data}","L");
    l->AddEntry(grModl,"CDF_{model}","L");
    l->AddEntry(grKS_statistics,"Abs(CDF_{data}- CDF_{model})","P");
    l->AddEntry((TObject*)0,Form("Max difference: %.2f",KS_statistics_max),"");
    l->SetTextSize(0.045);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();
    // print the canvas
    c->Print("Results/" + str_subfolder + "PhotoCrossSec/KSTest/" + "KS_" + str + ".pdf");
    delete c;

    return;
}