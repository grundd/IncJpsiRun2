// PhotoCrossSec_KSTest.C
// David Grund, July 14, 2022

// my headers
#include "PhotoCrossSec_Utilities.h"

// function to perform Kilmogorov-Smirnov test
Double_t IntegrateData(Double_t t_min, Double_t t_max);
Double_t DoKSTest(Int_t n_data, Double_t *t_val, Double_t *sigma_val);

void PhotoCrossSec_KSTest(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/KSTest/");

    ReadInput_HSModel();
    DoKSTest(nData_HS,abs_t_HS,sig_HS_inc_hs);

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

Double_t DoKSTest(Int_t n_data, Double_t *t_val, Double_t *sigma_val)
{
    Double_t t_low, t_upp;
    Double_t t_step = 0.04;
    ReadInput_Measurement();
    
    Double_t integral_data, integral_model;
    for(Int_t i = 1; i < 25; i++)
    {
        t_low = 0.04;
        t_upp = 0.04 + i * t_step;
        TCanvas *c = new TCanvas(Form("cBin%i",i),Form("cBin%i",i),900,600);
        TString str = "";
        integral_data = IntegrateData(t_upp);
        integral_model = GraphIntegral_Calculate(c,str,n_data,t_val,sigma_val,t_low,t_upp);
    }

    return -1;
}