// STARlight_OptimalRA.C
// David Grund, Mar 20, 2022

// cpp headers
#include <fstream>
// root headers
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

Double_t R_A_val[13] = { 0 };
Double_t chi2NDF[13] = { 0 };

void STARlight_OptimalRA_FindMinimum();

void STARlight_OptimalRA(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    STARlight_OptimalRA_FindMinimum();

    return;
}

void STARlight_OptimalRA_FindMinimum()
{
    // read the input files
    for(Int_t iRA = 0; iRA < 13; iRA++){
        Double_t R_A = 6.6 + iRA*0.1;
        TString str_in = "Results/" + str_subfolder + Form("PtFit_NoBkg/OptimalRA/Fits/modRA_%.2f_chi2.txt", R_A);
        
        ifstream ifs;
        ifs.open(str_in.Data());
        if(!ifs.fail()){
            ifs >> R_A_val[iRA] >> chi2NDF[iRA];
        } else Printf("File %s not found.", str_in.Data());
        ifs.close();
        Printf("R_A = %.2f\t chi2/NDF = %.3f", R_A_val[iRA], chi2NDF[iRA]);
    }

    TGraph *gr_RA = new TGraph(13, R_A_val, chi2NDF); 
    gr_RA->SetMarkerStyle(kFullSquare);
    gr_RA->SetMarkerColor(kBlack);
    gr_RA->SetMarkerSize(1);

    // Fit using a parabola
    TF1 *f = new TF1("f", "[2] * x * x + [1] * x + [0]"); 
    f->SetLineStyle(1);
    f->SetLineColor(215);
    f->SetLineWidth(3);
    gr_RA->Fit(f);

    Double_t par_val[3];
    Double_t par_err[3];
    for(Int_t ipar = 0; ipar < 3; ipar++){
        par_val[ipar] = f->GetParameter(ipar);
        par_err[ipar] = f->GetParError(ipar);
        Printf("p%i = (%.3f pm %.3f)", ipar, par_val[ipar], par_err[ipar]);
    }
    Double_t RA_min_val = -par_val[1] / 2 / par_val[2];
    Double_t RA_min_err = RA_min_val * TMath::Sqrt(TMath::Power(par_err[1]/par_val[1], 2) + TMath::Power(par_err[2]/par_val[2], 2));
    Printf("Minimum in R_A = (%.3f pm %.3f) fm", RA_min_val, RA_min_err);

    // Prepare dashed line at y = 1
    TLine *line_min = new TLine(RA_min_val,0.6,RA_min_val,1.45);
    line_min->SetLineColor(kRed+1);
    line_min->SetLineWidth(3);
    line_min->SetLineStyle(2);

    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Canvas
    TCanvas *c = new TCanvas("c","c",900,600);
    // Margins
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);
    //Plot the graphs
    TH1 *h = (TH1*) gr_RA->GetHistogram();
    h->SetTitle(";#it{R}_{A} (fm);#chi^{2}/NDF");
    // Vertical axis
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetRangeUser(0.6,4.2);
    // Horizontal axis
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetRangeUser(6.50,7.90);
    // Draw
    gr_RA->Draw("AP");
    f->Draw("AL SAME");
    line_min->Draw("SAME");

    TLegend *l2 = new TLegend(0.2,0.75,0.8,0.9);
    l2->AddEntry((TObject*)0,"STARlight value: #it{R}_{A} = 6.624 fm","");
    l2->AddEntry((TObject*)0,Form("Minimum at #it{R}_{A} = %.2f fm", RA_min_val),"");
    //l2->AddEntry((TObject*)0,Form("Minimum at #it{R}_{A} = (%.3f #pm %.3f) fm", RA_min_val, RA_min_err),"");
    l2->SetTextSize(0.06);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    // Print the results
    TString str_out = "Results/" + str_subfolder + "PtFit_NoBkg/OptimalRA/fit_RA";
    c->Print((str_out + ".pdf").Data());
    c->Print((str_out + ".png").Data());

    return;
}