// CrossSec_ExpFits.C
// David Grund, Sep 04, 2022

// my headers
#include "CrossSec_Utilities.h"

TGraphErrors *gr_data = NULL;
Double_t par_pure_val[2];
Double_t par_pure_err[2];
Double_t par_quad_val[3];
Double_t par_quad_err[3];

void DoExpFit(Int_t iM);
// iM ... index of a model
// iM == 0 ... 6 => models, == 7 => data

void CrossSec_ExpFits(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/ExpFits/");


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

    // prepare graphs for data without x errors
    gr_data = new TGraphErrors(nPtBins);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        gr_data->SetPoint(i,gr_data_uncr->GetPointX(i),gr_data_uncr->GetPointY(i));
        gr_data->SetPointError(i,0.,gr_data_uncr->GetErrorY(i));
    }

    // do the fits and print the results
    ofstream os;
    os.open("Results/" + str_subfolder + "CrossSec/ExpFits/#parameters.txt");
    os << std::fixed << std::setprecision(2)
       << "\t\tpure exponential \t\texp with a quadratic term\n"
       << "model \t\tA \terr \tb \terr \tA \terr \tb \terr \tc \terr \n";
    DoExpFit(7);
    os << "data\t\t" << par_pure_val[0] * 1e3 << "\t" << par_pure_err[0] * 1e3 << "\t"
                     << par_pure_val[1] << "\t" << par_pure_err[1] << "\t"
                     << par_quad_val[0] * 1e3 << "\t" << par_quad_err[0] * 1e3 << "\t"
                     << par_quad_val[1] << "\t" << par_quad_err[1] << "\t"
                     << par_quad_val[2] << "\t" << par_quad_err[2] << "\n";
    // do exp fits
    for(Int_t i = 0; i < 7; i++){
        DoExpFit(i);
        os << Form("%s\t", str_models[i].Data());
        if(i == 3 || i == 4) os << "\t";
        os << par_pure_val[0] * 1e3 << "\t" << par_pure_err[0] * 1e3 << "\t"
           << par_pure_val[1] << "\t" << par_pure_err[1] << "\t"
           << par_quad_val[0] * 1e3 << "\t" << par_quad_err[0] * 1e3 << "\t"
           << par_quad_val[1] << "\t" << par_quad_err[1] << "\t"
           << par_quad_val[2] << "\t" << par_quad_err[2] << "\n";
    } 
    os.close();

    return;
}

void DoExpFit(Int_t iM)
{
    TGraph *gr = NULL;
    if(iM == 7) gr = gr_data;
    else        gr = gr_models[iM];

    TF1 *f_exp_pure = new TF1("f_exp_pure", "[0] * exp(-[1] * x)", 0.04, 1.00);
    TF1 *f_exp_quad = new TF1("f_exp_quad", "[0] * exp(-[1] * x + [2] * x * x)", 0.04, 1.00);
    // set function properties
    f_exp_pure->SetLineStyle(1);
    f_exp_pure->SetLineColor(kRed);
    f_exp_pure->SetLineWidth(2);
    f_exp_quad->SetLineStyle(1);
    f_exp_quad->SetLineColor(215);
    f_exp_quad->SetLineWidth(2);
    
    // https://root.cern.ch/doc/master/classTH1.html#HFitRange
    gr->Fit(f_exp_pure, "R");
    gr->Fit(f_exp_quad, "R");
    for(Int_t ipar = 0; ipar < 3; ipar++){
        if(ipar != 2){
            par_pure_val[ipar] = f_exp_pure->GetParameter(ipar);
            par_pure_err[ipar] = f_exp_pure->GetParError(ipar);
        }
        par_quad_val[ipar] = f_exp_quad->GetParameter(ipar);
        par_quad_err[ipar] = f_exp_quad->GetParError(ipar);
    }
    // set data properties 
    gStyle->SetEndErrorSize(4);         
    gr->SetMarkerStyle(kFullCross);
    gr->SetMarkerSize(0.7);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kBlack);
    // canvas
    TCanvas *cFit = new TCanvas("cFit","cFit",900,600);
    // margins
    cFit->SetTopMargin(0.03);
    cFit->SetBottomMargin(0.14);
    cFit->SetRightMargin(0.03);
    cFit->SetLeftMargin(0.11);
    cFit->SetLogy();
    // plot the graph
    TH1 *h = (TH1*) gr->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    // vertical axis
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.15);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetDecimals(1);
    //h->SetMinimum(0.00002);
    //h->SetMaximum(0.04);
    // horizontal axis
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.045);
    TAxis *ax = gr->GetXaxis();
    ax->SetLimits(0.04,1.00);
    // draw graph and curves
    gr->Draw("AP SAME");
    f_exp_pure->Draw("AL SAME");
    f_exp_quad->Draw("AL SAME");

    TString name;
    if(iM == 7) name = "data";
    else        name = str_models[iM];

    TLegend *l = new TLegend(0.12,0.18,0.40,0.50);
    l->AddEntry((TObject*)0,name,"");
    l->AddEntry(f_exp_quad,"A exp(-bt + ct^{2})","L");
    l->AddEntry((TObject*)0,Form("A = (%.2f #pm %.2f) #mub #it{c}^{2} GeV^{-2}", par_quad_val[0] * 1e3, par_quad_err[0] * 1e3),"");
    l->AddEntry((TObject*)0,Form("b = (%.2f #pm %.2f) #it{c}^{2} GeV^{-2}", par_quad_val[1], par_quad_err[1]),"");
    l->AddEntry((TObject*)0,Form("c = (%.2f #pm %.2f) #it{c}^{4} GeV^{-4}", par_quad_val[2], par_quad_err[2]),"");
    l->AddEntry((TObject*)0,Form("integral 0.04 < |t| < 1.00 GeV^{2} #it{c}^{-2}: %.2f #mub", f_exp_quad->Integral(0.04,1.00) * 1e3),"");
    l->AddEntry(f_exp_pure,"A exp(-bt)","L");
    l->AddEntry((TObject*)0,Form("A = (%.2f #pm %.2f) #mub #it{c}^{2} GeV^{-2}", par_pure_val[0] * 1e3, par_pure_err[0] * 1e3),"");
    l->AddEntry((TObject*)0,Form("b = (%.2f #pm %.2f) #it{c}^{2} GeV^{-2}", par_pure_val[1], par_pure_err[1]),"");
    l->AddEntry((TObject*)0,Form("integral 0.04 < |t| < 1.00 GeV^{2} #it{c}^{-2}: %.2f #mub", f_exp_pure->Integral(0.04,1.00) * 1e3),"");
    l->SetTextSize(0.032);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    if(iM == 7) cFit->Print("Results/" + str_subfolder + "CrossSec/ExpFits/" + name + ".pdf");
    else cFit->Print("Results/" + str_subfolder + "CrossSec/ExpFits/" + name + ".pdf");

    return;
}