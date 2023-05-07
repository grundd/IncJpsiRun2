// CrossSec_ExpFits.C
// David Grund, Sep 04, 2022

// my headers
#include "CrossSec_Utilities.h"

TGraphErrors *gr_data = NULL;
TGraph *gr_binned[7] = { NULL };
Double_t par_val[3];
Double_t par_err[3];

void DoExpFit(Int_t iM, TString fit);
// iM ... index of a model
// iM = 0 ... 6 => models, == 7 => data
// fit = exp_pure => pure exponential
//     = exp_quad => exponential with a quadratic term
//     = log_line => fitting log points with a straight line

void CrossSec_ExpFits(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/ExpFits/exp_pure/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/ExpFits/exp_quad/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/ExpFits/log_line/");

    // open the file with histograms and graphs
    TFile *f = TFile::Open("Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/histograms_and_graphs.root","read");
    if(f) Printf("Input file %s loaded.", f->GetName()); 
    TList *l = (TList*) f->Get("graphs");
    if(l) Printf("List %s loaded.", l->GetName());
    // load the graphs
    gr_data_uncr = (TGraphAsymmErrors*)l->FindObject("gr_data_uncr");
    gr_data_corr = (TGraphAsymmErrors*)l->FindObject("gr_data_corr");
    for(Int_t i = 0; i < 9; i++) {
        gr_models[i] = (TGraph*)l->FindObject("gr_" + str_models[i]);
        gr_binned[i] = (TGraph*)l->FindObject("grBinned_" + str_models[i]);
    }

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
       << "model\t\ta [mub GeV^2] \tb [GeV^2]\n";
    // fit data
    if(true) {
        os << "data/model\t"; 
        DoExpFit(7,"log_line");
        os << par_val[0] * 1e3 << "\t" << par_err[0] * 1e3 << "\t"
           << par_val[1] << "\t" << par_err[1] << "\n";
        DoExpFit(7,"exp_pure");
        DoExpFit(7,"exp_quad");
    }
    // models
    if(true) {
        for(Int_t i = 0; i < 7; i++) {
            os << Form("%s\t", str_models[i].Data());
            if(i != 0 && i != 5) os << "\t";
            DoExpFit(i,"log_line");
            os << par_val[0] * 1e3 << "\t" << par_err[0] * 1e3 << "\t"
               << par_val[1] << "\t" << par_err[1] << "\n";
            DoExpFit(i,"exp_pure");
            DoExpFit(i,"exp_quad");
        } 
    }
    os.close();
    return;
}

void DoExpFit(Int_t iM, TString fit)
{
    TGraph *gr = NULL;
    TGraph *gr_log = new TGraph(nPtBins);
    if(iM == 7) {
        gr = gr_data;
        for(Int_t i = 0; i < nPtBins; i++) gr_log->SetPoint(i,gr_data->GetPointX(i),TMath::Log(gr_data->GetPointY(i)));
    } else {
        gr = gr_binned[iM];
        for(Int_t i = 0; i < nPtBins; i++) gr_log->SetPoint(i,gr_binned[iM]->GetPointX(i),TMath::Log(gr_binned[iM]->GetPointY(i)));
    }

    TF1 *f = NULL;
    TF1 *f_exp_pure = new TF1("f_exp_pure", "[0] * exp(-[1] * x)", 0.04, 1.0);
    TF1 *f_exp_quad = new TF1("f_exp_quad", "[0] * exp(-[1] * x + [2] * x * x)", 0.04, 1.0);
    TF1 *f_log = new TF1("f_log", "[0] - x * [1]", 0.04, 1.0);
    if(fit == "exp_pure" || fit == "log_line") {
        f = f_exp_pure;
        f->SetLineColor(kRed);
    } else if(fit == "exp_quad") {
        f = f_exp_quad;
        f->SetLineColor(kBlue);
    }
    f->SetLineStyle(9);
    f->SetLineWidth(2);
    
    // https://root.cern.ch/doc/master/classTH1.html#HFitRange
    if(fit != "log_line") gr->Fit(f, "R");
    else {
        gr_log->Fit(f_log, "R");
        f->SetParameter(0, TMath::Exp(f_log->GetParameter(0)));
        f->SetParameter(1, f_log->GetParameter(1));
    }
    int npar = 0;
    if(fit == "exp_pure") npar = 2;
    else if(fit == "exp_quad") npar = 3;
    else if(fit == "log_line") npar = 2;
    for(Int_t ipar = 0; ipar < npar; ipar++) {
        par_val[ipar] = f->GetParameter(ipar);
        par_err[ipar] = f->GetParError(ipar);
    }
    // set data properties 
    gStyle->SetEndErrorSize(4);         
    gr->SetMarkerStyle(kFullCross);
    gr->SetMarkerSize(1.0);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kBlack);
    // canvas
    TCanvas *cFit = new TCanvas("cFit","cFit",900,800);
    // margins
    cFit->SetTopMargin(0.03);
    cFit->SetBottomMargin(0.11);
    cFit->SetRightMargin(0.03);
    cFit->SetLeftMargin(0.14);
    cFit->SetLogy();
    // plot the graph
    TH1F* h = gPad->DrawFrame(0.04, 0.0004, 1.0, 0.06);
    h->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb GeV^{-2})");
    // x-axis
    h->GetXaxis()->SetTickLength(0.025); 
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetDecimals(1);   
    // y-axis
    h->GetYaxis()->SetTickLength(0.025); 
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.48);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetDecimals(1);
    
    h->Draw("AXIS");
    // draw graph and curves
    gr->Draw("P SAME");
    f->Draw("L SAME");

    TString name;
    if(iM == 7) name = "ALICE measurement";
    else        name = str_models[iM];

    TLegend *l = new TLegend(0.48,0.70,0.95,0.96);
    l->AddEntry((TObject*)0,Form("#bf{%s}",name.Data()),"");
    if(fit == "exp_pure") {
        l->AddEntry(f,"#it{a} exp(#minus#it{bt})","L");
        l->AddEntry((TObject*)0,Form("#it{a} = (%.2f #pm %.2f) #mub GeV^{-2}", par_val[0] * 1e3, par_err[0] * 1e3),"");
        l->AddEntry((TObject*)0,Form("#it{b} = (%.2f #pm %.2f) GeV^{-2}", par_val[1], par_err[1]),"");
        l->AddEntry((TObject*)0,Form("integral: %.2f #mub", f->Integral(0.04,1.0) * 1e3),""); // over 0.04 < |t| < 1.0 GeV^{2} #it{c}^{-2}
    } else if(fit == "exp_quad") {
        l->AddEntry(f,"#it{a} exp(#minus#it{bt} + #it{ct}^{2})","L");
        l->AddEntry((TObject*)0,Form("#it{a} = (%.2f #pm %.2f) #mub GeV^{-2}", par_val[0] * 1e3, par_err[0] * 1e3),"");
        l->AddEntry((TObject*)0,Form("#it{b} = (%.2f #pm %.2f) GeV^{-2}", par_val[1], par_err[1]),"");
        l->AddEntry((TObject*)0,Form("#it{c} = (%.2f #pm %.2f) GeV^{-4}", par_val[2], par_err[2]),"");
        l->AddEntry((TObject*)0,Form("integral: %.2f #mub", f->Integral(0.04,1.0) * 1e3),""); // over 0.04 < |t| < 1.0 GeV^{2} #it{c}^{-2}
    } else if(fit == "log_line") {
        l->AddEntry(f,"#it{a} exp(#minus#it{bt})","L");
        l->AddEntry((TObject*)0,Form("#it{a} = (%.2f #pm %.2f) #mub GeV^{-2}", par_val[0] * 1e3, par_err[0] * 1e3),"");
        l->AddEntry((TObject*)0,Form("#it{b} = (%.2f #pm %.2f) GeV^{-2}", par_val[1], par_err[1]),"");
        l->AddEntry((TObject*)0,Form("integral: %.2f #mub", f->Integral(0.04,1.0) * 1e3),""); // over 0.04 < |t| < 1.0 GeV^{2} #it{c}^{-2}
    }
    l->SetTextSize(0.039);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetMargin(0.16);
    l->Draw();
    
    if(iM == 7) cFit->Print("Results/" + str_subfolder + "CrossSec/ExpFits/" + fit + "/" + name + ".pdf");
    else cFit->Print("Results/" + str_subfolder + "CrossSec/ExpFits/" + fit + "/" + name + ".pdf");

    delete cFit;
    delete gr_log;
    delete f_exp_pure;
    delete f_exp_quad;
    delete f_log;
    return;
}