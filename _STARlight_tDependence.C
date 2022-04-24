// _STARlight_tDependence.C
// David Grund, Apr 24, 2022

// cpp headers
#include <fstream>
#include <vector>
// root headers
#include "TSystem.h"
#include "TH1.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "_STARlight_Utilities.h"

TString str_in = "Trees/STARlight/IncJ_tDep/";
Double_t sig_gPb = 0.015499; // see output.txt, line y = +0.000
Double_t rap_cut = 0.01; // rapidity cut (we need y around 0)
Double_t t_step = 0.02; // GeV^2

void Calculate_tDep(Double_t t_low, Double_t t_upp);

void _STARlight_tDependence()
{
    Calculate_tDep(0.00,1.10);
    
    Calculate_tDep(0.00,2.50);

    return;
}

void Calculate_tDep(Double_t t_low, Double_t t_upp)
{
    Int_t t_nBins = (t_upp - t_low) / t_step;
    Printf("%i bins will be defined.", t_nBins);

    Double_t pT_low = TMath::Sqrt(t_low);
    Double_t pT_upp = TMath::Sqrt(t_upp);

    TFile *fSL = TFile::Open((str_in + "tree_STARlight.root").Data(), "read");
    if(fSL) Printf("Input file SL loaded.");

    TTree *tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
    if(tSL) Printf("Tree %s loaded.", tSL->GetName());
    ConnectTreeVariables_tSL(tSL);

    TFile *fPt = TFile::Open((str_in + "trees_tPtGammaVMPom.root").Data(), "read");
    if(fPt) Printf("Input file Pt loaded.");

    TList *lPt = (TList*) fPt->Get("TreeList");
    if(lPt) Printf("List %s loaded.", lPt->GetName()); 

    TTree *tPtGammaVMPom = (TTree*)lPt->FindObject("tPtGammaVMPom");
    if(tPtGammaVMPom) Printf("Tree %s loaded.", tPtGammaVMPom->GetName());
    ConnectTreeVariables_tPtGammaVMPom(tPtGammaVMPom);

    Printf("STARlight tree contains %lli entries.", tSL->GetEntries());
    Printf("tPtGammaVMPom tree contains %lli entries.", tPtGammaVMPom->GetEntries());

    TH1D *hSigmaPhotoNuc = new TH1D("hSigmaPhotoNuc","hSigmaPhotoNuc", t_nBins, t_low, t_upp);

    Double_t nEv_tot = 0;
    for(Int_t iEntry = 0; iEntry < tSL->GetEntries(); iEntry++)
    {
        tSL->GetEntry(iEntry);
        tPtGammaVMPom->GetEntry(iEntry);

        // if the values differ by more than 1% => something wrong
        if(TMath::Abs(fPtVM - parent->Pt())/fPtVM > 0.01){
            Printf("(!) Entry no. %i: fPtVM = %.6f, from SL = %.6f", iEntry+1, fPtVM, parent->Pt());
        }

        if(TMath::Abs(parent->Rapidity()) < rap_cut){
            nEv_tot++;
            if(parent->Pt() > pT_low && parent->Pt() < pT_upp){
                hSigmaPhotoNuc->Fill(fPtPm*fPtPm);
            }
        }
    }

    Printf("No. of events with |y| < %.2f: %.0f", rap_cut, nEv_tot);
    Printf("No. of entries in histogram: %.0f", hSigmaPhotoNuc->GetEntries());

    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Plot the histogram with events
    TCanvas *c1 = new TCanvas("c1","c1",900,600);
    c1->SetLogy(); 
    // Margins
    c1->SetTopMargin(0.03);
    c1->SetBottomMargin(0.14);
    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.12);    
    hSigmaPhotoNuc->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); Events per bin");
    hSigmaPhotoNuc->SetLineStyle(1);
    hSigmaPhotoNuc->SetLineColor(kRed);
    hSigmaPhotoNuc->SetLineWidth(2);

    hSigmaPhotoNuc->Draw("HIST");
    gSystem->Exec("mkdir -p Trees/PhotoCrossSec/STARlight/");
    c1->Print(Form("Trees/PhotoCrossSec/STARlight/gen_events_%.2f-%.2f.pdf", t_low, t_upp));
    c1->Print(Form("Trees/PhotoCrossSec/STARlight/gen_events_%.2f-%.2f.png", t_low, t_upp));

    // Normalize the histogram
    Bool_t debug = kFALSE;
    hSigmaPhotoNuc->Scale(1.0/nEv_tot);
    if(debug) Printf("Integral is %.5f", hSigmaPhotoNuc->Integral());
    if(debug) Printf("Integral is %.5f", hSigmaPhotoNuc->Integral("width"));
    hSigmaPhotoNuc->Scale(sig_gPb);
    if(debug) Printf("Integral is %.5f", hSigmaPhotoNuc->Integral());
    if(debug) Printf("Integral is %.5f", hSigmaPhotoNuc->Integral("width"));
    hSigmaPhotoNuc->Scale(1.0, "width");
    if(debug) Printf("Integral is %.5f", hSigmaPhotoNuc->Integral());
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral("width"));

    // Check the value of the integral
    Double_t integral = 0;
    for(Int_t iBin = 1; iBin <= t_nBins; iBin++){
        integral += hSigmaPhotoNuc->GetBinContent(iBin) * (hSigmaPhotoNuc->GetBinLowEdge(iBin+1) - hSigmaPhotoNuc->GetBinLowEdge(iBin));
    }
    Printf("Manually calculated integral is %.5f", integral);

    // Fill the arrays
    vector<Double_t> t_abs; 
    vector<Double_t> sigma; 
    for(Int_t iBin = 1; iBin <= t_nBins; iBin++){
        t_abs.push_back(hSigmaPhotoNuc->GetBinCenter(iBin));
        sigma.push_back(hSigmaPhotoNuc->GetBinContent(iBin));
    }
    // Define TGraph
    Double_t *sigma_ptr, *t_abs_ptr;
    sigma_ptr = &sigma[0];
    t_abs_ptr = &t_abs[0];
    TGraph *grSL = new TGraph(t_nBins, t_abs_ptr, sigma_ptr);
    grSL->SetLineStyle(1);
    grSL->SetLineColor(kBlue);
    grSL->SetLineWidth(2);
    grSL->GetXaxis()->SetRangeUser(t_low,t_upp);
    grSL->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2}); d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");

    //grSL->Print();
    
    // Plots
    TCanvas *c2 = new TCanvas("c2","c2",900,600);
    c2->SetLogy(); 
    // Margins
    c2->SetTopMargin(0.03);
    c2->SetBottomMargin(0.14);
    c2->SetRightMargin(0.03);
    c2->SetLeftMargin(0.12);

    TLegend *l1 = new TLegend(0.55,0.75,0.80,0.95);
    l1->SetMargin(0.);
    l1->AddEntry((TObject*)0,Form("Total #sigma_{#gammaPb} = %.6f #mub", hSigmaPhotoNuc->Integral("width")), "");
    l1->SetTextSize(0.045);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent

    grSL->Draw("AL");
    l1->Draw();
    c2->Print(Form("Trees/PhotoCrossSec/STARlight/sigma_gPb_%.2f-%.2f.pdf", t_low, t_upp));
    c2->Print(Form("Trees/PhotoCrossSec/STARlight/sigma_gPb_%.2f-%.2f.png", t_low, t_upp));

    // Print the results to text file
    TString str_out = Form("Trees/PhotoCrossSec/STARlight/IncJ_tDep_%.2f-%.2f.txt", t_low, t_upp);
    ofstream outfile(str_out.Data());
    for(Int_t iBin = 0; iBin < t_nBins; iBin++)
    {
        outfile << Form("%.2f \t%.6f \n", t_abs[iBin], sigma[iBin]);
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str_out.Data());
    
    return;
}