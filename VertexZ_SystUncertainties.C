// VertexZ_SystUncertainties.C
// David Grund, June 17, 2022

// cpp headers
#include <fstream>
// root headers
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "AxE_PtBins_Utilities.h"

void NewCutZ_CompareCounts();
void NewCutZ_FillHistograms(TTree *t, TH1D *h, Double_t fCutZ);
void NewCutZ_PrepareTree(Double_t fCutZ);
void NewCutZ_AxE_PtBins(Double_t fCutZ);
void ConnectTreeVariables_tData(TTree *t);
Double_t CalculateErrorBinomial(Double_t k, Double_t n);

Bool_t debug = kTRUE;

void VertexZ_SystUncertainties(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VertexZ_SystUncertainties/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VertexZ_SystUncertainties/");

    NewCutZ_CompareCounts();

    return;
}

void NewCutZ_CompareCounts()
{
    NewCutZ_PrepareTree(15.0);
    NewCutZ_PrepareTree(10.0);

    // load events with Z_cut < 15 cm
    TFile *f15 = TFile::Open(("Trees/" + str_subfolder + "VertexZ_SystUncertainties/Zcut15.0_DataTree.root").Data(), "read");
    if(f15) Printf("Input data with Z_vtx < 15 cm loaded.");
    TTree *t15 = dynamic_cast<TTree*> (f15->Get("tData"));
    if(t15) Printf("Input tree loaded.");
    // load events with Z_cut < 10 cm
    TFile *f10 = TFile::Open(("Trees/" + str_subfolder + "VertexZ_SystUncertainties/Zcut10.0_DataTree.root").Data(), "read");
    if(f10) Printf("Input data with Z_vtx < 10 cm loaded.");
    TTree *t10 = dynamic_cast<TTree*> (f10->Get("tData"));
    if(t10) Printf("Input tree loaded.");

    TH1D *hEv15 = new TH1D("hEv15","hEv15",nPtBins,ptBoundaries);
    TH1D *hEv10 = new TH1D("hEv10","hEv10",nPtBins,ptBoundaries);

    NewCutZ_FillHistograms(t15,hEv15,15.);
    NewCutZ_FillHistograms(t10,hEv10,10.);
    
    // histogram of ratios with sumw2
    TH1D *hEvRat = (TH1D*)hEv10->Clone("hEvRat");
    hEvRat->SetTitle("hEvRat");
    hEvRat->Sumw2();
    hEvRat->Divide(hEv15);
    // histogram of ratios w/o sumw2
    TH1D *hEvRat_NoSumw2 = (TH1D*)hEv10->Clone("hEvRat_NoSumw2");
    hEvRat_NoSumw2->SetTitle("hEvRat_NoSumw2");
    hEvRat_NoSumw2->Divide(hEv15);

    TCanvas *cEv10 = new TCanvas("cEv10","cEv10",900,600);
    hEv10->Draw("E0");
    TCanvas *cEv15 = new TCanvas("cEv15","cEv15",900,600);
    hEv15->Draw("E0");
    TCanvas *cEvRat = new TCanvas("cEvRat","cEvRat",900,600);
    hEvRat->Draw("E0");

    // calculate the ratios of yields manually using various methods to compute errors
    // index 0 -> the 'allbins' range
    // indices 1 to 5 -> 5 pT bins
    Double_t nEv10_val[6] = { 0 };
    Double_t nEv10_err[6] = { 0 };
    Double_t nEv15_val[6] = { 0 };
    Double_t nEv15_err[6] = { 0 };
    Double_t nEvRat_val[6] = { 0 };
    Double_t nEvRat_err1[6] = { 0 };
    Double_t nEvRat_err2[6] = { 0 };
    Double_t nEvRat_err3[6] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        if(iBin == 0) {
            for(Int_t iBin = 1; iBin <= nPtBins; iBin++) {
                nEv10_val[0] += hEv10->GetBinContent(iBin);
                nEv15_val[0] += hEv15->GetBinContent(iBin);
            }
            nEv10_err[0] = TMath::Sqrt(nEv10_val[0]);
            nEv15_err[0] = TMath::Sqrt(nEv15_val[0]);
        } else {
            nEv10_val[iBin] = hEv10->GetBinContent(iBin);
            nEv10_err[iBin] = hEv10->GetBinError(iBin); // poisson error
            nEv15_val[iBin] = hEv15->GetBinContent(iBin);
            nEv15_err[iBin] = hEv15->GetBinError(iBin); // poisson error
        }
        // calculate the ratio
        nEvRat_val[iBin] = nEv10_val[iBin] / nEv15_val[iBin];
        // calculate the error using error propagation formula
        nEvRat_err1[iBin] = nEvRat_val[iBin] * TMath::Sqrt(TMath::Power(nEv10_err[iBin] / nEv10_val[iBin], 2) 
            + TMath::Power(nEv15_err[iBin] / nEv15_val[iBin], 2));
        // calculate the error using binomial distribution
        // N = number of events with Z_vtx < 15 cm in a given bin
        // p = we estimate it as N(Z_vtx < 15 cm) / N(Z_vtx < 10 cm)
        // q = 1 - p
        Double_t N = nEv15_val[iBin];
        Double_t p = nEv10_val[iBin] / nEv15_val[iBin];
        Double_t q = 1 - p;
        Double_t variance = N * p * q;
        Double_t sigma = TMath::Sqrt(variance); // this is the error of the mean (Np, hEv10->GetBinContent(iBin+1))
        // to get the error of the ratio, we need to divide the previous value by N again
        nEvRat_err2[iBin] = sigma / N;
        // we would get the same result by using the function CalculateErrorBinomial(.,.), see the paper by Ullrich
        // calculate the error using a Bayesian formula
        nEvRat_err3[iBin] = CalculateErrorBayes(nEv10_val[iBin], nEv15_val[iBin]);
    }
    // print the results
    TString str_out = "Results/" + str_subfolder + "VertexZ_SystUncertainties/";
    ofstream outfile(Form("%snEv.txt",str_out.Data()));
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        outfile << "***\n";
        if(iBin == 0) outfile << "Bin: allbins\n";
        else          outfile << "Bin: " << iBin << "\n";
        outfile << Form("nEv15: %.0f pm %.3f (Poisson)\n", nEv15_val[iBin], nEv15_err[iBin])
                << Form("nEv10: %.0f pm %.3f (Poisson)\n", nEv10_val[iBin], nEv10_err[iBin])
                << Form("ratio: %.4f\n", nEvRat_val[iBin]);
        if(iBin != 0) {
                // error calculated as sqrt(bin content):
        outfile << Form("err: %.4f (ROOT w/o Sumw2 = Poisson)\n", hEvRat_NoSumw2->GetBinError(iBin))
                // error calculated as sqrt(sum of squares of weights)
                // https://www-zeuthen.desy.de/~wischnew/amanda/discussion/wgterror/working.html 
                << Form("err: %.4f (ROOT Sumw2)\n", hEvRat->GetBinError(iBin));
        }    
                // error calculated from the error propagation formula: 
        outfile << Form("err: %.4f (err propagation formula)\n", nEvRat_err1[iBin])
                // error calculated from a binomial distribution
                << Form("err: %.4f (binomial)\n", nEvRat_err2[iBin])
                // error calculated using Bayesian formula
                << Form("err: %.4f (Bayes)\n", nEvRat_err3[iBin]);         
    }
    outfile << "***\n";
    outfile.close();
    // calculate the values of AxE with the cut on Z_vtx < 15 cm and < 10 cm
    NewCutZ_AxE_PtBins(15.0);
    NewCutZ_AxE_PtBins(10.0);
    // to calculate the ratios of AxE: we will compare NRec (NGen are the same in all bins)
    // the errors of the ratios will again be calculated from binomial distribution
    Double_t nNRec10_val[6] = { 0 };
    Double_t nNRec15_val[6] = { 0 };
    TString str_AxE15 = "";
    if(cut_fVertexZ == 15.0) str_AxE15 = "Results/" + str_subfolder + Form("AxE_PtBins/NRec_%ibins.txt", nPtBins);
    else                     str_AxE15 = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut15.0_AxE_PtBins/NRec_%ibins.txt", nPtBins);
    TString str_AxE10 = "";
    if(cut_fVertexZ == 10.0) str_AxE10 = "Results/" + str_subfolder + Form("AxE_PtBins/NRec_%ibins.txt", nPtBins);
    else                     str_AxE10 = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut10.0_AxE_PtBins/NRec_%ibins.txt", nPtBins);
    // load the values of NRec
    ifstream ifs;
    Int_t bin;
    ifs.open(str_AxE15);    
    for(Int_t iBin = 1; iBin < nPtBins+1; iBin++) ifs >> bin >> nNRec15_val[iBin];
    Printf("Values of AxE for Z_vtx < 15 cm calculated.");
    ifs.close();
    ifs.open(str_AxE10);    
    for(Int_t iBin = 1; iBin < nPtBins+1; iBin++) ifs >> bin >> nNRec10_val[iBin];
    Printf("Values of AxE for Z_vtx < 10 cm calculated.");
    ifs.close();
    // calculate the ratios of yields manually using various methods to compute errors
    // index 0 -> the 'allbins' range
    // indices 1 to 5 -> 5 pT bins
    Double_t nNRec10_err[6] = { 0 };
    Double_t nNRec15_err[6] = { 0 };
    Double_t nNRecRat_val[6] = { 0 };
    Double_t nNRecRat_err1[6] = { 0 };
    Double_t nNRecRat_err2[6] = { 0 };
    Double_t nNRecRat_err3[6] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        if(iBin == 0) {
            for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
                nNRec10_val[0] += nNRec10_val[iBin];
                nNRec15_val[0] += nNRec15_val[iBin];
            }
        }
        nNRec10_err[iBin] = TMath::Sqrt(nNRec10_val[iBin]);
        nNRec15_err[iBin] = TMath::Sqrt(nNRec15_val[iBin]);
        // calculate the ratio
        nNRecRat_val[iBin] = nNRec10_val[iBin] / nNRec15_val[iBin];
        // calculate the error using error propagation formula
        nNRecRat_err1[iBin] = nNRecRat_val[iBin] * TMath::Sqrt(TMath::Power(nNRec10_err[iBin] / nNRec10_val[iBin], 2) 
            + TMath::Power(nNRec15_err[iBin] / nNRec15_val[iBin], 2));
        // calculate the error using binomial distribution
        nNRecRat_err2[iBin] = CalculateErrorBinomial(nNRec10_val[iBin], nNRec15_val[iBin]);
        // calculate the error using Bayesian formula
        nNRecRat_err3[iBin] = CalculateErrorBayes(nNRec10_val[iBin], nNRec15_val[iBin]);
    }
    // print the results
    outfile.open(Form("%snNrec.txt",str_out.Data()));
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        outfile << "***\n";
        if(iBin == 0) outfile << "Bin: allbins\n";
        else          outfile << "Bin: " << iBin << "\n";
        outfile << Form("nNRec15: %.0f pm %.0f (Poisson)\n", nNRec15_val[iBin], nNRec15_err[iBin])
                << Form("nNRec10: %.0f pm %.0f (Poisson)\n", nNRec10_val[iBin], nNRec10_err[iBin])
                << Form("ratio: %.4f\n", nNRecRat_val[iBin])
                // error calculated from the error propagation formula:
                << Form("err: %.4f (err propagation formula)\n", nNRecRat_err1[iBin])
                // error calculated from a binomial distribution
                << Form("err: %.4f (binomial)\n", nNRecRat_err2[iBin])
                // error calculated using Bayesian formula
                << Form("err: %.4f (Bayes)\n", nNRecRat_err3[iBin]);
    }
    outfile << "***\n";
    outfile.close();
    // calculate the systematic uncertainty
    Double_t syst_uncr[6] = { 0 };
    for(Int_t i = 0; i < nPtBins+1; i++) syst_uncr[i] = (1 - nEvRat_val[i] / nNRecRat_val[i]) * 100;
    // print its values
    outfile.open(Form("%ssyst_uncertainties_%ibins.txt",str_out.Data(),nPtBins));
    outfile << std::fixed << std::setprecision(1);
    for(Int_t i = 0; i < nPtBins+1; i++) outfile << i << "\t" << TMath::Abs(syst_uncr[i]) << "\n";
    outfile.close();

    // plot the differences
    Double_t x_len = 1200;
    if(nPtBins == 5) x_len += 300;
    TCanvas *c = new TCanvas("c","c",x_len,300);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.02);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.13);
    c->Divide(5,1,0,0);
    TGraphErrors *gr[5] = { NULL }; 
    TLegend *lg[5] = { NULL };
    TLine *l = new TLine(0.925,0.925,1.005,1.005);
    l->SetLineColor(kRed);
    l->SetLineStyle(kDashed);
    l->SetLineWidth(2);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        c->cd(i+1);
        gr[i] = new TGraphErrors(1,&nNRecRat_val[i+1],&nEvRat_val[i+1],&nNRecRat_err3[i+1],&nEvRat_err3[i+1]);
        gr[i]->SetMarkerColor(kBlue);
        gr[i]->SetMarkerStyle(5);
        gr[i]->SetMarkerSize(2);
        TAxis *axis = gr[i]->GetXaxis();
        // set range on x-axis
        axis->SetLimits(0.925,1.005);
        TH1 *h = (TH1*) gr[i]->GetHistogram();
        // set range on y-axis
        h->SetMinimum(0.925);
        h->SetMaximum(1.005);
        // vertical axis format
        //h->GetYaxis()->SetTitle("N_{ev}(|z| < 10 cm)/N_{ev}(|z| < 15 cm)");
        h->GetYaxis()->SetTitle("R_{N}^{10/15} (-)");
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleOffset(1.4);
        h->GetYaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetDecimals(2);
        // horizontal axis format
        //h->GetXaxis()->SetTitle("(A#times#varepsilon)_{MC}(|z| < 10 cm)/(A#times#varepsilon)_{MC}(|z| < 15 cm)");
        h->GetXaxis()->SetTitle("R_{A#times#varepsilon}^{10/15} (-)");
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetXaxis()->SetDecimals(2);
        // draw graph and line
        gr[i]->Draw("AP");
        l->Draw("SAME");
        // draw legend
        Double_t x_low(0), x_upp(0);
        if(i == 0){ x_low = 0.19; x_upp = 0.29;}
        else      { x_low = 0.04; x_upp = 0.16;}
        lg[i] = new TLegend(x_low,0.72,x_upp,0.95);
        lg[i]->AddEntry((TObject*)0,Form("R_{N}^{10/15} = (%.1f pm %.1f)%%", nEvRat_val[i+1]*100, nEvRat_err3[i+1]*100),"");
        lg[i]->AddEntry((TObject*)0,Form("R_{A#times#varepsilon}^{10/15} = (%.1f pm %.1f)%%", nNRecRat_val[i+1]*100, nNRecRat_err3[i+1]*100),"");
        lg[i]->AddEntry((TObject*)0,Form("change = %.1f%%", syst_uncr[i+1]),"");
        lg[i]->AddEntry((TObject*)0,Form("syst. uncr. = %.1f%%", TMath::Abs(syst_uncr[i+1])),"");
        lg[i]->SetTextSize(0.05);
        lg[i]->SetBorderSize(0);
        lg[i]->SetFillStyle(0);
        lg[i]->Draw();
    }
    c->Print("Results/" + str_subfolder + "VertexZ_SystUncertainties/ratios.pdf");
    c->Print("Results/" + str_subfolder + "VertexZ_SystUncertainties/ratios.png");

    return;
}

void NewCutZ_FillHistograms(TTree *t, TH1D *h, Double_t fCutZ)
{
    ConnectTreeVariables_tData(t);

    Printf("%lli entries found in the tree with Z_vtx < %.0f cm.", t->GetEntries(), fCutZ);
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++)
    {
        t->GetEntry(iEntry);

        // go over bins in pT and save the number of surviving events with 3.0 < m < 3.2 GeV
        if(fM > 3.0 && fM < 3.2) h->Fill(fPt);

        if((iEntry+1) % 1000 == 0){
            nEntriesAnalysed += 1000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    return;
}

void NewCutZ_PrepareTree(Double_t fCutZ)
// to compare the difference in surviving number of events after applying all the cuts and setting 
// the cut on Z_vtx to either 10 or 15 cm
// mass within 3.0 and 3.2 GeV, pT in defined bins
// see Guilermo's email from June 16, 2022
{
    TString name = "Trees/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut%.1f_DataTree.root", fCutZ);

    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("Tree already created.");
        return;

    } else { 

        Printf("Tree will be created.");

        // data
        TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
        if(f_in) Printf("Input data loaded.");

        TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
        if(t_in) Printf("Input tree loaded.");

        ConnectTreeVariables(t_in);  

        // Create new data tree with applied cuts
        file = new TFile(name.Data(),"RECREATE");

        TTree *tData = new TTree("tData", "tData");
        tData->Branch("fPt", &fPt, "fPt/D");
        tData->Branch("fM", &fM, "fM/D");

        Printf("%lli entries found in the tree.", t_in->GetEntries());
        Int_t nEntriesAnalysed = 0;

        // save the original value of cut_fVertexZ
        Printf("Original cut on vertex Z: %.1f", cut_fVertexZ);
        Double_t fCutZ_orig = cut_fVertexZ;
        // set the new value of cut_fVertexZ
        cut_fVertexZ = fCutZ;
        Printf("New cut on vertex Z: %.1f", cut_fVertexZ);

        for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
        {
            t_in->GetEntry(iEntry);
            // inv mass cut: 2.2 < m < 4.5, pT cut: all (pT < 2.0)
            if(EventPassed(0, 2)) tData->Fill();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // set back the original value of cut_fVertexZ
        cut_fVertexZ = fCutZ_orig;
        Printf("Restoring the original cut on vertex Z: %.1f", cut_fVertexZ);        

        file->Write("",TObject::kWriteDelete);

        return;
    }
}

void NewCutZ_AxE_PtBins(Double_t fCutZ)
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut%.1f_AxE_PtBins/", fCutZ));
    AxE_PtBins_Calculate(fCutZ);

    return;
}

void ConnectTreeVariables_tData(TTree *t)
{
    t->SetBranchAddress("fPt", &fPt);
    t->SetBranchAddress("fM", &fM);

    Printf("Variables from %s connected.", t->GetName());

    return;
}

Double_t CalculateErrorBinomial(Double_t k, Double_t n)
{
    Double_t var = k * (n - k) / n / n / n;
    Double_t err = TMath::Sqrt(var);

    return err;
}