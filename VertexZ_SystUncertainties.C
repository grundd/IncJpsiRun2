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

    if(debug)
    {
        TCanvas *cx = new TCanvas("cx","cx",900,600);
        hEv10->Draw("E0");
        TCanvas *cy = new TCanvas("cy","cy",900,600);
        hEv15->Draw("E0");
        TCanvas *cz = new TCanvas("cz","cz",900,600);
        hEvRat->Draw("E0");
    }
    // Calculate the ratios of yields manually using various methods to compute errors
    Double_t nEvRat_val[5] = { 0 };
    Double_t nEvRat_err1[5] = { 0 };
    Double_t nEvRat_err2[5] = { 0 };
    Double_t nEvRat_err3[5] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        // Calculate the ratio
        nEvRat_val[iBin] = hEv10->GetBinContent(iBin+1) / hEv15->GetBinContent(iBin+1);
        // Calculate the error using error propagation formula
        nEvRat_err1[iBin] = nEvRat_val[iBin] * TMath::Sqrt(TMath::Power(hEv10->GetBinError(iBin+1) / hEv10->GetBinContent(iBin+1), 2) 
            + TMath::Power(hEv15->GetBinError(iBin+1) / hEv15->GetBinContent(iBin+1), 2));
        // Calculate the error using binomial distribution
        // N = number of events with Z_vtx < 15 cm in a given bin
        // p = we estimate it as N(Z_vtx < 15 cm) / N(Z_vtx < 10 cm)
        // q = 1 - p
        Double_t N = hEv15->GetBinContent(iBin+1);
        Double_t p = hEv10->GetBinContent(iBin+1) / hEv15->GetBinContent(iBin+1);
        Double_t q = 1 - p;
        Double_t variance = N * p * q;
        Double_t sigma = TMath::Sqrt(variance); // this is the error of the mean (Np, hEv10->GetBinContent(iBin+1))
        // to get the error of the ratio, we need to divide the previous value by N again
        nEvRat_err2[iBin] = sigma / N;
        // we would get the same result by using the function CalculateErrorBinomial(.,.), see the paper by Ullrich
        // Calculate the error using Bayesian formula
        nEvRat_err3[iBin] = CalculateErrorBayes(hEv10->GetBinContent(iBin+1), hEv15->GetBinContent(iBin+1));
    }
    // Print the results
    if(debug)
    {
        Printf("***");
        for(Int_t iBin = 0; iBin < nPtBins; iBin++)
        {
            Printf("BIN %i:", iBin+1);
            Printf("nEv15: %.0f pm %.3f (ROOT) pm %.3f (Poisson)", hEv15->GetBinContent(iBin+1), hEv15->GetBinError(iBin+1), TMath::Sqrt(hEv15->GetBinContent(iBin+1)));
            Printf("nEv10: %.0f pm %.3f (ROOT) pm %.3f (Poisson)", hEv10->GetBinContent(iBin+1), hEv10->GetBinError(iBin+1), TMath::Sqrt(hEv10->GetBinContent(iBin+1)));
            Printf("ratio: %.4f", nEvRat_val[iBin]);
            // error calculated as sqrt(bin content):
            Printf("err: %.4f (ROOT w/o Sumw2 = Poisson)", hEvRat_NoSumw2->GetBinError(iBin+1));
            // error calculated as sqrt(sum of squares of weights)
            // https://www-zeuthen.desy.de/~wischnew/amanda/discussion/wgterror/working.html 
            Printf("err: %.4f (ROOT Sumw2)", hEvRat->GetBinError(iBin+1));
            // error calculated from the error propagation formula:
            Printf("err: %.4f (err propagation formula)", nEvRat_err1[iBin]);
            // error calculated from a binomial distribution
            Printf("err: %.4f (binomial)", nEvRat_err2[iBin]);
            // error calculated using Bayesian formula
            Printf("err: %.4f (Bayes)", nEvRat_err3[iBin]);
            Printf("***");
        }
    }
    // Calculate the values of AxE with the cut on Z_vtx < 15 cm and < 10 cm
    NewCutZ_AxE_PtBins(15.0);
    NewCutZ_AxE_PtBins(10.0);
    // To calculate the ratios of AxE: we will compare NRec (NGen are the same in all bins)
    // the errors of the ratios will again be calculated from binomial distribution
    Double_t nNRec15[5] = { 0 };
    Double_t nNRec10[5] = { 0 };
    TString str_AxE15 = "";
    if(cut_fVertexZ == 15.0) str_AxE15 = "Results/" + str_subfolder + Form("AxE_PtBins/NRec_%ibins.txt", nPtBins);
    else                     str_AxE15 = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut15.0_AxE_PtBins/NRec_%ibins.txt", nPtBins);
    TString str_AxE10 = "";
    if(cut_fVertexZ == 10.0) str_AxE10 = "Results/" + str_subfolder + Form("AxE_PtBins/NRec_%ibins.txt", nPtBins);
    else                     str_AxE10 = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut10.0_AxE_PtBins/NRec_%ibins.txt", nPtBins);
    // Load the values of NRec
    ifstream ifs;
    Int_t bin;
    ifs.open(str_AxE15);    
    for(Int_t i = 0; i < nPtBins; i++)
    {
        ifs >> bin >> nNRec15[i];
    }
    Printf("Values of AxE for Z_vtx < 15 cm calculated.");
    ifs.close();
    ifs.open(str_AxE10);    
    for(Int_t i = 0; i < nPtBins; i++)
    {
        ifs >> bin >> nNRec10[i];
    }
    Printf("Values of AxE for Z_vtx < 10 cm calculated.");
    ifs.close();
    // Calculate the ratios of yields manually using various methods to compute errors
    Double_t nNRecRat_val[5] = { 0 };
    Double_t nNRecRat_err1[5] = { 0 };
    Double_t nNRecRat_err2[5] = { 0 };
    Double_t nNRecRat_err3[5] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        // Calculate the ratio
        nNRecRat_val[iBin] = nNRec10[iBin] / nNRec15[iBin];
        // Calculate the error using error propagation formula
        nNRecRat_err1[iBin] = nNRecRat_val[iBin] * TMath::Sqrt(TMath::Power(TMath::Sqrt(nNRec10[iBin]) / nNRec10[iBin], 2) 
            + TMath::Power(TMath::Sqrt(nNRec15[iBin]) / nNRec15[iBin], 2));
        // Calculate the error using binomial distribution
        nNRecRat_err2[iBin] = CalculateErrorBinomial(nNRec10[iBin], nNRec15[iBin]);
        // Calculate the error using Bayesian formula
        nNRecRat_err3[iBin] = CalculateErrorBayes(nNRec10[iBin], nNRec15[iBin]);
    }
    // Print the results
    if(debug)
    {
        Printf("***");
        for(Int_t iBin = 0; iBin < nPtBins; iBin++)
        {
            Printf("BIN %i:", iBin+1);
            Printf("nNRec15: %.0f pm %.0f (Poisson)", nNRec15[iBin], TMath::Sqrt(nNRec15[iBin]));
            Printf("nNRec10: %.0f pm %.0f (Poisson)", nNRec10[iBin], TMath::Sqrt(nNRec10[iBin]));
            Printf("ratio: %.4f", nNRecRat_val[iBin]);
            // error calculated from the error propagation formula:
            Printf("err: %.4f (err propagation formula)", nNRecRat_err1[iBin]);
            // error calculated from a binomial distribution
            Printf("err: %.4f (binomial)", nNRecRat_err2[iBin]);
            // error calculated using Bayesian formula
            Printf("err: %.4f (Bayes)", nNRecRat_err3[iBin]);
            Printf("***");
        }
    }
    // Calculate the systematic uncertainty
    Double_t syst_uncr[5] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++)
    {
        syst_uncr[i] = (1 - nEvRat_val[i] / nNRecRat_val[i]) * 100;
    }
    TString str_out = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/syst_uncertainties_%ibins.txt", nPtBins);
    ofstream outfile(str_out.Data());
    //outfile << "bin\tsyst\n";
    outfile << std::fixed << std::setprecision(1);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        outfile << i+1 << "\t" << TMath::Abs(syst_uncr[i]) << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());

    // Plot the differences
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
        gr[i] = new TGraphErrors(1,&nNRecRat_val[i],&nEvRat_val[i],&nNRecRat_err3[i],&nEvRat_err3[i]);
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
        lg[i]->AddEntry((TObject*)0,Form("R_{N}^{10/15} = (%.1f pm %.1f)%%", nEvRat_val[i]*100, nEvRat_err3[i]*100),"");
        lg[i]->AddEntry((TObject*)0,Form("R_{A#times#varepsilon}^{10/15} = (%.1f pm %.1f)%%", nNRecRat_val[i]*100, nNRecRat_err3[i]*100),"");
        lg[i]->AddEntry((TObject*)0,Form("change = %.1f%%", syst_uncr[i]),"");
        lg[i]->AddEntry((TObject*)0,Form("syst. uncr. = %.1f%%", TMath::Abs(syst_uncr[i])),"");
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