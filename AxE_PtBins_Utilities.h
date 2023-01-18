// AxE_PtBins_Utilities.h
// David Grund, Jun 07, 2022

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"

TH1F* hRec = NULL; 
TH1F* hGen = NULL; 
TH1F* hAxE = NULL;
Float_t NRec_tot_val, NRec_tot_err;
Float_t NGen_tot_val, NGen_tot_err;

void AxE_PtBins_SaveToFile(TString path, Float_t N_tot_val, Float_t N_tot_err, TH1F* h, Int_t prec = 0, Float_t fact = 1.) 
{
    ofstream of(path.Data());
    of << std::fixed << std::setprecision(prec)
       << "0\t" << fact * N_tot_val 
       << "\t" << fact * N_tot_err 
       << "\n";
    for(Int_t iBin = 1; iBin <= h->GetNbinsX(); iBin++) {
        of << iBin << "\t" << fact * h->GetBinContent(iBin) 
           << "\t" << fact * h->GetBinError(iBin)
           << "\n";
    }
    of.close();
    Printf("*** File saved in %s.***", path.Data());
}

Float_t CalculateErrorBayes(Float_t k, Float_t n){ // k = NRec, n = NGen

    Float_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Float_t err = TMath::Sqrt(var);

    return err;
}

void AxE_PtBins_FillHistNRec(Float_t fCutZ)
{
    // check if the corresponding text file already exists
    TString file;
    if(fCutZ == cut_fVertexZ) file = "Results/" + str_subfolder + "AxE_PtBins/";
    else                      file = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut%.1f_AxE_PtBins/", fCutZ);
    file.Append(Form("NRec_%ibins.txt", nPtBins));

    ifstream ifs;
    ifs.open(file);
    if(!ifs.fail())
    {
        // this configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // fill hRec with data from the file
        Int_t i_bin; Float_t NRec_val, NRec_err;
        // fiducial
        ifs >> i_bin >> NRec_tot_val >> NRec_tot_err;
        // pT bins
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++) {
            ifs >> i_bin >> NRec_val >> NRec_err;
            hRec->SetBinContent(iBin, NRec_val);
        }
        ifs.close(); 
        return;
    } 
    else 
    {
        // this configuration is yet to be calculated
        Printf("*** Calculating N rec per bin for %s... ***", file.Data());
        TFile *fRec = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fRec) Printf("MC rec file loaded.");
        TTree *tRec = dynamic_cast<TTree*> (fRec->Get(str_in_MC_tree_rec.Data()));
        if(tRec) Printf("MC rec tree loaded.");
        ConnectTreeVariablesMCRec(tRec);

        // |> *********** for VertexZ_SystUncertainties.C ***********
        // save the original value of cut_fVertexZ
        Printf("Original cut on vertex Z: %.1f", cut_fVertexZ);
        Float_t fCutZ_orig = cut_fVertexZ;
        if(fCutZ != cut_fVertexZ) {
            // set the new value of cut_fVertexZ
            cut_fVertexZ = fCutZ;
            Printf("New cut on vertex Z: %.1f", cut_fVertexZ);
        }
        // <| *****************************************************

        // go over tree entries and calculate NRec in the total range and in bins
        Int_t NRec[6] = { 0 };
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) {
            tRec->GetEntry(iEntry);
            // m between 2.2 and 4.5 GeV/c^2
            // & pT from 0.2 to 1.0 GeV/c
            if(EventPassedMCRec(0, 3)) NRec[0]++;
            // & pT in range of a given bin 
            for(Int_t iBin = 1; iBin <= nPtBins; iBin++) if(EventPassedMCRec(0, 4, iBin)) NRec[iBin]++;
        }
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++) hRec->SetBinContent(iBin, NRec[iBin]);
        Printf("*** Finished! ***");

        // |> *********** for VertexZ_SystUncertainties.C ***********
        if(cut_fVertexZ != fCutZ_orig)
        {
            // set back the original value of cut_fVertexZ
            cut_fVertexZ = fCutZ_orig;
            Printf("Restoring the original cut on vertex Z: %.1f", cut_fVertexZ);  
        }
        // <| *****************************************************
        
        NRec_tot_val = NRec[0];
        NRec_tot_err = TMath::Sqrt(NRec[0]);
        AxE_PtBins_SaveToFile(file,NRec[0],TMath::Sqrt(NRec[0]),hRec);
        return;
    }
}

void AxE_PtBins_FillHistNGen()
{
    // check if the corresponding text file already exists
    TString file = "Results/" + str_subfolder + "AxE_PtBins/";
    file.Append(Form("NGen_%ibins.txt", nPtBins));

    ifstream ifs;
    ifs.open(file);
    if(!ifs.fail())
    {
        // this configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // fill hGen with data from the text file
        Int_t i_bin; Float_t NGen_val, NGen_err;
        // fiducial
        ifs >> i_bin >> NGen_tot_val >> NGen_tot_err;
        // pT bins
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++) {
            ifs >> i_bin >> NGen_val >> NGen_err;
            hGen->SetBinContent(iBin, NGen_val);
        }
        ifs.close(); 
        return;
    } 
    else 
    {
        // this configuration is yet to be calculated
        Printf("*** Calculating N gen per bin for %s... ***", file.Data());
        TFile *fGen = TFile::Open((str_in_MC_fldr_gen + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fGen) Printf("MC gen file loaded.");
        TTree *tGen = dynamic_cast<TTree*> (fGen->Get(str_in_MC_tree_gen.Data()));
        if(tGen) Printf("MC gen tree loaded.");
        ConnectTreeVariablesMCGen(tGen);

        // go over tree entries and calculate NRec in the total range and in bins
        Int_t NGen[6] = { 0 };
        for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++) {
            tGen->GetEntry(iEntry);
            // m between 2.2 and 4.5 GeV/c^2
            // & pT from 0.2 to 1.0 GeV/c
            if(EventPassedMCGen(3)) NGen[0]++;
            // & pT in range of a given bin 
            for(Int_t iBin = 1; iBin <= nPtBins; iBin++) if(EventPassedMCGen(4, iBin)) NGen[iBin]++;
        }
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++) hGen->SetBinContent(iBin, NGen[iBin]);
        Printf("*** Finished! ***");
        
        NGen_tot_val = NGen[0];
        NGen_tot_err = TMath::Sqrt(NGen[0]);
        AxE_PtBins_SaveToFile(file,NGen[0],TMath::Sqrt(NGen[0]),hGen);
        return;
    }
    return;
}

void AxE_PtBins_Calculate(Float_t fCutZ)
{
    hRec = new TH1F("hRec","N_{rec} per bin",nPtBins,ptBoundaries);
    hGen = new TH1F("hGen","N_{gen} per bin",nPtBins,ptBoundaries);

    AxE_PtBins_FillHistNRec(fCutZ);
    AxE_PtBins_FillHistNGen();

    hAxE = (TH1F*)hRec->Clone("hAxE");
    hAxE->SetTitle("AxE per bin");
    hAxE->Sumw2();
    hAxE->Divide(hGen);

    // Draw the histogram:
    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.145);
    // gStyle
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");
    // Marker and line
    hAxE->SetLineColor(kBlue);
    hAxE->SetLineWidth(1.0);
    // Vertical axis
    hAxE->GetYaxis()->SetTitle("#it{N}_{rec}^{MC}/#it{N}_{gen}^{MC}");
    hAxE->GetYaxis()->SetTitleSize(0.056);
    hAxE->GetYaxis()->SetTitleOffset(1.3);
    hAxE->GetYaxis()->SetLabelSize(0.056);
    hAxE->GetYaxis()->SetDecimals(3);
    hAxE->GetYaxis()->SetRangeUser(0.0,hAxE->GetBinContent(1)*1.1);
    // Horizontal axis
    hAxE->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAxE->GetXaxis()->SetTitleSize(0.056);
    hAxE->GetXaxis()->SetTitleOffset(1.2);
    hAxE->GetXaxis()->SetLabelSize(0.056);
    hAxE->GetXaxis()->SetLabelOffset(0.015);
    hAxE->GetXaxis()->SetDecimals(1);
    // Eventually draw it
    hAxE->Draw("P E1");
    // legend
    TLegend *l = new TLegend(0.52,0.77,0.85,0.97);
    l->AddEntry((TObject*)0,Form("ALICE Simulation"),"");
    l->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l->SetTextSize(0.056);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();
    // legend 2
    TLegend *l2 = new TLegend(0.15,0.17,0.35,0.32);
    l2->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),"");
    l2->AddEntry((TObject*)0,Form("2.2 < #it{m} < 4.5 GeV/#it{c}^{2}"),"");
    l2->SetTextSize(0.056);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    // save the figures and print the results to txt file
    TString str;
    if(fCutZ == cut_fVertexZ) str = "Results/" + str_subfolder + Form("AxE_PtBins/AxE_%ibins", nPtBins);
    else                      str = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut%.1f_AxE_PtBins/AxE_%ibins", fCutZ, nPtBins);
    c->Print((str + ".pdf").Data());

    // compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kFALSE;
    if(DebugErrors){
        Float_t ErrRoot = 0;
        Float_t ErrBayes = 0;    
        for(Int_t i = 1; i <= nPtBins; i++) {
            ErrRoot = hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(hRec->GetBinContent(i),hGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // calculate the total value of AxE
    Float_t AxE_tot_val = NRec_tot_val / NGen_tot_val;
    Float_t AxE_tot_err = CalculateErrorBayes(NRec_tot_val, NGen_tot_val);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxE_tot_val*100., AxE_tot_err*100.);

    // print the results to a text file
    AxE_PtBins_SaveToFile(str + ".txt",AxE_tot_val,AxE_tot_err,hAxE,3,100.);
    return;
}