// _STARlight_NewPtShapes.C
// David Grund, Mar 12, 2022

// cpp headers
#include <fstream>
#include <stdio.h>
// root headers
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h" // to be able to use SetReducedRunList()

Bool_t drawCheck(kFALSE);

TString strMCArr[4] = {"CohJ","IncJ","CohP","IncP"};
Double_t fPtCutLowArr[4] = {0.0, 0.0, 0.0, 0.0}; // GeV/c
Double_t fPtCutUppArr[4] = {0.4, 1.2, 0.0, 0.0}; // GeV/c
Int_t nBinsArr[4] = {80, 240, 0, 0};

TString strMC;
Double_t fPtCutLow;
Double_t fPtCutUpp;
Int_t nBins;

Double_t fPtGenerated;
TLorentzVector *parent;
TClonesArray *daughters;

TH1D *hRecOld = NULL;
TH1D *hRatios = NULL;
TH1D *hRecNew = NULL;

void InitAnalysis(Int_t iMC, Bool_t pass3);
void FillTreeGen(const char* folder_in, Double_t R_A);
void FillHistRec();
void CalcAndPlotRatios(const char* subfolder_out, Double_t R_A);
void ConnectTreeVariables_tSL(TTree *tSL);


void _STARlight_NewPtShapes()
{
    // IncJ, R_A = 6.624 vs. 7.350 fm, 6.000.000 gen events
    InitAnalysis(1,kTRUE);
    FillTreeGen("Trees/STARlight/IncJ_6.624/",6.624);
    FillTreeGen("Trees/STARlight/IncJ_7.350/",7.350);
    FillHistRec();
    CalcAndPlotRatios("",7.350);

    return;
}

// #############################################################################################

void InitAnalysis(Int_t iMC, Bool_t pass3)
{
    // iMC == 0 => CohJ
    // iMC == 1 => IncJ
    // iMC == 2 => CohP
    // iMC == 3 => IncP

    isPass3 = pass3;
    SetReducedRunList(isPass3);
    strMC = strMCArr[iMC];
    fPtCutLow = fPtCutLowArr[iMC];
    fPtCutUpp = fPtCutUppArr[iMC];
    nBins = nBinsArr[iMC];
    hRecOld = new TH1D("hRecOld","hRecOld",nBins,fPtCutLow,fPtCutUpp);

    return;
}

// #############################################################################################

void FillTreeGen(const char* folder_in, Double_t R_A)
{
    Printf("*****");
    Printf("Process: %s", strMC.Data());
    Printf("Input folder: %s", folder_in);
    Printf("Filling the tree and histogram (generator level).");
    Printf("*****");

    Double_t nEvOld = 0;
    Double_t nEvNew = 0;

    // open the starlight file and starlight tree
    TFile *fSL = TFile::Open(Form("%stree_STARlight.root", folder_in), "read");
    if(fSL) Printf("File %s loaded.", fSL->GetName());

    // get the SL tree
    TTree *tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
    if(tSL) Printf("Tree %s loaded.", tSL->GetName());
    ConnectTreeVariables_tSL(tSL);

    // check if the output trees already created, if not, create them
    TString str_f_out = Form("Trees/STARlight/tGen_%s_RA_%.3f.root", strMC.Data(), R_A);
    TFile *fGen = TFile::Open(str_f_out.Data(),"read");

    if(fGen){

        Printf("File %s already created.", str_f_out.Data());

    } else {

        TH1D *hGen = new TH1D("hGen","hGen",nBins,fPtCutLow,fPtCutUpp);

        gROOT->cd();
        TTree *tGen = new TTree("tGen","tGen");
        tGen->Branch("fPtGen",&fPtGenerated,"fPtGen/D");

        Printf("Filling tGen and hGen.");
        Printf("tSL contains %lli entries.", tSL->GetEntries());

        // Loop over entries in tSL
        Int_t nEntriesAnalysed = 0;
        Int_t nEntriesProgress = (Double_t)tSL->GetEntries() / 20.;
        Int_t nPercent = 0;

        for(Int_t iEntry = 0; iEntry < tSL->GetEntries(); iEntry++){
            tSL->GetEntry(iEntry);
            if(TMath::Abs(fYGen) < 1.0){
                nEvOld++;
                fPtGenerated = parent->Pt();
                hGen->Fill(fPtGenerated);
                tGen->Fill();
            }
            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }
        Printf("No. of events with |y| < 1.0: %.0f", nEvOld);

        fGen = new TFile(str_f_out.Data(),"RECREATE");
        // open the file
        fGen->cd();        
        // write the histogram and tree to this directory
        hGen->Write("hGen", TObject::kSingleKey);
        tGen->Write("tGen",TObject::kSingleKey);
        // list the contents of the file
        fGen->ls();
        // close the file
        fGen->Close();
    }

    Printf("*****");
    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return;
}

// #############################################################################################

void CalcAndPlotRatios(const char* subfolder_out, Double_t R_A)
{
    Printf("*****");
    Printf("Process: %s", strMC.Data());
    Printf("Calculating ratios of nGenNew/nGenOld.");
    Printf("*****");

    // Load tGen with R_A = 6.624
    TString str_fGenOld = Form("Trees/STARlight/tGen_%s_RA_6.624.root", strMC.Data());
    TFile *fGenOld = fGenOld = TFile::Open(str_fGenOld.Data(), "read");
    if(fGenOld) Printf("File %s loaded.", fGenOld->GetName());

    TTree *tGenOld = (TTree*) fGenOld->Get("tGen");
    if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName()); 
    tGenOld->SetBranchAddress("fPtGen", &fPtGenerated);

    // Load tGen with new R_A
    TString str_fGenNew = Form("Trees/STARlight/tGen_%s_RA_%.3f.root", strMC.Data(), R_A);
    TFile *fGenNew = TFile::Open(str_fGenNew.Data(),"read");
    if(fGenNew) Printf("File %s loaded.", fGenNew->GetName());

    TTree *tGenNew = (TTree*) fGenNew->Get("tGen");
    if(tGenNew) Printf("Tree %s loaded.", tGenNew->GetName()); 
    tGenNew->SetBranchAddress("fPtGen", &fPtGenerated);

    // Load histograms with generated events
    TH1D *hGenOld = (TH1D*)fGenOld->Get("hGen");
    if(hGenOld) Printf("Histogram %sOld loaded.", hGenOld->GetName());
    // draw hGenOld
    if(drawCheck){
        TCanvas *c1 = new TCanvas("c1","c1",900,600);
        hGenOld->Draw();
    } 
    TH1D *hGenNew = (TH1D*)fGenNew->Get("hGen");
    if(hGenNew) Printf("Histogram %sNew loaded.", hGenNew->GetName());
    // draw hGenNew
    if(drawCheck){
        TCanvas *c2 = new TCanvas("c2","c2",900,600);
        hGenNew->Draw();
    } 
    hRatios = (TH1D*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("hRatios");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);
    // draw hRatios
    if(drawCheck){
        TCanvas *c3 = new TCanvas("c3","c3",900,600);
        hRatios->Draw();
    } 
    // draw hRecOld
    if(drawCheck){
        TCanvas *c4 = new TCanvas("c4","c4",900,600);
        hRecOld->Draw();
    }    
    hRecNew = (TH1D*)hRecOld->Clone("hRecNew");
    hRecNew->SetTitle("hRecNew");
    hRecNew->Multiply(hRatios);
    // draw hRecNew
    if(drawCheck){
        TCanvas *c5 = new TCanvas("c5","c5",900,600);
        hRecNew->Draw();
    } 

    // Print the results to the text file
    TString str_folder_out = Form("Results/_STARlight_NewPtShapes/%s/", subfolder_out);
    gSystem->Exec("mkdir -p " + str_folder_out);
    TString str_file_out = Form("%s%s_RA_%.3f.txt", str_folder_out.Data(), strMC.Data(), R_A);
    ofstream outfile(str_file_out.Data());
    outfile << "pT_low\tpT_upp\tnEvOld\tnEvNew\tratio\n";
    for(Int_t iBin = 1; iBin <= nBins; iBin++){
        outfile << Form("%.3f\t%.3f\t%.0f\t%.0f\t%.3f\n",
            hRatios->GetBinLowEdge(iBin), hRatios->GetBinLowEdge(iBin+1), 
            hGenOld->GetBinContent(iBin), hGenNew->GetBinContent(iBin), hRatios->GetBinContent(iBin));
    }
    outfile.close();

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    // Plot the results
    TCanvas *cRatios = new TCanvas("cRatios", "cRatios", 900, 600);
    // Canvas settings
    cRatios->SetTopMargin(0.04);
    cRatios->SetBottomMargin(0.15);
    cRatios->SetRightMargin(0.04);
    cRatios->SetLeftMargin(0.1);
    // Vertical axis
    hRatios->GetYaxis()->SetTitle(Form("#it{N}_{gen}(#it{R}_{A} = %.3f fm) / #it{N}_{gen}(#it{R}_{A} = 6.624 fm)", R_A));
    hRatios->GetYaxis()->SetTitleSize(0.05);
    hRatios->GetYaxis()->SetTitleOffset(0.9);
    hRatios->GetYaxis()->SetLabelSize(0.05);
    hRatios->GetYaxis()->SetDecimals(1);
    hRatios->GetYaxis()->SetRangeUser(0.0,5.0);
    // Horizontal axis
    hRatios->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRatios->GetXaxis()->SetTitleSize(0.05);
    hRatios->GetXaxis()->SetTitleOffset(1.2);
    hRatios->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRatios->SetLineColor(215);
    hRatios->Draw("E1");

    TCanvas *cRec = new TCanvas("cRec", "cRec", 900, 600);
    cRec->SetTopMargin(0.04);
    cRec->SetBottomMargin(0.15);
    cRec->SetRightMargin(0.04);
    cRec->SetLeftMargin(0.1);
    cRec->SetLogy();
    // Vertical axis
    hRecOld->GetYaxis()->SetTitle("#it{N}_{rec} (after selections)");
    hRecOld->GetYaxis()->SetTitleSize(0.05);
    hRecOld->GetYaxis()->SetTitleOffset(0.9);
    hRecOld->GetYaxis()->SetLabelSize(0.05);
    hRecOld->GetYaxis()->SetRangeUser(0.1,1e5);
    // Horizontal axis
    hRecOld->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRecOld->GetXaxis()->SetTitleSize(0.05);
    hRecOld->GetXaxis()->SetTitleOffset(1.2);
    hRecOld->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRecOld->SetLineColor(210);
    hRecOld->Draw("E1");
    hRecNew->SetLineColor(215);
    hRecNew->Draw("E1 SAME");

    // Print plots
    cRatios->Print(Form("%s%s_RA_%.3f_ratios.pdf", str_folder_out.Data(), strMC.Data(), R_A));
    cRatios->Print(Form("%s%s_RA_%.3f_ratios.png", str_folder_out.Data(), strMC.Data(), R_A));
    cRec->Print(Form("%s%s_RA_%.3f_recSpectra.pdf", str_folder_out.Data(), strMC.Data(), R_A));
    cRec->Print(Form("%s%s_RA_%.3f_recSpectra.png", str_folder_out.Data(), strMC.Data(), R_A));

    return;
}

// #############################################################################################

void FillHistRec()
{
    Printf("*****");
    Printf("Process: %s", strMC.Data());
    Printf("Filling the histogram (reconstructed level).");
    Printf("*****");

    TString str_f_in = "Trees/";
    if(!isPass3) str_f_in += "AnalysisDataMC_pass1/";
    else         str_f_in += "AnalysisDataMC_pass3/";
    // choose MC dataset
    // https://www.cplusplus.com/reference/cstring/strncmp/
    if     (strncmp(strMC.Data(),"CohJ",4) == 0) str_f_in += "AnalysisResults_MC_kCohJpsiToMu.root";
    else if(strncmp(strMC.Data(),"IncJ",4) == 0) str_f_in += "AnalysisResults_MC_kIncohJpsiToMu.root";
    else if(strncmp(strMC.Data(),"CohP",4) == 0) str_f_in += "AnalysisResults_MC_kCohPsi2sToMuPi.root";
    else if(strncmp(strMC.Data(),"IncP",4) == 0) str_f_in += "AnalysisResults_MC_kIncohPsi2sToMuPi.root";
    // open the input file
    TFile *fSL = NULL; 
    fSL = TFile::Open(str_f_in.Data(), "read");
    if(fSL) Printf("File %s loaded.", fSL->GetName());
    // get the MCRec tree
    TString str_t_in = "";
    if(!isPass3) str_t_in = "AnalysisOutput/fTreeJPsiMCRec";
    else         str_t_in = "AnalysisOutput/fTreeJpsi";
    TTree *tRec = dynamic_cast<TTree*> (fSL->Get(str_t_in.Data()));
    if(tRec) Printf("Tree %s loaded.", tRec->GetName());
    // connect tree varibles, first set if pass3
    ConnectTreeVariablesMCRec(tRec);

    Int_t nEntriesAnalysed = 0;
    Int_t nEntriesProgress = (Double_t)tRec->GetEntries() / 20.;
    Int_t nPercent = 0;

    // run over reconstructed events and fill the histogram hRecOld
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        // m between 3.0 and 3.2 GeV/c^2, pT cut: all
        if(EventPassedMCRec(1, 2)) hRecOld->Fill(fPt);

        if((iEntry+1) % nEntriesProgress == 0){
            nPercent += 5;
            nEntriesAnalysed += nEntriesProgress;
            Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
        }
    } 

    if(drawCheck){
        TCanvas *c0 = new TCanvas("c0","c0",900,600);
        hRecOld->Draw();
    }

    Printf("*****");
    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return;
}

// #############################################################################################

void ConnectTreeVariables_tSL(TTree *tSL){

    tSL->SetBranchAddress("parent", &parent);
    tSL->SetBranchAddress("daughters", &daughters);

    Printf("Variables from %s connected.", tSL->GetName());

    return;
}

// #############################################################################################