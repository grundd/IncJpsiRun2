// _STARlight_NewPtShapes.C
// David Grund, Mar 07, 2022

// root headers
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
// my headers
//#include "AnalysisManager.h"
//#include "STARlight_Utilities.h"

Double_t fPtGenerated;
Double_t fPtCutLow_CohJ = 0.0; // GeV/c
Double_t fPtCutUpp_CohJ = 0.4; // GeV/c
Int_t nBins_CohJ = 80;

TH1D *hRecOld = new TH1D("hRecOld","hRecOld",nBins_CohJ,fPtCutLow_CohJ,fPtCutUpp_CohJ);
TH1D *hRatios = NULL;
TH1D *hRecNew = NULL;

void FillHistRec();
void FillTreeGen(Bool_t b6000k, TString sIn, TString sOut_subfolder, Double_t R_A);
void CalculateAndPlotRatios(Bool_t b6000k, TString sOut_subfolder, Double_t R_A);

TString sIn = "";
TString sOut_subfolder = "";

void _STARlight_NewPtShapes()
{
    Bool_t b4800k = kTRUE;
    if(b4800k){
        sIn = "";
        sOut_subfolder = "CoherentShape";
        FillHistRec();
        FillTreeGen(kFALSE,sIn,sOut_subfolder,7.530);
        CalculateAndPlotRatios(kFALSE,sOut_subfolder,7.530);
    }

    Bool_t b6000k = kTRUE;
    if(b6000k){
        sIn = "coh_6000000_modRA";
        sOut_subfolder = "CoherentShape";
        FillHistRec();
        FillTreeGen(kTRUE,sIn,sOut_subfolder,7.530);
        CalculateAndPlotRatios(kTRUE,sOut_subfolder,7.530);
    }

    Bool_t bOptimalRA = kTRUE;
    if(bOptimalRA){
        FillHistRec();

        sIn = "OptimalRA/coh_modRA_0_6.624";
        sOut_subfolder = "OptimalRA";
        FillTreeGen(kTRUE,sIn,sOut_subfolder,6.624);
        CalculateAndPlotRatios(kTRUE,sOut_subfolder,6.624);

        Int_t iUpTo = 13;
        Double_t R_A = 0.;
        for(Int_t i = 0; i < iUpTo; i++){
            R_A = 6.60 + i * 0.10;  
            sIn = Form("OptimalRA/coh_modRA_%i_%.3f", i+1, R_A);
            sOut_subfolder = "OptimalRA";
            FillTreeGen(kTRUE,sIn,sOut_subfolder,R_A);
            CalculateAndPlotRatios(kTRUE,sOut_subfolder,R_A);
        }
    }

    return;
}

void FillTreeGen(Bool_t b6000k, TString sIn, TString sOut_subfolder, Double_t R_A)
{

    Printf("*****");
    if(b6000k)  Printf("Filling trees with nGen = 6.000.000.");
    else        Printf("Filling trees with nGen = 4.880.000.");
    Printf("*****");

    Double_t nEvOld = 0;
    Double_t nEvNew = 0;
    
    // #############################################################################################
    
    TFile *fSLOld = NULL;
    TTree *tSLOld = NULL;
    if(b6000k){
        // NGen = 6.000.000
        // Load generated coherent MC events with R_A = 6.624 fm
        TFile *fSLOld = TFile::Open("Trees/STARlight/coh_6000000_stdRA/trees_starlight.root", "read");
        if(fSLOld) Printf("File %s loaded.", fSLOld->GetName());
        // Get the SL tree
        tSLOld = dynamic_cast<TTree*> (fSLOld->Get("starlightTree"));
        if(tSLOld) Printf("Tree %s loaded.", tSLOld->GetName());
        ConnectTreeVariables_tSL(tSLOld);
    } else {
        // NGen = 4.880.000
        // Load official root file: kCohJpsiToMu (R_A = 6.624 fm)
        fSLOld = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
        if(fSLOld) Printf("File %s loaded.", fSLOld->GetName());
        // Get the MCGen tree
        tSLOld = dynamic_cast<TTree*> (fSLOld->Get("AnalysisOutput/fTreeJPsiMCGen"));
        if(tSLOld) Printf("Tree %s loaded.", tSLOld->GetName());
        ConnectTreeVariablesMCGen(tSLOld);
    }
    TFile *fSLNew = NULL;
    TTree *tSLNew = NULL;
    if(b6000k){
        // NGen = 6.000.000
        // Load generated coherent MC events with new R_A
        fSLNew = TFile::Open(Form("Trees/STARlight/%s/trees_starlight.root", sIn.Data()), "read");
        if(fSLNew) Printf("File %s loaded.", fSLNew->GetName());
    } else {
        // NGen = 4.880.000
        // Load generated coherent MC events (with R_A = 7.53 fm)
        fSLNew = TFile::Open("Trees/STARlight/coh_4880000_modRA/trees_starlight.root", "read");
        if(fSLNew) Printf("File %s loaded.", fSLNew->GetName());
    }
    // Get the SL tree
    tSLNew = dynamic_cast<TTree*> (fSLNew->Get("starlightTree"));
    if(tSLNew) Printf("Tree %s loaded.", tSLNew->GetName());
    ConnectTreeVariables_tSL(tSLNew);

    // #############################################################################################
    // Check if output already created, if not, create it

    TString str_fGenOld = "";
    if(b6000k)  str_fGenOld = Form("Trees/STARlight/%s/tGenOld_6000k.root", sOut_subfolder.Data());
    else        str_fGenOld = Form("Trees/STARlight/%s/tGenOld_4880k.root", sOut_subfolder.Data());
    TFile *fGenOld = TFile::Open(str_fGenOld.Data(),"read");

    if(fGenOld){

        Printf("%s already created.", str_fGenOld.Data());

    } else {

        fGenOld = new TFile(str_fGenOld.Data(),"RECREATE");

        TH1D *hGenOld = new TH1D("hGenOld","hGenOld",nBins_CohJ,fPtCutLow_CohJ,fPtCutUpp_CohJ);

        TTree *tGenOld = new TTree("tGenOld","tGenOld");
        tGenOld->Branch("fPtGen",&fPtGenerated,"fPtGen/D");

        Printf("Filling tGenOld and hGenOld.");
        Printf("Old tree contains %lli entries.", tSLOld->GetEntries());

        // Loop over entries in tGen
        Int_t nEntriesAnalysed = 0;
        Int_t nEntriesProgress = (Double_t)tSLOld->GetEntries() / 20.;
        Int_t nPercent = 0;

        for(Int_t iEntry = 0; iEntry < tSLOld->GetEntries(); iEntry++){
            tSLOld->GetEntry(iEntry);
            if(TMath::Abs(fYGen) < 1.0){
                nEvOld++;
                if(b6000k)  fPtGenerated = parent->Pt();
                else        fPtGenerated = fPtGen;
                hGenOld->Fill(fPtGenerated);
                tGenOld->Fill();
            }
            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }
        Printf("No. of events with |y| < 1.0: %.0f", nEvOld);

        fGenOld->Write("",TObject::kWriteDelete);
    }

    TString str_fGenNew = "";
    if(b6000k)  str_fGenNew = Form("Trees/STARlight/%s/tGenNew_6000k_%.3f.root", sOut_subfolder.Data(), R_A);
    else        str_fGenNew = Form("Trees/STARlight/%s/tGenNew_4880k_%.3f.root", sOut_subfolder.Data(), R_A);
    TFile *fGenNew = TFile::Open(str_fGenNew.Data(),"read");

    if(fGenNew){

        Printf("%s already created.", str_fGenNew.Data());

    } else {

        fGenNew = new TFile(str_fGenNew.Data(),"RECREATE");

        TH1D *hGenNew = new TH1D("hGenNew","hGenNew",nBins_CohJ,fPtCutLow_CohJ,fPtCutUpp_CohJ);

        TTree *tGenNew = new TTree("tGenNew","tGenNew");
        tGenNew->Branch("fPtGen",&fPtGenerated,"fPtGen/D");

        Printf("Filling tGenNew and hGenNew.");
        Printf("New tree contains %lli entries.", tSLNew->GetEntries());

        // Loop over entries in new SL tree
        Int_t nEntriesAnalysed = 0;
        Int_t nEntriesProgress = (Double_t)tSLNew->GetEntries() / 20.;
        Int_t nPercent = 0;


        for(Int_t iEntry = 0; iEntry < tSLNew->GetEntries(); iEntry++){
            tSLNew->GetEntry(iEntry);
            if(TMath::Abs(parent->Rapidity()) < 1.0){
                nEvNew++;
                fPtGenerated = parent->Pt();
                hGenNew->Fill(fPtGenerated);
                tGenNew->Fill();
            }
            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }
        Printf("No. of events with |y| < 1.0: %.0f", nEvNew);

        fGenNew->Write("",TObject::kWriteDelete);
    }

    // #############################################################################################

    Printf("*****");
    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return;
}

void CalculateAndPlotRatios(Bool_t b6000k, TString sOut_subfolder, Double_t R_A)
{
    Printf("*****");
    if(b6000k)  Printf("Calculating ratios nGen = 6.000.000.");
    else        Printf("Calculating ratios nGen = 4.880.000.");
    Printf("*****");

    // Load tGenOld
    TString str_fGenOld = "";
    if(b6000k)  str_fGenOld = Form("Trees/STARlight/%s/tGenOld_6000k.root", sOut_subfolder.Data());
    else        str_fGenOld = Form("Trees/STARlight/%s/tGenOld_4880k.root", sOut_subfolder.Data());
    TFile *fGenOld = fGenOld = TFile::Open(str_fGenOld.Data(), "read");
    if(fGenOld) Printf("File %s loaded.", fGenOld->GetName());

    TTree *tGenOld = (TTree*) fGenOld->Get("tGenOld");
    if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName()); 
    tGenOld->SetBranchAddress("fPtGen", &fPtGenerated);

    // Load tGenNew
    TString str_fGenNew = "";
    if(b6000k)  str_fGenNew = Form("Trees/STARlight/%s/tGenNew_6000k_%.3f.root", sOut_subfolder.Data(), R_A);
    else        str_fGenNew = Form("Trees/STARlight/%s/tGenNew_4880k_%.3f.root", sOut_subfolder.Data(), R_A);
    TFile *fGenNew = TFile::Open(str_fGenNew.Data(),"read");
    if(fGenNew) Printf("File %s loaded.", fGenNew->GetName());

    TTree *tGenNew = (TTree*) fGenNew->Get("tGenNew");
    if(tGenNew) Printf("Tree %s loaded.", tGenNew->GetName()); 
    tGenNew->SetBranchAddress("fPtGen", &fPtGenerated);

    // Load histograms with generated events
    TH1D *hGenOld = (TH1D*)fGenOld->Get("hGenOld");
    if(hGenOld) Printf("Histogram %s loaded.", hGenOld->GetName());
    TH1D *hGenNew = (TH1D*)fGenNew->Get("hGenNew");
    if(hGenNew) Printf("Histogram %s loaded.", hGenNew->GetName());

    hRatios = (TH1D*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("hRatios");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);

    // Stop weighting at 0.26 GeV
    Int_t iStopWeight = (Int_t)(nBins_CohJ * 0.26/fPtCutUpp_CohJ + 1);
    Printf("\n");
    Printf("*****");
    Printf("iStopWeight = %i", iStopWeight);
    Printf("*****");
    for(Int_t iBin = iStopWeight; iBin <= nBins_CohJ; iBin++){
        hRatios->SetBinContent(iBin, 1.0);
        hRatios->SetBinError(iBin, 0.0);
    }

    hRecNew = (TH1D*)hRecOld->Clone("hRecNew");
    hRecNew->SetTitle("hRecNew");
    hRecNew->Multiply(hRatios);

    // Print the results to console
    Printf("\n");
    Printf("*****");

    Printf("Output:");
    Printf("fPtCutLow_CohJ\tfPtCutUpp_CohJ\tnEvOld\tnEvNew\tRatio");
    for(Int_t iBin = 1; iBin <= nBins_CohJ; iBin++){
        Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
            hRatios->GetBinLowEdge(iBin), hRatios->GetBinLowEdge(iBin+1), 
            hGenOld->GetBinContent(iBin), hGenNew->GetBinContent(iBin), hRatios->GetBinContent(iBin));
    }

    Printf("*****");
    Printf("\n");  

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
    cRatios->Print(Form("Results/STARlight/%s/Ratios/coh_modRA_%.3f.pdf", sOut_subfolder.Data(), R_A));
    cRatios->Print(Form("Results/STARlight/%s/Ratios/coh_modRA_%.3f.png", sOut_subfolder.Data(), R_A));
    cRec->Print(Form("Results/STARlight/%s/RecSpectra/coh_modRA_%.3f.pdf", sOut_subfolder.Data(), R_A));
    cRec->Print(Form("Results/STARlight/%s/RecSpectra/coh_modRA_%.3f.png", sOut_subfolder.Data(), R_A));

    delete hGenOld;
    delete hGenNew;

    return;
}

void FillHistRec()
{
    Printf("*****");
    Printf("Calculating histogram of original reconstructed events.");
    Printf("*****");

    TFile *fSL = NULL;
    fSL = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
    if(fSL) Printf("File %s loaded.", fSL->GetName());

    // Get the MCRec tree
    TTree *tRec = dynamic_cast<TTree*> (fSL->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("Tree %s loaded.", tRec->GetName());
    ConnectTreeVariablesMCRec(tRec);

    Int_t nEntriesAnalysed = 0;
    Int_t nEntriesProgress = (Double_t)tRec->GetEntries() / 20.;
    Int_t nPercent = 0;

    // Prepare reconstructed events
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(1, 2)){
            hRecOld->Fill(fPt);
        } 

        if((iEntry+1) % nEntriesProgress == 0){
            nPercent += 5;
            nEntriesAnalysed += nEntriesProgress;
            Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
        }
    }

    Printf("*****");
    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return;
}