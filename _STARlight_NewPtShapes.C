// _STARlight_NewPtShapes.C
// David Grund, Mar 12, 2022

// cpp headers
#include <fstream>
#include <stdio.h>
// root headers
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TAxis.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h" // to be able to use SetReducedRunList()

Bool_t drawCheck(kFALSE);

TString strMCArr[4] = {"CohJ","IncJ","CohP","IncP"};
TString strMCFilesArr[4] = {"kCohJpsiToMu","kIncohJpsiToMu","kCohPsi2sToMuPi","kIncohPsi2sToMuPi"};
Double_t fPtCutLowArr[4] = {0.0, 0.0, 0.0, 0.0}; // GeV/c
Double_t fPtCutUppArr[4] = {0.4, 1.4, 0.6, 1.4}; // GeV/c
Int_t nBinsArr[4] = {80, 280, 120, 280};

TString strMC;
TString strMCFiles;
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
TTree* GetTreeRec();
void CalcAndPlotRatios(const char* subfolder_out, Double_t R_A);
void ConnectTreeVariables_tSL(TTree *tSL);


void _STARlight_NewPtShapes()
{
    // no title in plots
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // to study which R_A is optimal for CohJ to describe the measured pT distribution
    if(kTRUE){
        // CohJ, 6.000.000 gen events,
        // R_A = 6.624 vs. 6.600-7.800 fm (step = 0.100 fm)
        for(Int_t i = 0; i < 13; i++){
            Double_t R_A = 6.600 + i*0.100; // fm
            InitAnalysis(0,kTRUE);
            FillTreeGen("Trees/STARlight/CohJ_6.624/",6.624);
            FillTreeGen(Form("Trees/STARlight/CohJ_%.3f/", R_A),R_A);
            CalcAndPlotRatios("",R_A);
        }
    }

    if(kTRUE){
        // CohJ, R_A = 6.624 vs. 7.350 fm, 6.000.000 gen events
        InitAnalysis(0,kTRUE);
        FillTreeGen("Trees/STARlight/CohJ_6.624/",6.624);
        FillTreeGen("Trees/STARlight/CohJ_7.350/",7.350);
        CalcAndPlotRatios("",7.350);
    }
    if(kTRUE){
        // IncJ, R_A = 6.624 vs. 7.350 fm, 6.000.000 gen events
        InitAnalysis(1,kTRUE);
        FillTreeGen("Trees/STARlight/IncJ_6.624/",6.624);
        FillTreeGen("Trees/STARlight/IncJ_7.350/",7.350);
        CalcAndPlotRatios("",7.350);
    }
    if(kTRUE){
        // CohP, R_A = 6.624 vs. 7.350 fm, 6.000.000 gen events
        InitAnalysis(2,kTRUE);
        FillTreeGen("Trees/STARlight/CohP_6.624/",6.624);
        FillTreeGen("Trees/STARlight/CohP_7.350/",7.350);
        CalcAndPlotRatios("",7.350);
    }
    if(kTRUE){
        // IncP, R_A = 6.624 vs. 7.350 fm, 6.000.000 gen events
        InitAnalysis(3,kTRUE);
        FillTreeGen("Trees/STARlight/IncP_6.624/",6.624);
        FillTreeGen("Trees/STARlight/IncP_7.350/",7.350);
        CalcAndPlotRatios("",7.350);
    }

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
    strMCFiles = strMCFilesArr[iMC];
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

    if(fGen)
    {
        Printf("File %s already created.", str_f_out.Data());
    } 
    else 
    {
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

        for(Int_t iEntry = 0; iEntry < tSL->GetEntries(); iEntry++)
        {
            tSL->GetEntry(iEntry);
            if(TMath::Abs(fYGen) < 1.0)
            {
                nEvOld++;
                fPtGenerated = parent->Pt();
                hGen->Fill(fPtGenerated);
                tGen->Fill();
            }
            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0)
            {
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

    TH1D *hGenNew = (TH1D*)fGenNew->Get("hGen");
    if(hGenNew) Printf("Histogram %sNew loaded.", hGenNew->GetName());

    // Draw loaded histograms
    if(drawCheck)
    {
        TCanvas *cGenOld = new TCanvas("cGenOld","cGenOld",900,600);
        hGenOld->Draw();
        TCanvas *cGenNew = new TCanvas("cGenNew","cGenNew",900,600);
        hGenNew->Draw();
    } 

    // Calculate ratios
    hRatios = (TH1D*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("hRatios");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);

    // Print the ratios to the text file
    TString str_folder_out = Form("Results/_STARlight_NewPtShapes/%s/", subfolder_out);
    gSystem->Exec("mkdir -p " + str_folder_out);
    TString str_file_out = Form("%s%s_RA_%.3f.txt", str_folder_out.Data(), strMC.Data(), R_A);
    ofstream outfile(str_file_out.Data());
    outfile << "pT_low\tpT_upp\tnEvOld\tnEvNew\tratio\n";
    for(Int_t iBin = 1; iBin <= nBins; iBin++)
    {
        outfile << Form("%.3f\t%.3f\t%.0f\t%.0f\t%.3f\n",
            hRatios->GetBinLowEdge(iBin), hRatios->GetBinLowEdge(iBin+1), 
            hGenOld->GetBinContent(iBin), hGenNew->GetBinContent(iBin), hRatios->GetBinContent(iBin));
    }
    outfile.close();

    // define the tRec tree
    TTree *tRec = GetTreeRec();
    // define the integers showing the progress bar
    Int_t nEntriesAnalysed = 0;
    Int_t nEntriesProgress = (Double_t)tRec->GetEntries() / 20.;
    Int_t nPercent = 0;
    // if J/psi datasets
    if (strncmp(strMC.Data(),"CohJ",4) == 0 || strncmp(strMC.Data(),"IncJ",4) == 0) 
    {
        // run over reconstructed events
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++)
        {
            tRec->GetEntry(iEntry);
            // m between 3.0 and 3.2 GeV/c^2, pT cut: all
            if(EventPassedMCRec(1, 2))
            {
                // fill the histogram hRecOld
                hRecOld->Fill(fPt);
            } 
            if((iEntry+1) % nEntriesProgress == 0)
            {
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        } 
        // draw hRecOld 
        if(drawCheck)
        {
            TCanvas *cRecOld2 = new TCanvas("cRecOld2","cRecOld2",900,600);
            hRecOld->Draw();
        }
        // create hRecNew by scaling hRecOld with ratios
        hRecNew = (TH1D*)hRecOld->Clone("hRecNew");
        hRecNew->SetTitle("hRecNew");
        hRecNew->Multiply(hRatios);
    }
    // if Psi(2s) datasets
    else if (strncmp(strMC.Data(),"CohP",4) == 0 || strncmp(strMC.Data(),"IncP",4) == 0) 
    {
        hRecNew = new TH1D("hRecNew","hRecNew",nBins,fPtCutLow,fPtCutUpp);
        TAxis *xAxis = hRatios->GetXaxis();
        // run over reconstructed events
        Int_t iBinP = 0;
        Int_t iBinJ = 0;
        Double_t fJpsi = 0;
        //Printf("iBinP\tPtGenP\tiBinJ\tPtGenJ\tfJpsi\toldVal\tnewVal");
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++)
        {
            tRec->GetEntry(iEntry);
            // m between 3.0 and 3.2 GeV/c^2, pT cut: all
            if(EventPassedMCRec(1, 2))
            {
                // find index of the bin to which the current fPtGen_Psi2s corresponds
                iBinP = xAxis->FindBin(fPtGen_Psi2s);
                // find index of the bin to which the current fPt corresponds
                iBinJ = xAxis->FindBin(fPt);
                // scale the J/psi entry by the ratio with the index iBinPsi2s
                fJpsi = hRatios->GetBinContent(iBinP);
                // add the entry to hRecNew
                hRecNew->SetBinContent(iBinJ,hRecNew->GetBinContent(iBinJ)+fJpsi);
                // print the values to the console
                //Printf("%i\t%.3f\t%i\t%.3f\t%.3f\t%.3f\t%.3f", iBinP, fPtGen_Psi2s, iBinJ, fPtGen, fJpsi, hRecNew->GetBinContent(iBinJ), hRecNew->GetBinContent(iBinJ)+fJpsi);
                // fill the histogram hRecOld
                hRecOld->Fill(fPt);
            }
            if((iEntry+1) % nEntriesProgress == 0)
            {
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        } 
        // draw hRecOld 
        if(drawCheck)
        {
            TCanvas *cRecOld2 = new TCanvas("cRecOld2","cRecOld2",900,600);
            hRecOld->Draw();
        }
    } 
    else 
    {
        Printf("This option is not supported. Terminating..."); 
        return;
    }

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
    if (strncmp(strMC.Data(),"CohP",4) == 0 || strncmp(strMC.Data(),"IncP",4) == 0) hRatios->GetXaxis()->SetTitle("#it{p}_{T,#psi(2S)}^{gen} (GeV/#it{c})");
    else                                                                            hRatios->GetXaxis()->SetTitle("#it{p}_{T,J/#psi}^{gen} (GeV/#it{c})");
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
    hRecOld->GetXaxis()->SetTitle("#it{p}_{T,J/#psi}^{gen} (GeV/#it{c})");   
    hRecOld->GetXaxis()->SetTitleSize(0.05);
    hRecOld->GetXaxis()->SetTitleOffset(1.2);
    hRecOld->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRecOld->SetLineColor(kRed);
    hRecOld->Draw("E1");
    hRecNew->SetLineColor(215);
    hRecNew->Draw("E1 SAME");

    TLegend *l = new TLegend(0.50,0.85,0.90,0.95);
    l->AddEntry(hRecOld,Form("rec. events with #it{R}_{A} = 6.624 fm"),"L");
    l->AddEntry(hRecNew,Form("rec. events with #it{R}_{A} = %.3f fm", R_A),"L");
    l->SetTextSize(0.04);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    // Print plots
    cRatios->Print(Form("%s%s_RA_%.3f_ratios.pdf", str_folder_out.Data(), strMC.Data(), R_A));
    cRatios->Print(Form("%s%s_RA_%.3f_ratios.png", str_folder_out.Data(), strMC.Data(), R_A));
    cRec->Print(Form("%s%s_RA_%.3f_recSpectra.pdf", str_folder_out.Data(), strMC.Data(), R_A));
    cRec->Print(Form("%s%s_RA_%.3f_recSpectra.png", str_folder_out.Data(), strMC.Data(), R_A));

    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return;
}

// #############################################################################################

TTree* GetTreeRec()
{
    Printf("*****");
    Printf("Process: %s", strMC.Data());
    Printf("Getting the tRec.");

    // define fRec and tRec
    TFile *fRec = NULL;
    TTree *tRec = NULL;
    // define the paths to the file and tRec
    TString str_f_in = "Trees/";
    TString str_t_in = "AnalysisOutput/";
    // if pass3 or not
    if(!isPass3) str_f_in += "AnalysisDataMC_pass1/";
    else         str_f_in += "AnalysisDataMC_pass3/";
    // choose MC dataset
    Bool_t isPsi2sDataset = kFALSE;
    // if CohJ or IncJ
    if(strncmp(strMC.Data(),"CohJ",4) == 0 || strncmp(strMC.Data(),"IncJ",4) == 0)
    {
        str_f_in += "AnalysisResults_MC_" + strMCFiles + ".root";
        if(!isPass3) str_t_in += "fTreeJPsiMCRec";
        else         str_t_in += "fTreeJpsi";
        // open the input file
        fRec = TFile::Open(str_f_in.Data(), "read");
        if(fRec) Printf("File %s loaded.", fRec->GetName());
        // get the MCRec tree
        tRec = dynamic_cast<TTree*> (fRec->Get(str_t_in.Data()));
        if(tRec) Printf("Tree %s loaded.", tRec->GetName());
    }
    // if CohP or IncP
    if(strncmp(strMC.Data(),"CohP",4) == 0 || strncmp(strMC.Data(),"IncP",4) == 0) 
    {
        isPsi2sDataset = kTRUE;
        str_f_in += "AnalysisResults_MC_" + strMCFiles + "_2.root";
        if(!isPass3)
        {
            Printf("This option is not supported. Terminating..."); 
            return NULL;
        }  
        // open the input file
        fRec = TFile::Open(str_f_in.Data(), "read");
        if(fRec) Printf("File %s loaded.", fRec->GetName());
        // get fOutputList
        TList *l = (TList*) fRec->Get("AnalysisOutput/fOutputListcharged");
        if(l) Printf("List %s loaded.", l->GetName()); 
        // get the MCRec tree
        tRec = (TTree*)l->FindObject("fTreeJpsi");
        if(tRec) Printf("Tree %s loaded.", tRec->GetName());
    }

    // connect tree varibles
    ConnectTreeVariablesMCRec(tRec, isPsi2sDataset);

    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return tRec;
}

// #############################################################################################

void ConnectTreeVariables_tSL(TTree *tSL)
{
    tSL->SetBranchAddress("parent", &parent);
    tSL->SetBranchAddress("daughters", &daughters);

    Printf("Variables from %s connected.", tSL->GetName());

    return;
}

// #############################################################################################