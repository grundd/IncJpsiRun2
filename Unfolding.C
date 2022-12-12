// Unfolding.C
// David Grund, Nov 28, 2022

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
// cpp headers
#include <iostream>
#include <fstream>
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TLatex.h"
#include "TMath.h"
#include "RooUnfoldResponse.h"
#endif
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

const Int_t nIter = 6;
Bool_t scale = kFALSE;
Float_t fTrain = 0.8; // fraction of MC events used to train the response matrix
// the rest will be used to test the matrix

void SetCanvas(TCanvas* c, Bool_t is2Dhist = kFALSE)
{
    c->SetTopMargin(0.065);
    c->SetLeftMargin(0.10);
    if(!is2Dhist) c->SetRightMargin(0.03);
    else          c->SetRightMargin(0.12);
    return;
}

void SetHistoStyle(TH1F* h, Color_t c)
{
    // style
    gStyle->SetEndErrorSize(1); 
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(c);
    h->SetLineColor(c);
    h->SetLineWidth(2);
    h->SetLineStyle(1);
    // x-axis
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetDecimals(1);
    // y-axis
    h->GetYaxis()->SetTitle("counts");
    h->GetYaxis()->SetMaxDigits(3);
    return;
}

void UnfoldAndPlotResults(TString subfolder, RooUnfoldResponse* response, TH1F* hToUnfold, TH1F* hGen = NULL)
{
    if(hGen && scale) hGen->Scale(1.,"width");
    // prepare the output path
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "Unfolding/" + subfolder + "/");
    // perform bayes unfolding per each iteration
    for(Int_t iIt = 1; iIt <= nIter; iIt++)
    {
        RooUnfoldBayes unfold(response,hToUnfold,iIt);
        // get the histograms from the result
        TH1F* hUnfo = (TH1F*)(unfold.Hreco()->Clone("hUnfo"));
        TH1F* hMeas = (TH1F*)(unfold.Hmeasured()->Clone("hMeas"));
        // prepare the latex title
        TLatex* ltx = new TLatex();
        ltx->SetTextSize(0.033);
        ltx->SetTextAlign(21);
        ltx->SetNDC();
        // ***********************
        // * plot pT distributions
        // ***********************
        TCanvas* cRes = new TCanvas("cRes","cRes",700,600);
        SetCanvas(cRes);
        cRes->cd();
        SetHistoStyle(hMeas,kBlue+1);
        SetHistoStyle(hUnfo,kGreen+1); 
        if(scale) {
            hMeas->Scale(1., "width");
            hUnfo->Scale(1., "width");
        }
        Float_t maxVal = hMeas->GetBinContent(hMeas->GetMaximumBin());
        Float_t minVal = hMeas->GetBinContent(hMeas->GetMinimumBin());
        hMeas->GetYaxis()->SetRangeUser(0.6*minVal, 1.15*maxVal);
        hMeas->Draw("E1");
        if(hGen) {
            SetHistoStyle(hGen,kRed+1); 
            hGen->Draw("HIST SAME");
        }
        hUnfo->Draw("E1 SAME");
        ltx->DrawLatex(0.55,0.96,Form("Unfolded #it{p}_{T} distribution: %i iterations",iIt));
        // plot legend
        Int_t nRows = 2;
        if(hGen) nRows++;
        TLegend* l = new TLegend(0.35,0.92-0.045*nRows,0.55,0.92);
        l->AddEntry(hUnfo, "unfolded", "EPL");
        l->AddEntry(hMeas, "measured", "EPL");
        if(hGen) l->AddEntry(hGen, "MC truth", "L");
        l->SetTextSize(0.032);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetMargin(0.30);
        l->Draw();
        // print the result
        cRes->Print("Results/" + str_subfolder + Form("Unfolding/%s/ptDistributions_it%02i.pdf",subfolder.Data(),iIt));
        delete cRes;
        // *****************
        // * plot unf matrix
        // *****************
        TCanvas* cUnfMtx = new TCanvas("cUnfMtx","cUnfMtx",700,600);
        SetCanvas(cUnfMtx,kTRUE);
        cUnfMtx->cd();
        TMatrixD UnfMtx = unfold.UnfoldingMatrix();
        UnfMtx.Draw("colzTEXT");
        ltx->DrawLatex(0.55,0.96,Form("Unfolding matrix: %i iterations", iIt));
        cUnfMtx->Print("Results/" + str_subfolder + Form("Unfolding/%s/unfMtx_it%02i.pdf",subfolder.Data(),iIt));
        delete cUnfMtx;
        // *****************
        // * plot cov matrix
        // *****************
        TCanvas* cCov = new TCanvas("cCov","cCov",700,600);
        SetCanvas(cCov,kTRUE);
        cCov->cd();
        TMatrixD CovMtx = unfold.Ereco();
        CovMtx.Draw("colzTEXT");
        ltx->DrawLatex(0.55,0.96,Form("Covariance matrix: %i iterations", iIt));
        cCov->Print("Results/" + str_subfolder + Form("Unfolding/%s/covMtx_it%02i.pdf",subfolder.Data(),iIt));
        delete cCov;
        // ******************
        // * errs of unf dist
        // ******************
        ofstream of("Results/" + str_subfolder + Form("Unfolding/%s/hUnfoErrs_it%02i.txt",subfolder.Data(),iIt));
        of << std::fixed << std::setprecision(1);
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++)
        {
            of << iBin << "\t" 
               << hUnfo->GetBinContent(iBin) << "\t" 
               << hUnfo->GetBinError(iBin) << "\t" 
               << TMath::Sqrt(CovMtx[iBin-1][iBin-1]) << "\n";
        }
        of.close();
        // delete the histograms
        delete hMeas;
        delete hUnfo;
    }
    return;
}

void Unfolding(Int_t iAnalysis)
{
    gSystem->Load("RooUnfold/libRooUnfold");
    InitAnalysis(iAnalysis);
    SetPtBinning();

    // create the output folder
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "Unfolding/");

    // histograms with binning at the gen and rec level
    // can be empty, needed only to specify the dimensions of the distributions in RooUnfoldResponse
    TH1F* hTrainGen = new TH1F("hTrainGen","binning at the gen level",nPtBins,ptBoundaries);
    TH1F* hTrainRec = new TH1F("hTrainRec","binning at the rec level",nPtBins,ptBoundaries);

    // create the response matrix
    RooUnfoldResponse response(hTrainRec, hTrainGen);

    // MC tree: kIncohJpsiToMu
    TFile *f = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
    if(f) Printf("Input file %s loaded.", f->GetName());
    TTree *t = dynamic_cast<TTree*> (f->Get(str_in_MC_tree_rec.Data()));
    if(t) Printf("Input tree %s loaded.", t->GetName());
    ConnectTreeVariablesMCRec(t);

    // go over the events and fill the response matrix
    Int_t nEn = t->GetEntries();
    Int_t nTrain = (Int_t)(nEn * fTrain);
    Printf("There is %i events in the dataset:", nEn);
    Printf(" - %i of those (%.0f%%) will be used to train the response matrix (RM)", nTrain, fTrain*100);
    Printf(" - the rest will be used to test the matrix");
    
    Printf("Training the RM:");
    // progress bar:
    Float_t progress = 0.; // perc
    for(Int_t iEn = 0; iEn < nTrain; iEn++)
    {
        t->GetEntry(iEn);
        // update the progress bar
        if((iEn+1) % (Int_t)(nTrain/10.) == 0) {
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }
        // event passes all the selection criteria
        if(EventPassedMCRec(0,-1)) {
            response.Fill(fPt,fPtGen);
        } 
        //else {
        //    response.Miss(fPtGen);
        //}
    }

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.4f");

    // plot the reponse matrix
    TMatrixD RespMtx = response.Mresponse();
    TPaveText* label = new TPaveText(0.35,0.9,0.65,1.0,"brNDC");
    label->AddText("Response matrix");
    TCanvas* cMtx = new TCanvas("cMtx","cMtx",1600,900);
    cMtx->cd();
    RespMtx.Draw("colzTEXT");
    label->Draw("same");
    cMtx->Print("Results/" + str_subfolder + "Unfolding/responseMatrix.pdf");

    // unfolding:

    // ***************************
    // testing the response matrix
    // ***************************
    TH1F* hTestGen = new TH1F("hTestGen","p_{T} spectrum at the gen level",nPtBins,ptBoundaries);
    TH1F* hTestRec = new TH1F("hTestRec","p_{T} spectrum at the rec level",nPtBins,ptBoundaries);
    // fill the testing histograms
    progress = 0.; // perc
    for(Int_t iEn = nTrain; iEn < nEn; iEn++)
    {
        t->GetEntry(iEn);
        // update the progress bar
        if(((iEn-nTrain)+1) % (Int_t)((nEn-nTrain)/10.) == 0) {
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }
        if(EventPassedMCRec(0,-1)) { 
            hTestGen->Fill(fPtGen);
            hTestRec->Fill(fPt);
        }
    }
    // unfold hTestRec:
    TString subf = "testMC";
    UnfoldAndPlotResults(subf,&response,hTestRec,hTestGen);

    // *************************
    // unfolding the measurement
    // *************************
    Float_t N_yield_val[5] = { 0 };
    Float_t N_yield_err[5] = { 0 };
    Float_t AxE_val[5] = { 0 };
    Float_t AxE_err[5] = { 0 };
    Float_t fC_val[5] = { 0 };
    Float_t fC_err[5] = { 0 };
    Float_t fD_val[5] = { 0 };
    Float_t fD_err[5] = { 0 };
    Float_t corr_val[5] = { 0 }; // 1. + fC + fD
    Float_t corr_err[5] = { 0 };
    ifstream ifs; TString sIn;
    // import the yields
    for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
        sIn = Form("Results/" + str_subfolder + "InvMassFit/%ibins/bin%i_signal.txt", nPtBins, iBin+1);
        ifs.open(sIn.Data());
        ifs >> N_yield_val[iBin] >> N_yield_err[iBin];
        ifs.close(); 
    } 
    // import the AxE
    sIn = Form("Results/" + str_subfolder + "AxE_PtBins/AxE_%ibins.txt", nPtBins);
    ifs.open(sIn.Data());
    for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
        Int_t i_bin;
        ifs >> i_bin >> AxE_val[iBin] >> AxE_err[iBin];
    }
    ifs.close();
    // import fCs
    sIn = "Results/" + str_subfolder + "PtFit_NoBkg/RecSh4_fD0_fC.txt";
    ifs.open(sIn.Data());
    Int_t i = 0;
    std::string str;
    while(std::getline(ifs,str)){
        istringstream iss(str);
        Int_t i_bin;
        if(i > 0) iss >> i_bin >> fC_val[i-1] >> fC_err[i-1]; // skip the first line
        i++;   
    }
    ifs.close();
    // import fDs
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/fD_syst_errors.txt";
    ifs.open(sIn.Data());
    for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
        ifs >> fD_val[iBin] >> fD_err[iBin];
    }
    ifs.close();
    // create histogram to unfold
    TH1F* hToUnfold = new TH1F("hToUnfold","",nPtBins,ptBoundaries);
    for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
        corr_err[iBin] = TMath::Sqrt(
              TMath::Power(fC_err[iBin],2)
            + TMath::Power(fD_err[iBin],2)) / 100.;
        corr_val[iBin] = 1. + fC_val[iBin] / 100. + fD_val[iBin] / 100.;
        Float_t value = N_yield_val[iBin] / corr_val[iBin] / (AxE_val[iBin] / 100.);
        Float_t error = value * TMath::Sqrt(
              TMath::Power(N_yield_err[iBin] / N_yield_val[iBin], 2) 
            + TMath::Power(AxE_err[iBin] / AxE_val[iBin], 2) 
            + TMath::Power(corr_err[iBin] / corr_val[iBin], 2));
        hToUnfold->SetBinContent(iBin+1,value);
        hToUnfold->SetBinError(iBin+1,error);
    }
    // print the values
    subf = "unfData";
    ofstream of("Results/" + str_subfolder + "Unfolding/" + subf + "/hToUnfold.txt");
    of << "bin \tN \terr \tAxE \terr \tfC \terr \tfD \terr \tcorr \terr \thisto \terr \n";
    for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
        of << iBin+1 << "\t" 
           << std::fixed << std::setprecision(0)
           << N_yield_val[iBin] << "\t" << N_yield_err[iBin] << "\t" 
           << std::fixed << std::setprecision(3)
           << AxE_val[iBin] << "\t" << AxE_err[iBin] << "\t" 
           << std::fixed << std::setprecision(2)
           << fC_val[iBin] << "\t" << fC_err[iBin] << "\t" 
           << fD_val[iBin] << "\t" << fD_err[iBin] << "\t" 
           << std::fixed << std::setprecision(4)
           << corr_val[iBin] << "\t" << corr_err[iBin] << "\t"
           << std::fixed << std::setprecision(0)
           << hToUnfold->GetBinContent(iBin+1) << "\t" << hToUnfold->GetBinError(iBin+1) << "\n";
    }
    of.close();
    // unfold hToUnfold:
    UnfoldAndPlotResults(subf,&response,hToUnfold);
    // how many iterations to use?
    // let's take a look at the absolute and relative errors of hUnfo after each iteration
    Float_t errsAbs[nIter][5] = { 0 };
    Float_t errsRel[nIter][5] = { 0 };
    // loop over iterations
    for(Int_t iIt = 1; iIt <= nIter; iIt++)
    {
        ifs.open("Results/" + str_subfolder + Form("Unfolding/%s/hUnfoErrs_it%02i.txt",subf.Data(),iIt));
        for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
            Int_t i; Float_t val, err1, err2;
            ifs >> i >> val >> err1 >> err2;
            errsAbs[iIt-1][iBin] = err1;
            errsRel[iIt-1][iBin] = err1 / val * 100.;
        }
        ifs.close();
    }
    // print the values
    of.open("Results/" + str_subfolder + Form("Unfolding/%s/errsAbsRel.txt",subf.Data()));
    // absolute errors
    of << "absolute errors vs iterations and pT bins:\n"
       << "it/bin\t1 \t2 \t3 \t4 \t5 \n"
       << std::fixed << std::setprecision(0);
    for(Int_t iIt = 0; iIt < nIter; iIt++) {
        of << iIt+1 << "\t";
        for(Int_t iBin = 0; iBin < nPtBins-1; iBin++) {
            of << errsAbs[iIt][iBin] << "\t";
        }
        of << errsAbs[iIt][nPtBins-1] << "\n";
    }
    of << "\n";
    // relative errors
    of << "relative errors vs iterations and pT bins:\n"
       << "it/bin\t1 \t2 \t3 \t4 \t5 \n"
       << std::fixed << std::setprecision(2);
    for(Int_t iIt = 0; iIt < nIter; iIt++) {
        of << iIt+1 << "\t";
        for(Int_t iBin = 0; iBin < nPtBins-1; iBin++) {
            of << errsRel[iIt][iBin] << "\t";
        }
        of << errsRel[iIt][nPtBins-1] << "\n";
    }
    of.close();

    return;
}