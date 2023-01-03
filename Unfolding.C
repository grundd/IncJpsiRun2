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
Float_t tBoundaries[6] = { 0 };
TLatex* ltx = new TLatex();

void SetCanvas(TCanvas* c, Bool_t is2Dhist = kFALSE)
{
    c->SetTopMargin(0.08);
    c->SetLeftMargin(0.11);
    if(!is2Dhist) c->SetRightMargin(0.03);
    else          c->SetRightMargin(0.12);
    c->SetBottomMargin(0.12);
    return;
}

void SetHistoStyle(TH1F* h, Color_t c)
{
    // style
    gStyle->SetEndErrorSize(3);
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(c);
    h->SetLineColor(c);
    h->SetLineWidth(3);
    h->SetLineStyle(1);
    // x-axis
    h->GetXaxis()->SetTitle("|#it{t}| (GeV^{2})");
    h->GetXaxis()->SetTitleOffset(1.15);
    h->GetXaxis()->SetTitleSize(0.048);
    h->GetXaxis()->SetLabelSize(0.048);
    h->GetXaxis()->SetDecimals(1);
    // y-axis
    h->GetYaxis()->SetTitle("counts");
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleSize(0.048);
    h->GetYaxis()->SetLabelSize(0.048);
    h->GetYaxis()->SetMaxDigits(3);
    return;
}

void UnfoldAndPlotResults(TString subf, TString unf, RooUnfoldResponse* resp, TH1F* hRec, TH1F* hGen = NULL)
{
    // prepare the output path
    TString subfolder = subf + unf;
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "Unfolding/" + subfolder + "/");
    // how many times to perform unfolding: depends on the selected algorithm
    Int_t N;
    if(unf=="Bayes") N = 6; // number of iterations
    else if(unf=="Svd") N = 5; // regularisation parameter
    else if(unf=="BinByBin") N = 1; // only once
    // perform unfolding:
    for(Int_t i = 1; i <= N; i++)
    {
        RooUnfold* unfold = NULL;
        if(unf=="Bayes") unfold = new RooUnfoldBayes(resp,hRec,i);
        else if(unf=="Svd") unfold = new RooUnfoldSvd(resp,hRec,i);
        else if(unf=="BinByBin") unfold = new RooUnfoldBinByBin(resp,hRec);
        // get the histograms from the result
        TH1F* hUnfo = (TH1F*)(unfold->Hreco()->Clone("hUnfo"));
        TH1F* hMeas = (TH1F*)(unfold->Hmeasured()->Clone("hMeas"));
        TH1F* hTrue = NULL;
        // prepare the labels
        TString label;
        if(unf=="Bayes") label = Form("Bayes (%i iters)", i);
        else if(unf=="Svd") label = Form("SVD (reg par %i)", i);
        else if(unf=="BinByBin") label = "Bin by bin";
        // **********************
        // plot |t| distributions
        // **********************
        TCanvas* cDist = new TCanvas("cDist","cDist",800,600);
        SetCanvas(cDist); 
        cDist->cd();
        SetHistoStyle(hMeas,kBlue+1);
        SetHistoStyle(hUnfo,kGreen+1); 
        Float_t maxVal = hMeas->GetBinContent(hMeas->GetMaximumBin());
        Float_t minVal = hMeas->GetBinContent(hMeas->GetMinimumBin());
        hMeas->GetYaxis()->SetRangeUser(0.6*minVal, 1.15*maxVal);
        hMeas->Draw("E1");
        if(hGen) {
            hTrue = (TH1F*)hGen->Clone("hTrue");
            SetHistoStyle(hTrue,kRed+1); 
            hTrue->Draw("HIST SAME");
        }
        hUnfo->Draw("E1 SAME");
        ltx->DrawLatex(0.55,0.95,Form("Unfolded |#it{t}| distribution: %s",label.Data()));
        // plot legend
        Int_t nRows = 8;
        TString labelMeas = "measured";
        if(hGen) { nRows++; labelMeas = "pseudo-data"; }
        TLegend* l = new TLegend(0.62,0.91-0.055*nRows,0.95,0.91);
        l->AddEntry(hUnfo, "unfolded", "EPL");
        l->AddEntry(hMeas, Form("%s",labelMeas.Data()), "EPL");
        if(hGen) l->AddEntry(hTrue, "MC truth", "L");
        l->AddEntry((TObject*)0,Form("diff unfo./%s:",labelMeas.Data()),"");
        for(Int_t i = 0; i < nPtBins; i++) {
            Float_t diff = (1. - hUnfo->GetBinContent(i+1) / hMeas->GetBinContent(i+1)) * 100.;
            l->AddEntry((TObject*)0,Form("bin %i: %.1f%%",i+1,diff),"");
        }
        l->SetTextSize(0.042);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetMargin(0.17);
        l->Draw();
        // print the result
        cDist->Print("Results/" + str_subfolder + Form("Unfolding/%s/tDists_%02i.pdf",subfolder.Data(),i));
        delete cDist;
        // ***************
        // plot unf matrix
        // ***************
        /*
        TCanvas* cUnfMtx = new TCanvas("cUnfMtx","cUnfMtx",800,600);
        SetCanvas(cUnfMtx,kTRUE);
        cUnfMtx->cd();
        TMatrixD UnfMtx = unfold->UnfoldingMatrix();
        UnfMtx.Draw("colzTEXT");
        ltx->DrawLatex(0.55,0.95,Form("Unfolding matrix: %i iterations", i));
        cUnfMtx->Print("Results/" + str_subfolder + Form("Unfolding/%s/unfMtx_02i.pdf",subfolder.Data(),i));
        delete cUnfMtx;
        */
        // ***************
        // plot cov matrix
        // ***************
        TCanvas* cCov = new TCanvas("cCov","cCov",800,600);
        SetCanvas(cCov,kTRUE);
        cCov->cd();
        TMatrixD CovMtx = unfold->Ereco();
        CovMtx.Draw("colzTEXT");
        ltx->DrawLatex(0.55,0.95,Form("Covariance matrix: %s",label.Data()));
        cCov->Print("Results/" + str_subfolder + Form("Unfolding/%s/covMtx_%02i.pdf",subfolder.Data(),i));
        delete cCov;
        // ****************
        // errs of unf dist
        // ****************
        ofstream of("Results/" + str_subfolder + Form("Unfolding/%s/hUnfoErrs_%02i.txt",subfolder.Data(),i));
        of << std::fixed << std::setprecision(1);
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++)
        {
            of << iBin << "\t" 
               << hUnfo->GetBinContent(iBin) << "\t" 
               << hUnfo->GetBinError(iBin) << "\t" 
               << TMath::Sqrt(CovMtx[iBin][iBin]) << "\n"; // first bin is underflow bin
        }
        of.close();
        // **********************
        // plot norm to bin width
        // **********************
        hMeas->Scale(1., "width");
        hUnfo->Scale(1., "width");
        if(hGen) hTrue->Scale(1., "width");
        TCanvas* cNorm = new TCanvas("cNorm","cNorm",800,600);
        SetCanvas(cNorm);
        cNorm->cd();
        hMeas->GetYaxis()->SetTitle("d#it{N}/|#it{t}|");
        hMeas->Draw("E1");
        if(hGen) hTrue->Draw("HIST SAME");
        hUnfo->Draw("E1 SAME");
        ltx->DrawLatex(0.55,0.95,Form("Unfolded |#it{t}| distribution: %s",label.Data()));
        l->Draw();
        cNorm->Print("Results/" + str_subfolder + Form("Unfolding/%s/normDists_%02i.pdf",subfolder.Data(),i));
        delete cNorm;
        // delete the histograms
        delete unfold;
        delete hMeas;
        delete hUnfo;
        delete hTrue;
    }
    // how many iterations to use?
    // let's take a look at the absolute errors of hUnfo after each iteration
    // and the relative differences in hUnfo values between subsequent iterations
    if(unf=="Bayes") {
        Float_t unfoVal[nIter][5] = { 0 };
        Float_t unfoErr[nIter][5] = { 0 };
        Float_t diffRel[nIter-1][5] = { 0 }; // in percent
        // loop over iterations
        for(Int_t iIt = 1; iIt <= nIter; iIt++)
        {
            ifstream ifs("Results/" + str_subfolder + Form("Unfolding/%s/hUnfoErrs_%02i.txt",subfolder.Data(),iIt));
            for(Int_t iBin = 0; iBin < nPtBins; iBin++) {
                Int_t i; Float_t val, err1, err2;
                ifs >> i >> val >> err1 >> err2;
                unfoErr[iIt-1][iBin] = err1;
                unfoVal[iIt-1][iBin] = val;
            }
            ifs.close();
        }
        // calculate the relative differences
        for(Int_t iIt = 1; iIt <= nIter-1; iIt++) 
            for(Int_t iBin = 0; iBin < nPtBins; iBin++) 
                diffRel[iIt-1][iBin] = TMath::Abs(1. - unfoVal[iIt][iBin] / unfoVal[iIt-1][iBin]) * 100.;
        // print the values
        ofstream of("Results/" + str_subfolder + Form("Unfolding/%s/iterationsErrsDiffs.txt",subfolder.Data()));
        // absolute errors
        of << "absolute errors vs iterations and |t| bins:\n"
           << "it/bin\t1 \t2 \t3 \t4 \t5 \n"
           << std::fixed << std::setprecision(0);
        for(Int_t iIt = 0; iIt < nIter; iIt++) {
            of << iIt+1 << "\t";
            for(Int_t iBin = 0; iBin < nPtBins-1; iBin++) of << unfoErr[iIt][iBin] << "\t";
            of << unfoErr[iIt][nPtBins-1] << "\n";
        }
        of << "\n";
        // relative errors
        of << "relative diff in consequent iterations:\n"
           << "it/bin\t1 \t2 \t3 \t4 \t5 \n"
           << std::fixed << std::setprecision(2);
        for(Int_t iIt = 1; iIt < nIter; iIt++) {
            of << iIt << "-" << iIt+1 << "\t";
            for(Int_t iBin = 0; iBin < nPtBins-1; iBin++) of << diffRel[iIt-1][iBin] << "\t";
            of << diffRel[iIt-1][nPtBins-1] << "\n";
        }
        of.close();
    }
    return;
}

void Unfolding(Int_t iAnalysis)
{
    gSystem->Load("RooUnfold/libRooUnfold");
    InitAnalysis(iAnalysis);
    SetPtBinning();
    for(Int_t i = 0; i < nPtBins+1; i++) tBoundaries[i] = TMath::Power(ptBoundaries[i], 2);
    ltx->SetTextSize(0.048);
    ltx->SetTextAlign(21);
    ltx->SetNDC();

    // create the output folder
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "Unfolding/");

    // histograms with binning at the gen and rec level
    // can be empty, needed only to specify the dimensions of the distributions in RooUnfoldResponse
    TH1F* hTrainGen = new TH1F("hTrainGen","binning at the gen level",nPtBins,tBoundaries);
    TH1F* hTrainRec = new TH1F("hTrainRec","binning at the rec level",nPtBins,tBoundaries);

    // create the response matrix
    RooUnfoldResponse response(hTrainRec, hTrainGen);
    // set under/overflow bins
    response.UseOverflow(kTRUE);
    Printf("Binning of the response matrix (RM) defined.");
    Printf(" - Use of under/overflow bins? %o", response.UseOverflowStatus());

    // MC tree: kIncohJpsiToMu
    TFile *f = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
    if(f) Printf("Input file %s loaded.", f->GetName());
    TTree *t = dynamic_cast<TTree*> (f->Get(str_in_MC_tree_rec.Data()));
    if(t) Printf("Input tree %s loaded.", t->GetName());
    ConnectTreeVariablesMCRec(t);

    // go over the events and fill the RM
    Int_t nEn = t->GetEntries();
    Int_t nTrain = (Int_t)(nEn * fTrain);
    Printf("There is %i events in the dataset:", nEn);
    Printf(" - %i of those (%.0f%%) will be used to train the RM", nTrain, fTrain*100);
    Printf(" - the rest will be used to test it");
    
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
            response.Fill(fPt*fPt,fPtGen*fPtGen);
        } 
        /*
        else {
            response.Miss(fPtGen*fPtGen);
        }
        */
    }

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.4f");

    // plot the reponse matrix
    TMatrixD RespMtx = response.Mresponse();
    TCanvas* cMtx = new TCanvas("cMtx","cMtx",800,600);
    cMtx->cd();
    TH2D* hRespMtx = new TH2D(*static_cast<TMatrixDBase*>(&RespMtx));
    hRespMtx->SetMarkerSize(2.0);
    hRespMtx->GetXaxis()->SetDecimals(1);
    hRespMtx->GetYaxis()->SetDecimals(1);
    hRespMtx->Draw("colzTEXT");
    ltx->DrawLatex(0.55,0.94,"Response matrix");
    cMtx->Print("Results/" + str_subfolder + "Unfolding/responseMatrix.pdf");
    delete hRespMtx;

    // unfolding:

    // ***************************
    // testing the response matrix
    // ***************************
    TH1F* hTestGen = new TH1F("hTestGen","|#it{t}| spectrum at the gen level",nPtBins,tBoundaries);
    TH1F* hTestRec = new TH1F("hTestRec","|#it{t}| spectrum at the rec level",nPtBins,tBoundaries);
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
            hTestGen->Fill(fPtGen*fPtGen);
            hTestRec->Fill(fPt*fPt);
        }
    }
    // unfold hTestRec:
    TString subf = "testMC";
    UnfoldAndPlotResults(subf,"Bayes",&response,hTestRec,hTestGen);
    //UnfoldAndPlotResults(subf,"Svd",&response,hTestRec,hTestGen);
    //UnfoldAndPlotResults(subf,"BinByBin",&response,hTestRec,hTestGen);

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
        Int_t bin;
        ifs >> bin >> AxE_val[iBin] >> AxE_err[iBin];
    }
    ifs.close();
    // import fCs
    sIn = "Results/" + str_subfolder + "PtFit_NoBkg/RecSh4_fD0_fC.txt";
    ifs.open(sIn.Data());
    Int_t i = 0;
    std::string str;
    while(std::getline(ifs,str)){
        istringstream iss(str);
        Int_t bin;
        if(i > 0) iss >> bin >> fC_val[i-1] >> fC_err[i-1]; // skip the first line
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
    TH1F* hToUnfold = new TH1F("hToUnfold","",nPtBins,tBoundaries);
    for(Int_t iBin = 0; iBin < nPtBins; iBin++) 
    {
        Float_t binWidth = tBoundaries[iBin+1] - tBoundaries[iBin];
        corr_err[iBin] = TMath::Sqrt(
              TMath::Power(fC_err[iBin],2)
            + TMath::Power(fD_err[iBin],2)) / 100.;
        corr_val[iBin] = 1. + fC_val[iBin] / 100. + fD_val[iBin] / 100.;
        Float_t value = N_yield_val[iBin] / corr_val[iBin] / (AxE_val[iBin] / 100.);
        Float_t error = value * TMath::Sqrt(TMath::Power(N_yield_err[iBin] / N_yield_val[iBin], 2));
            //+ TMath::Power(AxE_err[iBin] / AxE_val[iBin], 2)
            //+ TMath::Power(corr_err[iBin] / corr_val[iBin], 2)
        hToUnfold->SetBinContent(iBin+1,value);
        hToUnfold->SetBinError(iBin+1,error);
    }
    // import the cross section
    TH1F* hToUnfoldCS = new TH1F("hToUnfoldCS","",nPtBins,tBoundaries);
    ifs.open("Results/" + str_subfolder + "CrossSec/CrossSec_photo.txt");
    for(Int_t i = 0; i < nPtBins; i++)
    {
        // cross section values in mub
        Int_t bin; Float_t tLow, tUpp, sig_val, sig_err_stat, sig_err_syst_uncr, sig_err_syst_corr;
        ifs >> bin >> tLow >> tUpp >> sig_val >> sig_err_stat >> sig_err_syst_uncr >> sig_err_syst_corr;
        hToUnfoldCS->SetBinContent(i+1, sig_val * (tUpp - tLow));
        hToUnfoldCS->SetBinError(i+1, sig_err_stat * (tUpp - tLow));
    }
    ifs.close();
    // print the values
    ofstream of("Results/" + str_subfolder + "Unfolding/hToUnfold.txt");
    of << "bin \tN \terr \tAxE \terr \tfC \terr \tfD \terr \tcorr \terr \ttoUnf \terr \ttoUnfCS\t stat \n";
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
           << hToUnfold->GetBinContent(iBin+1) << "\t" << hToUnfold->GetBinError(iBin+1) << "\t"
           << std::fixed << std::setprecision(2)
           << hToUnfoldCS->GetBinContent(iBin+1) << "\t" << hToUnfoldCS->GetBinError(iBin+1) << "\n";
    }
    of.close();

    // unfold hToUnfold:
    subf = "unfMeas";
    ///*
    UnfoldAndPlotResults(subf,"Bayes",&response,hToUnfold);
    //UnfoldAndPlotResults(subf,"Svd",&response,hToUnfold);
    //*/
    //UnfoldAndPlotResults(subf,"BinByBin",&response,hToUnfold);

    // unfold the full cross section with the statistic error
    subf = "unfCS";
    ///*
    UnfoldAndPlotResults(subf,"Bayes",&response,hToUnfoldCS);
    //UnfoldAndPlotResults(subf,"Svd",&response,hToUnfoldCS);
    //*/
    //UnfoldAndPlotResults(subf,"BinByBin",&response,hToUnfoldCS);

    return;
}