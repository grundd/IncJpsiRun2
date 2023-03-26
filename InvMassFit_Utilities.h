// InvMassFit_Utilities.h
// David Grund, Apr 04, 2022

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

using namespace RooFit;

Double_t N_Jpsi_all[2]; // number of J/psi events with arbitrary mass
Double_t N_bkgr_peak[2];// number of J/psi events with mass from 3.0 to 3.2 GeV (around the J/psi peak)
Double_t N_Jpsi_peak[2];// number of bckgr events with mass from 3.0 to 3.2 GeV (around the J/psi peak)

void InvMassFit_DrawCorrMatrix(TCanvas *cCorrMat, RooFitResult* fResFit)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    cCorrMat->SetTopMargin(0.03);
    cCorrMat->SetBottomMargin(0.12);
    cCorrMat->SetRightMargin(0.19);
    cCorrMat->SetLeftMargin(0.13);

    TH2* hCorr = fResFit->correlationHist();
    hCorr->SetMarkerSize(3.);

    hCorr->GetXaxis()->SetBinLabel(1,"#it{M}_{J/#psi}");
    hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{bkg}");
    hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{J/#psi}");
    hCorr->GetXaxis()->SetBinLabel(4,"#lambda");
    hCorr->GetXaxis()->SetBinLabel(5,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(1,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(2,"#lambda");
    hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{J/#psi}");
    hCorr->GetYaxis()->SetBinLabel(4,"#it{N}_{bkg}");
    hCorr->GetYaxis()->SetBinLabel(5,"#it{M}_{J/#psi}");

    hCorr->GetXaxis()->SetLabelSize(0.1);
    hCorr->GetXaxis()->SetLabelOffset(0.012);
    hCorr->GetYaxis()->SetLabelSize(0.1);
    hCorr->GetZaxis()->SetLabelSize(0.07);
    // https://root-forum.cern.ch/t/colz-color-palette-font-and-size/15263
    hCorr->Draw("colz,text");

    return;
}

void InvMassFit_SetCanvas(TCanvas *c, Bool_t bLogScale)
{
    if(bLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    return;
}

void InvMassFit_PrepareData(Int_t iMassCut)
{
    TString name;
    if(iMassCut == 0) name = "Trees/" + str_subfolder + "InvMassFit/InvMassFit.root";
    if(iMassCut == 2) name = "Trees/" + str_subfolder + "InvMassFit/InvMassFit_SystUncertainties.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("Data trees already created.");
        return;

    } else { 

        Printf("Data trees will be created.");

        TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
        if(f_in) Printf("Input data loaded.");

        TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
        if(t_in) Printf("Input tree loaded.");

        ConnectTreeVariables(t_in);

        // Create new data tree with applied cuts
        file = new TFile(name.Data(),"RECREATE");

        TTree *tIncEnrSample = new TTree("tIncEnrSample", "tIncEnrSample");
        tIncEnrSample->Branch("fPt", &fPt, "fPt/D");
        tIncEnrSample->Branch("fM", &fM, "fM/D");
        tIncEnrSample->Branch("fY", &fY, "fY/D");

        TTree *tCohEnrSample = new TTree("tCohEnrSample", "tCohEnrSample");
        tCohEnrSample->Branch("fPt", &fPt, "fPt/D");
        tCohEnrSample->Branch("fM", &fM, "fM/D");
        tCohEnrSample->Branch("fY", &fY, "fY/D");

        TTree *tMixedSample = new TTree("tMixedSample", "tMixedSample");
        tMixedSample->Branch("fPt", &fPt, "fPt/D");
        tMixedSample->Branch("fM", &fM, "fM/D");
        tMixedSample->Branch("fY", &fY, "fY/D");

        Printf("%lli entries found in the tree.", t_in->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
            t_in->GetEntry(iEntry);
            // iMassCut (for syst uncertainties = 2, otherwise = 0), pT cut: inc, coh, all
            if(EventPassed(iMassCut, 0)) tIncEnrSample->Fill();
            if(EventPassed(iMassCut, 1)) tCohEnrSample->Fill();
            if(EventPassed(iMassCut, 2)) tMixedSample->Fill();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        file->Write("",TObject::kWriteDelete);

        return;
    }
}

void InvMassFit_DoFit(Int_t opt, Double_t fMCutLow, Double_t fMCutUpp, Double_t fAlpha_L, Double_t fAlpha_R, Double_t fN_L, Double_t fN_R, TString str_out, Bool_t isSystUncr = kFALSE, Double_t fCutZ = -1)
{
    // Fit the invariant mass distribution using Double-sided CB function
    // Fix the values of the tail parameters to MC values
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fPtCut     = -999;
    Double_t fPtCutLow  = -999;
    Double_t fPtCutUpp  = -999;
    Double_t fYCut      = 0.80;

    switch(opt){
        case 0: // 'inc': incoherent-enriched sample
            fPtCut = 0.20;
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 1: // 'coh': coherent-enriched sample
            fPtCut = 0.20;
            sprintf(fStrReduce,"abs(fY)<%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 2: // 'all': total sample (pT < 2.0 GeV/c)
            fPtCut = 2.00;
            sprintf(fStrReduce,"abs(fY)<%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 3: // 'allbins': sample with pT from 0.2 to 1 GeV/c 
            fPtCutLow = 0.20;
            fPtCutUpp = 1.00;
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 4: // pT bin 1
            fPtCutLow = ptBoundaries[0];
            fPtCutUpp = ptBoundaries[1];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 5: // pT bin 2
            fPtCutLow = ptBoundaries[1];
            fPtCutUpp = ptBoundaries[2];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 6: // pT bin 3
            fPtCutLow = ptBoundaries[2];
            fPtCutUpp = ptBoundaries[3];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 7: // pT bin 4
            fPtCutLow = ptBoundaries[3];
            fPtCutUpp = ptBoundaries[4];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 8: // pT bin 5
            fPtCutLow = ptBoundaries[4];
            fPtCutUpp = ptBoundaries[5];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
    }

    // Binning:
    Int_t nBins = 115; // so that each bin between 2.2 and 4.5 GeV is 20 MeV wide
    RooBinning binM(nBins,fMCutLow,fMCutUpp);
    Double_t BinSizeDouble = (fMCutUpp - fMCutLow) * 1000 / nBins; // in MeV
    BinSizeDouble = BinSizeDouble + 0.5;
    // https://stackoverflow.com/questions/9695329/c-how-to-round-a-double-to-an-int
    Int_t BinSize = (Int_t)BinSizeDouble;

    Printf("\n");
    Printf("*** Bin size (double): %.3f ***", BinSizeDouble);
    Printf("*** Bin size (int): %i ***\n", BinSize);
  
    // Roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0,10.);
    RooRealVar fY("fY","fY",-0.8,0.8);

    //fM.setBinning(binM);

    // Get the data trees
    TFile *f_in = NULL;
    // ordinary fits:
    if(isSystUncr == kFALSE) f_in = new TFile("Trees/" + str_subfolder + "InvMassFit/InvMassFit.root"); 
    // systematic uncertainties 
    else 
    {
        // related to signal extraction
        if(fCutZ == -1) f_in = new TFile("Trees/" + str_subfolder + "InvMassFit/InvMassFit_SystUncertainties.root"); 
        // related to modifications of Z vertex cut 
        else f_in = new TFile("Trees/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit.root", fCutZ)); 
    }

    TTree *t_in = NULL;
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        f_in->GetObject("tIncEnrSample",t_in);
    } else if(opt == 1){
        f_in->GetObject("tCohEnrSample",t_in);
    } else if(opt == 2){
        f_in->GetObject("tMixedSample",t_in);
    }

    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*t_in));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Crystal Ball parameters from MC (to be fixed)
    // loaded in InvMassFit_SetFit()

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,30.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,30.);
    alpha_L.setConstant(kTRUE);
    alpha_R.setConstant(kTRUE);
    n_L.setConstant(kTRUE);
    n_R.setConstant(kTRUE);

    // Crystal Ball for J/Psi
    RooRealVar mass_Jpsi("mass_Jpsi","J/psi mass",3.097,3.00,3.20); 
    //mass_Jpsi.setConstant(kTRUE);
    RooRealVar sigma_Jpsi("sigma_Jpsi","J/psi resolution",0.08,0.01,0.1);
    RooGenericPdf mean_R("mean_R","J/psi mass","mass_Jpsi",RooArgSet(mass_Jpsi));
    RooGenericPdf sigma_R("sigma_R","J/psi resolution","sigma_Jpsi",RooArgSet(sigma_Jpsi));
    RooRealVar N_Jpsi("N_Jpsi","number of J/psi events",0.4*nEvents,0,nEvents);

    // Background
    RooRealVar lambda("lambda","background exp",-1.2,-10.,0.);
    RooRealVar N_bkg("N_bkg","number of background events",0.6*nEvents,0,nEvents);

    // Functions for fitting
    // J/psi:
    RooCBShape CB_left("CB_left","CB_left",fM,mass_Jpsi,sigma_Jpsi,alpha_L,n_L);
    RooCBShape CB_right("CB_right","CB_right",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_left,CB_right),RooArgList(frac));
    // Background:
    RooGenericPdf BkgPdf("BkgPdf","exp(fM*lambda)",RooArgSet(fM,lambda));

    // Create model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // Calculate the number of all J/psi events
    fM.setRange("WholeMassRange",fMCutLow,fMCutUpp);
    RooAbsReal *iDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("WholeMassRange"));
    // Integral of the normalized PDF, DSCB => will range from 0 to 1

    N_Jpsi_all[0] = iDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_all[1] = iDSCB->getVal()*N_Jpsi.getError();

    // Calculate the number of J/psi and bkg events with mass from 3.0 to 3.2 GeV/c^2 (around the J/psi peak)
    fM.setRange("JpsiMassRange",3.0,3.2);
    RooAbsReal *iBkg = BkgPdf.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_bkgr_peak[0] = iBkg->getVal()*N_bkg.getVal();
    N_bkgr_peak[1] = iBkg->getVal()*N_bkg.getError();
    RooAbsReal *iDSCB2 = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_Jpsi_peak[0] = iDSCB2->getVal()*N_Jpsi.getVal();
    N_Jpsi_peak[1] = iDSCB2->getVal()*N_Jpsi.getError();

    // ##########################################################
    // Plot the results
    // Draw Correlation Matrix
    TCanvas *cCorrMat = new TCanvas("cCorrMat","cCorrMat",700,600);
    InvMassFit_DrawCorrMatrix(cCorrMat,fResFit);

    // Draw histogram with fit results
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    InvMassFit_SetCanvas(c1,kFALSE);

    gStyle->SetEndErrorSize(0.);

    RooPlot* fFrameM = fM.frame(Title("Mass fit")); 
    fDataSet->plotOn(fFrameM,Name("fDataSet"),Binning(binM),MarkerStyle(kFullCircle),MarkerSize(1.),LineWidth(2));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DoubleSidedCB"),Components(DoubleSidedCB),LineColor(kBlack),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("BkgPdf"),Components(BkgPdf),LineColor(kRed),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DSCBAndBkgPdf"),LineColor(kBlue),LineWidth(3));
    // Vertical axis
    fFrameM->GetYaxis()->SetTitle(Form("Counts per %i MeV/#it{c}^{2}", BinSize));
    fFrameM->GetYaxis()->SetTitleSize(0.05);
    fFrameM->GetYaxis()->SetTitleOffset(1.1);
    fFrameM->GetYaxis()->SetLabelSize(0.05);
    fFrameM->GetYaxis()->SetLabelOffset(0.01);
    fFrameM->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    fFrameM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    fFrameM->GetXaxis()->SetTitleSize(0.05);
    fFrameM->GetXaxis()->SetLabelSize(0.05);
    fFrameM->GetXaxis()->SetDecimals(1);
    fFrameM->Draw();

    // Get chi2 
    Double_t chi2 = fFrameM->chiSquare("DSCBAndBkgPdf","fDataSet",fResFit->floatParsFinal().getSize()); // last argument = number of parameters
    Printf("********************");
    Printf("chi2/NDF = %.3f", chi2);
    Printf("NDF = %i", fResFit->floatParsFinal().getSize());
    Printf("chi2/NDF = %.3f/%i", chi2*fResFit->floatParsFinal().getSize(), fResFit->floatParsFinal().getSize());
    Printf("********************");   

    // -------------------------------------------------------------------------------- 
    // Legend1
    TLegend *l1 = new TLegend(0.09,0.76,0.3,0.935);
    //l1->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < %.1f", fYCut),"");
    // Print the pt cut
    if(opt == 0){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} > %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 1 || opt == 2){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} < %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    }
    l1->SetTextSize(0.040);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    TLegend *lTitle = new TLegend(0.325,0.88,0.95,0.935);
    lTitle->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    lTitle->SetTextSize(0.05);
    lTitle->SetBorderSize(0);
    lTitle->SetFillStyle(0);
    lTitle->Draw();

    // Legend2
    TLegend *l2 = new TLegend(0.52,0.29,0.95,0.87);
    l2->SetMargin(0.14);
    l2->AddEntry("DSCBAndBkgPdf","sum","L");
    l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l2->AddEntry("DoubleSidedCB","J/#psi signal","L");
    l2->AddEntry((TObject*)0,Form("#it{N}_{J/#psi} = %.0f #pm %.0f",N_Jpsi_all[0],N_Jpsi_all[1]),"");
    // Incoherent: lower precision:
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    // No incoherent: higher precision:
    } else if(opt == 1 || opt == 2){
        l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.4f #pm %.4f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#sigma = %.4f #pm %.4f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    }
    l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.2f", alpha_L.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{R} = %.2f", (-1)*(alpha_R.getVal())),"");
    l2->AddEntry("BkgPdf","background","L");
    l2->AddEntry((TObject*)0,Form("#lambda = %.2f #pm %.2f GeV^{-1}#it{c}^{2}",lambda.getVal(), lambda.getError()),"");
    l2->AddEntry((TObject*)0,"with #it{m}_{#mu#mu} #in (3.0,3.2) GeV/#it{c}^{2}:","");
    l2->AddEntry((TObject*)0,Form("#it{N}_{bkg} = %.0f #pm %.0f",N_bkgr_peak[0],N_bkgr_peak[1]),"");
    l2->SetTextSize(0.040); // was 0.042
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    TLegend *l3 = NULL;
    //if(!isNParInDSCBFixed)
    //{
        l3 = new TLegend(0.74,0.48,0.85,0.58);
        l3->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f", n_L.getVal()),"");
        l3->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f", n_R.getVal()),"");
        l3->SetTextSize(0.040); // was 0.042
        l3->SetBorderSize(0);
        l3->SetFillStyle(0);
        l3->Draw();
    //}

    // Print the numbers of events to text file
    ofstream outfile((str_out + ".txt").Data());
    outfile << "Signal in whole mass region 2.2 < m < 4.5 GeV:" << endl;
    outfile << "N_J/psi:\t" << N_Jpsi_all[0] << " pm " << N_Jpsi_all[1] << endl;
    outfile << "Mass region 3.0 < m < 3.2 GeV:" << endl;
    outfile << "N_J/psi:\t" << N_Jpsi_peak[0] << " pm " << N_Jpsi_peak[1] << endl;    
    outfile << "N_bkg:  \t" << N_bkgr_peak[0] << " pm " << N_bkgr_peak[1] << endl;
    outfile.close();
    Printf("*** Results printed to %s.***", (str_out + ".txt").Data());
    // Print the signal to text file
    ofstream outfile2((str_out + "_signal.txt").Data());
    outfile2 << N_Jpsi_all[0] << "\t" << N_Jpsi_all[1] << endl;
    outfile2.close();
    Printf("*** Results printed to %s.***", (str_out + "_signal.txt").Data());
    // Print the background to text file
    ofstream outfile3((str_out + "_bkg.txt").Data());
    outfile3 << N_bkgr_peak[0] << "\t" << N_bkgr_peak[1] << endl;
    outfile3.close();
    Printf("*** Results printed to %s.***", (str_out + "_bkg.txt").Data());

    // Print the plots
    c1->Print((str_out + ".pdf").Data());
    c1->Print((str_out + ".png").Data());
    cCorrMat->Print((str_out + "_cm.pdf").Data());
    cCorrMat->Print((str_out + "_cm.png").Data());    

    delete c1;
    delete cCorrMat;

    // ****************************************************************
    // Draw the result: paper figure
    if(opt == 3)
    {
        TCanvas *c2 = new TCanvas("c2","c2",900,800);
        c2->SetTopMargin(0.03);
        c2->SetBottomMargin(0.12);
        c2->SetRightMargin(0.03);
        c2->SetLeftMargin(0.14);
        c2->cd();
        fFrameM->GetYaxis()->SetTitleOffset(1.35);
        fFrameM->GetYaxis()->SetRangeUser(0.,235.);
        //fFrameM->GetYaxis()->SetNdivisions(505);
        fFrameM->Draw();

        Bool_t preliminary = kTRUE;
        Double_t xMin = 0.26;
        if(preliminary) xMin = 0.18;
        TLegend *lx = new TLegend(xMin,0.90,0.90,0.96);
        if(preliminary) lx->AddEntry((TObject*)0,"ALICE Preliminary, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
        else lx->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
        lx->SetMargin(0.);
        lx->SetTextSize(0.05);
        lx->SetBorderSize(0);
        lx->SetFillStyle(0);
        lx->Draw();

        TLegend *ly = new TLegend(0.56,0.57,0.92,0.85);
        ly->AddEntry((TObject*)0,"J/#psi #rightarrow #mu^{+} #mu^{-}","");
        ly->AddEntry((TObject*)0,"UPC, L_{int} = 232 #pm 7 #mub^{-1}","");
        ly->AddEntry((TObject*)0,"0.2 < #it{p}_{T} < 1.0 GeV/#it{c}","");
        ly->AddEntry((TObject*)0,"|#it{y}| < 0.8","");
        ly->AddEntry((TObject*)0,"#it{N}_{J/#psi} = 512 #pm 26","");
        //ly->AddEntry((TObject*)0,Form("#chi^{2}/dof = %.2f", chi2),"");
        ly->SetMargin(0.);
        ly->SetTextSize(0.042);
        ly->SetBorderSize(0);
        ly->SetFillStyle(0);
        ly->Draw();

        if(preliminary) {
            c2->Print("Results/" + str_subfolder + "_PreliminaryFigures/massFit.pdf");
            c2->Print("Results/" + str_subfolder + "_PreliminaryFigures/massFit.eps");
        } else c2->Print("Results/" + str_subfolder + "_PaperFigures/massFit.pdf");
        delete c2;
    }
    // ****************************************************************

    return;
}