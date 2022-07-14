// BinsThroughMassFit.C
// David Grund, Mar 20, 2022

// cpp headers
#include <fstream>
#include <chrono>  // sleep_for, sleep_until
#include <thread>  // nanoseconds, system_clock, seconds
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

using namespace RooFit;
using namespace std::this_thread;
using namespace std::chrono; 

// Main function
void BinsThroughMassFit_DoFit(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save = kFALSE, Int_t bin = -1);
// Support functions
void BinsThroughMassFit_SetCanvas(TCanvas *c, Bool_t bLogScale);

Double_t YieldJpsi_val = 0;
Double_t YieldJpsi_err = 0;

Double_t *PtBinsNew = NULL;
Double_t PtBinsNew_4bins[5] = {0.2, 0., 0., 0., 1.0};
Double_t PtBinsNew_5bins[6] = {0.2, 0., 0., 0., 0., 1.0};
Double_t YieldPerBin_val[5] = { 0 };
Double_t YieldPerBin_err[5] = { 0 };

Bool_t UniformBinYields = kFALSE;

void BinsThroughMassFit(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    // Adding ptStep = 0.001 GeV/c until a bin with sufficient signal (EvPerBin) is found

    if(nPtBins == 4) PtBinsNew = PtBinsNew_4bins;
    if(nPtBins == 5) PtBinsNew = PtBinsNew_5bins;

    Double_t ptStep = 0.001;
    Double_t CurrPtCutUpp = 0.2;

    Double_t EvTotal, EvTotalErr;
    // Load the total number of signal events
    TString str_ifs = "Results/" + str_subfolder + "InvMassFit/allbins/allbins_signal.txt";
    ifstream ifs;
    ifs.open(str_ifs.Data());
    while(!ifs.eof()){
        ifs >> EvTotal >> EvTotalErr;
    }
    ifs.close();
    Printf("Total number of events loaded: %.3f pn %.3f", EvTotal, EvTotalErr);

    Double_t EvPerBin_arr[5] = { 0 };
    Double_t correction(0.);
    if(!isPass3) correction = -1.0;
    else         correction = -0.8;
    if(UniformBinYields) // if uniform bin yields
    {
        Double_t EvPerBin = (EvTotal / (Double_t)nPtBins) + correction;
        for(Int_t i = 0; i < nPtBins; i++) EvPerBin_arr[i] = EvPerBin;
        Printf("Optimal yields per bin: %.2f", EvPerBin);
    }
    else // non-uniform bin yields: aprox 125, 125, 80, 80, 80
    {
        if(nPtBins == 5)
        {
            Double_t nEvOpt4bins = (EvTotal / (Double_t)(nPtBins-1));
            EvPerBin_arr[0] = nEvOpt4bins + correction;
            EvPerBin_arr[1] = nEvOpt4bins + correction;
            EvPerBin_arr[2] = (EvTotal - 2 * nEvOpt4bins) / 3 - 0.3;
            EvPerBin_arr[3] = (EvTotal - 2 * nEvOpt4bins) / 3 - 0.3;
            EvPerBin_arr[4] = (EvTotal - 2 * nEvOpt4bins) / 3 - 0.3;
            Printf("Optimal yields per bins: %.2f, %.2f, %.2f, %.2f, %.2f", 
                EvPerBin_arr[0], EvPerBin_arr[1], EvPerBin_arr[2], EvPerBin_arr[3], EvPerBin_arr[4]);
        }
        else return;
    }
    
    // Small delay to be able to read the console
    sleep_until(system_clock::now() + seconds(3));

    // Print the results
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "BinsThroughMassFit/");
    TString name_1 = Form("Results/" + str_subfolder + "BinsThroughMassFit/log_%ibins.txt", nPtBins);
    ofstream out_1(name_1.Data());
    out_1 << Form("Using pT step %.3f GeV/c.\n", ptStep);

    for(Int_t i = 0; i < nPtBins-1; i++){
        // While the yield of J/psi candidates in the current bin is smaller than then optimal one
        while(YieldJpsi_val <= EvPerBin_arr[i]){
            CurrPtCutUpp += ptStep;
            BinsThroughMassFit_DoFit(PtBinsNew[i], CurrPtCutUpp);
            out_1 << Form("(%.3f, %.3f): %.3f\n", PtBinsNew[i], CurrPtCutUpp, YieldJpsi_val);
        }
        PtBinsNew[i+1] = CurrPtCutUpp;
        out_1 << Form("Bin %i defined as (%.3f, %.3f)\n", (i+1), PtBinsNew[i], PtBinsNew[i+1]);
        out_1 << Form("Going to next bin...\n");
        YieldJpsi_val = 0;
    }
    out_1 << Form("Bin %i defined as (%.3f, %.3f)\n", nPtBins, PtBinsNew[nPtBins-1], PtBinsNew[nPtBins]);

    out_1.close();
    Printf("*** Results printed to %s.***", name_1.Data());

    // Do fits in the four/five calculated bins
    for(Int_t i = 0; i < nPtBins; i++){
        BinsThroughMassFit_DoFit(PtBinsNew[i], PtBinsNew[i+1], kTRUE, i+1);
        YieldPerBin_val[i] = YieldJpsi_val;
        YieldPerBin_err[i] = YieldJpsi_err;
    } 

    // Print the bin boundaries separately to a file from which they will be read
    TString name_2 = Form("Results/" + str_subfolder + "BinsThroughMassFit/%ibins_defined.txt", nPtBins);
    ofstream out_2(name_2.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        out_2 << Form("%.3f \t%.3f \t%.3f\n", PtBinsNew[i], YieldPerBin_val[i], YieldPerBin_err[i]);
    }
    out_2 << Form("%.3f\n", PtBinsNew[nPtBins]);
    out_2.close();
    Printf("*** Results printed to %s.***", name_2.Data());

    return;
}

void BinsThroughMassFit_DoFit(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save, Int_t bin){
    // Fit the invariant mass distribution using Double-sided CB function
    // Fix the values of the tail parameters to MC values
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fYCut      = 0.80;
    Double_t fMCutLow   = 2.2;
    Double_t fMCutUpp   = 4.5;

    sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);

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
    TFile *f_in = new TFile("Trees/" + str_subfolder + "InvMassFit/InvMassFit.root"); // created in InvMassFit.c
    TTree *t_in = NULL;
    f_in->GetObject("tIncEnrSample",t_in);
        
    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*t_in));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Crystal Ball parameters from MC (to be fixed)
    Double_t fAlpha_L;
    Double_t fAlpha_R;
    Double_t fN_L;
    Double_t fN_R;

    char name[20];
    Double_t values[4];
    Double_t errors[4];

    TString path = "Results/" + str_subfolder + "InvMassFit_MC/inc/inc.txt";

    ifstream f_text_in;
    f_text_in.open(path.Data());
    if(f_text_in.fail()){
        Printf("\n");
        Printf("*** Warning! ***");
        Printf("*** MC values for tail parameters not found. Terminating... *** \n");
        return;
    } else {
        Int_t i_line = 0;
        while(!f_text_in.eof()){
            f_text_in >> name >> values[i_line] >> errors[i_line];
            i_line++;
        }
        f_text_in.close();
    }
    fAlpha_L = values[0];
    fAlpha_R = values[1];
    fN_L = values[2];
    fN_R = values[3];

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,20.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,20.);
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

    // Create Model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // Calculate the number of J/psi events
    Double_t N_Jpsi_out[2];
    fM.setRange("WholeMassRange",fMCutLow,fMCutUpp);
    RooAbsReal *intDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("WholeMassRange"));
    // Integral of the normalized PDF, DSCB => will range from 0 to 1

    N_Jpsi_out[0] = intDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_out[1] = intDSCB->getVal()*N_Jpsi.getError();
    YieldJpsi_val = N_Jpsi_out[0];
    YieldJpsi_err = N_Jpsi_out[1];

    // Calculate the number of bkg events with mass in 3.0 to 3.2 GeV/c^2
    Double_t N_bkg_out[2];
    fM.setRange("JpsiMassRange",3.0,3.2);
    RooAbsReal *iBkg = BkgPdf.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_bkg_out[0] = iBkg->getVal()*N_bkg.getVal();
    N_bkg_out[1] = iBkg->getVal()*N_bkg.getError();

    // ##########################################################
    // Plot the results
    // Draw histogram with fit results
    TCanvas *cHist = NULL;
    
    if(save){
        cHist = new TCanvas("cHist","cHist",800,600);
        BinsThroughMassFit_SetCanvas(cHist,kFALSE);

        RooPlot* fFrameM = fM.frame(Title("Mass fit")); 
        fDataSet->plotOn(fFrameM,Name("fDataSet"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
        DSCBAndBkgPdf.plotOn(fFrameM,Name("DoubleSidedCB"),Components(DoubleSidedCB),LineColor(kBlack),LineStyle(kDashed),LineWidth(3));
        DSCBAndBkgPdf.plotOn(fFrameM,Name("BkgPdf"),Components(BkgPdf),LineColor(kRed),LineStyle(kDashed),LineWidth(3));
        DSCBAndBkgPdf.plotOn(fFrameM,Name("DSCBAndBkgPdf"),LineColor(215),LineWidth(3));
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
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
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
        l2->AddEntry((TObject*)0,Form("#it{N}_{J/#psi} = %.0f #pm %.0f",N_Jpsi_out[0],N_Jpsi_out[1]),"");
        l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.2f", alpha_L.getVal()),"");
        l2->AddEntry((TObject*)0,Form("#alpha_{R} = %.2f", (-1)*(alpha_R.getVal())),"");
        l2->AddEntry("BkgPdf","background","L");
        l2->AddEntry((TObject*)0,Form("#lambda = %.2f #pm %.2f GeV^{-1}#it{c}^{2}",lambda.getVal(), lambda.getError()),"");
        l2->AddEntry((TObject*)0,"with #it{m}_{#mu#mu} #in (3.0,3.2) GeV/#it{c}^{2}:","");
        l2->AddEntry((TObject*)0,Form("#it{N}_{bkg} = %.0f #pm %.0f",N_bkg_out[0],N_bkg_out[1]),"");
        l2->SetTextSize(0.040); // was 0.042
        l2->SetBorderSize(0);
        l2->SetFillStyle(0);
        l2->Draw();

        TLegend *l3 = NULL;
        if(!isNParInDSCBFixed)
        {
            l3 = new TLegend(0.74,0.48,0.85,0.58);
            l3->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f", n_L.getVal()),"");
            l3->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f", n_R.getVal()),"");
            l3->SetTextSize(0.040); // was 0.042
            l3->SetBorderSize(0);
            l3->SetFillStyle(0);
            l3->Draw();
        }

        // Prepare path
        TString str = Form("Results/" + str_subfolder + "BinsThroughMassFit/bin%i", bin);
        // Print the plots
        cHist->Print((str + ".pdf").Data());
        cHist->Print((str + ".png").Data()); 
    }

    f_in->Close();

    return;
}

void BinsThroughMassFit_SetCanvas(TCanvas *c, Bool_t bLogScale){

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    if(bLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    return;
}