// PtFit_WithoutBkg.c
// David Grund, Sep 27, 2021

// roofit headers
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooTFnBinding.h"
#include "RooRealSumPdf.h"
// my headers
#include "AnalysisManager.h"
#include "PtFit_Utilities.h"

//#######################################
// Options to set:
// *** iCohJShape ***
// 0 => classic histogram from STARlight (R_A = 6.624 fm)
// 1 => histogram from STARlight generated with R_A = 7.53 fm
// 2 => fit using "Gaussian shape" pT * exp(-b * pT^2)
// 3 => fit using STARlight formfactor, R_A left free
// *** To study which value of R_A is optimal: ***
// 1001 => R_A = 6.60 fm
// 1002 => R_A = 6.70 fm
// 1003 => R_A = 6.80 fm
// 1004 => R_A = 6.90 fm
// 1005 => R_A = 7.00 fm
// 1006 => R_A = 7.10 fm
// 1007 => R_A = 7.20 fm
// 1008 => R_A = 7.30 fm
// 1009 => R_A = 7.40 fm
// 1010 => R_A = 7.50 fm
// 1011 => R_A = 7.60 fm
// 1012 => R_A = 7.70 fm
// 1013 => R_A = 7.80 fm
Int_t iNormFD = 0;
// 0 => taken from STARlight (f_D^coh ~ 7%), the same way as Roman
// 1 => ratio of coherent cross sections fixed to the value from Michal's paper
//#######################################

// Main function
void PrepareDataTree();
void SubtractBackground();
void DoPtFitNoBkg(Int_t iCohJShape);
void DoInvMassFitMain(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save = kFALSE, Int_t bin = -1);
void DrawCorrelationMatrixModified(TCanvas *cCM, RooFitResult* ResFit, Int_t iCohJShape);

TString OutputPtFitWithoutBkg = "Results/PtFitWithoutBkg/";
TString OutputTrees = "Trees/PtFit/";

Double_t N_Jpsi_val = 0;
Double_t N_Jpsi_err = 0;
Double_t N_Bkgr_val = 0;
Double_t N_Bkgr_err = 0;

const Double_t hbarC = 0.197; // GeV*fm
const Double_t PbAtomicNumber = 208.; // Pb

Double_t VMD_model(Double_t *pT, Double_t *par)
{
    // input and parameters
    Double_t q = pT[0];     // q = sqrt of transferred momentum |t|, unit is GeV 
    if (q == 0.) return 0.; // protection against floating point exception
    //Double_t norm = par[2]; // normalization
    Double_t R_A = par[0];  // fm (radius of the lead nucleus, SL: R_A = 6.62 fm)
    Double_t a = par[1];    // fm (SL: a = 0.7 fm)
    // computation of the form factor
    Double_t K = 4.0 * TMath::Pi() * hbarC * hbarC * hbarC / PbAtomicNumber; // Gev^3*fm^3
    Double_t qR_A = q * R_A / hbarC;                        // unit dimension is zero
    Double_t F1 = K / (q * q * q);                          // fm^3
    Double_t F2 = TMath::Sin(qR_A) - qR_A*TMath::Cos(qR_A); // unit dimension is zero
    Double_t F3 = 1 + (a * a * q * q / (hbarC*hbarC));      // unit dimension is zero
    if (F3 == 0.) return 0.; // protection against floating point exception 
    Double_t FF = F1 * F2 / F3; 

    return q * FF * FF; // q (= pT) from the Jacobian ! (see a photo from Dec 2)
}
Double_t Zero_func(Double_t pT)
{
    return 0.;
}

void PtFit_WithoutBkg()
{
    //PrepareDataTree();

    SetPtBins(2);

    SubtractBackground();

    PreparePDFs_MC();

    PreparePDFs_modRA_CohJ();

    DoPtFitNoBkg(0);
    DoPtFitNoBkg(1);
    DoPtFitNoBkg(2);
    DoPtFitNoBkg(3);

    for(Int_t i = 1001; i < 1014; i++) DoPtFitNoBkg(i);

    return;
}

void DoPtFitNoBkg(Int_t iCohJShape)
{
    SetPtBinning();

    // Load the file with PDFs
    TFile *file = TFile::Open(Form("%sPDFs_MC_Binning%i.root", OutputPDFs.Data(), BinningOpt),"read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    TList *list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    // Load histograms
    // 1) kCohJpsiToMu
    TH1D *hCohJ = NULL;
    if(iCohJShape == 0 || iCohJShape == 2 || iCohJShape == 3){

        hCohJ = (TH1D*)list->FindObject(NamesPDFs[0].Data());
        if(hCohJ) Printf("Histogram %s loaded.", hCohJ->GetName());

    } else if(iCohJShape == 1 || iCohJShape > 1000){

        TString str_name = "";
        if(bStopWeight) str_name = Form("%sPDFs_MC_modRA_Binning%i_StopWeight.root", OutputPDFs.Data(), BinningOpt);
        else            str_name = Form("%sPDFs_MC_modRA_Binning%i.root", OutputPDFs.Data(), BinningOpt);

        TFile *f_modRA = TFile::Open(str_name.Data(),"read");
        if(f_modRA) Printf("Input file %s loaded.", f_modRA->GetName()); 

        TList *l_modRA = (TList*) f_modRA->Get("HistList");
        if(l_modRA) Printf("List %s loaded.", l_modRA->GetName()); 

        Double_t R_A = 0;
        if(iCohJShape == 1) R_A = 7.53;
        if(iCohJShape > 1000) R_A = 6.6 + (Double_t)(iCohJShape-1001) * 0.1;

        hCohJ = (TH1D*)l_modRA->FindObject(Form("hCohJ_modRA_%.2f", R_A));
        if(hCohJ) Printf("Histogram %s loaded.", hCohJ->GetName());

        f_modRA->Close();
    } 

    // 2) kIncohJpsiToMu
    TH1D *hIncJ = (TH1D*)list->FindObject(NamesPDFs[1].Data());
    if(hIncJ) Printf("Histogram %s loaded.", hIncJ->GetName());

    // 3) kCohPsi2sToMuPi
    TH1D *hCohP = (TH1D*)list->FindObject(NamesPDFs[2].Data());
    if(hCohP) Printf("Histogram %s loaded.", hCohP->GetName());

    // 4) kincohPsi2sToMuPi
    TH1D *hIncP = (TH1D*)list->FindObject(NamesPDFs[3].Data());
    if(hIncP) Printf("Histogram %s loaded.", hIncP->GetName());

    // 6) Dissociative
    TH1D *hDiss = (TH1D*)list->FindObject(NamesPDFs[5].Data());
    if(hDiss) Printf("Histogram %s loaded.", hDiss->GetName());

    // Definition of roofit variables
    RooRealVar fPt("fPt", "fPt", fPtLow, fPtUpp);
    RooArgSet fSetOfVariables(fPt);

    // Create PDFs
    // 1) kCohJpsiToMu
    RooDataHist DHisCohJ("DHisCohJ","DHisCohJ",fPt,hCohJ);
    RooHistPdf  hPDFCohJ("hPDFCohJ","hPDFCohJ",fSetOfVariables,DHisCohJ,0);

    // iCohJShape == 2 => Fit with pT * exp(-b * pT^2)
    RooRealVar par_b("par_b","b",100.,1.,500.);
    RooGenericPdf *gPDFCohJ = new RooGenericPdf("gPDFCohJ","gPDFCohJ","fPt*exp(-par_b*pow(fPt,2))",RooArgSet(fPt,par_b));
    
    // iCohJShape == 3 => Fit with STARlight form factor and R_A left free
    // https://root.cern.ch/download/doc/RooFit_Users_Manual_2.91-33.pdf
    // https://root-forum.cern.ch/t/bind-tf1-into-roofit-pdf-using-bindpdf/26623
    // https://root-forum.cern.ch/t/defining-a-roogenericpdf-from-a-rooaddition-or-any-rooabsreal/20483 
    // https://root-forum.cern.ch/t/roofit-fitting-a-tf1-binded-function/7656

    // Create RooAbsReal from a TF1 function defined using C++ function VMD_model()
    TF1 *fFormFactorSL = new TF1("fFormFactorSL",VMD_model,fPtLow,fPtUpp,2); // 2 = number of parameters
    RooRealVar R_A("R_A","R_A",6.62,1.,12.); // 6.624 fm = original SL value
    //R_A.setConstant(kTRUE);
    RooRealVar a("a","a", 0.7, 0.7, 0.7); // 0.7 fm = SL value
    a.setConstant(kTRUE);
    //RooRealVar norm("norm","norm",1e3,0,1e5);
    RooAbsReal *absRealCohJ = bindFunction(fFormFactorSL, fPt, RooArgList(R_A, a));
    // Create zero function
    TF1 *fZero = new TF1("fZero","Zero_func(x)",fPtLow,fPtUpp);
    RooAbsReal *absRealZero = bindFunction(fZero, fPt);
    // Convert RooAbsReal to PDF: use RooRealSumPdf
    RooRealVar coef("coef","coef",1.,1.,1.);
    coef.setConstant(kTRUE);
    RooRealSumPdf *sumPdfCohJ = new RooRealSumPdf("sumPdfCohJ","sumPdfCohJ",RooArgList(*absRealCohJ,*absRealZero),RooArgList(coef));

    // 2) kIncohJpsiToMu
    RooDataHist DHisIncJ("DHisIncJ","DHisIncJ",fPt,hIncJ);
    RooHistPdf  hPDFIncJ("hPDFIncJ","hPDFIncJ",fSetOfVariables,DHisIncJ,0);

    // 3) kCohPsi2sToMuPi
    RooDataHist DHisCohP("DHisCohP","DHisCohP",fPt,hCohP);
    RooHistPdf  hPDFCohP("hPDFCohP","hPDFCohP",fSetOfVariables,DHisCohP,0);

    // 4) kincohPsi2sToMuPi
    RooDataHist DHisIncP("DHisIncP","DHisIncP",fPt,hIncP);
    RooHistPdf  hPDFIncP("hPDFIncP","hPDFIncP",fSetOfVariables,DHisIncP,0);

    // 6) Dissociative
    RooDataHist DHisDiss("DHisDiss","DHisDiss",fPt,hDiss);
    RooHistPdf  hPDFDiss("hPDFDiss","hPDFDiss",fSetOfVariables,DHisDiss,0);

    // Close the file with PDFs
    file->Close();

    // Get the binned dataset
    file = TFile::Open(Form("%sJpsiSignalNoBkg_Binning%i.root", OutputTrees.Data(), BinningOpt), "read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    TH1D *hData = (TH1D*)list->FindObject("hNJpsiBins");
    if(hData) Printf("Histogram %s loaded.", hData->GetName());

    RooDataHist DHisData("DHisData","DHisData",fPt,hData);
    Printf("Binned data with background subtracted loaded.");
    // Calculate the number of entries
    Double_t N_all = 0;
    for(Int_t i = 1; i <= hData->GetNbinsX(); i++){
        N_all += hData->GetBinContent(i);
    }
    Printf("Data contain %.0f entries in %i bins.", N_all, DHisData.numEntries());
    file->Close();

    // 8) Create the model for fitting
    // 8.1) Normalizations:
    RooRealVar NCohJ("NCohJ","Number of coh J/psi events", 0.90*N_all,0.50*N_all,1.0*N_all);
    RooRealVar NIncJ("NIncJ","Number of inc J/psi events", 0.05*N_all,0.01*N_all,0.3*N_all);
    RooRealVar NDiss("NDiss","Number of dis J/psi events", 0.05*N_all,0.01*N_all,0.3*N_all);

    // Load the values of fD coefficients
    Double_t fDCohCh, fDCohChErr, fDIncCh, fDIncChErr, fDCohNe, fDCohNeErr, fDIncNe, fDIncNeErr;
    char ch[8];
    ifstream file_in;
    file_in.open("Results/FeedDown/FeedDown_PtFit_ratMeas.txt");
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> fDCohCh >> fDCohChErr >> fDIncCh >> fDIncChErr >> fDCohNe >> fDCohNeErr >> fDIncNe >> fDIncNeErr;
            i++;   
        } 
        file_in.close();
        Printf("fD coefficients loaded.");
    } else {
        Printf("fD coefficients missing. Terminating...");
        return;
    }
    Double_t fDCoh = (fDCohCh + fDCohNe) / 100;
    Double_t fDInc = (fDIncCh + fDIncNe) / 100;
    Printf("********************");
    Printf("fD_coh = %.4f", fDCoh);
    Printf("fD_inc = %.4f", fDInc);
    Printf("********************");
    
    RooGenericPdf NCohP("NCohP","Number of coh FD events",Form("NCohJ*%.4f", fDCoh),RooArgSet(NCohJ));
    RooGenericPdf NIncP("NIncP","Number of inc FD events",Form("NIncJ*%.4f", fDInc),RooArgSet(NIncJ));

    // 8.2) The model:
    RooAddPdf *Mod = NULL;
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000){
        Mod = new RooAddPdf("Mod","Sum of all PDFs",
            RooArgList(hPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
            RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
        );
    } else if(iCohJShape == 2){
        Mod = new RooAddPdf("Mod","Sum of all PDFs",
            RooArgList(*gPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
            RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
        );        
    } else if(iCohJShape == 3){
        Mod = new RooAddPdf("Mod","Sum of all PDFs",
            RooArgList(*sumPdfCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
            RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
        );        
    }

    // 9) Perform fitting
    RooFitResult *ResFit = Mod->fitTo(DHisData,Extended(kTRUE),Save(),Range(fPtLow,fPtUpp));

    // ###############################################################################################################
    // ###############################################################################################################
    // 10) Output to text file
    TString *str = NULL;
    str = new TString(Form("%s%ibins_Binn%i_CohSh%i", OutputPtFitWithoutBkg.Data(), nPtBins, BinningOpt, iCohJShape));
    if(iCohJShape > 1000 && !bStopWeight){
        Double_t R_A = 6.6 + (Double_t)(iCohJShape-1001) * 0.1;        
        str = new TString(Form("%sOptimalRA/WeightOverAll/%ibins_coh_modRA_%.2f", OutputPtFitWithoutBkg.Data(), nPtBins, R_A));
    } 
    if(iCohJShape > 1000 && bStopWeight){
        Double_t R_A = 6.6 + (Double_t)(iCohJShape-1001) * 0.1;        
        str = new TString(Form("%sOptimalRA/StopWeight/%ibins_coh_modRA_%.2f", OutputPtFitWithoutBkg.Data(), nPtBins, R_A));
    } 
    // Integrals of the PDFs in the whole pt range 
    fPt.setRange("fPtAll",0.0,2.0);
    RooAbsReal *fN_CohJ_all = NULL;
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000) fN_CohJ_all = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    else if(iCohJShape == 2) fN_CohJ_all = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    else if(iCohJShape == 3) fN_CohJ_all = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_IncJ_all = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_Diss_all = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));  
    RooAbsReal *fN_CohP_all = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtAll")); 
    RooAbsReal *fN_IncP_all = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    // Integrals of the PDFs in the incoherent-enriched sample (IES)
    fPt.setRange("fPtIES",0.2,2.0);
    RooAbsReal *fN_CohJ_ies = NULL;
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000) fN_CohJ_ies = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    else if(iCohJShape == 2) fN_CohJ_ies = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    else if(iCohJShape == 3) fN_CohJ_ies = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_IncJ_ies = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_Diss_ies = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    RooAbsReal *fN_CohP_ies = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtIES")); 
    RooAbsReal *fN_IncP_ies = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    // Integrals of the PDFs with 0.2 < pt < 1.0 GeV/c
    fPt.setRange("fPtTo1",0.2,1.0);
    RooAbsReal *fN_CohJ_to1 = NULL;
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000) fN_CohJ_to1 = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    else if(iCohJShape == 2) fN_CohJ_to1 = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    else if(iCohJShape == 3) fN_CohJ_to1 = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    RooAbsReal *fN_IncJ_to1 = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    RooAbsReal *fN_Diss_to1 = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    RooAbsReal *fN_CohP_to1 = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtTo1")); 
    RooAbsReal *fN_IncP_to1 = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    // Number of events in the whole pt range
        // values
        Double_t N_CohJ_all_val = fN_CohJ_all->getVal()*NCohJ.getVal();
        Double_t N_IncJ_all_val = fN_IncJ_all->getVal()*NIncJ.getVal();
        Double_t N_Diss_all_val = fN_Diss_all->getVal()*NDiss.getVal();
        Double_t N_CohP_all_val = fN_CohP_all->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_all_val = fN_IncP_all->getVal()*fDInc*NIncJ.getVal();
        // errors
        Double_t N_CohJ_all_err = fN_CohJ_all->getVal()*NCohJ.getError();
        Double_t N_IncJ_all_err = fN_IncJ_all->getVal()*NIncJ.getError();
        Double_t N_Diss_all_err = fN_Diss_all->getVal()*NDiss.getError();
        Double_t N_CohP_all_err = fN_CohP_all->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_all_err = fN_IncP_all->getVal()*fDInc*NIncJ.getError();
    // Number of events with 0.2 < pt < 2.0 GeV/c
        // values
        Double_t N_CohJ_ies_val = fN_CohJ_ies->getVal()*NCohJ.getVal();
        Double_t N_IncJ_ies_val = fN_IncJ_ies->getVal()*NIncJ.getVal();
        Double_t N_Diss_ies_val = fN_Diss_ies->getVal()*NDiss.getVal();
        Double_t N_CohP_ies_val = fN_CohP_ies->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_ies_val = fN_IncP_ies->getVal()*fDInc*NIncJ.getVal();
        // errors
        Double_t N_CohJ_ies_err = fN_CohJ_ies->getVal()*NCohJ.getError();
        Double_t N_IncJ_ies_err = fN_IncJ_ies->getVal()*NIncJ.getError();
        Double_t N_Diss_ies_err = fN_Diss_ies->getVal()*NDiss.getError();
        Double_t N_CohP_ies_err = fN_CohP_ies->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_ies_err = fN_IncP_ies->getVal()*fDInc*NIncJ.getError();  
    // Number of events with 0.2 < pt < 1.0 GeV/c
        // values
        Double_t N_CohJ_to1_val = fN_CohJ_to1->getVal()*NCohJ.getVal();
        Double_t N_IncJ_to1_val = fN_IncJ_to1->getVal()*NIncJ.getVal();
        Double_t N_Diss_to1_val = fN_Diss_to1->getVal()*NDiss.getVal();
        Double_t N_CohP_to1_val = fN_CohP_to1->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_to1_val = fN_IncP_to1->getVal()*fDInc*NIncJ.getVal();
        // errors
        Double_t N_CohJ_to1_err = fN_CohJ_to1->getVal()*NCohJ.getError();
        Double_t N_IncJ_to1_err = fN_IncJ_to1->getVal()*NIncJ.getError();
        Double_t N_Diss_to1_err = fN_Diss_to1->getVal()*NDiss.getError();
        Double_t N_CohP_to1_err = fN_CohP_to1->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_to1_err = fN_IncP_to1->getVal()*fDInc*NIncJ.getError();
    // Total fC correction (for the whole IES, with 0.2 < pt < 2.0 GeV/c)
    Double_t fC = N_CohJ_ies_val / (N_IncJ_ies_val + N_Diss_ies_val) * 100;
    Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_ies_err,2) + TMath::Power(N_Diss_ies_err,2));
    Double_t fC_err = fC * TMath::Sqrt(TMath::Power((N_CohJ_ies_err/N_CohJ_ies_val),2) + TMath::Power((denominator_err/(N_IncJ_ies_val + N_Diss_ies_val)),2));
    // Total fD correction (for the whole IES, with 0.2 < pt < 2.0 GeV/c)
    Double_t fDCoh_val = N_CohP_ies_val / (N_IncJ_ies_val + N_Diss_ies_val) * 100;
    Double_t fDInc_val = N_IncP_ies_val / (N_IncJ_ies_val + N_Diss_ies_val) * 100;
    Double_t fDCoh_err = fDCoh_val * TMath::Sqrt(TMath::Power((N_CohP_ies_err/N_CohP_ies_val),2) + TMath::Power((denominator_err/(N_IncJ_ies_val + N_Diss_ies_val)),2));
    Double_t fDInc_err = fDInc_val * TMath::Sqrt(TMath::Power((N_IncP_ies_err/N_IncP_ies_val),2) + TMath::Power((denominator_err/(N_IncJ_ies_val + N_Diss_ies_val)),2));
    Double_t fD_val = fDCoh_val + fDInc_val;
    Double_t fD_err = TMath::Sqrt(TMath::Power(fDCoh_err, 2) + TMath::Power(fDInc_err, 2));
    // Integrals of the PDFs in the pt bins
    RooAbsReal *fN_CohJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_Diss_bins[nPtBins] = { NULL };
    RooAbsReal *fN_CohP_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncP_bins[nPtBins] = { NULL };
    // values
    Double_t N_CohJ_bins_val[nPtBins] = { 0 };
    Double_t N_IncJ_bins_val[nPtBins] = { 0 };
    Double_t N_Diss_bins_val[nPtBins] = { 0 };
    Double_t N_CohP_bins_val[nPtBins] = { 0 };
    Double_t N_IncP_bins_val[nPtBins] = { 0 };
    // errors
    Double_t N_CohJ_bins_err[nPtBins] = { 0 };
    Double_t N_IncJ_bins_err[nPtBins] = { 0 };
    Double_t N_Diss_bins_err[nPtBins] = { 0 };
    Double_t N_CohP_bins_err[nPtBins] = { 0 };
    Double_t N_IncP_bins_err[nPtBins] = { 0 };
    // coefficients
    Double_t fC_bins_val[nPtBins] = { 0 };
    Double_t fC_bins_err[nPtBins] = { 0 };
    Double_t fDCoh_bins_val[nPtBins] = { 0 };
    Double_t fDCoh_bins_err[nPtBins] = { 0 };
    Double_t fDInc_bins_val[nPtBins] = { 0 };
    Double_t fDInc_bins_err[nPtBins] = { 0 };
    Double_t fD_bins_val[nPtBins] = { 0 };
    Double_t fD_bins_err[nPtBins] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        fPt.setRange(Form("fPtBin%i",i+1), ptBoundaries[i], ptBoundaries[i+1]);
        Printf("Now calculating for bin %i, (%.3f, %.3f) GeV", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000) fN_CohJ_bins[i] = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        else if(iCohJShape == 2) fN_CohJ_bins[i] = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        else if(iCohJShape == 3) fN_CohJ_bins[i] = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_IncJ_bins[i] = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_Diss_bins[i] = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_CohP_bins[i] = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_IncP_bins[i] = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));  
        // values  
        N_CohJ_bins_val[i] = fN_CohJ_bins[i]->getVal()*NCohJ.getVal();
        N_IncJ_bins_val[i] = fN_IncJ_bins[i]->getVal()*NIncJ.getVal();
        N_Diss_bins_val[i] = fN_Diss_bins[i]->getVal()*NDiss.getVal();
        N_CohP_bins_val[i] = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getVal();
        N_IncP_bins_val[i] = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getVal();
        // errors
        N_CohJ_bins_err[i] = fN_CohJ_bins[i]->getVal()*NCohJ.getError();
        N_IncJ_bins_err[i] = fN_IncJ_bins[i]->getVal()*NIncJ.getError();
        N_Diss_bins_err[i] = fN_Diss_bins[i]->getVal()*NDiss.getError();
        N_CohP_bins_err[i] = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getError();
        N_IncP_bins_err[i] = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getError();
        // fC correction        
        fC_bins_val[i] = N_CohJ_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]) * 100;
        Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_bins_err[i],2) + TMath::Power(N_Diss_bins_err[i],2));
        if(N_CohJ_bins_val[i] != 0) fC_bins_err[i] = fC_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohJ_bins_err[i]/N_CohJ_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));    
        else fC_bins_err[i] = 0.;
        fDCoh_bins_val[i] = N_CohP_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]) * 100;
        fDInc_bins_val[i] = N_IncP_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]) * 100;
        if(N_CohP_bins_val[i] != 0) fDCoh_bins_err[i] = fDCoh_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohP_bins_err[i]/N_CohP_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));
        else fDCoh_bins_err[i] = 0.;
        if(N_IncP_bins_val[i] != 0) fDInc_bins_err[i] = fDInc_bins_val[i] * TMath::Sqrt(TMath::Power((N_IncP_bins_err[i]/N_IncP_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));
        else fDInc_bins_val[i] = 0.;
        fD_bins_val[i] = fDCoh_bins_val[i] + fDInc_bins_val[i];
        fD_bins_err[i] = TMath::Sqrt(TMath::Power(fDCoh_bins_err[i], 2) + TMath::Power(fDInc_bins_err[i], 2));
    }
    Double_t sum_all = N_CohJ_all_val + N_IncJ_all_val + N_CohP_all_val + N_IncP_all_val + N_Diss_all_val;
    Double_t sum_ies = N_CohJ_ies_val + N_IncJ_ies_val + N_CohP_ies_val + N_IncP_ies_val + N_Diss_ies_val;
    Double_t sum_to1 = N_CohJ_to1_val + N_IncJ_to1_val + N_CohP_to1_val + N_IncP_to1_val + N_Diss_to1_val;
    // Print to text file
    ofstream outfile((*str + ".txt").Data());
    outfile << std::fixed << std::setprecision(2);
    outfile << Form("Dataset contains %.0f events.\n***\n", N_all);
    outfile << "pT range\t\tNCohJ \terr \tNIncJ \terr \tNDiss \terr \tNCohP \terr \tNIncP \terr \tfC \terr \tfD coh\terr \tfD inc\terr \tfD \terr\n";
    outfile << "(0.000, 2.000) GeV/c\t"
            << N_CohJ_all_val << "\t" << N_CohJ_all_err << "\t" 
            << N_IncJ_all_val << "\t" << N_IncJ_all_err << "\t" 
            << N_Diss_all_val << "\t" << N_Diss_all_err << "\t" 
            << N_CohP_all_val << "\t" << N_CohP_all_err << "\t" 
            << N_IncP_all_val << "\t" << N_IncP_all_err << "\n";
    outfile << "(0.200, 2.000) GeV/c\t"
            << N_CohJ_ies_val << "\t" << N_CohJ_ies_err << "\t" 
            << N_IncJ_ies_val << "\t" << N_IncJ_ies_err << "\t" 
            << N_Diss_ies_val << "\t" << N_Diss_ies_err << "\t" 
            << N_CohP_ies_val << "\t" << N_CohP_ies_err << "\t" 
            << N_IncP_ies_val << "\t" << N_IncP_ies_err << "\t"
            << fC << "\t" << fC_err << "\t" 
            << fDCoh_val << "\t" << fDCoh_err << "\t" 
            << fDInc_val << "\t" << fDInc_err << "\t"
            << fD_val << "\t" << fD_err << "\n";
    outfile << "(0.200, 1.000) GeV/c\t"
            << N_CohJ_to1_val << "\t" << N_CohJ_to1_err << "\t" 
            << N_IncJ_to1_val << "\t" << N_IncJ_to1_err << "\t" 
            << N_Diss_to1_val << "\t" << N_Diss_to1_err << "\t"  
            << N_CohP_to1_val << "\t" << N_CohP_to1_err << "\t" 
            << N_IncP_to1_val << "\t" << N_IncP_to1_err << "\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << Form("(%.3f, %.3f) GeV/c\t", ptBoundaries[i], ptBoundaries[i+1])
                << N_CohJ_bins_val[i] << "\t" << N_CohJ_bins_err[i] << "\t"
                << N_IncJ_bins_val[i] << "\t" << N_IncJ_bins_err[i] << "\t"
                << N_Diss_bins_val[i] << "\t" << N_Diss_bins_err[i] << "\t"
                << N_CohP_bins_val[i] << "\t" << N_CohP_bins_err[i] << "\t" 
                << N_IncP_bins_val[i] << "\t" << N_IncP_bins_err[i] << "\t"
                << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\t"
                << fDCoh_bins_val[i] << "\t" << fDCoh_bins_err[i] << "\t"
                << fDInc_bins_val[i] << "\t" << fDInc_bins_err[i] << "\t"
                << fD_bins_val[i] << "\t" << fD_bins_err[i] << "\n";
    }
    /*
    outfile << "Sum over bins:\n";
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss\n";
    outfile << std::fixed << std::setprecision(2);
    Double_t sum_bins[5] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        sum_bins[0] += N_CohJ_bins_val[i];
        sum_bins[1] += N_IncJ_bins_val[i];
        sum_bins[2] += N_CohP_bins_val[i];
        sum_bins[3] += N_IncP_bins_val[i];
        sum_bins[4] += N_Diss_bins_val[i];
    }
    for(Int_t i = 0; i < 4; i++){
        outfile << sum_bins[i] << "\t";
    }
    outfile << sum_bins[4] << "\n***\n";
    */
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + ".txt").Data());

    // Print the TeX table for fC
    outfile.open((*str + "_fC_TeX.txt").Data());
    outfile << std::fixed << std::setprecision(1);
    outfile << "$(0.2,1.0)$\t& $"
            << N_CohJ_to1_val << R"( \pm )" << N_CohJ_to1_err << "$\t& $"
            << N_IncJ_to1_val << R"( \pm )" << N_IncJ_to1_err << "$\t& $"
            << N_Diss_to1_val << R"( \pm )" << N_Diss_to1_err << "$\t& -" << R"( \\)" 
            //<< N_CohP_to1_val << R"( \pm )" << N_CohP_to1_err << "$ \t& $"
            //<< N_IncP_to1_val << R"( \pm )" << N_IncP_to1_err << "$ \t& $"
            << "\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(3)
                << "$(" << ptBoundaries[i] << "," << ptBoundaries[i+1] << ")$\t& $"
                << std::fixed << std::setprecision(1)
                << N_CohJ_bins_val[i] << R"( \pm )" << N_CohJ_bins_err[i] << "$\t& $"
                << N_IncJ_bins_val[i] << R"( \pm )" << N_IncJ_bins_err[i] << "$\t& $"
                << N_Diss_bins_val[i] << R"( \pm )" << N_Diss_bins_err[i] << "$\t& $"
                //<< N_CohP_bins_val[i] << R"( \pm )" << N_CohP_bins_err[i] << "$ \t& $"
                //<< N_IncP_bins_val[i] << R"( \pm )" << N_IncP_bins_err[i] << "$ \t& $"
                << std::fixed << std::setprecision(3)
                << fC_bins_val[i] << R"( \pm )" << fC_bins_err[i] << R"($ \\)"
                << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + "_fC_TeX.txt").Data());

    // Print the TeX table for fD
    outfile.open((*str + "_fD_TeX.txt").Data());
    outfile << std::fixed << std::setprecision(1);
    outfile << "$(0.2,1.0)$\t& $"
            //<< N_CohJ_to1_val << R"( \pm )" << N_CohJ_to1_err << "$ \t& $"
            << N_IncJ_to1_val << R"( \pm )" << N_IncJ_to1_err << "$\t& $"
            << N_Diss_to1_val << R"( \pm )" << N_Diss_to1_err << "$\t& $"
            << N_CohP_to1_val << R"( \pm )" << N_CohP_to1_err << "$\t& $" 
            << N_IncP_to1_val << R"( \pm )" << N_IncP_to1_err << "$\t& -" << R"( \\)" 
            << "\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(3)
                << "$(" << ptBoundaries[i] << "," << ptBoundaries[i+1] << ")$ & $"
                << std::fixed << std::setprecision(1)
                //<< N_CohJ_bins_val[i] << R"( \pm )" << N_CohJ_bins_err[i] << "$ \t& $"
                << N_IncJ_bins_val[i] << R"( \pm )" << N_IncJ_bins_err[i] << "$\t& $"
                << N_Diss_bins_val[i] << R"( \pm )" << N_Diss_bins_err[i] << "$\t& $"
                << N_CohP_bins_val[i] << R"( \pm )" << N_CohP_bins_err[i] << "$\t& $"
                << N_IncP_bins_val[i] << R"( \pm )" << N_IncP_bins_err[i] << "$\t& $"
                << std::fixed << std::setprecision(3)
                << fDCoh_bins_val[i] << R"( \pm )" << fDCoh_bins_err[i] << "$\t& $"
                << fDInc_bins_val[i] << R"( \pm )" << fDInc_bins_err[i] << "$\t& $"
                << fD_bins_val[i] << R"( \pm )" << fD_bins_err[i] << R"($ \\)"
                << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + "_fD_TeX.txt").Data());

    // Print to another file from which the values of fC for CalculateCrossSection.c will be loaded
    str = new TString(Form("%sfC_%ibins_Binn%i_CohSh%i", OutputPtFitWithoutBkg.Data(), nPtBins, BinningOpt, iCohJShape));
    ofstream outfile_fC((*str + ".txt").Data());
    outfile_fC << std::fixed << std::setprecision(3);
    outfile_fC << Form("Bin \tfC [%%]\terr \n");
    for(Int_t i = 0; i < nPtBins; i++){
        outfile_fC << i+1 << "\t" << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\n";
    }
    outfile_fC.close();
    Printf("*** Results printed to %s. ***", (*str + ".txt").Data());

    // Print to another file from which the values of fD for CalculateCrossSection.c will be loaded
    str = new TString(Form("%sfD_%ibins_Binn%i_CohSh%i", OutputPtFitWithoutBkg.Data(), nPtBins, BinningOpt, iCohJShape));
    ofstream outfile_fD((*str + ".txt").Data());
    outfile_fD << std::fixed << std::setprecision(1);
    outfile_fD << Form("Bin \tfD [%%]\terr \n");
    for(Int_t i = 0; i < nPtBins; i++){
        outfile_fD << i+1 << "\t" << fD_bins_val[i] << "\t" << fD_bins_err[i] << "\n";
    }
    outfile_fD.close();
    Printf("*** Results printed to %s. ***", (*str + ".txt").Data());

    // ###############################################################################################################
    // ###############################################################################################################

    // 11) Plot the results
    SetStyle();

    // 11.1) Draw the Correlation Matrix
    TCanvas *cCM = new TCanvas("cCM","cCM",600,500);
    DrawCorrelationMatrixModified(cCM,ResFit,iCohJShape);

    // 11.2) Draw the pt fit
    // Without log scale
    RooPlot* PtFrame = fPt.frame(Title("Pt fit"));
    DHisData.plotOn(PtFrame,Name("DSetData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins));
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000){
        Mod->plotOn(PtFrame,Name("hPDFCohJ"),Components(hPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    } else if(iCohJShape == 2){
        Mod->plotOn(PtFrame,Name("gPDFCohJ"),Components(*gPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    } else if(iCohJShape == 3){
        Mod->plotOn(PtFrame,Name("sumPdfCohJ"),Components(*sumPdfCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    }
    Mod->plotOn(PtFrame,Name("hPDFIncJ"),Components(hPDFIncJ), LineColor(kRed),  LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("hPDFCohP"),Components(hPDFCohP), LineColor(222),   LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("hPDFIncP"),Components(hPDFIncP), LineColor(kRed),  LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("hPDFDiss"),Components(hPDFDiss), LineColor(15),    LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("Mod"),                           LineColor(215),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));

    PtFrame->SetAxisRange(0,1,"X");
    // Set X axis
    PtFrame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    PtFrame->GetXaxis()->SetTitleSize(0.05);
    PtFrame->GetXaxis()->SetLabelSize(0.05);
    // Set Y axis
    //PtFrame->GetYaxis()->SetTitle(Form("Counts per bin"));
    PtFrame->GetYaxis()->SetTitle(Form("d#it{N}/d#it{p}_{T}"));
    PtFrame->GetYaxis()->SetTitleSize(0.05);
    PtFrame->GetYaxis()->SetTitleOffset(0.95);
    PtFrame->GetYaxis()->SetLabelSize(0.05);
    PtFrame->GetYaxis()->SetLabelOffset(0.01);
    PtFrame->GetYaxis()->SetMaxDigits(2);
    PtFrame->GetYaxis()->SetDecimals(1);
    // Draw
    TCanvas *cPt = new TCanvas("cPt","cPt",900,600);
    SetCanvas(cPt, kFALSE); 
    cPt->SetTopMargin(0.06);
    cPt->SetLeftMargin(0.10);
    PtFrame->Draw("]["); 

    // With log scale
    RooPlot* PtFrameLog = fPt.frame(Title("Pt fit log"));
    DHisData.plotOn(PtFrameLog,Name("DSetData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins));
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000){
        Mod->plotOn(PtFrameLog,Name("hPDFCohJ"),Components(hPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    } else if(iCohJShape == 2){
        Mod->plotOn(PtFrameLog,Name("gPDFCohJ"),Components(*gPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    } else if(iCohJShape == 3){
        Mod->plotOn(PtFrameLog,Name("sumPdfCohJ"),Components(*sumPdfCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    }
    Mod->plotOn(PtFrameLog,Name("hPDFIncJ"),Components(hPDFIncJ), LineColor(kRed),  LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("hPDFCohP"),Components(hPDFCohP), LineColor(222),   LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("hPDFIncP"),Components(hPDFIncP), LineColor(kRed),  LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("hPDFDiss"),Components(hPDFDiss), LineColor(15),    LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("Mod"),                           LineColor(215),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));

    PtFrameLog->SetAxisRange(0,2,"X");
    // Set X axis
    PtFrameLog->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    PtFrameLog->GetXaxis()->SetTitleSize(0.05);
    PtFrameLog->GetXaxis()->SetLabelSize(0.05);
    // Set Y axis
    //PtFrame->GetYaxis()->SetTitle(Form("Counts per bin"));
    PtFrameLog->GetYaxis()->SetTitle(Form("d#it{N}/d#it{p}_{T}"));
    PtFrameLog->GetYaxis()->SetTitleSize(0.05);
    PtFrameLog->GetYaxis()->SetTitleOffset(0.95);
    PtFrameLog->GetYaxis()->SetLabelSize(0.05);
    PtFrameLog->GetYaxis()->SetLabelOffset(0.01);
    PtFrameLog->GetYaxis()->SetMaxDigits(2);
    PtFrameLog->GetYaxis()->SetDecimals(1);
    // Draw
    TCanvas *cPtLog = new TCanvas("cPtLog","cPtLog",900,600);
    SetCanvas(cPtLog, kTRUE);
    PtFrameLog->Draw("][");

    // 11.3) Get chi2 
    Double_t chi2 = PtFrameLog->chiSquare("Mod","DSetData",ResFit->floatParsFinal().getSize()); // last argument = number of parameters
    Printf("********************");
    Printf("chi2/ndof = %.3f", chi2);
    Printf("ndof = %i", ResFit->floatParsFinal().getSize());
    Printf("chi2/ndof = %.3f/%i", chi2*ResFit->floatParsFinal().getSize(), ResFit->floatParsFinal().getSize());
    Printf("********************");    

    // 12) Draw the legends
    // Legend 1
    TLegend *l1 = new TLegend(0.16,0.705,0.52,0.93);
    l1->SetHeader("ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),"");
    l1->AddEntry((TObject*)0,Form("#it{m}_{#mu#mu} #in (3.0,3.2) GeV/#it{c}^{2}"),"");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent

    // Legend 2
    TLegend *l2 = new TLegend(0.59,0.555,0.9,0.93);
    //leg2->SetTextSize(0.027);
    l2->AddEntry("DSetData","Data", "EP");
    l2->AddEntry("Mod","sum","L");
    l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f = %.3f/%i",chi2,chi2*ResFit->floatParsFinal().getSize(), ResFit->floatParsFinal().getSize()),"");
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000) l2->AddEntry("hPDFCohJ","coherent J/#psi", "L");
    else if(iCohJShape == 2) l2->AddEntry("gPDFCohJ","coherent J/#psi", "L");
    else if(iCohJShape == 3) l2->AddEntry("sumPdfCohJ","coherent J/#psi", "L");
    l2->AddEntry("hPDFIncJ","incoherent J/#psi", "L");
    l2->AddEntry("hPDFDiss","inc. J/#psi with nucl. diss.", "L");
    l2->AddEntry("hPDFCohP","J/#psi from coh. #psi(2#it{S}) decay", "L");
    l2->AddEntry("hPDFIncP","J/#psi from inc. #psi(2S) decay", "L");
    //l2->AddEntry((TObject*)0,Form("f_{C} = %.4f #pm %.4f", f_C, f_C_err),""); // (#it{p}_T #in (0.2,2.0) GeV/#it{c})
    l2->SetTextSize(0.043);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);

    cPt->cd();
    l1->Draw();
    l2->Draw();

    cPtLog->cd();
    l1->Draw();
    l2->Draw();

    // 13) Print the results to pdf and png
    str = new TString(Form("%sBinn%i_CohSh%i", OutputPtFitWithoutBkg.Data(), BinningOpt, iCohJShape));
    if(iCohJShape > 1000 && !bStopWeight){
        Double_t R_A = 6.6 + (Double_t)(iCohJShape-1001) * 0.1;        
        str = new TString(Form("%sOptimalRA/WeightOverAll/coh_modRA_%.2f", OutputPtFitWithoutBkg.Data(), R_A));
    } 
    if(iCohJShape > 1000 && bStopWeight){
        Double_t R_A = 6.6 + (Double_t)(iCohJShape-1001) * 0.1;        
        str = new TString(Form("%sOptimalRA/StopWeight/coh_modRA_%.2f", OutputPtFitWithoutBkg.Data(), R_A));
    } 
    cCM->Print((*str + "_CM.pdf").Data());
    cCM->Print((*str + "_CM.png").Data());
    cPt->Print((*str + ".pdf").Data());
    cPt->Print((*str + ".png").Data());
    cPtLog->Print((*str + "_log.pdf").Data());
    cPtLog->Print((*str + "_log.png").Data());
    // Print chi2 vs. R_A to text file
    if(iCohJShape > 1000){
        outfile.open((*str + "_chi2.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        Double_t R_A = 6.6 + (Double_t)(iCohJShape-1001) * 0.1;
        outfile << R_A << "\t" << chi2;
        outfile.close();
    }

    return;
}

void PrepareDataTree()
{
    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    // Create new data tree with applied cuts
    TFile fFileOut("Trees/PtFit/PtFitWithoutBkgTree.root","RECREATE");

    TTree *Tree = new TTree("Tree", "Tree");
    Tree->Branch("fPt", &fPt, "fPt/D");
    Tree->Branch("fM", &fM, "fM/D");
    Tree->Branch("fY", &fY, "fY/D");

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(0,2)) Tree->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}

void DrawCorrelationMatrixModified(TCanvas *cCM, RooFitResult* ResFit, Int_t iCohJShape)
{
    // Set margins
    cCM->SetTopMargin(0.03);
    cCM->SetBottomMargin(0.11);
    cCM->SetRightMargin(0.17);
    cCM->SetLeftMargin(0.15);
    // Get 2D corr hist
    TH2* hCorr = ResFit->correlationHist();
    // Set X and Y axis
    if(iCohJShape == 0 || iCohJShape == 1 || iCohJShape > 1000){
        hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
        hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
        //
        hCorr->GetYaxis()->SetBinLabel(1,"#it{N}_{inc}");
        hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{coh}");
    } else if(iCohJShape == 2){
        hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
        hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
        hCorr->GetXaxis()->SetBinLabel(4,"#it{b}");
        //
        hCorr->GetYaxis()->SetBinLabel(1,"#it{b}");
        hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{inc}");
        hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{diss}");
        hCorr->GetYaxis()->SetBinLabel(4,"#it{N}_{coh}");        
    } else if(iCohJShape == 3){
        hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
        hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
        hCorr->GetXaxis()->SetBinLabel(4,"#it{R}_{A}");
        //
        hCorr->GetYaxis()->SetBinLabel(1,"#it{R}_{A}");
        hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{inc}");
        hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{diss}");
        hCorr->GetYaxis()->SetBinLabel(4,"#it{N}_{coh}");            
    }

    // Set corr hist and draw it
    hCorr->SetMarkerSize(3.6);
    hCorr->GetXaxis()->SetLabelSize(0.13);
    hCorr->GetYaxis()->SetLabelSize(0.13);
    hCorr->GetZaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}