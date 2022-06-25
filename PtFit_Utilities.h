// PtFit_Utilities.h
// David Grund, June 21, 2022

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// roofit headers
#include "RooTFnBinding.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRealSumPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning_PtFit.h"
#include "SetPtBinning.h"

using namespace RooFit;

TString NamesPDFs[6] = {"CohJ","IncJ","CohP","IncP","Bkgr","Diss"};

// #############################################################################################
// Options of pT fit without background:
// iRecShape:
// different shapes of CohJ:
//  = 0 => classic histogram from STARlight (R_A = 6.624 fm)
//  = 1 => histogram from STARlight data generated with R_A = 7.53 fm (~ Roman)
//  = 2 => fit using a "Gaussian shape" pT * exp(-b * pT^2)
//  = 3 => fit using a pure STARlight formfactor, R_A left free
// to study which value of R_A is optimal (finding a minimum chi2):
//  = 1001 => R_A = 6.60 fm
//  = 1002 => R_A = 6.70 fm
//  = 1003 => R_A = 6.80 fm
//  = 1004 => R_A = 6.90 fm
//  = 1005 => R_A = 7.00 fm
//  = 1006 => R_A = 7.10 fm
//  = 1007 => R_A = 7.20 fm
//  = 1008 => R_A = 7.30 fm
//  = 1009 => R_A = 7.40 fm
//  = 1010 => R_A = 7.50 fm
//  = 1011 => R_A = 7.60 fm
//  = 1012 => R_A = 7.70 fm
//  = 1013 => R_A = 7.80 fm
// from the above study of various CohJ shapes, it was found that the value of R_A = 7.330 fm is optimal
// all shapes modified accordingly: CohJ, IncJ, CohP, IncP (~ no effect at the last three)
//  = 4 => R_A = 7.330 fm for all shapes
// #############################################################################################

// Functions and variables needed to define STARlight formfactor
const Double_t hbarC = 0.197; // GeV*fm
const Double_t PbAtomicNumber = 208.; // Pb
Double_t VMD_model(Double_t *pT, Double_t *par)
{
    // input and parameters
    Double_t q = pT[0];     // q = sqrt of transferred momentum |t|, unit is GeV 
    if (q <= 0.) return 0.; // protection against floating point exception
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

    return q * FF * FF; // q (= pT) is from the Jacobian! (see the photo from Dec 02, 2021)
}
Double_t Zero_func(Double_t pT){ return 0.; }

// #############################################################################################

void PtFit_NoBkg_DrawCorrMatrix(TCanvas *cCM, RooFitResult* ResFit, Int_t iRecShape)
{
    // Set margins
    cCM->SetTopMargin(0.03);
    cCM->SetBottomMargin(0.11);
    cCM->SetRightMargin(0.17);
    cCM->SetLeftMargin(0.15);
    // Get 2D corr hist
    TH2* hCorr = ResFit->correlationHist();

    // Set X and Y axes
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000){
        // CohJ from STARlight
        // x axis
        hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
        hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
        // y axis
        hCorr->GetYaxis()->SetBinLabel(1,"#it{N}_{inc}");
        hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{coh}");
    } else if(iRecShape == 2){
        // Fit with a Gaussian
        // x axis
        hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
        hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
        hCorr->GetXaxis()->SetBinLabel(4,"#it{b}");
        // y axis
        hCorr->GetYaxis()->SetBinLabel(1,"#it{b}");
        hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{inc}");
        hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{diss}");
        hCorr->GetYaxis()->SetBinLabel(4,"#it{N}_{coh}");        
    } else if(iRecShape == 3){
        // Fit with the STARlight formfactor
        // x axis
        hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
        hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
        hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
        hCorr->GetXaxis()->SetBinLabel(4,"#it{R}_{A}");
        // y axis
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

// #############################################################################################

void PtFit_SetCanvas(TCanvas *c, Bool_t isLogScale)
{
    if(isLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.02);

    return;
}

// #############################################################################################

void PtFit_NoBkg_DoFit(Int_t iRecShape, Int_t ifD = 0)
// ifD = 0 => R_coh = R_inc = R = 0.18 (Michal's measured value)
// systematic uncertainties:
//     = -2 => R = 0.16
//     = -1 => R = 0.17
//     = 1 => R = 0.19
//     = 2 => R = 0.20
{
    Printf("###########################################");
    // ratio of the coherent psi(2S) and J/psi cross sections
    // needed to fix the normalizations of feed-down curves
    // (the ratio of incoherent cross sections is fixed to the same value)
    Double_t R = 0.18 + ifD * 0.01;
    Printf("Ratio of the cross sections: R = %.2f", R);

    // Load the values of fD coefficients
    Double_t fDCoh, fDInc;
    ifstream ifs;
    ifs.open(Form("Results/%sPtFit_FeedDownNormalization/fD_only_R_coh%.3f_R_inc%.3f.txt", str_subfolder.Data(), R, R));
    if(!ifs.fail())
    {
        // Read data from the file
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> fDCoh >> fDInc;
            i++;   
        } 
        ifs.close();
    } else {
        Printf("fD coefficients missing. Terminating...");
        return;
    }
    // from percent to decimal number
    fDCoh = fDCoh / 100.;
    fDInc = fDInc / 100.;
    // print loaded values
    Printf("fD_coh = %.4f", fDCoh);
    Printf("fD_inc = %.4f", fDInc);

    // nuclear radius (Pb)
    Double_t fR_A = 0;
    if(iRecShape == 0) fR_A = 6.624;
    if(iRecShape == 1) fR_A = 7.53;
    if(iRecShape == 4) fR_A = 7.33;
    if(iRecShape > 1000) fR_A = 6.6 + (Double_t)(iRecShape-1001) * 0.1;
    Printf("Pb radius used for CohJ: %.3f", fR_A);
    Printf("###########################################");

    // Load the file with PDFs
    TFile *file = TFile::Open(("Trees/" + str_subfolder + "PtFit/MCTemplates.root").Data(),"read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    TList *list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    //###################################################################################
    // Load histograms
    TH1D *hCohJ = NULL;
    TH1D *hIncJ = NULL;
    TH1D *hCohP = NULL;
    TH1D *hIncP = NULL;

    // if hIncJ, hCohP and hIncP taken from SL with R_A = 6.624 fm
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 2 || iRecShape == 3 || iRecShape > 1000)
    {
        // if CohJ from SL with regular R_A = 6.624 fm
        if(iRecShape == 0 || iRecShape == 2 || iRecShape == 3)
        {
            // 1) kCohJpsiToMu
            hCohJ = (TH1D*)list->FindObject(("h" + NamesPDFs[0]).Data());
            if(hCohJ) Printf("Histogram %s loaded.", hCohJ->GetName());
        }
        // if CohJ from SL with modified RA
        else if(iRecShape == 1 || iRecShape > 1000)
        {
            TString name_file = "Trees/" + str_subfolder + "PtFit/MCTemplates_modRA_CohJ.root";
            TString name_hist = Form("hCohJ_modRA_%.2f", fR_A);     

            TFile *f_modRA = TFile::Open(name_file.Data(),"read");
            if(f_modRA) Printf("Input file %s loaded.", f_modRA->GetName()); 

            TList *l_modRA = (TList*) f_modRA->Get("HistList");
            if(l_modRA) Printf("List %s loaded.", l_modRA->GetName()); 
            // 1) kCohJpsiToMu
            hCohJ = (TH1D*)l_modRA->FindObject(name_hist.Data());
            if(hCohJ) Printf("Histogram %s loaded.", hCohJ->GetName());

            f_modRA->Close();
        }
        // 2) kIncohJpsiToMu
        hIncJ = (TH1D*)list->FindObject(("h" + NamesPDFs[1]).Data());
        if(hIncJ) Printf("Histogram %s loaded.", hIncJ->GetName());
        // 3) kCohPsi2sToMuPi
        hCohP = (TH1D*)list->FindObject(("h" + NamesPDFs[2]).Data());
        if(hCohP) Printf("Histogram %s loaded.", hCohP->GetName());
        // 4) kincohPsi2sToMuPi
        hIncP = (TH1D*)list->FindObject(("h" + NamesPDFs[3]).Data());
        if(hIncP) Printf("Histogram %s loaded.", hIncP->GetName());
    } 
    // if all CohJ hIncJ, hCohP and hIncP taken with R_A = 7.330 fm
    else if(iRecShape == 4)
    {
        if(!isPass3)
        {
            Printf("This option is not supported. Terminating..."); 
            return;
        }

        TString name_file = "Trees/" + str_subfolder + "PtFit/MCTemplates_modRA_all.root";

        TFile *f_modRA = TFile::Open(name_file.Data(),"read");
        if(f_modRA) Printf("Input file %s loaded.", f_modRA->GetName()); 

        TList *l_modRA = (TList*) f_modRA->Get("HistList");
        if(l_modRA) Printf("List %s loaded.", l_modRA->GetName()); 

        l_modRA->ls();

        TString hNames_modRA[4] = {"hCohJ_modRA_7.330",
                                   "hIncJ_modRA_7.330",
                                   "hCohP_modRA_7.330",
                                   "hIncP_modRA_7.330"};

        // 1) kCohJpsiToMu
        hCohJ = (TH1D*)l_modRA->FindObject(hNames_modRA[0].Data());
        if(hCohJ) Printf("Histogram %s loaded.", hCohJ->GetName());
        // 2) kIncohJpsiToMu
        hIncJ = (TH1D*)l_modRA->FindObject(hNames_modRA[1].Data());
        if(hIncJ) Printf("Histogram %s loaded.", hIncJ->GetName());
        // 3) kCohPsi2sToMuPi
        hCohP = (TH1D*)l_modRA->FindObject(hNames_modRA[2].Data());
        if(hCohP) Printf("Histogram %s loaded.", hCohP->GetName());
        // 4) kincohPsi2sToMuPi
        hIncP = (TH1D*)l_modRA->FindObject(hNames_modRA[3].Data());
        if(hIncP) Printf("Histogram %s loaded.", hIncP->GetName());

        f_modRA->Close();
    } 
    // 5) Dissociative
    TH1D *hDiss = (TH1D*)list->FindObject(("h" + NamesPDFs[5]).Data());
    if(hDiss) Printf("Histogram %s loaded.", hDiss->GetName());

    //###################################################################################
    // Create PDFs
    // Definition of roofit variables
    RooRealVar fPt("fPt", "fPt", fPtCutLow_PtFit, fPtCutUpp_PtFit);
    RooArgSet fSetOfVariables(fPt);

    // 1) kCohJpsiToMu
    // iRecShape == 0, 1, 4, >1000
    RooDataHist DHisCohJ("DHisCohJ","DHisCohJ",fPt,hCohJ);
    RooHistPdf  hPDFCohJ("hPDFCohJ","hPDFCohJ",fSetOfVariables,DHisCohJ,0);

    // iRecShape == 2 => fit using the formula pT * exp(-b * pT^2)
    RooRealVar par_b("par_b","b",100.,1.,500.);
    RooGenericPdf *gPDFCohJ = new RooGenericPdf("gPDFCohJ","gPDFCohJ","abs(fPt)*exp(-par_b*pow(fPt,2))",RooArgSet(fPt,par_b));
    // use of abs(.) => so that the PDF is never negative (without abs, we get a WARNING)
    // other option would be to define it through a C++ function Gaussian(), inside of which we put
    // if (pT <= 0.) return 0.;
    
    // iRecShape == 3 => fit using the STARlight form factor and R_A left free
    // https://root.cern.ch/download/doc/RooFit_Users_Manual_2.91-33.pdf
    // https://root-forum.cern.ch/t/bind-tf1-into-roofit-pdf-using-bindpdf/26623
    // https://root-forum.cern.ch/t/defining-a-roogenericpdf-from-a-rooaddition-or-any-rooabsreal/20483 
    // https://root-forum.cern.ch/t/roofit-fitting-a-tf1-binded-function/7656

    // Create RooAbsReal from a TF1 function defined using C++ function VMD_model()
    TF1 *fFormFactorSL = new TF1("fFormFactorSL",VMD_model,fPtCutLow_PtFit,fPtCutUpp_PtFit,2); // 2 = number of parameters
    RooRealVar R_A("R_A","R_A",6.62,1.,12.); // 6.624 fm = the original SL value, will be used as a starting value
    //R_A.setConstant(kTRUE);
    RooRealVar a("a","a", 0.7, 0.7, 0.7); // 0.7 fm = SL value
    a.setConstant(kTRUE);
    RooAbsReal *absRealCohJ = bindFunction(fFormFactorSL, fPt, RooArgList(R_A, a));
    // Create zero function
    TF1 *fZero = new TF1("fZero","Zero_func(x)",fPtCutLow_PtFit,fPtCutUpp_PtFit);
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

    // 5) Dissociative
    RooDataHist DHisDiss("DHisDiss","DHisDiss",fPt,hDiss);
    RooHistPdf  hPDFDiss("hPDFDiss","hPDFDiss",fSetOfVariables,DHisDiss,0);

    // Close the file with MC templates (PDFs)
    file->Close();
    //###################################################################################

    // Get the binned dataset
    file = TFile::Open(("Trees/" + str_subfolder + "PtFit/SignalWithBkgSubtracted.root").Data(), "read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    TH1D *hData = (TH1D*)list->FindObject("hNSigPerBin");
    if(hData) Printf("Histogram %s loaded.", hData->GetName());

    RooDataHist DHisData("DHisData","DHisData",fPt,hData);
    Printf("Binned data with subtracted background loaded.");

    // Calculate the number of entries
    Double_t N_all = 0;
    for(Int_t i = 1; i <= hData->GetNbinsX(); i++){
        N_all += hData->GetBinContent(i);
    }
    Printf("Data contain %.0f entries in %i bins.", N_all, DHisData.numEntries());
    file->Close();

    // Create the model for fitting
    // Normalizations:
    RooRealVar NCohJ("NCohJ","Number of coh J/psi events", 0.90*N_all,0.50*N_all,1.0*N_all);
    RooRealVar NIncJ("NIncJ","Number of inc J/psi events", 0.05*N_all,0.01*N_all,0.3*N_all);
    RooRealVar NDiss("NDiss","Number of dis J/psi events", 0.05*N_all,0.01*N_all,0.3*N_all);
    // Feed-down datasets: fix to primary normalizations through fDs
    RooGenericPdf NCohP("NCohP","Number of coh FD events",Form("NCohJ*%.4f", fDCoh),RooArgSet(NCohJ));
    RooGenericPdf NIncP("NIncP","Number of inc FD events",Form("NIncJ*%.4f", fDInc),RooArgSet(NIncJ));

    // The model:
    RooAddPdf *Mod = NULL;
    // if CohJ from SL => use hPDFCohJ
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000){
        Mod = new RooAddPdf("Mod","Sum of all PDFs",
            RooArgList(hPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
            RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
        );
    // if CohJ as a gaussian => use gPDFCohJ
    } else if(iRecShape == 2){
        Mod = new RooAddPdf("Mod","Sum of all PDFs",
            RooArgList(*gPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
            RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
        );    
    // if CohJ as the SL formfactor => use sumPdfCohJ
    } else if(iRecShape == 3){
        Mod = new RooAddPdf("Mod","Sum of all PDFs",
            RooArgList(*sumPdfCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
            RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
        );        
    }

    // Perform fitting
    RooFitResult *ResFit = Mod->fitTo(DHisData,Extended(kTRUE),Range(""),SumW2Error(kTRUE),Save()); // Range("") => full range

    // ###############################################################################################################
    // Output to text file
    // different types of ifD are only distinguished in iRecShape == 4 (in other cases, we use ifD = 0 by default)
    TString name = "Results/" + str_subfolder;
    // if CohJ
    if(iRecShape == 0 || 
       iRecShape == 1 || 
       iRecShape == 2 || 
       iRecShape == 3) name += Form("PtFit_NoBkg/CohJ%i", iRecShape);
    // if RecShape == 4
    // if ifD == 0 => optimal pT fit
    else if(iRecShape == 4 && ifD == 0){
        name += Form("PtFit_NoBkg/RecSh%i_fD%i", iRecShape, ifD);
    } 
    // if ifD != 0 => systematic uncertainties
    else if(iRecShape == 4 && ifD != 0){
        gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_SystUncertainties/");
        name += Form("PtFit_SystUncertainties/RecSh%i_fD%i", iRecShape, ifD);
    }
    // if study of the optimal value of R_A
    else {
        gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_NoBkg/OptimalRA/Fits/");
        name += Form("PtFit_NoBkg/OptimalRA/Fits/modRA_%.2f", fR_A);
    }
    
    // Integrals of the PDFs in the whole pT range 
    fPt.setRange("fPt_all",fPtCutLow_PtFit,fPtCutUpp_PtFit);
    RooAbsReal *fN_CohJ_all = NULL;
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000) fN_CohJ_all = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPt_all"));
    else if(iRecShape == 2) fN_CohJ_all = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range("fPt_all"));
    else if(iRecShape == 3) fN_CohJ_all = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range("fPt_all"));
    RooAbsReal *fN_IncJ_all = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPt_all"));
    RooAbsReal *fN_Diss_all = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPt_all"));  
    RooAbsReal *fN_CohP_all = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPt_all")); 
    RooAbsReal *fN_IncP_all = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPt_all"));
    // Integrals of the PDFs in the incoherent-enriched sample (IES, pT > 0.2 GeV/c)
    fPt.setRange("fPt_inc",0.2,2.0);
    RooAbsReal *fN_CohJ_inc = NULL;
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000) fN_CohJ_inc = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPt_inc"));
    else if(iRecShape == 2) fN_CohJ_inc = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range("fPt_inc"));
    else if(iRecShape == 3) fN_CohJ_inc = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range("fPt_inc"));
    RooAbsReal *fN_IncJ_inc = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPt_inc"));
    RooAbsReal *fN_Diss_inc = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPt_inc"));  
    RooAbsReal *fN_CohP_inc = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPt_inc")); 
    RooAbsReal *fN_IncP_inc = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPt_inc"));  
    // Integrals of the PDFs with 0.2 < pT < 1.0 GeV/c (allbins)
    fPt.setRange("fPt_to1",0.2,1.0);
    RooAbsReal *fN_CohJ_to1 = NULL;
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000) fN_CohJ_to1 = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPt_to1"));
    else if(iRecShape == 2) fN_CohJ_to1 = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range("fPt_to1"));
    else if(iRecShape == 3) fN_CohJ_to1 = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range("fPt_to1"));
    RooAbsReal *fN_IncJ_to1 = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPt_to1"));
    RooAbsReal *fN_Diss_to1 = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPt_to1"));  
    RooAbsReal *fN_CohP_to1 = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPt_to1")); 
    RooAbsReal *fN_IncP_to1 = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPt_to1"));  
    // Number of events in the whole pT range
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
    // Number of events with 0.2 < pT < 2.0 GeV/c
        // values
        Double_t N_CohJ_inc_val = fN_CohJ_inc->getVal()*NCohJ.getVal();
        Double_t N_IncJ_inc_val = fN_IncJ_inc->getVal()*NIncJ.getVal();
        Double_t N_Diss_inc_val = fN_Diss_inc->getVal()*NDiss.getVal();
        Double_t N_CohP_inc_val = fN_CohP_inc->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_inc_val = fN_IncP_inc->getVal()*fDInc*NIncJ.getVal();
        // errors
        Double_t N_CohJ_inc_err = fN_CohJ_inc->getVal()*NCohJ.getError();
        Double_t N_IncJ_inc_err = fN_IncJ_inc->getVal()*NIncJ.getError();
        Double_t N_Diss_inc_err = fN_Diss_inc->getVal()*NDiss.getError();
        Double_t N_CohP_inc_err = fN_CohP_inc->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_inc_err = fN_IncP_inc->getVal()*fDInc*NIncJ.getError();  
    // Number of events with 0.2 < pT < 1.0 GeV/c
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
    // Total fC and fD corrections
    // a) for the whole IES (0.2 < pT < 2.0 GeV/c)
    // fC
    Double_t fC_inc_val = N_CohJ_inc_val / (N_IncJ_inc_val + N_Diss_inc_val) * 100;
    Double_t denominator_inc_err = TMath::Sqrt(TMath::Power(N_IncJ_inc_err,2) + TMath::Power(N_Diss_inc_err,2));
    Double_t fC_inc_err = fC_inc_val * TMath::Sqrt(TMath::Power((N_CohJ_inc_err/N_CohJ_inc_val),2) + TMath::Power((denominator_inc_err/(N_IncJ_inc_val + N_Diss_inc_val)),2));
    // fD
    Double_t fDCoh_inc_val = N_CohP_inc_val / (N_IncJ_inc_val + N_Diss_inc_val) * 100;
    Double_t fDInc_inc_val = N_IncP_inc_val / (N_IncJ_inc_val + N_Diss_inc_val) * 100;
    Double_t fDCoh_inc_err = fDCoh_inc_val * TMath::Sqrt(TMath::Power((N_CohP_inc_err/N_CohP_inc_val),2) + TMath::Power((denominator_inc_err/(N_IncJ_inc_val + N_Diss_inc_val)),2));
    Double_t fDInc_inc_err = fDInc_inc_val * TMath::Sqrt(TMath::Power((N_IncP_inc_err/N_IncP_inc_val),2) + TMath::Power((denominator_inc_err/(N_IncJ_inc_val + N_Diss_inc_val)),2));
    Double_t fD_inc_val = fDCoh_inc_val + fDInc_inc_val;
    Double_t fD_inc_err = TMath::Sqrt(TMath::Power(fDCoh_inc_err, 2) + TMath::Power(fDInc_inc_err, 2));
    // b) for the allbins range (0.2 < pT < 1.0 GeV/c)
    // fC
    Double_t fC_to1_val = N_CohJ_to1_val / (N_IncJ_to1_val + N_Diss_to1_val) * 100;
    Double_t denominator_to1_err = TMath::Sqrt(TMath::Power(N_IncJ_to1_err,2) + TMath::Power(N_Diss_to1_err,2));
    Double_t fC_to1_err = fC_to1_val * TMath::Sqrt(TMath::Power((N_CohJ_to1_err/N_CohJ_to1_val),2) + TMath::Power((denominator_to1_err/(N_IncJ_to1_val + N_Diss_to1_val)),2));
    // fD
    Double_t fDCoh_to1_val = N_CohP_to1_val / (N_IncJ_to1_val + N_Diss_to1_val) * 100;
    Double_t fDInc_to1_val = N_IncP_to1_val / (N_IncJ_to1_val + N_Diss_to1_val) * 100;
    Double_t fDCoh_to1_err = fDCoh_to1_val * TMath::Sqrt(TMath::Power((N_CohP_to1_err/N_CohP_to1_val),2) + TMath::Power((denominator_to1_err/(N_IncJ_to1_val + N_Diss_to1_val)),2));
    Double_t fDInc_to1_err = fDInc_to1_val * TMath::Sqrt(TMath::Power((N_IncP_to1_err/N_IncP_to1_val),2) + TMath::Power((denominator_to1_err/(N_IncJ_to1_val + N_Diss_to1_val)),2));
    Double_t fD_to1_val = fDCoh_to1_val + fDInc_to1_val;
    Double_t fD_to1_err = TMath::Sqrt(TMath::Power(fDCoh_to1_err, 2) + TMath::Power(fDInc_to1_err, 2));
    // fC and fD corrections in the pT bins
    // Integrals of the PDFs in the pT bins:
    vector<RooAbsReal*> fN_CohJ_bins(nPtBins);
    vector<RooAbsReal*> fN_IncJ_bins(nPtBins);
    vector<RooAbsReal*> fN_Diss_bins(nPtBins);
    vector<RooAbsReal*> fN_CohP_bins(nPtBins);
    vector<RooAbsReal*> fN_IncP_bins(nPtBins);
    // values of N
    vector<Double_t> N_CohJ_bins_val(nPtBins);
    vector<Double_t> N_IncJ_bins_val(nPtBins);
    vector<Double_t> N_Diss_bins_val(nPtBins);
    vector<Double_t> N_CohP_bins_val(nPtBins);
    vector<Double_t> N_IncP_bins_val(nPtBins);
    // errors of N
    vector<Double_t> N_CohJ_bins_err(nPtBins);
    vector<Double_t> N_IncJ_bins_err(nPtBins);
    vector<Double_t> N_Diss_bins_err(nPtBins);
    vector<Double_t> N_CohP_bins_err(nPtBins);
    vector<Double_t> N_IncP_bins_err(nPtBins);
    // coefficients
    vector<Double_t> fC_bins_val(nPtBins);
    vector<Double_t> fC_bins_err(nPtBins);
    vector<Double_t> fDCoh_bins_val(nPtBins);
    vector<Double_t> fDCoh_bins_err(nPtBins);
    vector<Double_t> fDInc_bins_val(nPtBins);
    vector<Double_t> fDInc_bins_err(nPtBins);
    vector<Double_t> fD_bins_val(nPtBins);
    vector<Double_t> fD_bins_err(nPtBins);
    // calculate the values of everything
    for(Int_t i = 0; i < nPtBins; i++)
    {
        fPt.setRange(Form("fPt_bin%i",i+1), ptBoundaries[i], ptBoundaries[i+1]);
        Printf("Now calculating for bin %i, (%.3f, %.3f) GeV", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000) fN_CohJ_bins[i] = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));
        else if(iRecShape == 2) fN_CohJ_bins[i] = gPDFCohJ->createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));
        else if(iRecShape == 3) fN_CohJ_bins[i] = sumPdfCohJ->createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));
        fN_IncJ_bins[i] = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));
        fN_Diss_bins[i] = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));
        fN_CohP_bins[i] = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));
        fN_IncP_bins[i] = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range(Form("fPt_bin%i",i+1)));  
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
        // fD correction
        fDCoh_bins_val[i] = N_CohP_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]) * 100;
        fDInc_bins_val[i] = N_IncP_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]) * 100;
        if(N_CohP_bins_val[i] != 0) fDCoh_bins_err[i] = fDCoh_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohP_bins_err[i]/N_CohP_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));
        else fDCoh_bins_err[i] = 0.;
        if(N_IncP_bins_val[i] != 0) fDInc_bins_err[i] = fDInc_bins_val[i] * TMath::Sqrt(TMath::Power((N_IncP_bins_err[i]/N_IncP_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));
        else fDInc_bins_val[i] = 0.;
        fD_bins_val[i] = fDCoh_bins_val[i] + fDInc_bins_val[i];
        fD_bins_err[i] = TMath::Sqrt(TMath::Power(fDCoh_bins_err[i], 2) + TMath::Power(fDInc_bins_err[i], 2));
    }
    // total sums of events
    Double_t sum_all_val = N_CohJ_all_val + N_IncJ_all_val + N_CohP_all_val + N_IncP_all_val + N_Diss_all_val;
    Double_t sum_inc_val = N_CohJ_inc_val + N_IncJ_inc_val + N_CohP_inc_val + N_IncP_inc_val + N_Diss_inc_val;
    Double_t sum_to1_val = N_CohJ_to1_val + N_IncJ_to1_val + N_CohP_to1_val + N_IncP_to1_val + N_Diss_to1_val;
    Double_t sum_all_err = TMath::Sqrt(TMath::Power(N_CohJ_all_err,2) + TMath::Power(N_IncJ_all_err,2) + 
                                       TMath::Power(N_CohP_all_err,2) + TMath::Power(N_IncP_all_err,2) +
                                       TMath::Power(N_Diss_all_err,2));
    Double_t sum_inc_err = TMath::Sqrt(TMath::Power(N_CohJ_inc_err,2) + TMath::Power(N_IncJ_inc_err,2) + 
                                       TMath::Power(N_CohP_inc_err,2) + TMath::Power(N_IncP_inc_err,2) +
                                       TMath::Power(N_Diss_inc_err,2));
    Double_t sum_to1_err = TMath::Sqrt(TMath::Power(N_CohJ_to1_err,2) + TMath::Power(N_IncJ_to1_err,2) + 
                                       TMath::Power(N_CohP_to1_err,2) + TMath::Power(N_IncP_to1_err,2) +
                                       TMath::Power(N_Diss_to1_err,2));
    // sums over bins
    vector<Double_t> sum_bins_val(nPtBins);
    vector<Double_t> sum_bins_err(nPtBins);
    for(Int_t i = 0; i < nPtBins; i++){
        sum_bins_val[i] = N_CohJ_bins_val[i] + N_IncJ_bins_val[i] + N_CohP_bins_val[i] + N_IncP_bins_val[i] + N_Diss_bins_val[i];
        sum_bins_err[i] = TMath::Sqrt(TMath::Power(N_CohJ_bins_err[i],2) + TMath::Power(N_IncJ_bins_err[i],2) + 
                                      TMath::Power(N_CohP_bins_err[i],2) + TMath::Power(N_IncP_bins_err[i],2) +
                                      TMath::Power(N_Diss_bins_err[i],2));
    }    

    // Print the numbers to the text file
    ofstream outfile((name + ".txt").Data());
    outfile << Form("Dataset contains %.0f events.\n***\n", N_all);
    outfile << "pT range\t\tNCohJ \terr \tNIncJ \terr \tNDiss \terr \tNCohP \terr \tNIncP \terr \tsum \terr \tfC \terr \tfD coh\terr \tfD inc\terr \tfD \terr\n";
    outfile << "(0.000, 2.000) GeV/c\t"
            << std::fixed << std::setprecision(1)
            << N_CohJ_all_val << "\t" << N_CohJ_all_err << "\t" 
            << N_IncJ_all_val << "\t" << N_IncJ_all_err << "\t" 
            << N_Diss_all_val << "\t" << N_Diss_all_err << "\t" 
            << N_CohP_all_val << "\t" << N_CohP_all_err << "\t" 
            << N_IncP_all_val << "\t" << N_IncP_all_err << "\t"
            << sum_all_val << "\t" << sum_all_err << "\n";
    outfile << "(0.200, 2.000) GeV/c\t"
            << std::fixed << std::setprecision(1)
            << N_CohJ_inc_val << "\t" << N_CohJ_inc_err << "\t" 
            << N_IncJ_inc_val << "\t" << N_IncJ_inc_err << "\t" 
            << N_Diss_inc_val << "\t" << N_Diss_inc_err << "\t" 
            << N_CohP_inc_val << "\t" << N_CohP_inc_err << "\t" 
            << N_IncP_inc_val << "\t" << N_IncP_inc_err << "\t"
            << sum_inc_val    << "\t" << sum_inc_err << "\t"
            << std::fixed << std::setprecision(2)
            << fC_inc_val << "\t" << fC_inc_err << "\t" 
            << fDCoh_inc_val << "\t" << fDCoh_inc_err << "\t" 
            << fDInc_inc_val << "\t" << fDInc_inc_err << "\t"
            << fD_inc_val << "\t" << fD_inc_err << "\n";
    outfile << "(0.200, 1.000) GeV/c\t"
            << std::fixed << std::setprecision(1)
            << N_CohJ_to1_val << "\t" << N_CohJ_to1_err << "\t" 
            << N_IncJ_to1_val << "\t" << N_IncJ_to1_err << "\t" 
            << N_Diss_to1_val << "\t" << N_Diss_to1_err << "\t"  
            << N_CohP_to1_val << "\t" << N_CohP_to1_err << "\t" 
            << N_IncP_to1_val << "\t" << N_IncP_to1_err << "\t"
            << sum_to1_val    << "\t" << sum_to1_err << "\t"
            << std::fixed << std::setprecision(2)
            << fC_to1_val << "\t" << fC_to1_err << "\t" 
            << fDCoh_to1_val << "\t" << fDCoh_to1_err << "\t" 
            << fDInc_to1_val << "\t" << fDInc_to1_err << "\t"
            << fD_to1_val << "\t" << fD_to1_err << "\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << Form("(%.3f, %.3f) GeV/c\t", ptBoundaries[i], ptBoundaries[i+1])
                << std::fixed << std::setprecision(1)
                << N_CohJ_bins_val[i] << "\t" << N_CohJ_bins_err[i] << "\t"
                << N_IncJ_bins_val[i] << "\t" << N_IncJ_bins_err[i] << "\t"
                << N_Diss_bins_val[i] << "\t" << N_Diss_bins_err[i] << "\t"
                << N_CohP_bins_val[i] << "\t" << N_CohP_bins_err[i] << "\t" 
                << N_IncP_bins_val[i] << "\t" << N_IncP_bins_err[i] << "\t"
                << sum_bins_val[i]    << "\t" << sum_bins_err[i] << "\t"
                << std::fixed << std::setprecision(2)
                << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\t"
                << fDCoh_bins_val[i] << "\t" << fDCoh_bins_err[i] << "\t"
                << fDInc_bins_val[i] << "\t" << fDInc_bins_err[i] << "\t"
                << fD_bins_val[i] << "\t" << fD_bins_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (name + ".txt").Data());

    // Print the TeX table for fC
    outfile.open((name + "_fC_TeX.txt").Data());
    outfile << std::fixed << std::setprecision(1);
    outfile << "$(0.200,1.000)$\t& $"
            << N_CohJ_to1_val << R"( \pm )" << N_CohJ_to1_err << "$\t& $"
            << N_IncJ_to1_val << R"( \pm )" << N_IncJ_to1_err << "$\t& $"
            << N_Diss_to1_val << R"( \pm )" << N_Diss_to1_err << "$\t& $"  
            //<< N_CohP_to1_val << R"( \pm )" << N_CohP_to1_err << "$ \t& $"
            //<< N_IncP_to1_val << R"( \pm )" << N_IncP_to1_err << "$ \t& $"
            << std::fixed << std::setprecision(2)
            << fC_to1_val << R"( \pm )" << fC_to1_err << R"($ \\)"
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
    Printf("*** Results printed to %s. ***", (name + "_fC_TeX.txt").Data());

    // Print the TeX table for fD
    outfile.open((name + "_fD_TeX.txt").Data());
    outfile << std::fixed << std::setprecision(1);
    outfile << "$(0.200,1.000)$ & $"
            //<< N_CohJ_to1_val << R"( \pm )" << N_CohJ_to1_err << "$ \t& $"
            << N_IncJ_to1_val << R"( \pm )" << N_IncJ_to1_err << "$\t& $"
            << N_Diss_to1_val << R"( \pm )" << N_Diss_to1_err << "$\t& $"
            << N_CohP_to1_val << R"( \pm )" << N_CohP_to1_err << "$\t& $" 
            << N_IncP_to1_val << R"( \pm )" << N_IncP_to1_err << "$\t& $"
            << fDCoh_to1_val  << R"( \pm )" << fDCoh_to1_err  << "$\t& $" 
            << fDInc_to1_val  << R"( \pm )" << fDInc_to1_err  << "$\t& $" 
            << fD_to1_val     << R"( \pm )" << fD_to1_err     << R"($ \\)"
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
    Printf("*** Results printed to %s. ***", (name + "_fD_TeX.txt").Data());

    // Print the values of fC to another file from which they will be loaded in CalculateCrossSection.h
    outfile.open((name + "_fC.txt").Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << Form("Bin \tfC [%%]\terr \n");
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << i+1 << "\t" << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (name + "_fC.txt").Data());

    // Print the values of fD to another file from which they will be loaded in CalculateCrossSection.h
    outfile.open((name + "_fD.txt").Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << Form("Bin \tfD [%%]\terr \n");
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << i+1 << "\t" << fD_bins_val[i] << "\t" << fD_bins_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (name + "_fD.txt").Data());

    // ###############################################################################################################

    // Plot the results
    gStyle->SetOptTitle(0); // suppress title
    gStyle->SetOptStat(0);  // the type of information printed in the histogram statistics box
                            // 0 = no information
    gStyle->SetPalette(1);  // set color map
    gStyle->SetPaintTextFormat("4.2f"); // precision if plotted with "TEXT"

    // Draw the Correlation Matrix
    TCanvas *cCM = new TCanvas("cCM","cCM",600,500);
    PtFit_NoBkg_DrawCorrMatrix(cCM,ResFit,iRecShape);

    // Draw the pT fit
    // Without log scale
    RooPlot* PtFrame = fPt.frame(Title("pT fit"));
    DHisData.plotOn(PtFrame,Name("DHisData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins_PtFit));
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000){
        Mod->plotOn(PtFrame,Name("hPDFCohJ"),Components(hPDFCohJ),     LineColor(222),LineStyle(1),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    } else if(iRecShape == 2){
        Mod->plotOn(PtFrame,Name("gPDFCohJ"),Components(*gPDFCohJ),    LineColor(222),LineStyle(1),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    } else if(iRecShape == 3){
        Mod->plotOn(PtFrame,Name("sumPdfCohJ"),Components(*sumPdfCohJ),LineColor(222),LineStyle(1),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    }
    Mod->plotOn(PtFrame,Name("hPDFIncJ"),Components(hPDFIncJ),LineColor(kRed),LineStyle(1),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("hPDFCohP"),Components(hPDFCohP),LineColor(222), LineStyle(7),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("hPDFIncP"),Components(hPDFIncP),LineColor(kRed),LineStyle(7),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("hPDFDiss"),Components(hPDFDiss),LineColor(15),  LineStyle(1),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrame,Name("Mod"),                          LineColor(215), LineStyle(1),LineWidth(3),Range(""),Normalization(sum_all_val,RooAbsReal::NumEvent));

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
    PtFit_SetCanvas(cPt, kFALSE); 
    cPt->SetTopMargin(0.06);
    cPt->SetLeftMargin(0.10);
    PtFrame->Draw("]["); 

    // With log scale
    RooPlot* PtFrameLog = fPt.frame(Title("pT fit log scale"));
    DHisData.plotOn(PtFrameLog,Name("DHisData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins_PtFit));
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000){
        Mod->plotOn(PtFrameLog,Name("hPDFCohJ"),Components(hPDFCohJ),     LineColor(222),LineStyle(1),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    } else if(iRecShape == 2){
        Mod->plotOn(PtFrameLog,Name("gPDFCohJ"),Components(*gPDFCohJ),    LineColor(222),LineStyle(1),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    } else if(iRecShape == 3){
        Mod->plotOn(PtFrameLog,Name("sumPdfCohJ"),Components(*sumPdfCohJ),LineColor(222),LineStyle(1),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    }
    Mod->plotOn(PtFrameLog,Name("hPDFIncJ"),Components(hPDFIncJ),LineColor(kRed),LineStyle(1),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("hPDFCohP"),Components(hPDFCohP),LineColor(222), LineStyle(7),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("hPDFIncP"),Components(hPDFIncP),LineColor(kRed),LineStyle(7),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("hPDFDiss"),Components(hPDFDiss),LineColor(15),  LineStyle(1),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));
    Mod->plotOn(PtFrameLog,Name("Mod"),                          LineColor(215), LineStyle(1),LineWidth(3),Normalization(sum_all_val,RooAbsReal::NumEvent));

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
    PtFit_SetCanvas(cPtLog, kTRUE);
    PtFrameLog->Draw("][");

    // Get chi2 
    Double_t chi2 = PtFrameLog->chiSquare("Mod","DHisData",ResFit->floatParsFinal().getSize()); // last argument = number of parameters
    Printf("********************");
    Printf("chi2/ndof = %.3f", chi2);
    Printf("ndof = %i", ResFit->floatParsFinal().getSize());
    Printf("chi2/ndof = %.3f/%i", chi2*ResFit->floatParsFinal().getSize(), ResFit->floatParsFinal().getSize());
    Printf("********************");    

    // Draw the legends
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
    l2->AddEntry("DHisData","Data", "EP");
    l2->AddEntry("Mod","sum","L");
    l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f = %.3f/%i",chi2,chi2*ResFit->floatParsFinal().getSize(), ResFit->floatParsFinal().getSize()),"");
    if(iRecShape == 0 || iRecShape == 1 || iRecShape == 4 || iRecShape > 1000) l2->AddEntry("hPDFCohJ","coherent J/#psi", "L");
    else if(iRecShape == 2) l2->AddEntry("gPDFCohJ","coherent J/#psi", "L");
    else if(iRecShape == 3) l2->AddEntry("sumPdfCohJ","coherent J/#psi", "L");
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

    // Print the results to pdf and png
    cCM->Print((name + "_CM.pdf").Data());
    cCM->Print((name + "_CM.png").Data());
    cPt->Print((name + ".pdf").Data());
    cPt->Print((name + ".png").Data());
    cPtLog->Print((name + "_log.pdf").Data());
    cPtLog->Print((name + "_log.png").Data());
    // If we study the optimal value of R_A, print chi2 vs. R_A to a text file
    if(iRecShape > 1000){
        outfile.open((name + "_chi2.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << fR_A << "\t" << chi2;
        outfile.close();
    }

    delete cCM;
    delete cPt;
    delete cPtLog;

    return;
}