// AnalysisConfig.h
// David Grund, Feb 26, 2022
// Configure values of the parameters

// cpp headers
#include <stdio.h> // printf
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
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

using namespace RooFit;

// ******************************
Bool_t START_FROM_CLEAN = kFALSE;
Int_t iWhichAnalysisToRun = 10;
// ******************************

// Options to set:
Int_t cut_fVertexContrib = 2;
Double_t cut_fVertexZ = 15.;
Double_t cut_fY = 0.8;
Double_t cut_fEta = 0.8;
Double_t cut_fZN_neutrons = 10.5;

// Options that will be set by the choice of iAnalysis in InitAnalysis:
Int_t nPtBins = -1;
Bool_t isPass3;
Bool_t isZNcut;
TString str_subfolder = "";
TString str_in_DT_fldr = "";
TString str_in_DT_tree = "";
TString str_in_MC_fldr = "";
TString str_in_MC_tree = "";

void InitAnalysis(Int_t iAnalysis){
    // iAnalysis = 0:  pass1, 4 bins, isZNcut OFF
    // iAnalysis = 1:  pass1, 5 bins, isZNcut OFF
    // iAnalysis = 10: pass3, 4 bins, isZNcut OFF
    // iAnalysis = 11: pass3, 5 bins, isZNcut OFF
    // iAnalysis = 20: pass3, 4 bins, isZNcut ON
    // iAnalysis = 21: pass3, 5 bins, isZNcut ON
    if(iAnalysis == 0){
        str_subfolder = "pass1_4bins/";
        nPtBins = 4;
        isPass3 = kFALSE;
        isZNcut = kFALSE;
    }
    if(iAnalysis == 1){
        str_subfolder = "pass1_5bins/";
        nPtBins = 5;
        isPass3 = kFALSE;
        isZNcut = kFALSE;
    }
    if(iAnalysis == 10){
        str_subfolder = "pass3_4bins/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isZNcut = kFALSE;
    }
    if(iAnalysis == 11){
        str_subfolder = "pass3_5bins/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isZNcut = kFALSE;
    }
    if(iAnalysis == 20){
        str_subfolder = "pass3_4bins_ZNcut/";
        nPtBins = 4;
        isPass3 = kTRUE;
        isZNcut = kTRUE;
    }
    if(iAnalysis == 21){
        str_subfolder = "pass3_5bins_ZNcut/";
        nPtBins = 5;
        isPass3 = kTRUE;
        isZNcut = kTRUE;
    }
    // set the path to the input data
    if(!isPass3){
        str_in_DT_fldr = "Trees/AnalysisData_pass1/";
        str_in_DT_tree = "AnalysisOutput/fTreeJPsi";
        str_in_MC_fldr = "Trees/AnalysisDataMC_pass1/";
        str_in_MC_tree = "AnalysisOutput/fTreeJPsiMCRec";
    } else {
        str_in_DT_fldr = "Trees/AnalysisData_pass3/";
        str_in_DT_tree = "AnalysisOutput/fTreeJpsi";
        str_in_MC_fldr = "Trees/AnalysisDataMC_pass3/";
        str_in_MC_tree = "AnalysisOutput/fTreeJpsi";
    }
    // if we want to start with a clean folder
    if(START_FROM_CLEAN){
        // delete previous subfolder inside the Results/ folder (if exists)
        gSystem->Exec("rm -r Results/" + str_subfolder);
        // create a new one of the same name
        gSystem->Exec("mkdir -p Results/" + str_subfolder);
    }

    Printf("*** Analysis initiated successfully. ***");
    Printf("*** Results will be stored in the folder: Results/%s ***", str_subfolder.Data());

    return;
}