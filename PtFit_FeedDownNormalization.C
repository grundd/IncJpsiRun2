// PtFit_FeedDownNormalization.C
// David Grund, May 29, 2022
// only implemented for pass3 input data

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include <string> // getline
// root headers
#include "TFile.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

// Arrays of AxEs
Double_t AxE_CohJ_val = 0;
Double_t AxE_CohJ_err = 0;
Double_t AxE_IncJ_val = 0;
Double_t AxE_IncJ_err = 0;
Double_t AxE_Psi2s_val[4] = { 0 }; // order: coh_ch, coh_ne, inc_ch, inc_ne
Double_t AxE_Psi2s_err[4] = { 0 };
// Results
Double_t fD_val[4] = { 0 };
Double_t fD_err[4] = { 0 };
// STARlight cross sections
Double_t sig_SL_j_inc = 5.247; //mb
Double_t sig_SL_j_coh = 12.504;//mb
Double_t sig_SL_p_inc = 0.92;  //mb
Double_t sig_SL_p_coh = 2.52;  //mb
Double_t ratio_coh = 0.18; // from Michal's paper
// Branching ratios
Double_t BR_ch = 0.3468;
Double_t BR_ch_err = 0.0030;
Double_t BR_ne = 0.2538;
Double_t BR_ne_err = 0.0032;

TString NamesMC[6] = {"kCohJpsiToMu","kIncohJpsiToMu","kCohPsi2sToMuPi_ch","kCohPsi2sToMuPi_ne","kIncohPsi2sToMuPi_ch","kIncohPsi2sToMuPi_ne"};

void Calc_AxE(Int_t iMC);
void Calc_fD(Bool_t ratioFromMeas);
Double_t Calc_ErrBayes(Double_t k, Double_t n);

void PtFit_FeedDownNormalization(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_FeedDownNormalization/");

    for(Int_t iMC = 0; iMC < 6; iMC++)
    {
        //Calc_AxE(iMC);
    }

    Calc_fD();

    return;
}

void Calc_AxE(Int_t iMC)
{
    Bool_t isPsi2sDataset = kFALSE;
    TFile *f = NULL;
    TTree *tRec = NULL;
    TTree *tGen = NULL;

    // kCohJpsiToMu or kIncohJpsiToMu
    if(iMC == 0 || iMC == 1)
    {
        if(iMC == 0) f = new TFile("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kCohJpsiToMu.root");
        if(iMC == 1) f = new TFile("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root");
        if(f) Printf("Input file loaded.");

        tRec = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJpsi"));
        tGen = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJpsiMCGen"));
        
    }
    // if feed-down datasets
    else
    {
        isPsi2sDataset = kTRUE;

        if(iMC == 2 || iMC == 3) f = new TFile("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kCohPsi2sToMuPi_2.root");
        if(iMC == 4 || iMC == 5) f = new TFile("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohPsi2sToMuPi_2.root");

        if(iMC == 2 || iMC == 4)
        {
            TList *l = (TList*) f->Get("AnalysisOutput/fOutputListcharged");
            if(l) Printf("List %s loaded.", l->GetName());
            tRec = (TTree*)l->FindObject("fTreeJpsi");
            tGen = (TTree*)l->FindObject("fTreeJpsiMCGen");
        }
        if(iMC == 3 || iMC == 5)
        {
            TList *l = (TList*) f->Get("AnalysisOutput/fOutputListneutral");
            if(l) Printf("List %s loaded.", l->GetName());
            tRec = (TTree*)l->FindObject("fTreeJpsi");
            tGen = (TTree*)l->FindObject("fTreeJpsiMCGen");
        }
    }
    if(tRec) Printf("Input tRec loaded.");
    if(tGen) Printf("Input tGen loaded.");

    ConnectTreeVariablesMCRec(tRec, isPsi2sDataset);
    ConnectTreeVariablesMCGen(tGen);

    Double_t NRec(0.), NGen(0.);

    // go over reconstructed events
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++)
    {
        tRec->GetEntry(iEntry);
        // m between 3.0 and 3.2 GeV/c^2, pT between 0.0 and 2.0 GeV/c
        if(EventPassedMCRec(1, 2)) NRec++;
    }
    // go over generated events
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++)
    {
        tGen->GetEntry(iEntry);
        // no pT cut
        if(EventPassedMCGen(-1)) NGen++;
    }

    Double_t AxE_val = NRec / NGen;
    Double_t AxE_err = Calc_ErrBayes(NRec, NGen);

    TString str_out = "Results/" + str_subfolder + "PtFit_FeedDownNormalization/AxE_" + NamesMC[iMC];
    ofstream outfile((str_out + ".txt").Data());
    outfile << std::fixed << std::setprecision(5);
    outfile << AxE_val << "\t" << AxE_err;
    outfile.close();
    Printf("*** AxE values saved to %s.txt ***", str_out.Data());

    return;
}

void Calc_fD(Bool_t ratioFromMeas)
{
    TString str_in = "Results/" + str_subfolder + "PtFit_FeedDownNormalization/AxE_";
    ifstream infile;
    // Load the values of AxE for CohJ
    infile.open((str_in + NamesMC[0] + ".txt").Data());
    infile >> AxE_CohJ_val >> AxE_CohJ_err;
    infile.close();
    // Load the values of AxE for IncJ
    infile.open((str_in + NamesMC[1] + ".txt").Data());
    infile >> AxE_IncJ_val >> AxE_IncJ_err;
    infile.close();
    // Load the values of AxE for feed-down processes
    for(Int_t i = 0; i < 4; i++)
    {
        infile.open((str_in + NamesMC[i] + ".txt").Data());
        infile >> AxE_Psi2s_val[i] >> AxE_Psi2s_err[i];
        infile.close();  
    }
    Printf("Input files loaded...");

    // 2) Define output file
    TString str_out = "Results/" + str_subfolder + "PtFit_FeedDownNormalization/fD_PtFit";
    if(Bool_t ratioFromMeas) str_out.Append("_ratioFromMeas.txt");
    str_out.Append(".txt");
    ofstream outfile(str_out.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("fD[%%] \tCohCh \tErr \tIncCh \tErr \tCohNe \tErr \tIncNe \tErr \n");
    /*
    // 3) Calculate fD correction
    for(Int_t i = 0; i < 1; i++){
        // all values of fD in percent already
        if(!bRatioMeasured) fD_coh_ch_val[i] = sig_SL_p_coh / sig_SL_j_coh * AxE_Psi2s_val[0][i] / AxE_CohJ_val * BR_ch * 100;
        else fD_coh_ch_val[i] = ratio_coh * AxE_Psi2s_val[0][i] / AxE_CohJ_val * BR_ch * 100;
        fD_inc_ch_val[i] = sig_SL_p_inc / sig_SL_j_inc * AxE_Psi2s_val[1][i] / AxE_IncJ_val[i] * BR_ch * 100;
        if(!bRatioMeasured) fD_coh_ne_val[i] = sig_SL_p_coh / sig_SL_j_coh * AxE_Psi2s_val[2][i] / AxE_CohJ_val * BR_ne * 100;
        else fD_coh_ne_val[i] = ratio_coh * AxE_Psi2s_val[2][i] / AxE_CohJ_val * BR_ne * 100;
        fD_inc_ne_val[i] = sig_SL_p_inc / sig_SL_j_inc * AxE_Psi2s_val[3][i] / AxE_IncJ_val[i] * BR_ne * 100;
        if(AxE_Psi2s_val[0][i] != 0){
            fD_coh_ch_err[i] = fD_coh_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[0][i]/AxE_Psi2s_val[0][i]),2) + 
                TMath::Power((AxE_CohJ_err/AxE_CohJ_val),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_coh_ch_err[i] = 0;
        if(AxE_Psi2s_val[1][i] != 0){
            fD_inc_ch_err[i] = fD_inc_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[1][i]/AxE_Psi2s_val[1][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_inc_ch_err[i] = 0;
        if(AxE_Psi2s_val[2][i] != 0){
            fD_coh_ne_err[i] = fD_coh_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[2][i]/AxE_Psi2s_val[2][i]),2) + 
                TMath::Power((AxE_CohJ_err/AxE_CohJ_val),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_coh_ne_err[i] = 0;
        if(AxE_Psi2s_val[3][i] != 0){
            fD_inc_ne_err[i] = fD_inc_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[3][i]/AxE_Psi2s_val[3][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_inc_ne_err[i] = 0; 
        // Print the results to text file
        outfile << "Total:";
        outfile << "\t" << fD_coh_ch_val[i] << "\t" << fD_coh_ch_err[i] 
                << "\t" << fD_inc_ch_val[i] << "\t" << fD_inc_ch_err[i] 
                << "\t" << fD_coh_ne_val[i] << "\t" << fD_coh_ne_err[i] 
                << "\t" << fD_inc_ne_val[i] << "\t" << fD_inc_ne_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());
    */
    
    return;
}


Double_t Calc_ErrBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}