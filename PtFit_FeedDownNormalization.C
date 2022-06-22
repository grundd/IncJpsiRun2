// PtFit_FeedDownNormalization.C
// David Grund, June 20, 2022
// only implemented for pass3 input data

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
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
// Branching ratios
Double_t BR_ch = 0.3468;
Double_t BR_ch_err = 0.0030;
Double_t BR_ne = 0.2538;
Double_t BR_ne_err = 0.0032;

TString NamesMC[6] = {"kCohJpsiToMu",
                      "kIncohJpsiToMu",
                      "kCohPsi2sToMuPi_charged",
                      "kCohPsi2sToMuPi_neutral",
                      "kIncohPsi2sToMuPi_charged",
                      "kIncohPsi2sToMuPi_neutral"};

void Calc_AxE(Int_t iMC);
void Calc_fD(Double_t R_coh, Double_t R_inc); // R = ratios of the corresponding cross sections
Double_t Calc_ErrBayes(Double_t k, Double_t n);

void PtFit_FeedDownNormalization(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_FeedDownNormalization/");

    for(Int_t iMC = 0; iMC < 6; iMC++) Calc_AxE(iMC);

    // to compare with AN_v1
    Calc_fD(0.202,0.176);

    // using the ratio measured by Michal for both values
    // R = 0.1800 pm 0.0185 (stat.) pm 0.0280 (syst.) pm 0.0050 (BR)
    // error when all contributions combined in quadrature: R = 0.18 pm 0.03
    for(Double_t R = 0.16; R < 0.21; R+=0.01) Calc_fD(R,R);

    return;
}

void Calc_AxE(Int_t iMC)
{
    TString str_AxE = "Results/" + str_subfolder + "PtFit_FeedDownNormalization/AxE_" + NamesMC[iMC] + ".txt";
    ifstream ifs;
    ifs.open(str_AxE);
    if(!ifs.fail())
    {
        Printf("This AxE has already been calculated.");
        ifs.close();
        return;
    }
    else
    {
        // kCohJpsiToMu or kIncohJpsiToMu
        Bool_t isPsi2sDataset = kFALSE;
        // if feed-down datasets
        if(iMC > 1) isPsi2sDataset = kTRUE;

        TFile *fRec = new TFile(Form("%sAnalysisResults_MC_%s.root", str_in_MC_fldr_rec.Data(), NamesMC[iMC].Data()));
        if(fRec) Printf("Input fRec loaded.");
        TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJpsi"));
        if(tRec) Printf("Input tRec loaded.");
        TFile *fGen = NULL;
        TList *lGen = NULL;
        TTree *tGen = NULL;
        TString str_f_gen = "";
        // kCohJpsiToMu or kIncohJpsiToMu
        if(iMC <= 1)
        {
            fGen = new TFile(Form("%sAnalysisResults_MC_%s.root", str_in_MC_fldr_gen.Data(), NamesMC[iMC].Data()));
            if(fGen) Printf("Input fGen loaded.");
            tGen = dynamic_cast<TTree*> (fGen->Get("AnalysisOutput/fTreeJpsiMCGen"));
            if(tGen) Printf("Input tGen loaded.");
        }
        // if feed-down datasets
        else
        {
            if(iMC == 2 || iMC == 3) str_f_gen = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kCohPsi2sToMuPi_2.root";
            if(iMC == 4 || iMC == 5) str_f_gen = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohPsi2sToMuPi_2.root";
            fGen = new TFile(Form("%s", str_f_gen.Data()));
            if(fGen) Printf("Input fGen loaded.");
            TString str_l = "";
            // charged feed-down
            if(iMC == 2 || iMC == 4) str_l = "AnalysisOutput/fOutputListcharged";
            // neutral feed-down
            else                     str_l = "AnalysisOutput/fOutputListneutral";
            lGen = (TList*)fGen->Get(str_l.Data());
            if(lGen) Printf("List %s loaded.", lGen->GetName()); 
            tGen = (TTree*)lGen->FindObject("fTreeJpsiMCGen");
            if(tGen) Printf("Tree %s loaded.", tGen->GetName());
        }
        
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

        Double_t AxE_val = NRec / NGen * 100; // in percent
        Double_t AxE_err = Calc_ErrBayes(NRec, NGen) * 100; // in percent
        
        ofstream outfile(str_AxE.Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << AxE_val << "\t" << AxE_err;
        outfile.close();
        Printf("*** AxE values saved to %s. ***", str_AxE.Data());

        return;
    }
}

void Calc_fD(Double_t R_coh, Double_t R_inc)
{
    // Load the values of AxE
    TString str_in = "Results/" + str_subfolder + "PtFit_FeedDownNormalization/AxE_";
    ifstream ifs;
    // CohJ
    ifs.open((str_in + NamesMC[0] + ".txt").Data());
    ifs >> AxE_CohJ_val >> AxE_CohJ_err;
    ifs.close();
    // IncJ
    ifs.open((str_in + NamesMC[1] + ".txt").Data());
    ifs >> AxE_IncJ_val >> AxE_IncJ_err;
    ifs.close();
    // Feed-down datasets
    for(Int_t i = 0; i < 4; i++)
    {
        ifs.open((str_in + NamesMC[i+2] + ".txt").Data());
        ifs >> AxE_Psi2s_val[i] >> AxE_Psi2s_err[i];
        ifs.close();  
    }
    Printf("Input values loaded...");
    // Define the output file
    TString str_out = "Results/" + str_subfolder + Form("PtFit_FeedDownNormalization/fD_R_coh%.3f_R_inc%.3f.txt", R_coh, R_inc);
    ofstream outfile(str_out.Data());
    // Print the values of AxE
    outfile << std::fixed << std::setprecision(3)
            << "\tCohJ \tIncJ \tCohCh \tCohNe \tIncCh \tIncNe \n"
            << Form("AxE [%%]\t")
            << AxE_CohJ_val << "\t" << AxE_IncJ_val << "\t" 
            << AxE_Psi2s_val[0] << "\t" << AxE_Psi2s_val[1] << "\t" 
            << AxE_Psi2s_val[2] << "\t" << AxE_Psi2s_val[3] << "\n"
            << Form("err [%%]\t")
            << AxE_CohJ_err << "\t" << AxE_IncJ_err << "\t" 
            << AxE_Psi2s_err[0] << "\t" << AxE_Psi2s_err[1] << "\t" 
            << AxE_Psi2s_err[2] << "\t" << AxE_Psi2s_err[3] << "\n\n";

    // Calculate the values of fD and then print them values of fD
    Double_t R;
    Double_t AxE_J_val, AxE_J_err;
    Double_t BR_val, BR_err;
    for(Int_t i = 0; i < 4; i++)
    {
        // coherent
        if(i == 0 || i == 1){R = R_coh; AxE_J_val = AxE_CohJ_val; AxE_J_err = AxE_CohJ_err;}
        // incoherent
        else                {R = R_inc; AxE_J_val = AxE_IncJ_val; AxE_J_err = AxE_IncJ_err;}
        // charged feed-down
        if(i == 0 || i == 2){BR_val = BR_ch; BR_err = BR_ch_err;}
        // neutral feed-down
        else                {BR_val = BR_ne; BR_err = BR_ne_err;}
        // calculate fD (in percent already)
        fD_val[i] = R * AxE_Psi2s_val[i] / AxE_J_val * BR_val * 100;
        // calculate error
        fD_err[i] = fD_val[i] * TMath::Sqrt(
            TMath::Power(AxE_Psi2s_err[i] / AxE_Psi2s_val[i], 2) +
            TMath::Power(AxE_CohJ_err / AxE_CohJ_val, 2) +
            TMath::Power(BR_err / BR_val, 2)
        );
    }
    // Print the results to text file
    outfile << std::fixed << std::setprecision(2)
            << "\tCohCh \tCohNe \tIncCh \tIncNe \n"
            << Form("fD  [%%]\t")
            << fD_val[0] << "\t" << fD_val[1] << "\t" << fD_val[2] << "\t" << fD_val[3] << "\n"
            << Form("err [%%]\t")
            << fD_err[0] << "\t" << fD_err[1] << "\t" << fD_err[2] << "\t" << fD_err[3] << "\n\n";
    outfile.close();
    Printf("*** Results printed to %s.***", str_out.Data());
    // Print only the values of fD_coh and fD_inc to be loaded by PtFit_NoBkg
    TString str_out2 = "Results/" + str_subfolder + Form("PtFit_FeedDownNormalization/fD_only_R_coh%.3f_R_inc%.3f.txt", R_coh, R_inc);
    outfile.open(str_out2.Data());
    outfile << std::fixed << std::setprecision(2)
            << "fD_coh \t fD_inc \n" 
            << fD_val[0] + fD_val[1] << "\t" << fD_val[2] + fD_val[3];
    outfile.close();
    Printf("*** Results printed to %s.***", str_out2.Data());
    
    return;
}

Double_t Calc_ErrBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}