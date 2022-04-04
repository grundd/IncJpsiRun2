// PhotoCrossSec_Calculate.c
// David Grund, Sep 21, 2021

// cpp headers
#include <fstream>  // print output to txt file
#include <iomanip>  // std::setprecision()
#include <string>   // getline
// root headers
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

//*************************************************
// Systematic uncertainties (in percent)
Double_t ErrSyst_SigExtr[nPtBins] = { 0 };
Double_t ErrSyst_AxE[nPtBins] = { 0 };
Double_t ErrSyst_fD[nPtBins] = { 0 };
Double_t ErrSyst_fC[nPtBins] = { 0 };
Double_t ErrSyst_lumi = 2.7;
Double_t ErrSyst_veto = 3.0;
Double_t ErrSyst_EMD = 2.0;
Double_t ErrSyst_tracks = 2.8;
Double_t ErrSyst_CCUP31 = 1.3;
Double_t ErrSyst_flux = 2.0;
//*************************************************
Double_t Lumi18q = 90.114;  // 1/(mu barn)
Double_t Lumi18r = 137.812; // 1/(mu barn)
Double_t LumiAll_val = Lumi18q + Lumi18r;
Double_t LumiAll_err = LumiAll_val * ErrSyst_lumi / 100.;
Double_t BR_val = 0.05961;
Double_t BR_err = 0.00033;
Double_t ErrSyst_BR = BR_err / BR_val * 100.;
Double_t RapWidth = 1.6;
Double_t Eff_veto_val = 92.0; // !!!!!!!!!!!!!!!!!!!!!!! 94.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Double_t Eff_veto_err = Eff_veto_val * ErrSyst_veto / 100.;
Double_t Eff_EMD_val = 92.0;
Double_t Eff_EMD_err = Eff_EMD_val * ErrSyst_EMD / 100.;
Double_t PhotonFlux_val = 84.9;
Double_t PhotonFlux_err = PhotonFlux_val * ErrSyst_flux / 100.;
// Old values
Double_t Lumi18q_aod = 91.431;  // 1/(mu barn)
Double_t Lumi18r_aod = 147.292; // 1/(mu barn)
Double_t LumiAll_aod_val = Lumi18q_aod + Lumi18r_aod;
Double_t LumiAll_aod_err = LumiAll_aod_val * ErrSyst_lumi / 100.;
// Cross section in pT^2 bins
Double_t N_yield_val[nPtBins] = { 0 };
Double_t N_yield_err[nPtBins] = { 0 };
Double_t Pt2Widths[nPtBins] = { 0 };
Double_t AxE_val[nPtBins] = { 0 };
Double_t AxE_err[nPtBins] = { 0 };
Double_t CorrFD_val[4][nPtBins] = { 0 };
Double_t CorrFD_err[4][nPtBins] = { 0 };
Double_t CorrFD_sum_val[nPtBins] = { 0 };
Double_t CorrFD_sum_err[nPtBins] = { 0 };
Double_t CorrFC_val[nPtBins] = { 0 };
Double_t CorrFC_err[nPtBins] = { 0 };
Double_t Sigma_UPC_val[nPtBins] = { 0 };
Double_t Sigma_UPC_err_stat[nPtBins] = { 0 };
Double_t Sigma_UPC_err_syst[nPtBins] = { 0 };
Double_t Sigma_photo_val[nPtBins] = { 0 };
Double_t Sigma_photo_err_stat[nPtBins] = { 0 };
Double_t Sigma_photo_err_syst[nPtBins] = { 0 };
Double_t t_avg_val[nPtBins] = { 0 };
// Total cross section for pt > 0.2 GeV/c
Double_t CorrFC_tot = 0;
Double_t CorrFC_tot_err = 0;
Double_t Sigma_tot = 0;
Double_t Sigma_tot_err = 0;

Int_t i_bin;

void CalculateCrossSectionTotal(Int_t iFeedDown, Bool_t bAOD = kFALSE);
void CalculateCrossSectionBins(Int_t iFeedDown);
void PrintErr(TString str);

void CalculateCrossSection(){

    // ESDs
    //CalculateCrossSectionTotal(kFALSE);

    // AODs
    //CalculateCrossSectionTotal(kTRUE);

    CalculateCrossSectionBins(0);
    CalculateCrossSectionBins(1);

    return;
}

void CalculateCrossSectionTotal(Int_t iFeedDown, Bool_t bAOD){
    // for pt > 0.2 GeV/c

    ifstream ifs;

    // 1) Load N_yield
    if(!bAOD){
        TString str_yield = "Results/InvMassFit/inc/inc_signal.txt";
        ifs.open(str_yield.Data());
        // Read data from the file
        if(!ifs.fail()){
            ifs >> N_yield_val[0] >> N_yield_err[0];
        } else {
            PrintErr(str_yield);
            return;
        }
        ifs.close(); 
    } else {
        N_yield_val[0] = 643;
        N_yield_err[0] = 31;         
    }
    Printf("1) N_yield loaded.");

    // 2) Load AxE
    TString str_AxE;
    if(!bAOD) str_AxE = "Results/AccAndEffMC/AxE_JInc_MassCut0_PtCut0.txt";
    else str_AxE = "Results/AccAndEffMC/AOD/AxE_AOD_JInc_MassCut1_PtCut0.txt";
    ifs.open(str_AxE.Data());
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> AxE_val[0] >> AxE_err[0];
            i++;   
        }
    } else {
        PrintErr(str_AxE);
        return;            
    }
    ifs.close();
    Printf("2) AxE for kIncohJpsiToMu loaded.");

    // 3) Load all FD corrections
    TString str_FD; 
    if(     !bAOD && iFeedDown == 0) return;
    else if( bAOD && iFeedDown == 0) return;
    else if(!bAOD && iFeedDown == 1) str_FD = "Results/FeedDown/FeedDown_Total.txt";
    else if( bAOD && iFeedDown == 1) str_FD = "Results/FeedDown/FeedDown_Total_AOD.txt";
    ifs.open(str_FD.Data());
    // Read data from the file
    if(!ifs.fail()){
        if(iFeedDown == 0){
            // 0 = feed-down from PtFitWithoutBkg.c
            Int_t i = 0;
            std::string str;
            while(std::getline(ifs,str)){
                istringstream in_stream(str);
                // skip first line
                if(i > 0) in_stream >> i_bin >> CorrFD_val[0][i-1] >> CorrFD_err[0][i-1]
                                             >> CorrFD_val[1][i-1] >> CorrFD_err[1][i-1];
                // Cross-check:
                // Printf("%.4f pm %.4f", CorrFD_val[0][i-1], CorrFD_err[0][i-1]);
                i++;   
            }
        } else {
            // 1 = feed-down from FeedDown.c (FeedDown_debug.c)
            Int_t i = 0;
            char ch[8];
            std::string str;
            while(std::getline(ifs,str)){
                istringstream in_stream(str);
                // skip first line
                if(i > 0) in_stream >> ch >> CorrFD_val[0][0] >> CorrFD_err[0][0] 
                                          >> CorrFD_val[1][0] >> CorrFD_err[1][0]
                                          >> CorrFD_val[2][0] >> CorrFD_err[2][0]
                                          >> CorrFD_val[3][0] >> CorrFD_err[3][0];
                // Cross-check:
                // Printf("%.4f pm %.4f", CorrFD_val[0][i-1], CorrFD_err[0][i-1]);
                i++;   
            }
        }
        Printf("3) FD corrections loaded.");
        ifs.close();
    } else {
        PrintErr(str_FD);
        return;
    }
    // Calculate total FD per bin
    Double_t CorrFD_sum_val = 0;
    Double_t CorrFD_sum_err = 0;
    Double_t SumOfSquares = 0;
    for(Int_t iFD = 0; iFD < 4; iFD++){
        CorrFD_sum_val += CorrFD_val[iFD][0];
        SumOfSquares += TMath::Power(CorrFD_err[iFD][0],2);
    }
    CorrFD_sum_err = TMath::Sqrt(SumOfSquares);

    // 4) Load FC corr
    // (...)

    /*
    // Calculate the total cross section
    Double_t Factors;
    Double_t Factors_err;
    Factors = 1.0 + CorrFD_sum / 100 + CorrFC_tot / 100;
    Factors_err = TMath::Sqrt(TMath::Power(CorrFD_sum_err / 100,2) + TMath::Power(CorrFC_tot_err / 100,2));
    Sigma_tot = N_yield_tot_val / (Factors * EffVetoes * AxE_tot/100 * BR * RapWidth * LumiAll);
    Sigma_tot_err = Sigma_tot * TMath::Sqrt(
        TMath::Power(N_yield_tot_val_err / N_yield_tot_val, 2) +
        TMath::Power(AxE_tot_err / AxE_tot, 2) +
        TMath::Power(Factors_err / Factors, 2) +
        TMath::Power(BR_err / BR, 2));

    // Define output text file to print results
    TString FilePath;
    if(!bAOD) FilePath = "Results/CrossSection/total_ESDs_output.txt";
    else FilePath = "Results/CrossSection/total_AODs_output.txt";
    ofstream fout_sigmaUPC(FilePath.Data());
    fout_sigmaUPC << std::fixed << std::setprecision(3);
    // Print results to the text file
    fout_sigmaUPC << Form("Lumi = %.3f 1/(mili barn)\n", LumiAll);
    fout_sigmaUPC << Form("BR(J/psi -> mu mu) = (%.3f pm %.3f)%%\n", BR*100, BR_err*100);
    fout_sigmaUPC << Form("Delta y = %.1f\n", RapWidth);
    fout_sigmaUPC << Form("EffVetoes = %.1f%%\n", EffVetoes*100);
    fout_sigmaUPC << Form("N \tN_er \tAxE \tAxE_er\tFD [%%]\tFD_err \tFC [%%]\tFC_er \tf [%%]\tf_er \tsig \tsig_er \n");
    fout_sigmaUPC << std::fixed << std::setprecision(2);
    fout_sigmaUPC << N_yield_tot_val << "\t";
    fout_sigmaUPC << N_yield_tot_val_err << "\t";
    fout_sigmaUPC << std::fixed << std::setprecision(3);
    fout_sigmaUPC << AxE_tot << "\t";
    fout_sigmaUPC << AxE_tot_err << "\t";
    fout_sigmaUPC << CorrFD_sum << "\t";
    fout_sigmaUPC << CorrFD_sum_err << "\t";
    fout_sigmaUPC << CorrFC_tot << "\t";
    fout_sigmaUPC << CorrFC_tot_err << "\t";
    fout_sigmaUPC << Factors << "\t";
    fout_sigmaUPC << Factors_err << "\t";
    fout_sigmaUPC << Sigma_tot << "\t";
    fout_sigmaUPC << Sigma_tot_err << "\n";
    fout_sigmaUPC.close();
    Printf("Results printed to %s.", FilePath.Data()); 
    */

    return;
}

void CalculateCrossSectionBins(Int_t iFeedDown){
    // 0 = feed-down from PtFitWithoutBkg.c
    // 1 = feed-down from FeedDown.c (FeedDown_debug.c)

    SetPtBinning();

    ifstream ifs;

    Printf("Calculating cross section in %i bins.", nPtBins);

    // 1) Load N_yield per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        TString str_yield = Form("Results/InvMassFit/%ibins/bin%i_signal.txt", nPtBins, iBin+1);
        ifs.open(str_yield.Data());
        // Read data from the file
        if(!ifs.fail()){
            ifs >> N_yield_val[iBin] >> N_yield_err[iBin];
        } else {
            PrintErr(str_yield);
            return;
        }
        ifs.close(); 
    } 
    Printf("1) N_yield for %ibins loaded.", nPtBins);

    // 2) Load AxE per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        TString str_AxE = Form("Results/AccAndEffMC/AxE_%ibins/AxE_bin%i_JInc.txt", nPtBins, iBin+1);
        ifs.open(str_AxE.Data());
        // Read data from the file
        if(!ifs.fail()){
            Int_t i = 0;
            std::string str;
            while(std::getline(ifs,str)){
                istringstream in_stream(str);
                // skip first line
                if(i == 1) in_stream >> AxE_val[iBin] >> AxE_err[iBin];
                i++;   
            }
        } else {
            PrintErr(str_AxE);
            return;            
        }
        ifs.close();
    }
    Printf("2) AxE from kIncohJpsiToMu for %ibins loaded.", nPtBins);

    // 3) Load all FD corrections per bin
    TString str_FD; 
    if(iFeedDown == 0) str_FD = Form("Results/PtFitWithoutBkg/fD_%ibins_Binn2_CohSh0.txt", nPtBins);
    else if(iFeedDown == 1) str_FD = Form("Results/FeedDown/FeedDown_%ibins.txt", nPtBins);
    ifs.open(str_FD.Data());
    // Read data from the file
    if(!ifs.fail()){
        if(iFeedDown == 0){
            // 0 = feed-down from PtFitWithoutBkg.c
            Int_t i = 0;
            std::string str;
            while(std::getline(ifs,str)){
                istringstream in_stream(str);
                // skip first line
                if(i > 0) in_stream >> i_bin >> CorrFD_val[0][i-1] >> CorrFD_err[0][i-1]
                                             >> CorrFD_val[1][i-1] >> CorrFD_err[1][i-1];
                // Cross-check:
                // Printf("%.4f pm %.4f", CorrFD_val[0][i-1], CorrFD_err[0][i-1]);
                i++;   
            }
        } else {
            // 1 = feed-down from FeedDown.c (FeedDown_debug.c)
            Int_t i = 0;
            std::string str;
            while(std::getline(ifs,str)){
                istringstream in_stream(str);
                // skip first line
                if(i > 0) in_stream >> i_bin >> CorrFD_val[0][i-1] >> CorrFD_err[0][i-1] 
                                    >> CorrFD_val[1][i-1] >> CorrFD_err[1][i-1]
                                    >> CorrFD_val[2][i-1] >> CorrFD_err[2][i-1]
                                    >> CorrFD_val[3][i-1] >> CorrFD_err[3][i-1];
                // Cross-check:
                // Printf("%.4f pm %.4f", CorrFD_val[0][i-1], CorrFD_err[0][i-1]);
                i++;   
            }
        }
        Printf("3) FD corrections for %ibins loaded.", nPtBins);
        ifs.close();
    } else {
        PrintErr(str_FD);
        return;
    }
    // Calculate total FD per bin
    Double_t CorrFD_sum_val[nPtBins] = { 0 };
    Double_t CorrFD_sum_err[nPtBins] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Double_t SumOfSquares = 0;
        for(Int_t iFD = 0; iFD < 4; iFD++){
            CorrFD_sum_val[iBin] += CorrFD_val[iFD][iBin];
            SumOfSquares += TMath::Power(CorrFD_err[iFD][iBin],2);
        }
        CorrFD_sum_err[iBin] = TMath::Sqrt(SumOfSquares);
    }

    // 4) Load FC corr per bin
    TString str_FC = Form("Results/PtFitWithoutBkg/fC_%ibins_Binn2_CohSh0.txt", nPtBins);
    ifs.open(str_FC.Data());
    // Read data from the file
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFC_val[i-1] >> CorrFC_err[i-1];
            // Cross-check:
            // Printf("%.4f pm %.4f", CorrFC_val[i-1], CorrFC_err[i-1]);
            i++;   
        }
        Printf("4) FC corrections for %ibins loaded.", nPtBins);
        ifs.close();
    } else {
        PrintErr(str_FC);
        return;
    }

    // 5) Widths of pt intervals [in GeV^2]
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Pt2Widths[iBin] = TMath::Power(ptBoundaries[iBin+1], 2) - TMath::Power(ptBoundaries[iBin], 2);
    }
    Printf("5) pt^2 widths calculated.");

    // Calculate the UPC cross section per bin
    Double_t Factors_val[nPtBins];
    Double_t Factors_err[nPtBins];
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Sigma_UPC_val[iBin] = N_yield_val[iBin] / (
            (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_val[iBin] / 100.) * 
            (AxE_val[iBin] / 100.) * 
            (Eff_veto_val / 100.) * 
            (Eff_EMD_val / 100.) * 
            (LumiAll_val * 1000) *
            BR_val * 
            RapWidth * Pt2Widths[iBin]);
        Sigma_UPC_err_stat[iBin] = Sigma_UPC_val[iBin] * N_yield_err[iBin] / N_yield_val[iBin];
    }

    // Systematic uncertainties
    // Signal extraction
    TString str_ErrSystSigExtr = Form("Results/InvMassFit_SystUncertainties/%ibins/ErrSystSignalExtraction_%ibins.txt", nPtBins, nPtBins);
    ifs.open(str_ErrSystSigExtr.Data());
    if(!ifs.fail()){
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            ifs >> i_bin >> ErrSyst_SigExtr[iBin];
        }
        ifs.close();
    } else {
        PrintErr(str_ErrSystSigExtr);
        return;        
    }
    // AxE MC
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        ErrSyst_AxE[iBin] = AxE_err[iBin] / AxE_val[iBin] * 100.;
    }
    // fC, fD and total systematic uncertainty of sigma UPC
    Double_t CorrFD_upp[nPtBins] = { 0 };
    Double_t CorrFD_low[nPtBins] = { 0 };
    Double_t CorrFC_upp[nPtBins] = { 0 };
    Double_t CorrFC_low[nPtBins] = { 0 };
    Double_t Sigma_FD_upp[nPtBins] = { 0 };
    Double_t Sigma_FD_low[nPtBins] = { 0 };
    Double_t Sigma_FC_upp[nPtBins] = { 0 };
    Double_t Sigma_FC_low[nPtBins] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        CorrFD_upp[iBin] = CorrFD_sum_val[iBin] + CorrFD_sum_err[iBin];
        CorrFD_low[iBin] = CorrFD_sum_val[iBin] - CorrFD_sum_err[iBin];
        CorrFC_upp[iBin] = CorrFC_val[iBin] + CorrFC_err[iBin];
        CorrFC_low[iBin] = CorrFC_val[iBin] - CorrFC_err[iBin];
        Sigma_FD_upp[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_upp[iBin] / 100. + CorrFC_val[iBin] / 100.);
        Sigma_FD_low[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_low[iBin] / 100. + CorrFC_val[iBin] / 100.);
        Sigma_FC_upp[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_upp[iBin] / 100.);
        Sigma_FC_low[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_sum_val[iBin] / 100. + CorrFC_low[iBin] / 100.);
        Double_t Sigma_fD_upp_diff, Sigma_fD_low_diff, Sigma_fC_upp_diff, Sigma_fC_low_diff;
        Sigma_fD_upp_diff = TMath::Abs(Sigma_FD_upp[iBin] - Sigma_UPC_val[iBin]);
        Sigma_fD_low_diff = TMath::Abs(Sigma_FD_low[iBin] - Sigma_UPC_val[iBin]);
        Sigma_fC_upp_diff = TMath::Abs(Sigma_FC_upp[iBin] - Sigma_UPC_val[iBin]);
        Sigma_fC_low_diff = TMath::Abs(Sigma_FC_low[iBin] - Sigma_UPC_val[iBin]);
        ErrSyst_fD[iBin] = TMath::Max(Sigma_fD_upp_diff / Sigma_UPC_val[iBin], Sigma_fD_low_diff / Sigma_UPC_val[iBin]) * 100.;
        ErrSyst_fC[iBin] = TMath::Max(Sigma_fC_upp_diff / Sigma_UPC_val[iBin], Sigma_fC_low_diff / Sigma_UPC_val[iBin]) * 100.;

        Sigma_UPC_err_syst[iBin] = Sigma_UPC_val[iBin] * TMath::Sqrt(
            TMath::Power(ErrSyst_SigExtr[iBin] / 100., 2) +
            TMath::Power(ErrSyst_AxE[iBin] / 100., 2) +
            TMath::Power(ErrSyst_fD[iBin] / 100., 2) +
            TMath::Power(ErrSyst_fC[iBin] / 100., 2) +
            TMath::Power(ErrSyst_lumi / 100., 2) + 
            TMath::Power(ErrSyst_veto / 100., 2) + 
            TMath::Power(ErrSyst_EMD / 100., 2) + 
            TMath::Power(ErrSyst_tracks / 100., 2) + 
            TMath::Power(ErrSyst_CCUP31 / 100., 2) + 
            TMath::Power(BR_err / BR_val, 2)
        );
    }

    // Calculate the photonuclear cross section per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Sigma_photo_val[iBin] = Sigma_UPC_val[iBin] / 2. / PhotonFlux_val * 1000;
        Sigma_photo_err_stat[iBin] = Sigma_UPC_err_stat[iBin] / 2. / PhotonFlux_val * 1000;
        Sigma_photo_err_syst[iBin] = Sigma_photo_val[iBin] * TMath::Sqrt(
            TMath::Power(Sigma_UPC_err_syst[iBin] / Sigma_UPC_val[iBin], 2) + 
            TMath::Power(ErrSyst_flux / 100., 2)
        );
    }

    // Load avg values of |t| per bin
    TString str_avgT = Form("DependenceOnT/output_%ibins.txt", nPtBins);
    ifs.open(str_avgT.Data()); 
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            istr >> i_bin >> t_avg_val[i];
            i++;   
        }
        ifs.close();
    } else {
        PrintErr(str_avgT);
        return;         
    }
    Printf("6) Values of avg |t| values per bin loaded.");

    // Print results to the text file
    // 1a) Print the UPC cross section 
    TString str_out_1a = Form("Results/CrossSection/%ibins_FeedDown%i.txt", nPtBins, iFeedDown);
    ofstream fout_sigmaUPC(str_out_1a.Data());
    fout_sigmaUPC << Form("Lumi\terr\tRapW\tBR\terr\te_veto\terr\te_EMD\terr\tflux\terr\n")
                  << Form("%.1f \t%.1f \t%.1f \t%.3f \t%.3f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n\n",
                            LumiAll_val, LumiAll_err, 
                            RapWidth, 
                            BR_val*100., BR_err*100., 
                            Eff_veto_val, Eff_veto_err, 
                            Eff_EMD_val, Eff_EMD_err,
                            PhotonFlux_val, PhotonFlux_err);
    fout_sigmaUPC << Form("Bin\tPt2Low\tPt2Upp\tPt2_W\tN\terr\tAxE\terr\tFD [%%]\terr\tFC [%%]\terr\tsig\terr_sta\terr_sys\n");
    for(Int_t i = 0; i < nPtBins; i++){
        fout_sigmaUPC << std::fixed << std::setprecision(3)
                << i+1 << "\t"
                << ptBoundaries[i] * ptBoundaries[i] << "\t"
                << ptBoundaries[i+1] * ptBoundaries[i+1] << "\t"
                << std::fixed << std::setprecision(4)
                << Pt2Widths[i] << "\t"
                << std::fixed << std::setprecision(2)
                << N_yield_val[i] << "\t" << N_yield_err[i] << "\t"
                << AxE_val[i] << "\t" << AxE_err[i] << "\t"
                << std::fixed << std::setprecision(3)
                << CorrFD_sum_val[i] << "\t" << CorrFD_sum_err[i] << "\t"
                << CorrFC_val[i] << "\t" << CorrFC_err[i] << "\t"
                << std::fixed << std::setprecision(2)
                << Sigma_UPC_val[i] << "\t" << Sigma_UPC_err_stat[i] << "\t" << Sigma_UPC_err_syst[i] << "\n";
    }
    fout_sigmaUPC.close();
    Printf("Results printed to %s.", str_out_1a.Data()); 
    // 1b) Print the UPC cross section: TeX file
    TString str_out_1b = Form("Results/CrossSection/%ibins_FeedDown%i_TeX.txt", nPtBins, iFeedDown);
    ofstream fout_sigmaUPC_TeX(str_out_1b.Data());
    fout_sigmaUPC_TeX << Form("$%.0f", LumiAll_val) << R"( \pm )" << Form("%.0f$", LumiAll_err) << " &\t"
                      << Form("%.1f", RapWidth) << "\t"
                      << Form("$%.3f", BR_val*100.) << R"( \pm )" << Form("%.3f$", BR_err*100.) << " &\t"
                      << Form("$%.1f",Eff_veto_val) << R"( \pm )" << Form("%.1f$",Eff_veto_err) << " &\t"
                      << Form("$%.1f", Eff_EMD_val) << R"( \pm )" << Form("%.1f$", Eff_EMD_err) << " &\t"
                      << Form("$%.1f", PhotonFlux_val) << R"( \pm )" << Form("%.1f$", PhotonFlux_err) << R"( \\)" << "\n\n";
    for(Int_t i = 0; i < nPtBins; i++){
        fout_sigmaUPC_TeX << std::fixed << std::setprecision(3)
                          << "$(" << ptBoundaries[i] << "," << ptBoundaries[i+1] << ")$ & "
                          << std::fixed << std::setprecision(4)
                          << Pt2Widths[i] << " &\t$"
                          << std::fixed << std::setprecision(0)
                          << N_yield_val[i] << R"( \pm )" << N_yield_err[i] << "$ &\t$"
                          << std::fixed << std::setprecision(2)
                          << AxE_val[i] << R"( \pm )" << AxE_err[i] << "$ &\t$"
                          << std::fixed << std::setprecision(2)
                          << CorrFD_sum_val[i] << R"( \pm )" << CorrFD_sum_err[i] << "$ &\t$"
                          << std::fixed << std::setprecision(3)
                          << CorrFC_val[i] << R"( \pm )" << CorrFC_err[i] << "$ &\t$"
                          << std::fixed << std::setprecision(2)
                          << Sigma_UPC_val[i] << R"( \pm )" << Sigma_UPC_err_stat[i] << R"((stat.) \pm )" << Sigma_UPC_err_syst[i] << R"((syst.)$ \\)" << "\n";
    }                  
    fout_sigmaUPC_TeX.close();
    Printf("Results printed to %s.", str_out_1b.Data());

    // 2a) Print the systematic uncertainties
    TString str_out2 = Form("Results/CrossSection/%ibins_FeedDown%i_systematics.txt", nPtBins, iFeedDown);
    ofstream fout_systErr(str_out2.Data());
    fout_systErr << "All in percent\n";
    fout_systErr << "lumi \tveto \tEMD \ttracks \tCCUP31 \tBR \n"
                 << Form("%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n\n",
                    ErrSyst_lumi, ErrSyst_veto, ErrSyst_EMD, ErrSyst_tracks, ErrSyst_CCUP31, ErrSyst_BR);
    fout_systErr << "Bin \tSigExt \tAxE MC \tfD \tfC \n";
    for(Int_t i = 0; i < nPtBins; i++){
        fout_systErr << i+1 << std::fixed << std::setprecision(1) << "\t"
                            << ErrSyst_SigExtr[i] << "\t"
                            << ErrSyst_AxE[i] << "\t"
                            << ErrSyst_fD[i] << "\t"
                            << ErrSyst_fC[i] << "\n";
    }    
    fout_systErr.close();
    Printf("Results printed to %s.", str_out2.Data());
    // 2b) Print the systematic uncertainties: TeX file
    TString str_out2b = Form("Results/CrossSection/%ibins_FeedDown%i_systematics_TeX.txt", nPtBins, iFeedDown);
    ofstream fout_systErr_TeX(str_out2b.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        fout_systErr_TeX << std::fixed << std::setprecision(3)
                        << "$(" << ptBoundaries[i] << "," << ptBoundaries[i+1] << ")$ & "
                        << std::fixed << std::setprecision(1)
                        << ErrSyst_SigExtr[i] << " & "
                        << ErrSyst_AxE[i] << " & "
                        << ErrSyst_fD[i] << " & "
                        << ErrSyst_fC[i] << R"( \\)" << "\n";
                            
    }
    fout_systErr_TeX.close();
    Printf("Results printed to %s.", str_out2b.Data()); 

    // 3a) Print the photonuclear cross section (also an input for PhenoPredictions.c)
    TString str_out3a = Form("Results/CrossSection/%ibins_FeedDown%i_photo.txt", nPtBins, iFeedDown);
    ofstream fout_sigmaPhoto(str_out3a.Data());
    fout_sigmaPhoto << std::fixed << std::setprecision(4);
    fout_sigmaPhoto << "Bin \tt_low \tt_upp \tsig \terr_sta\terr_syst\n";
    for(Int_t i = 0; i < nPtBins; i++){
        fout_sigmaPhoto << i+1 << std::fixed << std::setprecision(4)
                               << "\t" << ptBoundaries[i] * ptBoundaries[i] << "\t" 
                               << ptBoundaries[i+1] * ptBoundaries[i+1] << "\t" 
                               << std::fixed << std::setprecision(1)
                               << Sigma_photo_val[i] << "\t"
                               << Sigma_photo_err_stat[i] << "\t"
                               << Sigma_photo_err_syst[i] << "\n";
    }
    fout_sigmaPhoto.close();
    Printf("Results printed to %s.", str_out3a.Data()); 
    // 3b) Print the photonuclear cross section: TeX file
    TString str_out3b = Form("Results/CrossSection/%ibins_FeedDown%i_photo_TeX.txt", nPtBins, iFeedDown);
    ofstream fout_sigmaPhoto_TeX(str_out3b.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        fout_sigmaPhoto_TeX << std::fixed << std::setprecision(3)
                            << "$(" << ptBoundaries[i] * ptBoundaries[i] << "," << ptBoundaries[i+1] * ptBoundaries[i+1] << ")$ & "
                            << t_avg_val[i] << " &\t$"
                            << std::fixed << std::setprecision(1)
                            << Sigma_photo_val[i] << R"( \pm )" << Sigma_photo_err_stat[i] << R"((stat.) \pm )" << Sigma_photo_err_syst[i] << R"((syst.)$ \\)" << "\n";
    }
    fout_sigmaPhoto_TeX.close();
    Printf("Results printed to %s.", str_out3b.Data()); 

    return;
}

void PrintErr(TString str){
    Printf("ERR: file %s missing. Terminating.", str.Data());
    return;
}