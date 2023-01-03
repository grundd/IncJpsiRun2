// PtFit_SystUncertainties.C
// David Grund, June 21, 2022

#include "PtFit_Utilities.h"

// values of fD returned from the fit:
// index 0 -> the 'allbins' range (=> fiducial cross section)
// remaining indices -> cross section in pT bins
Double_t fD_fromFit_val[6] = { 0 };
Double_t fD_fromFit_err[6] = { 0 };
// values of fD when the value of R is varied:
Double_t fD_16_val[6] = { 0 };
Double_t fD_16_err[6] = { 0 };
Double_t fD_20_val[6] = { 0 };
Double_t fD_20_err[6] = { 0 };
// values of fD when the dissociative shape is varied:
Double_t fD_DissLL_val[6] = { 0 }; // diss low low
Double_t fD_DissLL_err[6] = { 0 };
Double_t fD_DissUL_val[6] = { 0 }; // diss upp low
Double_t fD_DissUL_err[6] = { 0 };
Double_t fD_DissLU_val[6] = { 0 }; // diss upp low
Double_t fD_DissLU_err[6] = { 0 };
Double_t fD_DissUU_val[6] = { 0 }; // diss upp low
Double_t fD_DissUU_err[6] = { 0 };
// arrays showing the final errors
// (errors returned by the pT fit and the differences added in quadrature)
Double_t fD_err[6] = { 0 };

void PtFit_ReadResultsFromFile(TString sIn, Double_t val[6], Double_t err[6])
{
    ifstream ifs;
    Int_t i_bin;
    ifs.open(sIn.Data());
    // read data from the file
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)) {
            istringstream inStream(str);
            // skip first two lines
            if(i > 1) inStream >> i_bin >> val[i-2] >> err[i-2];
            i++;   
        }
    }
    ifs.close();
    return;
}

void PtFit_SystUncertainties(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning_PtFit();
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_SystUncertainties/");
    TString sIn;

    // load values of fD returned from the fit
    sIn = "Results/" + str_subfolder + "PtFit_NoBkg/RecSh4_fD0_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_fromFit_val,fD_fromFit_err);

    // *******************************************************************************************
    // Try various values of R
    // *******************************************************************************************

    for(Int_t ifD = -2; ifD <= 2; ifD++) PtFit_NoBkg_DoFit(4,5,ifD);

    // load values of fD with R = 0.16
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_fD-2_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_16_val,fD_16_err);

    // load values of fD with R = 0.20
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_fD2_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_20_val,fD_20_err);

    // calculate the differences between 0.16/0.18 and 0.18/0.20
    Double_t diff_low[6] = { 0 };
    Double_t diff_upp[6] = { 0 };
    Double_t diff_mean[6] = { 0 };
    for(Int_t i = 0; i < nPtBins+1; i++)
    {
        diff_low[i] = fD_fromFit_val[i] - fD_16_val[i];
        diff_upp[i] = fD_20_val[i] - fD_fromFit_val[i];
        diff_mean[i] = (diff_low[i] + diff_upp[i]) / 2;
    }

    // print the results
    ofstream outfile;
    TString str_out = "Results/" + str_subfolder + "PtFit_SystUncertainties/differences_valueOfR.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(2)
            << "18_val\t18_err\tdiff_l\tdiff_u\tdiff_m\n";
    for(Int_t i = 0; i < nPtBins+1; i++)
    {
        outfile << fD_fromFit_val[i] << "\t" << fD_fromFit_err[i] << "\t" << diff_low[i] << "\t" << diff_upp[i] << "\t" 
                << diff_mean[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());

    // *******************************************************************************************
    // Try various dissociative shapes
    // *******************************************************************************************

    for(Int_t iDiss = 6; iDiss < 10; iDiss++) PtFit_NoBkg_DoFit(4,iDiss);

    // load values of fD with diss low low
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_Diss6_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_DissLL_val,fD_DissLL_err);

    // load values of fD with diss upp low
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_Diss7_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_DissUL_val,fD_DissUL_err);

    // load values of fD with diss low upp
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_Diss8_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_DissLU_val,fD_DissLU_err);

    // load values of fD with diss upp upp
    sIn = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_Diss9_fD.txt";
    PtFit_ReadResultsFromFile(sIn,fD_DissUU_val,fD_DissUU_err);

    // calculate the differences between the values from the original fit and LL, UL, LU, UU
    Double_t diff_LL[6] = { 0 };
    Double_t diff_UL[6] = { 0 };
    Double_t diff_LU[6] = { 0 };
    Double_t diff_UU[6] = { 0 };
    Double_t diff_max[6] = { 0 };
    for(Int_t i = 0; i < nPtBins+1; i++)
    {
        diff_LL[i] = fD_fromFit_val[i] - fD_DissLL_val[i];
        diff_UL[i] = fD_fromFit_val[i] - fD_DissUL_val[i];
        diff_LU[i] = fD_fromFit_val[i] - fD_DissLU_val[i];
        diff_UU[i] = fD_fromFit_val[i] - fD_DissUU_val[i];
        // find the maximum difference for each bin
        Double_t maxDiff = 0;
        if(TMath::Abs(diff_LL[i]) > TMath::Abs(maxDiff)) maxDiff = diff_LL[i];
        if(TMath::Abs(diff_UL[i]) > TMath::Abs(maxDiff)) maxDiff = diff_UL[i];
        if(TMath::Abs(diff_LU[i]) > TMath::Abs(maxDiff)) maxDiff = diff_LU[i];
        if(TMath::Abs(diff_UU[i]) > TMath::Abs(maxDiff)) maxDiff = diff_UU[i];
        diff_max[i] = maxDiff;
    }

    // print the results
    str_out = "Results/" + str_subfolder + "PtFit_SystUncertainties/differences_dissShape.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(2)
            << "18_val\t18_err\tdiff_LL\tdiff_UL\tdiff_LU\tdiff_UU\tmaxDiff\n";
    for(Int_t i = 0; i < nPtBins+1; i++)
    {
        outfile << fD_fromFit_val[i] << "\t" << fD_fromFit_err[i] << "\t" << diff_LL[i] << "\t" << diff_UL[i] << "\t" 
                << diff_LU[i] << "\t" << diff_UU[i] << "\t" << diff_max[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());

    // *******************************************************************************************
    // Summarize the results
    // *******************************************************************************************

    for(Int_t i = 0; i < nPtBins+1; i++) 
    {
        fD_err[i] = TMath::Sqrt(TMath::Power(diff_mean[i], 2)); 
        //+ TMath::Power(diff_max[i], 2)
        //+ TMath::Power(fD_fromFit_err[i], 2)
    }

    // print the results to be read by CrossSec_Calculate()
    str_out = "Results/" + str_subfolder + "PtFit_SystUncertainties/fD_syst_errors.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(1);
    for(Int_t i = 0; i < nPtBins+1; i++)
    {
        outfile << fD_fromFit_val[i] << "\t" << fD_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());

    return;
}