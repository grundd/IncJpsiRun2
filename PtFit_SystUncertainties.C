// PtFit_SystUncertainties.C
// David Grund, June 21, 2022

#include "PtFit_Utilities.h"

// #############################################################################################

Double_t CorrFD_16_val[5] = { 0 };
Double_t CorrFD_16_err[5] = { 0 };
Double_t CorrFD_18_val[5] = { 0 };
Double_t CorrFD_18_err[5] = { 0 };
Double_t CorrFD_20_val[5] = { 0 };
Double_t CorrFD_20_err[5] = { 0 };
// arrays showing the final errors
// (errors returned by the pT fit and the differences added in quadrature)
Double_t CorrFD_err[5] = { 0 };

void PtFit_SystUncertainties(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning_PtFit();
    SetPtBinning();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_SystUncertainties/");

    // Try various values of R
    for(Int_t ifD = -2; ifD <= 2; ifD++) PtFit_NoBkg_DoFit(4,ifD);

    ifstream ifs;
    TString str_in;
    Int_t i_bin;
    // load values of fD with R = 0.16
    str_in = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_fD-2_fD.txt";
    ifs.open(str_in.Data());
    // Read data from the file
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFD_16_val[i-1] >> CorrFD_16_err[i-1];
            i++;   
        }
    }
    ifs.close();  
    // load values of fD with R = 0.18
    str_in = "Results/" + str_subfolder + "PtFit_NoBkg/RecSh4_fD0_fD.txt";
    ifs.open(str_in.Data());
    // Read data from the file
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFD_18_val[i-1] >> CorrFD_18_err[i-1];
            i++;   
        }
    }
    ifs.close();  
    // load values of fD with R = 0.20
    str_in = "Results/" + str_subfolder + "PtFit_SystUncertainties/RecSh4_fD2_fD.txt";
    ifs.open(str_in.Data());
    // Read data from the file
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFD_20_val[i-1] >> CorrFD_20_err[i-1];
            i++;   
        }
    }
    ifs.close();
    // calculate the differences between 0.16/0.18 and 0.18/0.20
    Double_t diff_low[5] = { 0 };
    Double_t diff_upp[5] = { 0 };
    Double_t diff_mean[5] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++)
    {
        diff_low[i] = CorrFD_18_val[i] - CorrFD_16_val[i];
        diff_upp[i] = CorrFD_20_val[i] - CorrFD_18_val[i];
        diff_mean[i] = (diff_low[i] + diff_upp[i]) / 2;
        CorrFD_err[i] = TMath::Sqrt(TMath::Power(diff_mean[i], 2)); //+ TMath::Power(CorrFD_18_err[i], 2)
    }
    // print the results
    ofstream outfile;
    TString str_out = "Results/" + str_subfolder + "PtFit_SystUncertainties/differences.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(2)
            << "18_val\t18_err\tdiff_l\tdiff_u\tdiff_m\tERROR\n";
    for(Int_t i = 0; i < nPtBins; i++)
    {
        outfile << CorrFD_18_val[i] << "\t" << CorrFD_18_err[i] << "\t" << diff_low[i] << "\t" << diff_upp[i] << "\t" 
                << diff_mean[i] << "\t" << CorrFD_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());
    // print the results to be read by PhotoCrossSec_Calculate()
    str_out = "Results/" + str_subfolder + "PtFit_SystUncertainties/fD_syst_errors.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(1);
    for(Int_t i = 0; i < nPtBins; i++)
    {
        outfile << CorrFD_18_val[i] << "\t" << CorrFD_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());
    return;
}