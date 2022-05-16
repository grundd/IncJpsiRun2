// _CompareCountsPass1Pass3.C
// David Grund, Mar 06, 2022

#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include "TSystem.h"

void CompareCounts(Bool_t isPass1Calibrated, Bool_t isPass3Calibrated);

void _CompareCountsPass1Pass3()
{
    CompareCounts(kFALSE,kTRUE);

    CompareCounts(kFALSE,kFALSE);

    return;
}

void CompareCounts(Bool_t isPass1Calibrated = kFALSE, Bool_t isPass3Calibrated = kFALSE)
{
    gSystem->Exec("mkdir -p Results/_CompareCountsPass1Pass3/");

    Double_t counterMC_pass1[19] = { 0 };
    Double_t counterMC_pass3[19] = { 0 };
    TString cutsNames[14] = {"CCUP31","runNum","ADAofl","ADCofl","V0Aofl","V0Cofl",
                             "SPDmat","PIDmuo","rap<.8","eta<.8","oppCha","ms2245",
                             "pT3032","ms3032"};

    ifstream ifs;
    // MC
    // pass1
    TString str_MC_pass1 = "";
    if(!isPass1Calibrated) str_MC_pass1 = "Results/pass1_4bins/CountEvents_MC/cuts_numbersOnly.txt";
    else                   str_MC_pass1 = "Results/pass1_4bins_calibPID/CountEvents_MC/cuts_numbersOnly.txt";
    ifs.open(str_MC_pass1.Data());
    for(Int_t i = 0; i < 19; i++){
        ifs >> counterMC_pass1[i];
        Printf("%i: %.0f", i, counterMC_pass1[i]);
    } 
    ifs.close();

    // pass3
    TString str_MC_pass3 = "";
    if(!isPass3Calibrated) str_MC_pass3 = "Results/pass3_4bins/CountEvents_MC/cuts_numbersOnly.txt";
    else                   str_MC_pass3 = "Results/pass3_4bins_calibPID/CountEvents_MC/cuts_numbersOnly.txt";
    ifs.open(str_MC_pass3.Data());
    for(Int_t i = 0; i < 19; i++){
        ifs >> counterMC_pass3[i];
        Printf("%i: %.0f", i, counterMC_pass3[i]);
    } 
    ifs.close();

    // print the ratios to a new file
    TString str_p1 = "No";
    TString str_p3 = "No";
    if(isPass1Calibrated) str_p1 = "Yes";
    if(isPass3Calibrated) str_p3 = "Yes";
    TString str_out_MC = Form("Results/_CompareCountsPass1Pass3/MC_calib_p1%s_p3%s.txt", str_p1.Data(), str_p3.Data());
    ofstream outfile;
    outfile.open(str_out_MC.Data());
    outfile << "cutName nEv_p1\tnEv_p3\tratio\n";
    for(Int_t i = 3; i < 17; i++){
        Double_t ratio = counterMC_pass1[i] / counterMC_pass3[i];
        outfile << std::fixed << std::setprecision(0)
                << Form("%s", cutsNames[i-3].Data()) << "\t"
                << counterMC_pass1[i] << "\t" << counterMC_pass3[i] << "\t" 
                << std::fixed << std::setprecision(3)
                << ratio << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str_out_MC.Data());

    Double_t counterData_pass1[18] = { 0 };
    Double_t counterData_pass3[18] = { 0 };

    // data
    // pass1
    TString str_data_pass1 = "";
    if(!isPass1Calibrated) str_data_pass1 = "Results/pass1_4bins/CountEvents/cuts_numbersOnly.txt";
    else                   str_data_pass1 = "Results/pass1_4bins_calibPID/CountEvents/cuts_numbersOnly.txt";
    ifs.open(str_data_pass1.Data());
    for(Int_t i = 0; i < 18; i++){
        ifs >> counterData_pass1[i];
        Printf("%i: %.0f", i, counterData_pass1[i]);
    } 
    ifs.close();

    // pass3
    TString str_data_pass3 = "";
    if(!isPass3Calibrated) str_data_pass3 = "Results/pass3_4bins/CountEvents/cuts_numbersOnly.txt";
    else                   str_data_pass3 = "Results/pass3_4bins_calibPID/CountEvents/cuts_numbersOnly.txt";
    ifs.open(str_data_pass3.Data());
    for(Int_t i = 0; i < 18; i++){
        ifs >> counterData_pass3[i];
        Printf("%i: %.0f", i, counterData_pass3[i]);
    } 
    ifs.close();

    // print the ratios to a new file
    TString str_out_data = Form("Results/_CompareCountsPass1Pass3/data_calib_p1%s_p3%s.txt", str_p1.Data(), str_p3.Data());

    outfile.open(str_out_data.Data());
    outfile << "cutName nEv_p1\tnEv_p3\tratio\n";
    for(Int_t i = 3; i < 16; i++){
        Double_t ratio = counterData_pass1[i] / counterData_pass3[i];
        outfile << std::fixed << std::setprecision(0)
                << Form("%s", cutsNames[i-2].Data()) << "\t"
                << counterData_pass1[i] << "\t" << counterData_pass3[i] << "\t" 
                << std::fixed << std::setprecision(3)
                << ratio << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str_out_data.Data());

    return;
}