// CrossSec_Calculate.C
// David Grund, Apr 05, 2022

// cpp headers
#include <fstream>  // print output to txt file
#include <iomanip>  // std::setprecision()
#include <string>   // getline
// root headers
#include "TMath.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

// values to calculate the UPC cross section:
// index 0 -> the 'allbins' range (=> fiducial cross section)
// remaining indices -> cross section in pT bins
Float_t N_yield_val[6] = { 0 };
Float_t N_yield_err[6] = { 0 };
Float_t pT2_widths[6] = { 0 };
Float_t AxE_MC_val[6] = { 0 };
Float_t AxE_MC_err[6] = { 0 };
Float_t fD_val[6] = { 0 };
Float_t fD_err[6] = { 0 };
Float_t fC_val[6] = { 0 };
Float_t fC_err[6] = { 0 };
Float_t avgt_val[6] = { 0 };
// UPC cross section:
Float_t sig_upc_val[6] = { 0 };
Float_t sig_upc_err_stat[6] = { 0 };
Float_t sig_upc_err_syst_uncr[6] = { 0 };
Float_t sig_upc_err_syst_corr[6] = { 0 };
// photonuclear cross section:
Float_t sig_gPb_val[6] = { 0 };
Float_t sig_gPb_err_stat[6] = { 0 };
Float_t sig_gPb_err_syst_uncr[6] = { 0 };
Float_t sig_gPb_err_syst_corr[6] = { 0 };
// systematic uncertainties (in percent):
Float_t syst_sig_extr[6] = { 0 };
Float_t syst_z_vtx[6] = { 0 };
Float_t syst_fD[6] = { 0 };
Float_t syst_fC[6] = { 0 };
Float_t syst_BR = 0.; // calculated later
Float_t syst_lumi = 2.7;
Float_t syst_eff_veto_pileup = 3.0;
Float_t syst_eff_veto_diss = 3.8;
Float_t syst_tracking = 2.8; // added quadratically (the PF committee suggested to add linearly)
Float_t syst_trig_eff = 1.3;
Float_t syst_flux = 2.0;
// global factors 
Float_t rap_width = 1.6;
Float_t lumi_val = 0; // 1/(mu barn), loaded later
Float_t lumi_err = 0; // 1/(mu barn)
Float_t BR_val = 0.05961;
Float_t BR_err = 0.00033;
Float_t flux_val = 84.9;
Float_t flux_err = flux_val * syst_flux / 100.;
Float_t eff_veto_pileup_val = 94.0;
Float_t eff_veto_pileup_err = eff_veto_pileup_val * syst_eff_veto_pileup / 100.;
Float_t eff_veto_diss_val = 63.7;
Float_t eff_veto_diss_err = eff_veto_diss_val * syst_eff_veto_diss / 100.;

// temporary variable used when loading data
Int_t i_bin;

void CalculateCrossSec_PtBins();
void PrintErr(TString str);

void CrossSec_Calculate(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/");

    CalculateCrossSec_PtBins();

    return;
}

void CalculateCrossSec_PtBins()
{
    SetPtBinning();

    // *******************************************************************************************
    // load values and calculate the UPC cross section and statistical errors
    // *******************************************************************************************
    Printf("Calculating the UPC cross section and statistical uncertainties...");
    ifstream ifs;

    // integrated lumi for both periods
    Float_t lumi_periods[2] = { 0 };
    TString str_period[2] = {"18q", "18r"};
    for(Int_t iPeriod = 0; iPeriod < 2; iPeriod++)
    {
        TString str_lumi = Form("Results/" + str_subfolder + "Lumi/lumi_%s.txt", str_period[iPeriod].Data());
        ifs.open(str_lumi.Data());
        if(!ifs.fail()) ifs >> lumi_periods[iPeriod];
        else {
            PrintErr(str_lumi);
            return;
        }
        ifs.close(); 
    }
    lumi_val = lumi_periods[0] + lumi_periods[1];
    lumi_err = lumi_val * syst_lumi / 100.;
    Printf("Loaded: integrated lumi (%.0f pm %.0f)", lumi_val, lumi_err);

    // yields 
    // total value + in pT bins
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        TString str_yield;
        if(iBin == 0) str_yield = "Results/" + str_subfolder + "InvMassFit/allbins/allbins_signal.txt"; 
        else          str_yield = "Results/" + str_subfolder + Form("InvMassFit/%ibins/bin%i_signal.txt",nPtBins,iBin);
        ifs.open(str_yield.Data());
        if(!ifs.fail()) ifs >> N_yield_val[iBin] >> N_yield_err[iBin];
        else {
            PrintErr(str_yield);
            return;
        }
        ifs.close(); 
    } 
    Printf("Loaded: yields");

    // AxE_MC
    // total value -> fiducial
    TString str_AxE = Form("Results/" + str_subfolder + "AxE_PtBins/AxE_%ibins_total.txt",nPtBins);
    ifs.open(str_AxE.Data());
    if(!ifs.fail()) ifs >> i_bin >> AxE_MC_val[0] >> AxE_MC_err[0];
    else {
        PrintErr(str_AxE);
        return;            
    }
    ifs.close();
    // in pT bins
    TString str_AxEs = Form("Results/" + str_subfolder + "AxE_PtBins/AxE_%ibins.txt",nPtBins);
    ifs.open(str_AxEs.Data());
    if(!ifs.fail()) for(Int_t iBin = 0; iBin < nPtBins; iBin++) ifs >> i_bin >> AxE_MC_val[iBin+1] >> AxE_MC_err[iBin+1];
    else {
        PrintErr(str_AxEs);
        return;            
    }
    ifs.close();
    Printf("Loaded: AxE_MC values");

    // fD corrections 
    // total value + in pT bins
    TString str_FDs = "Results/" + str_subfolder + "PtFit_SystUncertainties/fD_syst_errors.txt";
    ifs.open(str_FDs.Data());
    if(!ifs.fail()) for(Int_t i = 0; i < nPtBins+1; i++) ifs >> fD_val[i] >> fD_err[i];
    else {
        PrintErr(str_FDs);
        return;
    }
    ifs.close();
    Printf("Loaded: fD correction factors");

    // fC corrections 
    // total value + in pT bins
    TString str_FCs = "Results/" + str_subfolder + "PtFit_NoBkg/RecSh4_fD0_fC.txt";
    ifs.open(str_FCs.Data());
    if(!ifs.fail()) {
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)) {
            istringstream in_stream(str);
            // skip the first two lines
            if(i > 1) in_stream >> i_bin >> fC_val[i-2] >> fC_err[i-2];
            i++;   
        }
    } else {
        PrintErr(str_FCs);
        return;
    }
    ifs.close();
    Printf("Loaded: fC correction factors");

    // widths of intervals in pT^2 [GeV^2]
    pT2_widths[0] = ptBoundaries[nPtBins] - ptBoundaries[0];
    for(Int_t iBin = 1; iBin < nPtBins+1; iBin++)
    {
        pT2_widths[iBin] = TMath::Power(ptBoundaries[iBin], 2) - TMath::Power(ptBoundaries[iBin-1], 2);
    }
    Printf("Calculated: pT^2 widths");

    // cross-check: print the loaded values
    Printf("The following values will be used:");
    Printf("pT_low\tpT_upp\tpT2_w\tN_val\tN_err\tAxE_val\tAxE_err\tfD_val\tfD_err\tfC_val\tfC_err");
    for(Int_t iBin = 1; iBin < nPtBins+1; iBin++)
    {
        Printf("%.3f\t%.3f\t%.4f\t%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
            ptBoundaries[iBin-1], ptBoundaries[iBin], pT2_widths[iBin], N_yield_val[iBin], N_yield_err[iBin], 
            AxE_MC_val[iBin], AxE_MC_err[iBin], fD_val[iBin], fD_err[iBin], fC_val[iBin], fC_err[iBin]);
    }
    
    // calculate the UPC cross section in pT bins and statistical errors
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        sig_upc_val[iBin] = N_yield_val[iBin] / (
            (1.0 + fD_val[iBin] / 100. + fC_val[iBin] / 100.) * 
            (AxE_MC_val[iBin] / 100.) * 
            (eff_veto_pileup_val / 100.) * 
            (eff_veto_diss_val / 100.) * 
            (lumi_val * 1000) *
            BR_val * 
            rap_width * pT2_widths[iBin] );
        sig_upc_err_stat[iBin] = sig_upc_val[iBin] * TMath::Sqrt(TMath::Power(N_yield_err[iBin] / N_yield_val[iBin], 2)
            + TMath::Power(AxE_MC_err[iBin] / AxE_MC_val[iBin], 2));
    }
    Printf("Calculated: UPC cross section and stat errors");

    // *******************************************************************************************
    // load values and calculate systematic uncertainties
    // *******************************************************************************************
    Printf("Calculating systematic uncertainties...");

    // signal extraction
    TString str_sig_extr = "Results/" + str_subfolder + Form("InvMassFit_SystUncertainties/ErrSystSignalExtraction_%ibins.txt",nPtBins);
    ifs.open(str_sig_extr.Data());
    if(!ifs.fail()) for(Int_t iBin = 0; iBin < nPtBins+1; iBin++) ifs >> i_bin >> syst_sig_extr[iBin];
    else {
        PrintErr(str_sig_extr);
        return;        
    }
    ifs.close();
    Printf("Loaded: signal extraction syst errs");

    // |z_vtx| selection
    TString str_SystZVtx = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/syst_uncertainties_%ibins.txt",nPtBins);
    ifs.open(str_SystZVtx.Data());
    if(!ifs.fail()) for(Int_t iBin = 0; iBin < nPtBins+1; iBin++) ifs >> i_bin >> syst_z_vtx[iBin];
    else {
        PrintErr(str_SystZVtx);
        return;        
    }
    ifs.close();
    Printf("Loaded: vtx selection syst errs");

    // fC and fD
    Float_t fD_upp[6] = { 0 };
    Float_t fD_low[6] = { 0 };
    Float_t fC_upp[6] = { 0 };
    Float_t fC_low[6] = { 0 };
    Float_t sig_fD_upp[6] = { 0 };
    Float_t sig_fD_low[6] = { 0 };
    Float_t sig_fC_upp[6] = { 0 };
    Float_t sig_fC_low[6] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        fD_upp[iBin] = fD_val[iBin] + fD_err[iBin];
        fD_low[iBin] = fD_val[iBin] - fD_err[iBin];
        fC_upp[iBin] = fC_val[iBin] + fC_err[iBin];
        fC_low[iBin] = fC_val[iBin] - fC_err[iBin];
        sig_fD_upp[iBin] = sig_upc_val[iBin] * (1.0 + fD_val[iBin] / 100. + fC_val[iBin] / 100.) / (1.0 + fD_upp[iBin] / 100. + fC_val[iBin] / 100.);
        sig_fD_low[iBin] = sig_upc_val[iBin] * (1.0 + fD_val[iBin] / 100. + fC_val[iBin] / 100.) / (1.0 + fD_low[iBin] / 100. + fC_val[iBin] / 100.);
        sig_fC_upp[iBin] = sig_upc_val[iBin] * (1.0 + fD_val[iBin] / 100. + fC_val[iBin] / 100.) / (1.0 + fD_val[iBin] / 100. + fC_upp[iBin] / 100.);
        sig_fC_low[iBin] = sig_upc_val[iBin] * (1.0 + fD_val[iBin] / 100. + fC_val[iBin] / 100.) / (1.0 + fD_val[iBin] / 100. + fC_low[iBin] / 100.);
        Float_t sig_fD_upp_diff, sig_fD_low_diff, sig_fC_upp_diff, sig_fC_low_diff;
        sig_fD_upp_diff = TMath::Abs(sig_fD_upp[iBin] - sig_upc_val[iBin]);
        sig_fD_low_diff = TMath::Abs(sig_fD_low[iBin] - sig_upc_val[iBin]);
        sig_fC_upp_diff = TMath::Abs(sig_fC_upp[iBin] - sig_upc_val[iBin]);
        sig_fC_low_diff = TMath::Abs(sig_fC_low[iBin] - sig_upc_val[iBin]);
        syst_fD[iBin] = TMath::Max(sig_fD_upp_diff / sig_upc_val[iBin], sig_fD_low_diff / sig_upc_val[iBin]) * 100.;
        syst_fC[iBin] = TMath::Max(sig_fC_upp_diff / sig_upc_val[iBin], sig_fC_low_diff / sig_upc_val[iBin]) * 100.;
    }
    Printf("Calculated: fC and fD syst errs");

    // calculate systematic errors of the UPC cross section
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        syst_BR = BR_err / BR_val * 100.;
        sig_upc_err_syst_corr[iBin] = sig_upc_val[iBin] * TMath::Sqrt(
            TMath::Power(syst_fD[iBin] / 100., 2) +
            TMath::Power(syst_fC[iBin] / 100., 2) +
            TMath::Power(syst_lumi / 100., 2) + 
            TMath::Power(syst_eff_veto_pileup / 100., 2) + 
            TMath::Power(syst_eff_veto_diss / 100., 2) + 
            TMath::Power(syst_tracking / 100., 2) + 
            TMath::Power(syst_trig_eff / 100., 2) + 
            TMath::Power(syst_BR / 100., 2)
        );
        sig_upc_err_syst_uncr[iBin] = sig_upc_val[iBin] * TMath::Sqrt(
            TMath::Power(syst_sig_extr[iBin] / 100., 2) +
            TMath::Power(syst_z_vtx[iBin] / 100., 2)
        );
    }
    Printf("Calculated: UPC cross section syst errs");

    // *******************************************************************************************
    // load values and calculate the photonuclear cross section
    // *******************************************************************************************
    Printf("Calculating the photonuclear cross section and its errors...");

    // calculate the photonuclear cross section in pT bins
    for(Int_t iBin = 0; iBin < nPtBins+1; iBin++)
    {
        sig_gPb_val[iBin] = sig_upc_val[iBin] / 2. / flux_val * 1000;
        sig_gPb_err_stat[iBin] = sig_upc_err_stat[iBin] / 2. / flux_val * 1000;
        sig_gPb_err_syst_uncr[iBin] = sig_upc_err_syst_uncr[iBin] / 2. / flux_val * 1000;
        sig_gPb_err_syst_corr[iBin] = sig_gPb_val[iBin] * TMath::Sqrt(
            TMath::Power(sig_upc_err_syst_corr[iBin] / sig_upc_val[iBin], 2) + 
            TMath::Power(syst_flux / 100., 2)
        );
    }
    Printf("Calculated: photonuclear cross section and errors");

    // avg values of |t| per bin
    TString str_t_avg = "Results/" + str_subfolder + "STARlight_tVsPt/AvgTPerBin.txt";
    ifs.open(str_t_avg.Data());
    if(!ifs.fail()) for(Int_t iBin = 0; iBin < nPtBins; iBin++) { 
        ifs >> i_bin >> avgt_val[iBin];  
        if(kFALSE) Printf("Reading: bin %i, |t| = %.4f", i_bin, avgt_val[iBin]);
    } else {
        PrintErr(str_t_avg);
        return;
    }
    ifs.close();
    Printf("Loaded: average |t| values");

    // print the results to text files
    // print the UPC cross section 
    TString str_out = "Results/" + str_subfolder + "CrossSec/CrossSec_UPC.txt";
    ofstream outfile(str_out.Data());
    outfile << Form("Lumi\terr\tRapW\tBR\terr\te_p-up\terr\te_diss\terr\tflux\terr\n")
                  << Form("%.1f \t%.1f \t%.1f \t%.3f \t%.3f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n\n",
                            lumi_val, lumi_err, 
                            rap_width, 
                            BR_val*100., BR_err*100., 
                            eff_veto_pileup_val, eff_veto_pileup_err, 
                            eff_veto_diss_val, eff_veto_diss_err,
                            flux_val, flux_err);
    outfile << Form("Bin\tPt2Low\tPt2Upp\tPt2_W\tN\terr\tAxE\terr\tFD [%%]\terr\tFC [%%]\terr\tsig\tstat\tsyst u\tsyst c\n");
    for(Int_t i = 1; i <= nPtBins; i++) {
        outfile << std::fixed << std::setprecision(3)
                << i << "\t"
                << ptBoundaries[i-1] * ptBoundaries[i-1] << "\t"
                << ptBoundaries[i] * ptBoundaries[i] << "\t"
                << std::fixed << std::setprecision(4) << pT2_widths[i] << "\t"
                << std::fixed << std::setprecision(1) << N_yield_val[i] << "\t" << N_yield_err[i] << "\t"
                << std::fixed << std::setprecision(2) << AxE_MC_val[i] << "\t" << AxE_MC_err[i] << "\t"
                << std::fixed << std::setprecision(1) << fD_val[i] << "\t" << fD_err[i] << "\t"
                << std::fixed << std::setprecision(3) << fC_val[i] << "\t" << fC_err[i] << "\t"
                << std::fixed << std::setprecision(2) 
                << sig_upc_val[i] << "\t" << sig_upc_err_stat[i] << "\t" << sig_upc_err_syst_uncr[i] << "\t" << sig_upc_err_syst_corr[i] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", str_out.Data()); 

    // print the UPC cross section: TeX table
    str_out = "Results/" + str_subfolder + "CrossSec/CrossSec_UPC_TeX.txt";
    outfile.open(str_out.Data());
    outfile << Form("$%.0f", lumi_val) << R"( \pm )" << Form("%.0f$", lumi_err) << " &\n" 
            << Form("%.1f", rap_width) << " &\n"
            << Form("$%.3f", BR_val*100.) << R"( \pm )" << Form("%.3f$", BR_err*100.) << " &\n"
            << Form("$%.1f", eff_veto_pileup_val) << R"( \pm )"<< Form("%.1f$", eff_veto_pileup_err)<< " &\n"
            << Form("$%.1f", eff_veto_diss_val) << R"( \pm )" << Form("%.1f$", eff_veto_diss_err) << " &\n"
            << Form("$%.1f", flux_val) << R"( \pm )" << Form("%.1f$", flux_err) << R"( \\)" << "\n\n";
    for(Int_t i = 1; i <= nPtBins; i++) {
        outfile << std::fixed << std::setprecision(3) << "$(" 
                << ptBoundaries[i-1] * ptBoundaries[i-1] << "," 
                << ptBoundaries[i] * ptBoundaries[i] << ")$ & "
                << std::fixed << std::setprecision(4) << pT2_widths[i] << " &\t$"
                << std::fixed << std::setprecision(0) << N_yield_val[i] << R"( \pm )" << N_yield_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(2) << AxE_MC_val[i] << R"( \pm )" << AxE_MC_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(1) << fD_val[i] << R"( \pm )" << fD_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(3) << fC_val[i] << R"( \pm )" << fC_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(2) 
                << sig_upc_val[i] << R"( \pm )" << sig_upc_err_stat[i] << R"( \pm )" << sig_upc_err_syst_uncr[i] << R"( \pm )" << sig_upc_err_syst_corr[i] << R"($ \\)" << "\n";
    }                  
    outfile.close();
    Printf("Results printed to %s.", str_out.Data());

    // print the systematic uncertainties
    str_out = "Results/" + str_subfolder + "CrossSec/Systematics.txt";
    outfile.open(str_out.Data());
    outfile << "[all in percent]\n"
            << "CORRELATED:\n"
            << "lumi\tveto\tEMD\ttracks\tCCUP31\tBR\n"
            << Form("%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n\n",
                syst_lumi, syst_eff_veto_pileup, syst_eff_veto_diss, syst_tracking, syst_trig_eff, syst_BR);
    outfile << "UNCORRELATED:\n"
            << "Bin\tSigExt\tZVtx\tfD\tfC\n";
    for(Int_t i = 1; i <= nPtBins; i++) {
        outfile << i << std::fixed << std::setprecision(1) << "\t"
                << syst_sig_extr[i] << "\t"
                << syst_z_vtx[i] << "\t"
                << syst_fD[i] << "\t"
                << syst_fC[i] << "\n";
    }    
    outfile.close();
    Printf("Results printed to %s.", str_out.Data());

    // print the systematic uncertainties: TeX table
    str_out = "Results/" + str_subfolder + "CrossSec/Systematics_TeX.txt";
    outfile.open(str_out.Data());
    for(Int_t i = 1; i <= nPtBins; i++) {
        outfile << std::fixed << std::setprecision(3)
                << "$(" << ptBoundaries[i-1] << "," << ptBoundaries[i] << ")$ & "
                << std::fixed << std::setprecision(1)
                << syst_sig_extr[i] << " & "
                << syst_z_vtx[i] << " & "
                << syst_fD[i] << " & "
                << syst_fC[i] << R"( \\)" << "\n";
                            
    }
    outfile.close();
    Printf("Results printed to %s.", str_out.Data()); 

    // print the photonuclear cross section
    str_out = "Results/" + str_subfolder + "CrossSec/CrossSec_photo.txt";
    outfile.open(str_out.Data());
    //outfile << "Bin \tt_low \tt_upp \tsig \tstat\tsyst u\tsyst c\n";
    for(Int_t i = 1; i <= nPtBins; i++) {
        outfile << i << std::fixed << std::setprecision(4) << "\t" 
                << ptBoundaries[i-1] * ptBoundaries[i-1] << "\t" 
                << ptBoundaries[i] * ptBoundaries[i] << "\t" 
                << std::fixed << std::setprecision(2)
                << sig_gPb_val[i] << "\t"
                << sig_gPb_err_stat[i] << "\t"
                << sig_gPb_err_syst_uncr[i] << "\t"
                << sig_gPb_err_syst_corr[i] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", str_out.Data()); 

    // print the photonuclear cross section: TeX table
    str_out = "Results/" + str_subfolder + "CrossSec/CrossSec_photo_TeX.txt";
    outfile.open(str_out.Data());
    for(Int_t i = 1; i <= nPtBins; i++) {
        outfile << std::fixed << std::setprecision(3) << "$(" 
                << ptBoundaries[i-1] * ptBoundaries[i-1] << "," 
                << ptBoundaries[i] * ptBoundaries[i] << ")$ & "
                << avgt_val[i] << " &\t$"
                << std::fixed << std::setprecision(2)
                << sig_gPb_val[i] << R"( \pm )" << sig_gPb_err_stat[i] << R"( \pm )" << sig_gPb_err_syst_uncr[i] << R"( \pm )" << sig_gPb_err_syst_corr[i] << R"($ \\)" << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", str_out.Data()); 

    return;
}

void PrintErr(TString str)
{
    Printf("ERR: file %s missing. Terminating.", str.Data());
    return;
}