// PhotoCrossSec_Calculate.C
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

//*************************************************
// Cross section in pT^2 bins
Double_t N_yield_val[5] = { 0 };
Double_t N_yield_err[5] = { 0 };
Double_t Pt2Widths[5] = { 0 };
Double_t AxE_val[5] = { 0 };
Double_t AxE_err[5] = { 0 };
Double_t CorrFD_val[5] = { 0 };
Double_t CorrFD_err[5] = { 0 };
Double_t CorrFC_val[5] = { 0 };
Double_t CorrFC_err[5] = { 0 };
Double_t Sigma_UPC_val[5] = { 0 };
Double_t Sigma_UPC_err_stat[5] = { 0 };
Double_t Sigma_UPC_err_syst[5] = { 0 };
Double_t Sigma_photo_val[5] = { 0 };
Double_t Sigma_photo_err_stat[5] = { 0 };
Double_t Sigma_photo_err_syst[5] = { 0 };
Double_t t_avg_val[5] = { 0 };
//*************************************************
// Systematic uncertainties (in percent)
Double_t ErrSyst_SigExtr[5] = { 0 };
Double_t ErrSyst_ZVertex[5] = { 0 };
Double_t ErrSyst_fD[5] = { 0 };
Double_t ErrSyst_fC[5] = { 0 };
Double_t ErrSyst_lumi = 2.7;
Double_t ErrSyst_veto = 3.0;
Double_t ErrSyst_EMD = 3.7;
Double_t ErrSyst_tracks = 2.8;
Double_t ErrSyst_CCUP31 = 1.3;
Double_t ErrSyst_flux = 2.0;
//*************************************************
Double_t Lumi_val = 0; // 1/(mu barn)
Double_t Lumi_err = 0; // 1/(mu barn)
Double_t BR_val = 0.05961;
Double_t BR_err = 0.00033;
Double_t ErrSyst_BR = BR_err / BR_val * 100.;
Double_t RapWidth = 1.6;
Double_t Eff_veto_val = 94.0;
Double_t Eff_veto_err = Eff_veto_val * ErrSyst_veto / 100.;
Double_t Eff_EMD_val = 63.7;
Double_t Eff_EMD_err = Eff_EMD_val * ErrSyst_EMD / 100.;
Double_t PhotonFlux_val = 84.9;
Double_t PhotonFlux_err = PhotonFlux_val * ErrSyst_flux / 100.;
//*************************************************

Int_t i_bin;

void CalculateCrossSec_PtBins();
void PrintErr(TString str);

void PhotoCrossSec_Calculate(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/");

    CalculateCrossSec_PtBins();

    return;
}

void CalculateCrossSec_PtBins()
{
    SetPtBinning();

    ifstream ifs;

    Printf("Calculating the (photonuclear) cross section in %i pT bins.", nPtBins);

    //#####################################################################################################
    // 1) Load integrated luminosity for both periods
    Double_t Lumi_periods[2] = { 0 };
    TString str_period[2] = {"18q", "18r"};
    for(Int_t iPeriod = 0; iPeriod < 2; iPeriod++)
    {
        TString str_lumi = Form("Results/" + str_subfolder + "Lumi/lumi_%s.txt", str_period[iPeriod].Data());
        ifs.open(str_lumi.Data());
        // Read data from the file
        if(!ifs.fail()){
            ifs >> Lumi_periods[iPeriod];
        } else {
            PrintErr(str_lumi);
            return;
        }
        ifs.close(); 
    }
    Lumi_val = Lumi_periods[0] + Lumi_periods[1];
    Lumi_err = Lumi_val * ErrSyst_lumi / 100.;

    Printf("1) Integrated lumi loaded.");

    //#####################################################################################################
    // 2) Load N_yield per pT bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        TString str_yield = Form("Results/" + str_subfolder + "InvMassFit/%ibins/bin%i_signal.txt", nPtBins, iBin+1);
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
    Printf("2) N_yield loaded.");

    //#####################################################################################################
    // 3) Load AxE per pT bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        TString str_AxE = Form("Results/" + str_subfolder + "AxE_PtBins/AxE_%ibins.txt", nPtBins);
        ifs.open(str_AxE.Data());
        // Read data from the file
        if(!ifs.fail()){
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                ifs >> i_bin >> AxE_val[iBin] >> AxE_err[iBin];
            }
        } else {
            PrintErr(str_AxE);
            return;            
        }
        ifs.close();
    }
    Printf("3) AxE loaded.");

    //#####################################################################################################
    // 4) Load FD corrections per pT bin
    TString str_FD = "Results/" + str_subfolder + "PtFit_SystUncertainties/fD_syst_errors.txt";
    ifs.open(str_FD.Data());
    // Read data from the file
    if(!ifs.fail()){
        for(Int_t i = 0; i < nPtBins; i++)
        {
            ifs >> CorrFD_val[i] >> CorrFD_err[i];
        }
    } else {
        PrintErr(str_FD);
        return;
    }
    ifs.close();
    Printf("4) FD corrections loaded.");

    //#####################################################################################################
    // 5) Load FC corrections per pT bin
    TString str_FC = "Results/" + str_subfolder + "PtFit_NoBkg/RecSh4_fD0_fC.txt";
    ifs.open(str_FC.Data());
    // Read data from the file
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFC_val[i-1] >> CorrFC_err[i-1];
            i++;   
        }
    } else {
        PrintErr(str_FC);
        return;
    }
    ifs.close();
    Printf("5) FC corrections loaded.");

    //#####################################################################################################
    // 6) Widths of intervals in pT^2 [GeV^2]
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        Pt2Widths[iBin] = TMath::Power(ptBoundaries[iBin+1], 2) - TMath::Power(ptBoundaries[iBin], 2);
    }
    Printf("6) Widths of pT^2 bins calculated.");

    //#####################################################################################################
    // Cross-check: print the loaded values
    Printf("Lumi: %.2f pm %.2f", Lumi_val, Lumi_err);
    Printf("pT_low\tpT_upp\tpT2_w\tN_val\tN_err\tAxE_val\tAxE_err\tfD_val\tfD_err\tfC_val\tfC_err");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        Printf("%.3f\t%.3f\t%.4f\t%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
            ptBoundaries[iBin], ptBoundaries[iBin+1], Pt2Widths[iBin], N_yield_val[iBin], N_yield_err[iBin], 
            AxE_val[iBin], AxE_err[iBin], CorrFD_val[iBin], CorrFD_err[iBin], CorrFC_val[iBin], CorrFC_err[iBin]);
    }

    //#####################################################################################################
    // Calculate the UPC cross section per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        Sigma_UPC_val[iBin] = N_yield_val[iBin] / (
            (1.0 + CorrFD_val[iBin] / 100. + CorrFC_val[iBin] / 100.) * 
            (AxE_val[iBin] / 100.) * 
            (Eff_veto_val / 100.) * 
            (Eff_EMD_val / 100.) * 
            (Lumi_val * 1000) *
            BR_val * 
            RapWidth * Pt2Widths[iBin] );
        Sigma_UPC_err_stat[iBin] = Sigma_UPC_val[iBin] * TMath::Sqrt(TMath::Power(N_yield_err[iBin] / N_yield_val[iBin], 2)
            + TMath::Power(AxE_err[iBin] / AxE_val[iBin], 2));
    }

    //#####################################################################################################
    // Systematic uncertainties: SIGNAL EXTRACTION
    TString str_SystUncr = "Results/" + str_subfolder + Form("InvMassFit_SystUncertainties/ErrSystSignalExtraction_%ibins.txt", nPtBins);
    ifs.open(str_SystUncr.Data());
    // Read data from the file
    if(!ifs.fail()){
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            ifs >> i_bin >> ErrSyst_SigExtr[iBin];
        }
    } else {
        PrintErr(str_SystUncr);
        return;        
    }
    ifs.close();

    //#####################################################################################################
    // Systematic uncertainties: Z_VERTEX SELECTION
    TString str_SystZVtx = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/syst_uncertainties_%ibins.txt", nPtBins);
    ifs.open(str_SystZVtx.Data());
    // Read data from the file
    if(!ifs.fail()){
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            ifs >> i_bin >> ErrSyst_ZVertex[iBin];
        }
    } else {
        PrintErr(str_SystZVtx);
        return;        
    }
    ifs.close();

    //#####################################################################################################
    // Systematic uncertainties: FC AND FD
    Double_t CorrFD_upp[5] = { 0 };
    Double_t CorrFD_low[5] = { 0 };
    Double_t CorrFC_upp[5] = { 0 };
    Double_t CorrFC_low[5] = { 0 };
    Double_t Sigma_FD_upp[5] = { 0 };
    Double_t Sigma_FD_low[5] = { 0 };
    Double_t Sigma_FC_upp[5] = { 0 };
    Double_t Sigma_FC_low[5] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        CorrFD_upp[iBin] = CorrFD_val[iBin] + CorrFD_err[iBin];
        CorrFD_low[iBin] = CorrFD_val[iBin] - CorrFD_err[iBin];
        CorrFC_upp[iBin] = CorrFC_val[iBin] + CorrFC_err[iBin];
        CorrFC_low[iBin] = CorrFC_val[iBin] - CorrFC_err[iBin];
        Sigma_FD_upp[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_upp[iBin] / 100. + CorrFC_val[iBin] / 100.);
        Sigma_FD_low[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_low[iBin] / 100. + CorrFC_val[iBin] / 100.);
        Sigma_FC_upp[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_val[iBin] / 100. + CorrFC_upp[iBin] / 100.);
        Sigma_FC_low[iBin] = Sigma_UPC_val[iBin] * (1.0 + CorrFD_val[iBin] / 100. + CorrFC_val[iBin] / 100.) / (1.0 + CorrFD_val[iBin] / 100. + CorrFC_low[iBin] / 100.);
        Double_t Sigma_fD_upp_diff, Sigma_fD_low_diff, Sigma_fC_upp_diff, Sigma_fC_low_diff;
        Sigma_fD_upp_diff = TMath::Abs(Sigma_FD_upp[iBin] - Sigma_UPC_val[iBin]);
        Sigma_fD_low_diff = TMath::Abs(Sigma_FD_low[iBin] - Sigma_UPC_val[iBin]);
        Sigma_fC_upp_diff = TMath::Abs(Sigma_FC_upp[iBin] - Sigma_UPC_val[iBin]);
        Sigma_fC_low_diff = TMath::Abs(Sigma_FC_low[iBin] - Sigma_UPC_val[iBin]);
        ErrSyst_fD[iBin] = TMath::Max(Sigma_fD_upp_diff / Sigma_UPC_val[iBin], Sigma_fD_low_diff / Sigma_UPC_val[iBin]) * 100.;
        ErrSyst_fC[iBin] = TMath::Max(Sigma_fC_upp_diff / Sigma_UPC_val[iBin], Sigma_fC_low_diff / Sigma_UPC_val[iBin]) * 100.;

        // Total systematic uncertainty of sigma UPC
        Sigma_UPC_err_syst[iBin] = Sigma_UPC_val[iBin] * TMath::Sqrt(
            TMath::Power(ErrSyst_SigExtr[iBin] / 100., 2) +
            TMath::Power(ErrSyst_ZVertex[iBin] / 100., 2) +
            TMath::Power(ErrSyst_fD[iBin] / 100., 2) +
            TMath::Power(ErrSyst_fC[iBin] / 100., 2) +
            TMath::Power(ErrSyst_lumi / 100., 2) + 
            TMath::Power(ErrSyst_veto / 100., 2) + 
            TMath::Power(ErrSyst_EMD / 100., 2) + 
            TMath::Power(ErrSyst_tracks / 100., 2) + 
            TMath::Power(ErrSyst_CCUP31 / 100., 2) + 
            TMath::Power(ErrSyst_BR / 100., 2)
        );
    }

    //#####################################################################################################
    // Calculate the photonuclear cross section per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        Sigma_photo_val[iBin] = Sigma_UPC_val[iBin] / 2. / PhotonFlux_val * 1000;
        Sigma_photo_err_stat[iBin] = Sigma_UPC_err_stat[iBin] / 2. / PhotonFlux_val * 1000;
        Sigma_photo_err_syst[iBin] = Sigma_photo_val[iBin] * TMath::Sqrt(
            TMath::Power(Sigma_UPC_err_syst[iBin] / Sigma_UPC_val[iBin], 2) + 
            TMath::Power(ErrSyst_flux / 100., 2)
        );
    }

    //#####################################################################################################
    // Load avg values of |t| per bin
    TString str_t_avg = "Results/" + str_subfolder + "STARlight_tVsPt/AvgTPerBin.txt";
    ifs.open(str_t_avg.Data()); 
    // Read data from the file
    if(!ifs.fail()){
        for(Int_t iBin = 0; iBin < nPtBins; iBin++)
        {
            Int_t bin;
            ifs >> bin >> t_avg_val[iBin];
            if(kFALSE) Printf("Reading: bin %i, |t| = %.4f", bin, t_avg_val[iBin]);
        }
    } else {
        PrintErr(str_t_avg);
        return;
    }
    ifs.close();
    Printf("Values of an avg |t| value per bin loaded.");

    //#####################################################################################################
    // Print the results to text files
    // 1a) Print the UPC cross section 
    TString str_out_1a = "Results/" + str_subfolder + "PhotoCrossSec/CrossSec_UPC.txt";
    ofstream outfile(str_out_1a.Data());
    outfile << Form("Lumi\terr\tRapW\tBR\terr\te_veto\terr\te_EMD\terr\tflux\terr\n")
                  << Form("%.1f \t%.1f \t%.1f \t%.3f \t%.3f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n\n",
                            Lumi_val, Lumi_err, 
                            RapWidth, 
                            BR_val*100., BR_err*100., 
                            Eff_veto_val, Eff_veto_err, 
                            Eff_EMD_val, Eff_EMD_err,
                            PhotonFlux_val, PhotonFlux_err);
    outfile << Form("Bin\tPt2Low\tPt2Upp\tPt2_W\tN\terr\tAxE\terr\tFD [%%]\terr\tFC [%%]\terr\tsig\tstat\tsyst\n");
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(3)
                << i+1 << "\t"
                << ptBoundaries[i] * ptBoundaries[i] << "\t"
                << ptBoundaries[i+1] * ptBoundaries[i+1] << "\t"
                << std::fixed << std::setprecision(4) << Pt2Widths[i] << "\t"
                << std::fixed << std::setprecision(1) << N_yield_val[i] << "\t" << N_yield_err[i] << "\t"
                << std::fixed << std::setprecision(2) << AxE_val[i] << "\t" << AxE_err[i] << "\t"
                << std::fixed << std::setprecision(1) << CorrFD_val[i] << "\t" << CorrFD_err[i] << "\t"
                << std::fixed << std::setprecision(3) << CorrFC_val[i] << "\t" << CorrFC_err[i] << "\t"
                << std::fixed << std::setprecision(2) << Sigma_UPC_val[i] << "\t" << Sigma_UPC_err_stat[i] << "\t" << Sigma_UPC_err_syst[i] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", str_out_1a.Data()); 

    // 1b) Print the UPC cross section: TeX table
    TString str_out_1b = "Results/" + str_subfolder + "PhotoCrossSec/CrossSec_UPC_TeX.txt";
    outfile.open(str_out_1b.Data());
    outfile << Form("$%.0f", Lumi_val) << R"( \pm )" << Form("%.0f$", Lumi_err) << " &\n" 
            << Form("%.1f", RapWidth) << " &\n"
            << Form("$%.3f", BR_val*100.) << R"( \pm )" << Form("%.3f$", BR_err*100.) << " &\n"
            << Form("$%.1f", Eff_veto_val) << R"( \pm )"<< Form("%.1f$", Eff_veto_err)<< " &\n"
            << Form("$%.1f", Eff_EMD_val) << R"( \pm )" << Form("%.1f$", Eff_EMD_err) << " &\n"
            << Form("$%.1f", PhotonFlux_val) << R"( \pm )" << Form("%.1f$", PhotonFlux_err) << R"( \\)" << "\n\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(3) << "$(" 
                << ptBoundaries[i] * ptBoundaries[i] << "," 
                << ptBoundaries[i+1] * ptBoundaries[i+1] << ")$ & "
                << std::fixed << std::setprecision(4) << Pt2Widths[i] << " &\t$"
                << std::fixed << std::setprecision(0) << N_yield_val[i] << R"( \pm )" << N_yield_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(2) << AxE_val[i] << R"( \pm )" << AxE_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(1) << CorrFD_val[i] << R"( \pm )" << CorrFD_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(3) << CorrFC_val[i] << R"( \pm )" << CorrFC_err[i] << "$ &\t$"
                << std::fixed << std::setprecision(2) << Sigma_UPC_val[i] << R"( \pm )" << Sigma_UPC_err_stat[i] << R"((\text{stat.}) \pm )" << Sigma_UPC_err_syst[i] << R"((\text{syst.})$ \\)" << "\n";
    }                  
    outfile.close();
    Printf("Results printed to %s.", str_out_1b.Data());

    // 2a) Print the systematic uncertainties
    TString str_out_2a = "Results/" + str_subfolder + "PhotoCrossSec/Systematics.txt";
    outfile.open(str_out_2a.Data());
    outfile << "[all in percent]\n"
            << "CORRELATED:\n"
            << "lumi\tveto\tEMD\ttracks\tCCUP31\tBR\n"
            << Form("%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n\n",
                ErrSyst_lumi, ErrSyst_veto, ErrSyst_EMD, ErrSyst_tracks, ErrSyst_CCUP31, ErrSyst_BR);
    outfile << "UNCORRELATED:\n"
            << "Bin\tSigExt\tZVtx\tfD\tfC\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << i+1 << std::fixed << std::setprecision(1) << "\t"
                << ErrSyst_SigExtr[i] << "\t"
                << ErrSyst_ZVertex[i] << "\t"
                << ErrSyst_fD[i] << "\t"
                << ErrSyst_fC[i] << "\n";
    }    
    outfile.close();
    Printf("Results printed to %s.", str_out_2a.Data());

    // 2b) Print the systematic uncertainties: TeX table
    TString str_out_2b = "Results/" + str_subfolder + "PhotoCrossSec/Systematics_TeX.txt";
    outfile.open(str_out_2b.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(3)
                << "$(" << ptBoundaries[i] << "," << ptBoundaries[i+1] << ")$ & "
                << std::fixed << std::setprecision(1)
                << ErrSyst_SigExtr[i] << " & "
                << ErrSyst_ZVertex[i] << " & "
                << ErrSyst_fD[i] << " & "
                << ErrSyst_fC[i] << R"( \\)" << "\n";
                            
    }
    outfile.close();
    Printf("Results printed to %s.", str_out_2b.Data()); 

    // 3a) Print the photonuclear cross section
    TString str_out_3a = "Results/" + str_subfolder + "PhotoCrossSec/CrossSec_photo.txt";
    outfile.open(str_out_3a.Data());
    outfile << std::fixed << std::setprecision(4)
            << "Bin \tt_low \tt_upp \tsig \tstat\tsyst\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << i+1 << std::fixed << std::setprecision(4) << "\t" 
                << ptBoundaries[i] * ptBoundaries[i] << "\t" 
                << ptBoundaries[i+1] * ptBoundaries[i+1] << "\t" 
                << std::fixed << std::setprecision(1)
                << Sigma_photo_val[i] << "\t"
                << Sigma_photo_err_stat[i] << "\t"
                << Sigma_photo_err_syst[i] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", str_out_3a.Data()); 

    // 3b) Print the photonuclear cross section: TeX table
    TString str_out_3b = "Results/" + str_subfolder + "PhotoCrossSec/CrossSec_photo_TeX.txt";
    outfile.open(str_out_3b.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(3) << "$(" 
                << ptBoundaries[i] * ptBoundaries[i] << "," 
                << ptBoundaries[i+1] * ptBoundaries[i+1] << ")$ & "
                << t_avg_val[i] << " &\t$"
                << std::fixed << std::setprecision(1)
                << Sigma_photo_val[i] << R"( \pm )" << Sigma_photo_err_stat[i] << R"((\text{stat.}) \pm )" << Sigma_photo_err_syst[i] << R"((\text{syst.})$ \\)" << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", str_out_3b.Data()); 

    return;
}

void PrintErr(TString str)
{
    Printf("ERR: file %s missing. Terminating.", str.Data());
    return;
}