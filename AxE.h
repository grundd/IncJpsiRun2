// AxE.h
// David Grund, Mar 05, 2022

TString path = "Results/AccAndEffMC/";

Double_t NRec;
Double_t NGen;

Double_t AxE;
Double_t AxE_err;

TString DatasetsMC[6] = {"JCoh","JInc","PCohCh","PIncCh", "PCohNe", "PIncNe"};

void CalculateAxE(Int_t iMC, Int_t iMassCut, Int_t iPtCut, Int_t iPtBin = -1);
// iMC == 0 => kCohJpsiToMu
//     == 1 => kIncohJpsiToMu
//     == 2 => kCohPsi2sToMuPi charged
//     == 3 => kIncohPsi2sToMuPi charged
//     == 4 => kCohPsi2sToMuPi neutral
//     == 5 => kIncohPsi2sToMuPi neutral
void CalculateAxE_AOD(Int_t iMC, Int_t iMassCut, Int_t iPtCut);
Double_t CalculateErrorBayes(Double_t k, Double_t n);

void AccAndEffMC(){

    // Calculate the total AxE from ESD data
    // IncJ, 2.2 < m < 4.5 GeV, pt > 0.2 GeV
    Bool_t bTotalESD = kFALSE;
    if(bTotalESD) CalculateAxE(1, 0, 0);

    // Calculate the total AxE from (old) AOD data
    // IncJ, pt > 0.2 GeV
    Bool_t bTotalAOD = kFALSE;
    if(bTotalAOD){
        CalculateAxE_AOD(1, 0, 0); // 2.2 < m < 4.5 GeV
        CalculateAxE_AOD(1, 1, 0); // 3.0 < m < 3.2 GeV
    }

    // Calculate AxE for the total FD correction
    Bool_t bCorrFDTotal = kFALSE;
    if(bCorrFDTotal){
        TString str = Form("%sAxE_FeedDown_Total.txt", path.Data());
        ofstream outfile(str.Data());
        outfile << std::fixed << std::setprecision(4);
        outfile << Form("\tAxE[%%] \tErr \n");
        // For all datasets (except JCoh):
        for(Int_t iMC = 1; iMC < 6; iMC++){
            CalculateAxE(iMC, 0, 0);
            outfile << DatasetsMC[iMC] << "\t" << AxE << "\t" << AxE_err << "\n";
        }
        outfile.close();
        Printf("*** Results printed to %s.***", str.Data());
    }

    // Calculate AxE for the total FD correction with AODs
    Bool_t bCorrFDTotalAOD = kFALSE;
    if(bCorrFDTotalAOD){
        TString str = Form("%sAOD/AxE_AOD_FeedDown_Total.txt", path.Data());
        ofstream outfile(str.Data());
        outfile << std::fixed << std::setprecision(4);
        outfile << Form("\tAxE[%%] \tErr \n");
        // For all datasets (except JCoh):
        for(Int_t iMC = 1; iMC < 6; iMC++){
            CalculateAxE_AOD(iMC, 1, 0);
            outfile << DatasetsMC[iMC] << "\t" << AxE << "\t" << AxE_err << "\n";
        }
        outfile.close();
        Printf("*** Results printed to %s.***", str.Data());
    }

    // Calculate AxE in pt bins for FD correction
    Bool_t bCorrFD = kFALSE;
    if(bCorrFD){
        SetPtBinning();
        TString str1 = Form("%sAxE_FeedDown_%ibins.txt", path.Data(), nPtBins);
        TString str2 = Form("%sAxE_FeedDown_%ibins_NGen_SL.txt", path.Data(), nPtBins);
        ofstream outfile1(str1.Data());
        ofstream outfile2(str2.Data());
        outfile1 << std::fixed << std::setprecision(4);
        outfile2 << std::fixed << std::setprecision(0);
        outfile1 << Form("AxE[%%] \tJInc \tErr \tPCohCh \tErr \tPIncCh \tErr \tPCohNe \tErr \tPIncNe \tErr \n");
        outfile2 << "NGen \tJInc \tPCohCh \tPIncCh \tPCohNe \tPIncNe \nTotal";
        // Go over datasets to get NGen_tot
        for(Int_t iMC = 1; iMC < 6; iMC++){
            CalculateAxE(iMC, 0, 0);
            outfile2 << "\t" << NGen;
        }
        outfile2 << "\n";
        // For all pt bins:
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++){
            outfile1 << iBin;
            outfile2 << iBin;
            // For all datasets (except JCoh):
            for(Int_t iMC = 1; iMC < 6; iMC++){
                CalculateAxE(iMC, 0, 4, iBin);
                outfile1 << "\t" << AxE << "\t" << AxE_err;
                outfile2 << "\t" << NGen;
            }
            outfile1 << "\n";
            outfile2 << "\n";
        }
        outfile1.close();
        Printf("*** Results printed to %s.***", str1.Data());
        Printf("*** Results printed to %s.***", str2.Data());
    }

    // Calculate AxE for PtFit (FC correction)
    // 3.0 < m < 3.2 GeV, pt < 2 GeV
    Bool_t bCorrFCTotal = kFALSE;
    if(bCorrFCTotal){
        TString str = Form("%sAxE_PtFit.txt", path.Data());
        ofstream outfile(str.Data());
        outfile << std::fixed << std::setprecision(4);
        outfile << Form("Dataset\tAxE [%%]\tAxE_err \n");
        // For all datasets:
        for(Int_t iMC = 0; iMC < 6; iMC++){ 
            CalculateAxE(iMC, 1, 2);
            outfile << DatasetsMC[iMC] << "\t" << AxE << "\t" << AxE_err << "\n";
        }
        outfile.close();
        Printf("*** Results printed to %s.***", str.Data());
    }

    // Calculate AxE for PtFit (FC correction) with AODs
    // 3.0 < m < 3.2 GeV, pt < 2 GeV
    Bool_t bCorrFCTotalAOD = kFALSE;
    if(bCorrFCTotalAOD){
        TString str = Form("%sAOD/AxE_AOD_PtFit.txt", path.Data());
        ofstream outfile(str.Data());
        outfile << std::fixed << std::setprecision(4);
        outfile << Form("Dataset\tAxE [%%]\tAxE_err \n");
        // For all datasets:
        for(Int_t iMC = 0; iMC < 6; iMC++){ 
            CalculateAxE_AOD(iMC, 1, 2);
            outfile << DatasetsMC[iMC] << "\t" << AxE << "\t" << AxE_err << "\n";
        }
        outfile.close();
        Printf("*** Results printed to %s.***", str.Data());
    }

    // Direct calculation of FC correction from AxEs: a total value
    Bool_t bCorrFCDirect = kFALSE;
    if(bCorrFCDirect){
        TString str = Form("%sAxE_CorrFC_Total.txt", path.Data());
        ofstream outfile(str.Data());
        outfile << std::fixed << std::setprecision(4);
        outfile << Form("AxE[%%] \tJCoh \tErr \tJInc \tErr \nTotal");
        // Only kCohJpsiToMu and kIncohJpsiToMu
        for(Int_t iMC = 0; iMC < 2; iMC++){ 
            CalculateAxE(iMC, 0, 0);
            outfile << "\t" << AxE << "\t" << AxE_err;
        }
        outfile.close();
        Printf("*** Results printed to %s.***", str.Data());
    }

    // Direct calculation of FC correction from AxEs: bins
    Bool_t bCorrFCDirectBins = kTRUE;
    if(bCorrFCDirectBins){
        SetPtBinning();
        TString str1 = Form("%sAxE_CorrFC_%ibins.txt", path.Data(), nPtBins);
        TString str2 = Form("%sAxE_CorrFC_%ibins_NGen_SL.txt", path.Data(), nPtBins);
        ofstream outfile1(str1.Data());
        ofstream outfile2(str2.Data());
        outfile1 << std::fixed << std::setprecision(4);
        outfile2 << std::fixed << std::setprecision(0);
        outfile1 << Form("AxE[%%] \tJCoh \tErr \tJInc \tErr\n");
        outfile2 << Form("NGen \tJCoh \tJInc \nTotal");
        // Go over kCohJpsiToMu and kIncohJpsiToMu to get NGen_tot
        for(Int_t iMC = 0; iMC < 2; iMC++){
            CalculateAxE(iMC, 0, 0);
            outfile2 << "\t" << NGen;
        }
        outfile2 << "\n";
        // For all pt bins:
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++){
            outfile1 << iBin;
            outfile2 << iBin;
            // For kCohJpsiToMu and kIncohJpsiToMu
            for(Int_t iMC = 0; iMC < 2; iMC++){
                CalculateAxE(iMC, 0, 4, iBin);
                outfile1 << "\t" << AxE << "\t" << AxE_err;
                outfile2 << "\t" << NGen;
            }
            outfile1 << "\n";
            outfile2 << "\n";
        }
        outfile1.close();
        Printf("*** Results printed to %s.***", str1.Data());
        Printf("*** Results printed to %s.***", str2.Data());
    }

    return;
}

void CalculateAxE(Int_t iMC, Int_t iMassCut, Int_t iPtCut, Int_t iPtBin){

    TFile *fRec = NULL;
    switch(iMC){
        case 0:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
            break;
        case 1:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
            break;        
        case 2:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            break;
        case 3:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
            break;
        case 4:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi_neutral.root", "read");
            break;
        case 5:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi_neutral.root", "read");
            break;
    }
    if(fRec) Printf("MC rec file for %s loaded.", DatasetsMC[iMC].Data());

    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    NRec = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(iMassCut, iPtCut, iPtBin)) NRec++;
    }

    TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    NGen = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(iPtCut, iPtBin)) NGen++;
    }

    Printf("N_rec = %.0f", NRec);
    Printf("N_gen = %.0f", NGen);

    if(NGen != 0){
        AxE = NRec / NGen * 100; // in percent already
        AxE_err = CalculateErrorBayes(NRec, NGen) * 100;
    } else {
        AxE = 0;
        AxE_err = 0;
    }

    Printf("AxE = (%.4f pm %.4f)%%", AxE, AxE_err);

    TString str;
    if(iPtBin <= 0) str = Form("%sAxE_%s_MassCut%i_PtCut%i.txt", path.Data(), DatasetsMC[iMC].Data(), iMassCut, iPtCut);
    if(iPtBin > 0) str = Form("%s/AxE_%ibins/AxE_bin%i_%s.txt", path.Data(), nPtBins, iPtBin, DatasetsMC[iMC].Data());
    ofstream outfile(str.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("AxE [%%]\tAxE_err \n");
    outfile << AxE << "\t" << AxE_err;
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());

    return;
}

void CalculateAxE_AOD(Int_t iMC, Int_t iMassCut, Int_t iPtCut){

    TFile *fRec = NULL;
    switch(iMC){
        case 0: 
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kCohJpsiToMu.root", "read");
            break;
        case 1:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohJpsiToMu_migr.root", "read");
            break;
        case 2:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kCohPsi2sToMuPi.root", "read");
            break;
        case 3:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohPsi2sToMuPi.root", "read");
            break;
        case 4:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kCohPsi2sToMuPi_NeutPi.root", "read");
            break;
        case 5:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohPsi2sToMuPi_NeutPi.root", "read");
            break;
    }
    if(fRec) Printf("MC rec file loaded.");

    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("analysisTree"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec_AOD(tRec);

    NRec = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec_AOD(iMassCut, iPtCut)) NRec++;
    }

    TFile *fGen = NULL;
    switch(iMC){
        case 0:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kCohJpsiToMu.root", "read");
            break;
        case 1:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kIncohJpsiToMu_migr.root", "read");
            break;
        case 2:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kCohPsi2sToMuPi.root", "read");
            break;
        case 3:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kIncohPsi2sToMuPi.root", "read");
            break;
        case 4:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kCohPsi2sToMuPi_NeutPi.root", "read");
            break;
        case 5:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kIncohPsi2sToMuPi_NeutPi.root", "read");
            break;
    }
    if(fGen) Printf("MC gen file loaded.");

    TTree *tGen = dynamic_cast<TTree*> (fGen->Get("MCgenTree"));
    if(tGen) Printf("MC gen tree loaded.");
    
    if(iMC == 1) ConnectTreeVariablesMCGen_AOD(tGen);
    else ConnectTreeVariablesMCGen_AOD_old(tGen);

    NGen = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(iPtCut)) NGen++;
    }

    Printf("N_rec = %.0f", NRec);
    Printf("N_gen = %.0f", NGen);

    if(NGen != 0){
        AxE = NRec / NGen * 100; // in percent already
        AxE_err = CalculateErrorBayes(NRec, NGen) * 100;
    } else {
        AxE = 0;
        AxE_err = 0;
    }
    

    Printf("AxE = (%.4f pm %.4f)%%", AxE, AxE_err);

    TString str = Form("%sAOD/AxE_AOD_%s_MassCut%i_PtCut%i", path.Data(), DatasetsMC[iMC].Data(), iMassCut, iPtCut);
    ofstream outfile((str + ".txt").Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("AxE [%%]\tAxE_err \n");
    outfile << AxE << "\t" << AxE_err;
    outfile.close();
    Printf("*** Results printed to %s.***", (str + ".txt").Data());

    return;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}