// _CheckPIDdEdx.C
// David Grund, Mar 06, 2022

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"

TString DatasetsMCNames[8] = {"kCohJpsiToMu",
                              "kCohPsi2sToMuPi",
                              "kIncohJpsiToMu",
                              "kIncohPsi2sToMuPi",
                              "kTwoGammaToMuMedium",
                              "kCohPsi2sToElPi",
                              "kIncohPsi2sToElPi",
                              "kTwoGammaToElMedium"};

void PlotAndFitHistograms(Int_t iDataset, Bool_t calibrated = kFALSE);
// iDataset = 0 => data
//          = 1 => kCohJpsiToMu
//          = 2 => kCohPsi2sToMuPi
//          = 3 => kIncohJpsiToMu
//          = 4 => kIncohPsi2sToMuPi
//          = 5 => kTwoGammaToMuMedium
//          = 6 => kCohPsi2sToElPi
//          = 7 => kIncohPsi2sToElPi
//          = 8 => kTwoGammaToElMedium
void CalculateAverageShift(Bool_t calibrated);
void ShiftPIDSignal(Int_t iDataset, Bool_t pass3);

void _CheckPIDdEdx()
{
    if(kTRUE){
        gSystem->Exec("mkdir -p Results/_CheckPIDdEdx/");

        // data
        PlotAndFitHistograms(0);
        
        // MC datasets without calibration of NSigmas
        for(Int_t i = 1; i <= 8; i++) PlotAndFitHistograms(i, kFALSE);

        // calculate average shifts before NSigmas are calibrated
        CalculateAverageShift(kFALSE);

        // pass3
        // repair NSigmas and kinematics in: kCohJpsiToMu, kCohPsi2sToMuPi, 
        // kIncohJpsiToMu, kIncohPsi2sToMuPi and kTwoGammaToMuMedium
        gSystem->Exec("mkdir -p Trees/AnalysisDataMC_pass3/PIDCalibrated/");
        for(Int_t i = 1; i <= 8; i++) ShiftPIDSignal(i, kTRUE);
    
        // MC datasets after NSigmas are calibrated
        for(Int_t i = 1; i <= 8; i++) PlotAndFitHistograms(i, kTRUE); 

        // calculate average shifts after NSigmas are calibrated
        CalculateAverageShift(kTRUE);          
    } 

    return;
}

void PlotAndFitHistograms(Int_t iDataset, Bool_t calibrated)
{
    // if MC, fit histograms with a gaussian

    // pass1
    TString str_file_pass1 = "";
    TString str_tree_pass1 = "";

    TFile *f_pass1 = NULL;
    TTree *t_pass1 = NULL;

    if(iDataset < 6 && !calibrated){
        if(iDataset == 0) str_file_pass1 = "Trees/AnalysisData_pass1/AnalysisResults.root";
        else if(!calibrated) str_file_pass1 = Form("Trees/AnalysisDataMC_pass1/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());
        else if(calibrated){
            Printf("Input files not prepared. Terminating...");
            return;
        }

        if(iDataset == 0) str_tree_pass1 = "AnalysisOutput/fTreeJPsi";
        else              str_tree_pass1 = "AnalysisOutput/fTreeJPsiMCRec";  

        f_pass1 = TFile::Open(str_file_pass1.Data(), "read");
        if(f_pass1) Printf("Input data loaded.");

        t_pass1 = dynamic_cast<TTree*> (f_pass1->Get(str_tree_pass1.Data()));
        if(t_pass1) Printf("Input tree loaded.");

        isPass3 = kFALSE;

        if(iDataset == 0) ConnectTreeVariables(t_pass1);
        else              ConnectTreeVariablesMCRec(t_pass1);
    }

    // pass3
    TString str_file_pass3 = "";
    TString str_tree_pass3 = "AnalysisOutput/fTreeJpsi";

    if(iDataset == 0) str_file_pass3 = "Trees/AnalysisData_pass3/AnalysisResults.root";
    else if(!calibrated) str_file_pass3 = Form("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());
    else if(calibrated)  str_file_pass3 = Form("Trees/AnalysisDataMC_pass3/PIDCalibrated/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());

    TFile *f_pass3 = TFile::Open(str_file_pass3.Data(), "read");
    if(f_pass3) Printf("Input data loaded.");

    TTree *t_pass3 = dynamic_cast<TTree*> (f_pass3->Get(str_tree_pass3.Data()));
    if(t_pass3) Printf("Input tree loaded.");

    isPass3 = kTRUE;

    if(iDataset == 0) ConnectTreeVariables(t_pass3);
    else              ConnectTreeVariablesMCRec(t_pass3);

    // histograms
    Int_t nBins = 100;
    Double_t range(15.);
    TH1D *hSigmaTPC[4] = { NULL };
    hSigmaTPC[0] = new TH1D("pass1_SigIfMu", "pass1_SigIfMu", nBins, -range, range);
    hSigmaTPC[1] = new TH1D("pass1_SigIfEl", "pass1_SigIfEl", nBins, -range, range);
    hSigmaTPC[2] = new TH1D("pass3_SigIfMu", "pass3_SigIfMu", nBins, -range, range);
    hSigmaTPC[3] = new TH1D("pass3_SigIfEl", "pass3_SigIfEl", nBins, -range, range);

    if(iDataset < 6 && !calibrated){
        Printf("%lli entries found in the tree.", t_pass1->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
            t_pass1->GetEntry(iEntry);

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }

            // Cuts 0-4 offline
            // if MC, trigger has to be replayed:
            if(iDataset != 0){
                Bool_t CCUP31 = kFALSE;
                if(
                    !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
                    !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
                    !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
                    !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
                    fTriggerInputsMC[10] &&  //  0STG (SPD topological)
                    fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
                ) CCUP31 = kTRUE;
                if(!CCUP31) continue;    
            }
            // Cut 5a: ADA offline veto (no effect on MC)
            if(!(fADA_dec == 0)) continue;
            // Cut 5b: ADC offline veto (no effect on MC)
            if(!(fADC_dec == 0)) continue;
            // Cut 6a) V0A offline veto (no effect on MC)
            if(!(fV0A_dec == 0)) continue;
            // Cut 6b) V0C offline veto (no effect on MC)
            if(!(fV0C_dec == 0)) continue;
            // Cut 7: SPD cluster matches FOhits
            if(fMatchingSPD == kFALSE) continue;
            // Cut 8: Dilepton rapidity |y| < 0.8
            if(!(abs(fY) < 0.8)) continue;
            // Cut 9: Pseudorapidity of both tracks |eta| < 0.8
            if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
            // Cut 10: Tracks have opposite charges
            if(!(fQ1 * fQ2 < 0)) continue;

            hSigmaTPC[0]->Fill(fTrk1SigIfMu);
            hSigmaTPC[0]->Fill(fTrk2SigIfMu);
            hSigmaTPC[1]->Fill(fTrk1SigIfEl);
            hSigmaTPC[1]->Fill(fTrk2SigIfEl);  
        }
    }

    Printf("%lli entries found in the tree.", t_pass3->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_pass3->GetEntries(); iEntry++){
        t_pass3->GetEntry(iEntry);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

        // Cuts 0-2 offline
        // if MC, trigger has to be replayed:
        if(iDataset != 0){
            Bool_t CCUP31 = kFALSE;
            if(
                !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
                !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
                !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
                !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
                fTriggerInputsMC[10] &&  //  0STG (SPD topological)
                fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
            ) CCUP31 = kTRUE;
            if(!CCUP31) continue;    
        }
        // Cut 3: more than 2 tracks associated with the primary vertex
        if(fVertexContrib < 2) continue;
        // Cut 4: z-distance from the IP lower than 15 cm
        if(fVertexZ > 15) continue;
        // Cut 5a: ADA offline veto (no effect on MC)
        if(!(fADA_dec == 0)) continue;
        // Cut 5b: ADC offline veto (no effect on MC)
        if(!(fADC_dec == 0)) continue;
        // Cut 6a) V0A offline veto (no effect on MC)
        if(!(fV0A_dec == 0)) continue;
        // Cut 6b) V0C offline veto (no effect on MC)
        if(!(fV0C_dec == 0)) continue;
        // Cut 7: SPD cluster matches FOhits
        if(fMatchingSPD == kFALSE) continue;
        // Cut 8: Dilepton rapidity |y| < 0.8
        if(!(abs(fY) < 0.8)) continue;
        // Cut 9: Pseudorapidity of both tracks |eta| < 0.8
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
        // Cut 10: Tracks have opposite charges
        if(!(fQ1 * fQ2 < 0)) continue;

        hSigmaTPC[2]->Fill(fTrk1SigIfMu);
        hSigmaTPC[2]->Fill(fTrk2SigIfMu);
        hSigmaTPC[3]->Fill(fTrk1SigIfEl);
        hSigmaTPC[3]->Fill(fTrk2SigIfEl);  
    }

    // plot histograms
    TCanvas *c = new TCanvas("c","c",900,700);
    c->Divide(2,2,0,0);
    for(Int_t i = 0; i < 4; i++){
        c->cd(i+1);
        hSigmaTPC[i]->Draw();
    }

    // fit histograms
    TF1 *fGauss[4] = { 0 };
    if(iDataset != 0){
        for(Int_t i = 0; i < 4; i++){
            fGauss[i] = new TF1(Form("fGauss%i", i+1), "gaus", -range, range);
            hSigmaTPC[i]->Fit(fGauss[i]);
        }
    }

    // print the plots and (if MC) print the values of the fitted means
    TString str_f_out = "";
    if(iDataset == 0) str_f_out = "Results/_CheckPIDdEdx/data_LHC18qr";
    else {
        if(!calibrated) str_f_out = "Results/_CheckPIDdEdx/no_calibration/";
        else            str_f_out = "Results/_CheckPIDdEdx/calibrated/";
        gSystem->Exec("mkdir -p " + str_f_out);

        TString str_dataset = Form("MC_%s", DatasetsMCNames[iDataset-1].Data());

        ofstream f_out((str_f_out + "means_" + str_dataset + ".txt").Data());
        f_out << std::fixed << std::setprecision(3);
        for(Int_t i = 0; i < 4; i++){
            f_out << fGauss[i]->GetParameter(1) << endl;
        }
        f_out.close();
        Printf("*** Results printed to %s.***", (str_f_out + "means_" + str_dataset + ".txt").Data());

        str_f_out += str_dataset;
    }

    c->Print((str_f_out + ".png").Data());
    c->Print((str_f_out + ".pdf").Data());

    delete c;

    return;
}

void CalculateAverageShift(Bool_t calibrated)
{
    Double_t fAvgShiftMu_pass1 = 0.;
    Double_t fAvgShiftEl_pass1 = 0.;
    Double_t fAvgShiftMu_pass3 = 0.;
    Double_t fAvgShiftEl_pass3 = 0.;
    Double_t fSigIfMu_pass1, fSigIfEl_pass1, fSigIfMu_pass3, fSigIfEl_pass3;

    // shift for muons
    for(Int_t iDataset = 0; iDataset < 5; iDataset++){
        ifstream ifs;
        TString str_in = "";
        if(!calibrated) str_in = Form("Results/_CheckPIDdEdx/no_calibration/means_MC_%s.txt", DatasetsMCNames[iDataset].Data());
        else            str_in = Form("Results/_CheckPIDdEdx/calibrated/means_MC_%s.txt", DatasetsMCNames[iDataset].Data());
        ifs.open(str_in.Data());
        ifs >> fSigIfMu_pass1 >> fSigIfEl_pass1 >> fSigIfMu_pass3 >> fSigIfEl_pass3;
        Printf("Loaded means for %s:", DatasetsMCNames[iDataset].Data());
        Printf("fSigIfMu_pass1 = %.3f", fSigIfMu_pass1);
        Printf("fSigIfEl_pass1 = %.3f", fSigIfEl_pass1);
        Printf("fSigIfMu_pass3 = %.3f", fSigIfMu_pass3);
        Printf("fSigIfEl_pass3 = %.3f", fSigIfEl_pass3);

        fAvgShiftMu_pass1 += fSigIfMu_pass1;
        fAvgShiftMu_pass3 += fSigIfMu_pass3;
    }
    fAvgShiftMu_pass1 = fAvgShiftMu_pass1 / 5.;
    fAvgShiftMu_pass3 = fAvgShiftMu_pass3 / 5.;

    // shift for electrons
    for(Int_t iDataset = 0; iDataset < 3; iDataset++){
        ifstream ifs;
        TString str_in = "";
        if(!calibrated) str_in = Form("Results/_CheckPIDdEdx/no_calibration/means_MC_%s.txt", DatasetsMCNames[5+iDataset].Data());
        else            str_in = Form("Results/_CheckPIDdEdx/calibrated/means_MC_%s.txt", DatasetsMCNames[5+iDataset].Data());
        ifs.open(str_in.Data());
        ifs >> fSigIfMu_pass1 >> fSigIfEl_pass1 >> fSigIfMu_pass3 >> fSigIfEl_pass3;
        Printf("Loaded means for %s:", DatasetsMCNames[5+iDataset].Data());
        Printf("fSigIfMu_pass1 = %.3f", fSigIfMu_pass1);
        Printf("fSigIfEl_pass1 = %.3f", fSigIfEl_pass1);
        Printf("fSigIfMu_pass3 = %.3f", fSigIfMu_pass3);
        Printf("fSigIfEl_pass3 = %.3f", fSigIfEl_pass3);

        fAvgShiftEl_pass1 += fSigIfEl_pass1;
        fAvgShiftEl_pass3 += fSigIfEl_pass3;
    }
    fAvgShiftEl_pass1 = fAvgShiftEl_pass1 / 3.;
    fAvgShiftEl_pass3 = fAvgShiftEl_pass3 / 3.; 

    // print the results
    TString str_out = "";
    if(!calibrated) str_out = "Results/_CheckPIDdEdx/shifts_no_calibration.txt";
    else            str_out = "Results/_CheckPIDdEdx/shifts_calibrated.txt";
    ofstream outfile(str_out.Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << "fAvgShiftMu_pass1" << "\t" << fAvgShiftMu_pass1 << "\n"
            << "fAvgShiftMu_pass3" << "\t" << fAvgShiftMu_pass3 << "\n"
            << "fAvgShiftEl_pass1" << "\t" << fAvgShiftEl_pass1 << "\n"
            << "fAvgShiftEl_pass3" << "\t" << fAvgShiftEl_pass3 << "\n";

    return;
}

void ShiftPIDSignal(Int_t iDataset, Bool_t pass3)
{
    TString str_f_in = "";
    TString str_t_in = "";
    if(!pass3){
        str_f_in = Form("Trees/AnalysisDataMC_pass1/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());
        str_t_in = "AnalysisOutput/fTreeJPsiMCRec";
    } else {
        str_f_in = Form("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());
        str_t_in = "AnalysisOutput/fTreeJpsi";
    }  
        
    TFile *f_in = TFile::Open(str_f_in.Data(), "read");
    if(f_in) Printf("Input file %s loaded.", str_f_in.Data());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_t_in.Data()));
    if(t_in) Printf("Input tree loaded.");

    TList *l_in = dynamic_cast<TList*> (f_in->Get("AnalysisOutput/fOutputList"));
    if(l_in) Printf("Input list loaded.");

    isPass3 = pass3;

    ConnectTreeVariablesMCRec(t_in);

    // define new histograms to which we copy the old ones
    TList *l_out = new TList();
    TH1D *hCounterCuts = (TH1D*)l_in->FindObject("hCounterCuts")->Clone("hCounterCuts");
    l_out->Add(hCounterCuts);

    // histograms where the differences will be stored
    Double_t range = 0.02;
    TH1D *hPt = new TH1D("hPt", "hPt", 100, -range, range);
    TH1D *hPhi = new TH1D("hPhi", "hPhi", 100, -range, range);
    TH1D *hY = new TH1D("hY", "hY", 100, -range, range);
    TH1D *hM = new TH1D("hM", "hM", 100, -range, range);

    // create paths to new files and trees
    TString str_f_out = "";
    TString str_t_out = "";
    if(!pass3){
        str_f_out = Form("Trees/AnalysisDataMC_pass1/PIDCalibrated/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());
        str_t_out = "AnalysisOutput/fTreeJPsiMCRec";
    } else {
        str_f_out = Form("Trees/AnalysisDataMC_pass3/PIDCalibrated/AnalysisResults_MC_%s.root", DatasetsMCNames[iDataset-1].Data());
        str_t_out = "AnalysisOutput/fTreeJpsi";
    }  

    // create new tree
    gROOT->cd();
    TTree *fTreeJpsi = new TTree(str_t_out.Data(),str_t_out.Data());
    // Basic things:
    fTreeJpsi->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    fTreeJpsi->Branch("fTriggerName", &fTriggerName);
    // PID, sigmas:
    if(pass3){
        fTreeJpsi->Branch("fTrk1dEdx", &fTrk1dEdx, "fTrk1dEdx/D");
        fTreeJpsi->Branch("fTrk2dEdx", &fTrk2dEdx, "fTrk2dEdx/D");
    }
    fTreeJpsi->Branch("fTrk1SigIfMu", &fTrk1SigIfMu, "fTrk1SigIfMu/D");
    fTreeJpsi->Branch("fTrk1SigIfEl", &fTrk1SigIfEl, "fTrk1SigIfEl/D");
    fTreeJpsi->Branch("fTrk2SigIfMu", &fTrk2SigIfMu, "fTrk2SigIfMu/D");
    fTreeJpsi->Branch("fTrk2SigIfEl", &fTrk2SigIfEl, "fTrk2SigIfEl/D");
    // Kinematics:
    fTreeJpsi->Branch("fPt", &fPt, "fPt/D");
    fTreeJpsi->Branch("fPhi", &fPhi, "fPhi/D");
    fTreeJpsi->Branch("fY", &fY, "fY/D");
    fTreeJpsi->Branch("fM", &fM, "fM/D");
    // Two tracks:
    fTreeJpsi->Branch("fPt1", &fPt1, "fPt1/D");
    fTreeJpsi->Branch("fPt2", &fPt2, "fPt2/D");
    fTreeJpsi->Branch("fEta1", &fEta1, "fEta1/D");
    fTreeJpsi->Branch("fEta2", &fEta2, "fEta2/D");
    fTreeJpsi->Branch("fPhi1", &fPhi1, "fPhi1/D");
    fTreeJpsi->Branch("fPhi2", &fPhi2, "fPhi2/D");
    fTreeJpsi->Branch("fQ1", &fQ1, "fQ1/D");
    fTreeJpsi->Branch("fQ2", &fQ2, "fQ2/D");
    // Vertex info:
    if(pass3){
        fTreeJpsi->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");
        fTreeJpsi->Branch("fVertexContrib", &fVertexContrib, "fVertexContrib/I"); 
    }
    // Info from the detectors:
    // ZDC:
    fTreeJpsi->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    fTreeJpsi->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    fTreeJpsi->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    fTreeJpsi->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    // V0, AD:
    fTreeJpsi->Branch("fV0A_dec", &fV0A_dec, "fV0A_dec/I");
    fTreeJpsi->Branch("fV0C_dec", &fV0C_dec, "fV0C_dec/I");
    fTreeJpsi->Branch("fADA_dec", &fADA_dec, "fADA_dec/I");
    fTreeJpsi->Branch("fADC_dec", &fADC_dec, "fADC_dec/I");
    if(!pass3){
        fTreeJpsi->Branch("fV0A_time", &fV0A_time, "fV0A_time/D");
        fTreeJpsi->Branch("fV0C_time", &fV0C_time, "fV0C_time/D");
        fTreeJpsi->Branch("fADA_time", &fADA_time, "fADA_time/D");
        fTreeJpsi->Branch("fADC_time", &fADC_time, "fADC_time/D");        
    }
    // Matching SPD clusters with FOhits:
    fTreeJpsi->Branch("fMatchingSPD", &fMatchingSPD, "fMatchingSPD/O");
    // Replayed trigger inputs:
    fTreeJpsi->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[11]/O"); 
    // Kinematics, MC gen:
    fTreeJpsi->Branch("fPtGen", &fPtGen, "fPtGen/D");
    fTreeJpsi->Branch("fMGen", &fMGen, "fMGen/D");
    fTreeJpsi->Branch("fYGen", &fYGen, "fYGen/D");
    fTreeJpsi->Branch("fPhiGen", &fPhiGen, "fPhiGen/D");

    Double_t fShiftMu = 0;
    Double_t fShiftEl = 0;
    // these values are taken from 'means.txt' created in the function CalculateAverageShift()
    if(!pass3){
        fShiftMu = 0.35; // 0.354
        fShiftEl = 0.00; // not calculated
    } else {
        fShiftMu = 1.41; // 1.410
        fShiftEl = 2.49; // 2.490
    }

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);

        fTrk1SigIfMu = fTrk1SigIfMu - fShiftMu;
        fTrk1SigIfEl = fTrk1SigIfEl - fShiftEl;
        fTrk2SigIfMu = fTrk2SigIfMu - fShiftMu;
        fTrk2SigIfEl = fTrk2SigIfEl - fShiftEl;

        // store the old values of J/psi kinematic variables
        Double_t fPt_old = fPt;
        Double_t fPhi_old = fPhi;
        Double_t fY_old = fY;
        Double_t fM_old = fM;

        // recalculate J/psi kinematics after assigning a proper mass to tracks
        Double_t isMuonPair = fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu;
        Double_t isElectronPair = fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl;
        Double_t massTracks = -1;
        if(isMuonPair < isElectronPair) massTracks = 0.105658; // GeV/c^2
        else                            massTracks = 0.000511; // GeV/c^2

        TLorentzVector vTrk1, vTrk2;
        vTrk1.SetPtEtaPhiM(fPt1, fEta1, fPhi1, massTracks);
        vTrk2.SetPtEtaPhiM(fPt2, fEta2, fPhi2, massTracks);
        TLorentzVector vTrkTrk = vTrk1 + vTrk2;

        // set tree variables
        fPt = vTrkTrk.Pt(); 
        fPhi = vTrkTrk.Phi();
        fY = vTrkTrk.Rapidity(); 
        fM = vTrkTrk.M();

        // calculate the differences
        Double_t Delta_pt = fPt_old - fPt;
        Double_t Delta_phi = fPhi_old - fPhi;
        Double_t Delta_y = fY_old - fY;
        Double_t Delta_m = fM_old - fM;
        // fill the histograms
        hPt->Fill(Delta_pt);
        hPhi->Fill(Delta_phi);
        hY->Fill(Delta_y);
        hM->Fill(Delta_m);

        fTreeJpsi->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

    }

    // create new file
    TFile *f_out = new TFile(str_f_out.Data(),"RECREATE");
    // open the file
    f_out->cd();
    // create a subdirectory
    gDirectory->mkdir("AnalysisOutput");
    // open it
    f_out->cd("AnalysisOutput");
    // write the list and tree to this directory
    l_out->Write("fOutputList", TObject::kSingleKey);
    fTreeJpsi->Write("fTreeJpsi",TObject::kSingleKey);
    // list the contents of the file
    f_out->ls();
    // close the file
    f_out->Close();

    Printf("New file created and saved.");

    // plot histograms
    TCanvas *c = new TCanvas("c","c",900,700);
    c->Divide(2,2,0,0);
    c->cd(1); hPt->Draw();
    c->cd(2); hPhi->Draw();
    c->cd(3); hY->Draw();
    c->cd(4); hM->Draw();

    TString str_out = "";
    if(!pass3){
        gSystem->Exec("mkdir -p Results/_CheckPIDdEdx/pass1_" + DatasetsMCNames[iDataset-1] + "/");
        str_out = Form("Results/_CheckPIDdEdx/pass1_%s/differences_pass1", DatasetsMCNames[iDataset-1].Data());
    } else {
        gSystem->Exec("mkdir -p Results/_CheckPIDdEdx/pass3_" + DatasetsMCNames[iDataset-1] + "/");
        str_out = Form("Results/_CheckPIDdEdx/pass3_%s/differences_pass3", DatasetsMCNames[iDataset-1].Data());
    }
    c->Print((str_out + ".png").Data());
    c->Print((str_out + ".pdf").Data());

    delete c;

    return;
}