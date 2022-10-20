// VertexZDistribution.C
// David Grund, Jun 06, 2022 

// root headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

void PrepareTree();
Bool_t EventPassedLocal(Bool_t isMC);
void ConnectTreeVariablesLocal(TTree *t);

void VertexZDistribution(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VertexZDistribution/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VertexZDistribution/");

    //PrepareTree();

    TString name = "Trees/" + str_subfolder + "VertexZDistribution/tree.root";

    // open the file with trees
    TFile *f_in = TFile::Open(name.Data(), "read");
    if(f_in) Printf("Input data loaded.");
    TTree *t_mc = dynamic_cast<TTree*> (f_in->Get("t_mc"));
    if(t_mc) Printf("Input MC tree loaded.");
    ConnectTreeVariablesLocal(t_mc); 
    TTree *t_dt = dynamic_cast<TTree*> (f_in->Get("t_dt"));
    if(t_dt) Printf("Input data tree loaded.");
    ConnectTreeVariablesLocal(t_dt);

    // define histograms
    Int_t nBins = 80; // so that each bin between -20 cm and +20 cm is 0.5 cm wide
    TH1D *hVertexZ_dt = new TH1D("hVertexZ_dt","hVertexZ_dt",80,-20.,20.);
    TH1D *hVertexZ_mc = new TH1D("hVertexZ_mc","hVertexZ_mc",80,-20.,20.);
    TH2D *hVertZvsPt = new TH2D("hVertZvsPt","hVertZvsPt",80,0.2,1.0,80,-20.,20.);
    Double_t Zboundaries[4] = {-15.,-10.,10.,15.};
    TH2D *hVertZvsPt_bins = new TH2D("hVertZvsPt_bins","hVertZvsPt_bins",nPtBins,ptBoundaries,3,Zboundaries);

    // go over data events
    Printf("%lli entries found in the data tree.", t_dt->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_dt->GetEntries(); iEntry++)
    {
        t_dt->GetEntry(iEntry);
        if(fM > 3.0 && fM < 3.2 && fPt > 0.2 && fPt < 1.0)
        {
            hVertexZ_dt->Fill(fVertexZ);
            hVertZvsPt->Fill(fPt,fVertexZ);
            hVertZvsPt_bins->Fill(fPt,fVertexZ);
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // go over MC events
    Printf("%lli entries found in the MC tree.", t_mc->GetEntries());
    nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_mc->GetEntries(); iEntry++)
    {
        t_mc->GetEntry(iEntry);
        if(fM > 3.0 && fM < 3.2 && fPt > 0.2 && fPt < 1.0) hVertexZ_mc->Fill(fVertexZ);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed); 
        }
    }

    // if there are 80 bins from -20 cm to +20 cm, each of them is 0.5 cm wide
    // indices of specific bins:
    // (-15.0,-14.5): 11
    // (-10.5,-10.0): 20
    // (+10.0,+10.5): 61
    // (+14.5,+15.0): 70
    Int_t iLeft_low = 11;
    Int_t iLeft_upp = 20;
    Int_t iRght_low = 61;
    Int_t iRght_upp = 70;

    Double_t nEv_dt_left = hVertexZ_dt->Integral(iLeft_low,iLeft_upp);
    Double_t nEv_dt_rght = hVertexZ_dt->Integral(iRght_low,iRght_upp);
    Printf("Total no. of events: %.2f", hVertexZ_dt->Integral());
    Printf("No. of data events between:");
    Printf("(-15.cm,-10.cm): %.2f", nEv_dt_left);
    Printf("(10.0cm,15.0cm): %.2f", nEv_dt_rght);

    // normalize both histograms to 1.0
    hVertexZ_dt->Scale(1.0/hVertexZ_dt->Integral());
    hVertexZ_mc->Scale(1.0/hVertexZ_mc->Integral());

    Double_t fEv_dt_left = hVertexZ_dt->Integral(iLeft_low,iLeft_upp);
    Double_t fEv_dt_rght = hVertexZ_dt->Integral(iRght_low,iRght_upp);
    Double_t fEv_mc_left = hVertexZ_mc->Integral(iLeft_low,iLeft_upp);
    Double_t fEv_mc_rght = hVertexZ_mc->Integral(iRght_low,iRght_upp);
    Printf("Integral of normalized vertex distributions within:");
    Printf("(-15.cm,-10.cm): MC: %.4f, data: %.4f", fEv_mc_left, fEv_dt_left);
    Printf("(10.0cm,15.0cm): MC: %.4f, data: %.4f", fEv_mc_rght, fEv_dt_rght);

    // Percentage of events lying in the tails:
    Double_t fPerc_dt = (fEv_dt_left + fEv_dt_rght) / hVertexZ_dt->Integral() * 100;
    Double_t fPerc_mc = (fEv_mc_left + fEv_mc_rght) / hVertexZ_mc->Integral() * 100;
    Printf("Percentage of events lying in the tails |Z| in (10,15)cm:");
    Printf("MC:   %.1f", fPerc_dt);
    Printf("data: %.1f", fPerc_mc);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.0f");

    TString str_out = "Results/" + str_subfolder + "VertexZDistribution/";

    TCanvas *c1 = new TCanvas("c1","c1",900,600);
    hVertexZ_dt->SetLineColor(kBlue);
    hVertexZ_mc->SetLineColor(kRed);
    hVertexZ_dt->Draw("HIST");
    hVertexZ_mc->Draw("HIST SAME");
    c1->Print((str_out + "data_vs_mc.pdf").Data());

    TCanvas *c2 = new TCanvas("c2","c2",900,600);
    hVertZvsPt->Draw("COLZ");
    c2->Print((str_out + "2d_ZvsPt.pdf").Data());

    TCanvas *c3 = new TCanvas("c3","c3",900,600);
    hVertZvsPt_bins->SetMarkerSize(2);
    hVertZvsPt_bins->Draw("COLZ TEXT");
    c3->Print((str_out + "2d_ZvsPt_bins.pdf").Data());

    return;
}

void PrepareTree()
{
    TString str_f_dt = str_in_DT_fldr + "AnalysisResults.root";
    TString str_f_mc = str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root";

    // open the data tree
    TFile *f_dt = TFile::Open(str_f_dt.Data(), "read");
    if(f_dt) Printf("Input data loaded.");
    TTree *t_DT = dynamic_cast<TTree*> (f_dt->Get(str_in_DT_tree.Data()));
    if(t_DT) Printf("Input tree loaded.");
    ConnectTreeVariables(t_DT); 
    // open the MC tree
    TFile *f_mc = TFile::Open(str_f_mc.Data(), "read");
    if(f_mc) Printf("Input data loaded.");
    TTree *t_MC = dynamic_cast<TTree*> (f_mc->Get(str_in_MC_tree_rec.Data()));
    if(t_MC) Printf("Input tree loaded.");
    ConnectTreeVariablesMCRec(t_MC);   

    // Create new tree
    TString name = "Trees/" + str_subfolder + "VertexZDistribution/tree.root";
    TFile *f_out = new TFile(name.Data(),"RECREATE");

    TTree *t_dt = new TTree("t_dt", "t_dt");
    t_dt->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");
    t_dt->Branch("fPt", &fPt, "fPt/D");
    t_dt->Branch("fM", &fM, "fM/D");

    TTree *t_mc = new TTree("t_mc", "t_mc");
    t_mc->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");
    t_mc->Branch("fPt", &fPt, "fPt/D");
    t_mc->Branch("fM", &fM, "fM/D");

    // go over data events
    Printf("%lli entries found in the data tree.", t_DT->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_DT->GetEntries(); iEntry++)
    {
        t_DT->GetEntry(iEntry);
        if(EventPassedLocal(kFALSE)) t_dt->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // go over MC events
    Printf("%lli entries found in the MC tree.", t_MC->GetEntries());
    nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_MC->GetEntries(); iEntry++)
    {
        t_MC->GetEntry(iEntry);
        if(EventPassedLocal(kTRUE)) t_mc->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    f_out->Write("",TObject::kWriteDelete);

    return;
}

Bool_t EventPassedLocal(Bool_t isMC)
{
        // Run number in the GoodHadronPID lists published by DPG
        if(!RunNumberInListOfGoodRuns()) return kFALSE;

        // if not MC:
        if(!isMC){
            // pass1:
            // 0) fEvent non-empty
            // 1) At least two tracks associated with the vertex
            // 2) Distance from the IP lower than 15 cm
            // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
            // 4) Central UPC trigger CCUP31
            // pass3:
            // 0) fEvent non-empty
            // 1) At least two tracks associated with the vertex
            // 2) Central UPC trigger CCUP31
            if(isPass3){
                // 3) At least two tracks associated with the vertex
                if(fVertexContrib < cut_fVertexContrib) return kFALSE;
                // 4) Distance from the IP lower than 15 cm
                // no cut
            }
        // if MC:
        } else {
            // pass1:
            // All selections applied on the GRID:
            // 0) fEvent non-empty
            // 1) At least two tracks associated with the vertex
            // 2) Distance from the IP lower than 15 cm
            // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
            // pass3:
            // 0) fEvent non-empty
            // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
            if(isPass3){
                // 3) At least two tracks associated with the vertex
                if(fVertexContrib < cut_fVertexContrib) return kFALSE;
                // 4) Distance from the IP lower than 15 cm
                // no cut
            }
            // 4) Central UPC trigger CCUP31:
            // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
            // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
            Bool_t CCUP31 = kFALSE;
            if(
                !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
                !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
                !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
                !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
                fTriggerInputsMC[10] &&  //  0STG (SPD topological)
                fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
            ) CCUP31 = kTRUE;
            if(!CCUP31) return kFALSE;
        }
    
    // 5a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;

    // 5b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;

    // 6a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;

    // 6b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;

    // 7) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;

    // 8) Muon pairs only
    if(!(TMath::Power(fTrk1SigIfMu,2) + TMath::Power(fTrk2SigIfMu,2) < TMath::Power(fTrk1SigIfEl,2) + TMath::Power(fTrk2SigIfEl,2))) return kFALSE;

    // 9) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;

    // 10) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;

    // 11) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;

    // 12) Invariant mass cut
    // no cut

    // 13) Transverse momentum cut
    // no cut

    // Event passed all the selections =>
    return kTRUE;
}

void ConnectTreeVariablesLocal(TTree *t)
{
    t->SetBranchAddress("fVertexZ", &fVertexZ);
    t->SetBranchAddress("fPt", &fPt);
    t->SetBranchAddress("fM", &fM);

    return;
}