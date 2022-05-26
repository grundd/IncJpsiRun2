// RaphaelleComments.C
// Macros to answer Raphaelle's questions

// root headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

void R02_AcceptanceDimuons();
void R02_FillHistogram(TTree *t, TH1 *h, Bool_t MC, Bool_t etaCut);
void R02_SetHistogram(TH1 *h);
void R15_DeltaPhiVsPt();
void R15_SetHistogram(TH1 *h);

void RaphaelleComments(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    ///*
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "RaphaelleComments/R02_AccDimuons/");
    R02_AcceptanceDimuons();
    //*/

    ///*
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "RaphaelleComments/R15_DeltaPhiVsPt/");
    R15_DeltaPhiVsPt();
    //*/

    return;
}

void R02_AcceptanceDimuons()
{
    TString str_f_data = str_in_DT_fldr + "AnalysisResults.root";
    TString str_f_MC = str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root";

    // open the data tree
    TFile *f_data = TFile::Open(str_f_data.Data(), "read");
    if(f_data) Printf("Input data loaded.");
    TTree *t_data = dynamic_cast<TTree*> (f_data->Get(str_in_DT_tree.Data()));
    if(t_data) Printf("Input tree loaded.");
    ConnectTreeVariables(t_data); 
    // open the MC tree
    TFile *f_MC = TFile::Open(str_f_MC.Data(), "read");
    if(f_MC) Printf("Input data loaded.");
    TTree *t_MC = dynamic_cast<TTree*> (f_MC->Get(str_in_MC_tree_rec.Data()));
    if(t_MC) Printf("Input tree loaded.");
    ConnectTreeVariablesMCRec(t_MC);     

    TH1D *h[4] = { NULL };
    TString titles[4] = {"h_data", "h_data_noEta", "h_MC", "h_MC_noEta"};
    for(Int_t i = 0; i < 4; i++){
        h[i] = new TH1D(titles[i].Data(), titles[i].Data(), 100, -1., 1.);
        h[i]->SetLineWidth(2);
    }

    R02_FillHistogram(t_data, h[0], kFALSE, kTRUE); // data
    R02_FillHistogram(t_data, h[1], kFALSE, kFALSE);// data without eta cut
    R02_FillHistogram(t_MC, h[2], kTRUE, kTRUE); // MC
    R02_FillHistogram(t_MC, h[3], kTRUE, kFALSE);// MC without eta cut

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);

    // with eta cut
    TCanvas *c_eta = new TCanvas("c_eta","c_eta",900,600);
    c_eta->SetTopMargin(0.03);
    c_eta->SetBottomMargin(0.13);
    c_eta->SetRightMargin(0.03);
    c_eta->SetLeftMargin(0.14);
    // without eta cut
    TCanvas *c_noEta = new TCanvas("c_noEta","c_noEta",900,600);
    c_noEta->SetTopMargin(0.03);
    c_noEta->SetBottomMargin(0.13);
    c_noEta->SetRightMargin(0.03);
    c_noEta->SetLeftMargin(0.14);
    // normalize all histograms to 1
    for(Int_t i = 0; i < 4; i++) h[i]->Scale(1./h[i]->Integral());
    // set histograms
    R02_SetHistogram(h[0]); // data
    R02_SetHistogram(h[1]); // data without eta cut
    h[2]->SetLineColor(kRed);
    h[3]->SetLineColor(kRed);

    c_eta->cd();
    h[0]->Draw();       // data
    h[2]->Draw("SAME"); // MC
    // legend
    TLegend *l = new TLegend(0.8,0.75,1.0,0.9);
    l->AddEntry(h[0],"data","L");
    l->AddEntry(h[2],"MC","L");
    l->SetTextSize(0.05);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    c_noEta->cd();
    h[1]->Draw();       // data without eta cut
    h[3]->Draw("SAME"); // MC without eta cut
    l->Draw();

    TString folder = "Results/" + str_subfolder + "RaphaelleComments/R02_AccDimuons/";

    c_eta->Print((folder + "etaCut.pdf").Data());
    c_eta->Print((folder + "etaCut.png").Data());
    c_noEta->Print((folder + "noEtaCut.pdf").Data());
    c_noEta->Print((folder + "noEtaCut.png").Data());

    return;    
}

void R02_FillHistogram(TTree *t, TH1 *h, Bool_t MC, Bool_t etaCut)
{
    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++)
    {
        t->GetEntry(iEntry);

        // Run number in the GoodHadronPID lists published by DPG
        if(!RunNumberInListOfGoodRuns()) continue;

        // if not MC:
        if(!MC){
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
                if(fVertexContrib < cut_fVertexContrib) continue;
                // 4) Distance from the IP lower than 15 cm
                if(fVertexZ > cut_fVertexZ) continue;
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
                if(fVertexContrib < cut_fVertexContrib) continue;
                // 4) Distance from the IP lower than 15 cm
                if(fVertexZ > cut_fVertexZ) continue;
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
            if(!CCUP31) continue;
        }
        // 5a) ADA offline veto (no effect on MC)
        if(!(fADA_dec == 0)) continue;
        // 5b) ADC offline veto (no effect on MC)
        if(!(fADC_dec == 0)) continue;
        // 6a) V0A offline veto (no effect on MC)
        if(!(fV0A_dec == 0)) continue;
        // 6b) V0C offline veto (no effect on MC)
        if(!(fV0C_dec == 0)) continue;
        // 7) SPD cluster matches FOhits
        if(!(fMatchingSPD == kTRUE)) continue;
        // 8) Muon pairs only
        if(!(TMath::Power(fTrk1SigIfMu,2) + TMath::Power(fTrk2SigIfMu,2) < TMath::Power(fTrk1SigIfEl,2) + TMath::Power(fTrk2SigIfEl,2))) continue;
        // 9) Dilepton rapidity |y| < 0.8
        // (skipped)
        // 10) Pseudorapidity of both tracks |eta| < 0.8
        if(etaCut && !(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
        // 11) Tracks have opposite charges
        if(!(fQ1 * fQ2 < 0)) continue;            
        // 12) invariant mass between 2.2 and 4.5 GeV/c^2
        if(!(fM > 2.2 && fM < 4.5)) continue;        
        // 13) transverse momentum cut
        // (skipped)            

        // event passed all the selections =>
        h->Fill(fY);
    }
    return;
}

void R02_SetHistogram(TH1 *h){

    h->SetLineColor(kBlue);
    // horizontal axis
    h->GetXaxis()->SetTitle("#it{y} (-)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.18);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetLabelOffset(0.015);
    h->GetXaxis()->SetDecimals(1);
    // vertical axis
    h->GetYaxis()->SetTitle("Counts per bin (normalized to 1.0)");
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetLabelSize(0.05);
    //h->GetYaxis()->SetDecimals(1);

    return;
}

void R15_DeltaPhiVsPt()
{
    TString str_file = str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root";

    TFile *f = TFile::Open(str_file.Data(), "read");
    if(f) Printf("MC rec file loaded.");
    TTree *t = dynamic_cast<TTree*> (f->Get(str_in_MC_tree_rec.Data()));
    if(t) Printf("MC rec tree loaded.");
    ConnectTreeVariablesMCRec(t);

    // horizontal axis = pT, vertical axis = DeltaPhi
    TH2D *hTriggerInputs[4] = { NULL };
    TProfile *hProfiles[4] = { NULL };
    TString titles[4] = {"0STG_passed", "0STG_rejected", "0OMU_passed", "0OMU_rejected"};
    for(Int_t i = 0; i < 4; i++){ 
        hTriggerInputs[i] = new TH2D(("h_" + titles[i]).Data(), ("h_" + titles[i]).Data(), 160,0.,1.8,68,1.8,3.16);
        hProfiles[i] = new TProfile(("hProf_" + titles[i]).Data(), ("hProf_" + titles[i]).Data(), 160,0.,1.8, 0., 3.16);
        hProfiles[i]->SetLineWidth(2);
    }
    hProfiles[0]->SetLineColor(7);  // 0STG passed
    hProfiles[1]->SetLineColor(4);  // 0STG rejected
    hProfiles[2]->SetLineColor(222);// 0OMU passed
    hProfiles[3]->SetLineColor(2);  // 0OMU rejected

    Double_t DeltaPhi;
    Int_t i0STG_passed(0), i0STG_rejected(0), i0OMU_passed(0), i0OMU_rejected(0);
    ///*
    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++)
    {
        t->GetEntry(iEntry);

        // Run number in the GoodHadronPID lists published by DPG
        if(!RunNumberInListOfGoodRuns()) continue;

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
            // 2) At least two tracks associated with the vertex
            if(fVertexContrib < 2) continue;
            // 3) Distance from the IP lower than 15 cm
            if(fVertexZ > cut_fVertexZ) continue;
        }

        // 4) Central UPC trigger CCUP31: V0 and AD vetoes
        Bool_t CCUP31_vetoes = kFALSE;
        if(
            !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
            !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
            !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
            !fTriggerInputsMC[3]     // !0UBC (no signal in the ADC)
        ) CCUP31_vetoes = kTRUE;
        if(!CCUP31_vetoes) continue;

        DeltaPhi = TMath::Abs(fPhi1 - fPhi2);
        if(DeltaPhi > TMath::Pi()) DeltaPhi = 2*TMath::Pi() - DeltaPhi;   

        if(fTriggerInputsMC[10]){ // 0STG passed
            hTriggerInputs[0]->Fill(fPt,DeltaPhi); 
            hProfiles[0]->Fill(fPt,DeltaPhi);
            i0STG_passed++;
        } else { // 0STG rejected
            hTriggerInputs[1]->Fill(fPt,DeltaPhi); 
            hProfiles[1]->Fill(fPt,DeltaPhi);
            i0STG_rejected++;
        }                     
        if(fTriggerInputsMC[4]){ // 0OMU passed
            hTriggerInputs[2]->Fill(fPt,DeltaPhi); 
            hProfiles[2]->Fill(fPt,DeltaPhi);
            i0OMU_passed++;
        } else { // 0OMU rejected
            hTriggerInputs[3]->Fill(fPt,DeltaPhi); 
            hProfiles[3]->Fill(fPt,DeltaPhi);
            i0OMU_rejected++;
        }                     
    }
    //*/

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);

    TCanvas *c[4] = { NULL };
    for(Int_t i = 0; i < 4; i++){
        c[i] = new TCanvas(Form("c%i", i+1),Form("c%i", i+1),900,600);
        c[i]->SetTopMargin(0.03);
        c[i]->SetBottomMargin(0.13);
        c[i]->SetRightMargin(0.11);
        c[i]->SetLeftMargin(0.95);
        c[i]->SetLogz();
    } 

    TString folder = "Results/" + str_subfolder + "RaphaelleComments/R15_DeltaPhiVsPt/";

    for(Int_t i = 0; i < 4; i++){
        R15_SetHistogram(hTriggerInputs[i]);
        R15_SetHistogram(hProfiles[i]);
        // plot 2d histograms
        c[i]->cd();
        hTriggerInputs[i]->Draw("COLZ");
        c[i]->Print((folder + titles[i] + ".pdf").Data());
        c[i]->Print((folder + titles[i] + ".png").Data());
    }
    // plot profile histograms
    TCanvas *c_profile = new TCanvas("c_profile","c_profile",900,600);
    c_profile->SetTopMargin(0.03);
    c_profile->SetBottomMargin(0.13);
    c_profile->SetRightMargin(0.03);
    c_profile->SetLeftMargin(0.11);
    c_profile->SetLogz();
    c_profile->cd();
    hProfiles[0]->Draw("");
    hProfiles[0]->GetYaxis()->SetRangeUser(2.0, 3.16);
    for(Int_t i = 1; i < 4; i++) hProfiles[i]->Draw("SAME");
    // legend
    TLegend *l = new TLegend(0.15,0.2,0.45,0.4);
    l->AddEntry(hProfiles[0],"0STG passed","L");
    l->AddEntry(hProfiles[1],"0STG rejected","L");
    l->AddEntry(hProfiles[2],"0OMU passed","L");
    l->AddEntry(hProfiles[3],"0OMU rejected","L");
    l->SetTextSize(0.05);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    c_profile->Print((folder + "profiles.pdf").Data());
    c_profile->Print((folder + "profiles.png").Data());

    Printf("0STG passed: %i", i0STG_passed);
    Printf("0STG rejected: %i", i0STG_rejected);
    Printf("0OMU passed: %i", i0OMU_passed);
    Printf("0OMU rejected: %i", i0OMU_rejected);

    return;    
}

void R15_SetHistogram(TH1 *h){
    // horizontal axis
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.18);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetLabelOffset(0.015);
    h->GetXaxis()->SetDecimals(1);
    // vertical axis
    h->GetYaxis()->SetTitle("#Delta#phi (-)");
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetDecimals(1);
}