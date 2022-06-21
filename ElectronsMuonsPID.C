// ElectronsMuonsPID.c
// David Grund, Dec 2, 2021

// root headers
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

void Make2dHistograms(Bool_t isMC);
Bool_t EventPassedLocal(Bool_t isMC);

void ElectronsMuonsPID(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "ElectronsMuonsPID/");

    // data
    Make2dHistograms(kFALSE);
    // MC, kIncohJpsiToMu
    Make2dHistograms(kTRUE);

    return;
}

void Make2dHistograms(Bool_t isMC)
{
    Int_t nBins1 = 100;
    Double_t sigma_low = 0.0;
    Double_t sigma_upp = 20.0;
    Int_t nBins2 = 200;
    Double_t dEdx_low = 0.0;
    Double_t dEdx_upp = 150.0;
    // horizontal axis = sigma if electron, vertical axis = sigma if muon
    TH2D *hSigmasTPC = new TH2D("hSigmasTPC","hSigmasTPC",nBins1,sigma_low,sigma_upp,nBins1,sigma_low,sigma_upp);
    // horizontal axis = dEdx for negative lepton, vertical axis = dEdx for positive lepton
    TH2D *hdEdxElec = new TH2D("hdEdxElec","hdEdxElec",nBins2,dEdx_low,dEdx_upp,nBins2,dEdx_low,dEdx_upp);
    TH2D *hdEdxMuon = new TH2D("hdEdxMuon","hdEdxMuon",nBins2,dEdx_low,dEdx_upp,nBins2,dEdx_low,dEdx_upp);

    // Load data
    TString str_f_in, str_t_in;
    if(isMC){
        str_f_in = str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root";
        str_t_in = str_in_MC_tree_rec;
    } else {
        str_f_in = str_in_DT_fldr + "AnalysisResults.root";
        str_t_in = str_in_DT_tree;
    }    
    TFile *f = TFile::Open(str_f_in.Data(), "read");
    if(f) Printf("Input data loaded.");

    TTree *t = dynamic_cast<TTree*> (f->Get(str_t_in.Data()));
    if(t) Printf("Input tree loaded.");

    if(isMC) ConnectTreeVariablesMCRec(t);
    else     ConnectTreeVariables(t);

    Printf("%lli entries found in the tree.", t->GetEntries());
    Int_t nEntriesAnalysed = 0;

    ///*
    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++)
    {
        t->GetEntry(iEntry);
        Double_t SigmaIfEls = TMath::Sqrt(fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl);
        Double_t SigmaIfMus = TMath::Sqrt(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu);
        if(EventPassedLocal(isMC))
        {
            // plot with sigmas from TPC:
            hSigmasTPC->Fill(SigmaIfEls,SigmaIfMus);

            // plot with dEdx for the negative and positive lepton:
            // if considered muons
            //if(SigmaIfMus < SigmaIfEls){ 
                if(fQ1 < 0) hdEdxMuon->Fill(fTrk1dEdx, fTrk2dEdx);
                else        hdEdxMuon->Fill(fTrk2dEdx, fTrk1dEdx);
            // if considered electrons
            //} else {
            //    if(fQ1 < 0) hdEdxElec->Fill(fTrk1dEdx, fTrk2dEdx);
            //    else        hdEdxElec->Fill(fTrk2dEdx, fTrk1dEdx);
            //}    
        }

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    //*/

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);

    // ******************************************************************************************
    // plot with sigmas from TPC:
    TCanvas *c1 = new TCanvas("c1","c1",900,800);
    c1->SetGrid();
    c1->SetTopMargin(0.03);
    c1->SetBottomMargin(0.145);
    c1->SetRightMargin(0.11);
    c1->SetLeftMargin(0.13);
    if(isMC) c1->SetLogz();

    // horizontal axis
    hSigmasTPC->GetXaxis()->SetTitle("#sqrt{#it{N}#sigma_{e}^{2}(+) + #it{N}#sigma_{e}^{2}(-)}");
    hSigmasTPC->GetXaxis()->SetTitleSize(0.05);
    hSigmasTPC->GetXaxis()->SetTitleOffset(1.15);
    hSigmasTPC->GetXaxis()->SetLabelSize(0.05);
    hSigmasTPC->GetXaxis()->SetDecimals(0);
    // vertical axis
    hSigmasTPC->GetYaxis()->SetTitle("#sqrt{#it{N}#sigma_{#mu}^{2}(+) + #it{N}#sigma_{#mu}^{2}(-)}");
    hSigmasTPC->GetYaxis()->SetTitleSize(0.05);
    hSigmasTPC->GetYaxis()->SetTitleOffset(1.15);
    hSigmasTPC->GetYaxis()->SetLabelSize(0.05);
    hSigmasTPC->GetYaxis()->SetDecimals(0);
    // Z-axis
    //hSigmasTPC->GetZaxis()->SetLabelSize(0.05);
    // Set ranges and draw
    hSigmasTPC->GetXaxis()->SetRangeUser(0.0,20.0);
    hSigmasTPC->GetYaxis()->SetRangeUser(0.0,20.0);
    hSigmasTPC->Draw("COLZ");

    // Draw dashed line y = x
    TLine *line = new TLine(0.0,0.0,20.0,20.0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(9);
    line->Draw("SAME");

    // Legends
    Double_t x_el = 0.55;
    Double_t y_el = 0.85;
    TLegend *leg_el = new TLegend(x_el,y_el,x_el+0.18,y_el+0.06);
    leg_el->AddEntry((TObject*)0,Form("electrons"),""); 
    leg_el->SetTextSize(0.05);
    leg_el->SetBorderSize(0);
    leg_el->SetMargin(0.05);
    //leg_el->SetFillColor(kBlue);
    leg_el->Draw();
    Double_t x_mu = 0.72;
    Double_t y_mu = 0.7;
    TLegend *leg_mu = new TLegend(x_mu,y_mu,x_mu+0.135,y_mu+0.06);
    leg_mu->AddEntry((TObject*)0,Form("muons"),""); 
    leg_mu->SetTextSize(0.05);
    leg_mu->SetBorderSize(0);
    leg_mu->SetMargin(0.05);
    //leg_mu->SetFillColor(kBlue);
    leg_mu->Draw();

    TString str_out1;
    if(isMC) str_out1 = "Results/" + str_subfolder + "ElectronsMuonsPID/h2D_SigmasTPC_MC";
    else     str_out1 = "Results/" + str_subfolder + "ElectronsMuonsPID/h2D_SigmasTPC_data";

    c1->Print((str_out1 + ".pdf").Data());
    c1->Print((str_out1 + ".png").Data());

    // ******************************************************************************************
    // plot with dEdx
    TCanvas *c2 = new TCanvas("c2","c2",900,800);
    c2->SetGrid();
    c2->SetTopMargin(0.03);
    c2->SetBottomMargin(0.145);
    c2->SetRightMargin(0.11);
    c2->SetLeftMargin(0.13);
    if(isMC) c2->SetLogz();

    // horizontal axis
    hdEdxMuon->GetXaxis()->SetTitle("dE/dx^{TPC}(l^{-}) (a.u.)");
    hdEdxMuon->GetXaxis()->SetTitleSize(0.05);
    hdEdxMuon->GetXaxis()->SetTitleOffset(1.15);
    hdEdxMuon->GetXaxis()->SetLabelSize(0.05);
    hdEdxMuon->GetXaxis()->SetDecimals(0);
    // vertical axis
    hdEdxMuon->GetYaxis()->SetTitle("dE/dx^{TPC}(l^{+}) (a.u.)");
    hdEdxMuon->GetYaxis()->SetTitleSize(0.05);
    hdEdxMuon->GetYaxis()->SetTitleOffset(1.15);
    hdEdxMuon->GetYaxis()->SetLabelSize(0.05);
    hdEdxMuon->GetYaxis()->SetDecimals(0);
    // Z-axis
    //hSigmasTPC->GetZaxis()->SetLabelSize(0.05);
    // Set ranges, colors and draw
    hdEdxMuon->GetYaxis()->SetRangeUser(30.,120.0);
    hdEdxMuon->GetXaxis()->SetRangeUser(30.,120.0);
    //hdEdxElec->SetMarkerStyle(20); 
    //hdEdxMuon->SetMarkerStyle(20);
    //hdEdxElec->SetMarkerColor(kBlue);  
    //hdEdxMuon->SetMarkerColor(kRed);
    //hdEdxElec->Draw("COLZ");
    hdEdxMuon->Draw("COLZ");

    TString str_out2;
    if(isMC) str_out2 = "Results/" + str_subfolder + "ElectronsMuonsPID/h2D_dEdx_MC";
    else     str_out2 = "Results/" + str_subfolder + "ElectronsMuonsPID/h2D_dEdx_data";

    c2->Print((str_out2 + ".pdf").Data());
    c2->Print((str_out2 + ".png").Data());

    delete hSigmasTPC;
    delete hdEdxElec;
    delete hdEdxMuon;

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
                // 4) Distance from the IP lower than cut_fVertexZ
                if(fVertexZ > cut_fVertexZ) return kFALSE;
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
                // 4) Distance from the IP lower than cut_fVertexZ
                if(fVertexZ > cut_fVertexZ) return kFALSE;
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
    // (no cut)

    // 9) Dilepton rapidity |y| < cut_fY
    if(!(abs(fY) < cut_fY)) return kFALSE;

    // 10) Pseudorapidity of both tracks |eta| < cut_fEta
    if(!(abs(fEta1) < cut_fEta && abs(fEta2) < cut_fEta)) return kFALSE;

    // 11) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;

    // 13) Transverse momentum cut
    if(!(fPt > 0.2 && fPt < 1.0)) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}