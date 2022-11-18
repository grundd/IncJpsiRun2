// MuonMomentumRange.C
// David Grund, Nov 18, 2022
// to answer a comment from the IRC

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

Bool_t EventPassedLocal();

void MuonMomentumRange(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "MuonMomentumRange/");

    TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
    if(t_in) Printf("Input tree loaded.");

    ConnectTreeVariables(t_in);

    TH2F* hMuonMomentumCorr = new TH2F("hMuonMomentumCorr","correlation of momentum of the two muons",40,0.,4.,40,0.,4.);
    TH1F* hMuonMomentum = new TH1F("hMuonMomentum","momentum of muons from J/psi decay",40,0.,4.);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;
    ///*
    for(Int_t i = 0; i < t_in->GetEntries(); i++)
    {
        t_in->GetEntry(i);
        if((i+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

        if(!EventPassedLocal()) continue;

        TLorentzVector muon1;
        muon1.SetPtEtaPhiM(fPt1,fEta1,fPhi1,0.105658);
        TLorentzVector muon2;
        muon2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,0.105658);

        if(fQ1 > 0) hMuonMomentumCorr->Fill(muon1.P(),muon2.P());
        else        hMuonMomentumCorr->Fill(muon2.P(),muon1.P());
        hMuonMomentum->Fill(muon1.P());
        hMuonMomentum->Fill(muon2.P());
    }
    //*/
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    /*
    TLegend *l = new TLegend(0.50,0.65,0.95,0.92);
    l->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    l->AddEntry((TObject*)0,"|#it{y}| < 0.8","");
    l->AddEntry((TObject*)0,"|#it{#eta}(#mu^{#pm})| < 0.8","");
    l->AddEntry((TObject*)0,"2.2 < |#it{m}_{#mu#mu}| < 4.5 GeV/#it{c}^{2}","");
    l->AddEntry((TObject*)0,"0.2 < |#it{p}_{T}| < 1.0 GeV/#it{c}","");
    l->SetMargin(0.);
    l->SetTextSize(0.042);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    */
    // plot the histograms
    // 2d correlation histogram
    TCanvas cCorr("cCorr","cCorr",900,800);
    //cCorr.SetLogz();
    cCorr.SetGrid();
    cCorr.SetTopMargin(0.03);
    cCorr.SetBottomMargin(0.12);
    cCorr.SetRightMargin(0.11);
    cCorr.SetLeftMargin(0.12);
    // horizontal axis
    hMuonMomentumCorr->GetXaxis()->SetTitle("#it{p}(#mu^{+}) [GeV/#it{c}]");
    hMuonMomentumCorr->GetXaxis()->SetTitleSize(0.05);
    hMuonMomentumCorr->GetXaxis()->SetTitleOffset(1.1);
    hMuonMomentumCorr->GetXaxis()->SetLabelSize(0.05);
    // vertical axis
    hMuonMomentumCorr->GetYaxis()->SetTitle("#it{p}(#mu^{-}) [GeV/#it{c}]");
    hMuonMomentumCorr->GetYaxis()->SetTitleSize(0.05);
    hMuonMomentumCorr->GetYaxis()->SetTitleOffset(1.12);
    hMuonMomentumCorr->GetYaxis()->SetLabelSize(0.05);
    hMuonMomentumCorr->Draw("COLZ");
    //l->Draw();
    cCorr.Print("Results/" + str_subfolder + "MuonMomentumRange/hMuonMomentumCorr.pdf");
    // 1d distribution histogram
    TCanvas c("c","c",900,800);
    //c.SetLogy();
    c.SetTopMargin(0.03);
    c.SetBottomMargin(0.12);
    c.SetRightMargin(0.03);
    c.SetLeftMargin(0.13);
    hMuonMomentum->SetLineWidth(2.);
    hMuonMomentum->SetLineColor(kBlue);
    // horizontal axis
    hMuonMomentum->GetXaxis()->SetTitle("#it{p}(#mu^{#pm}) [GeV/#it{c}]");
    hMuonMomentum->GetXaxis()->SetTitleSize(0.05);
    hMuonMomentum->GetXaxis()->SetTitleOffset(1.1);
    hMuonMomentum->GetXaxis()->SetLabelSize(0.05);
    // vertical axis
    hMuonMomentum->GetYaxis()->SetTitle("counts [-]");
    hMuonMomentum->GetYaxis()->SetTitleSize(0.05);
    hMuonMomentum->GetYaxis()->SetTitleOffset(1.25);
    hMuonMomentum->GetYaxis()->SetLabelSize(0.05);
    hMuonMomentum->Draw("HIST");
    //l->Draw();
    c.Print("Results/" + str_subfolder + "MuonMomentumRange/hMuonMomentum.pdf");

    return;
}

Bool_t EventPassedLocal()
{
    // Run number in the GoodHadronPID lists published by DPG
    if(!RunNumberInListOfGoodRuns()) return kFALSE;
    
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
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;

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