// AnalysisManager.h
// David Grund, Feb 26, 2022
// Functions to set branch addresses and to check if events pass selection criteria

// cpp headers
#include <stdio.h> // printf
// root headers
#include "TTree.h"
// my headers
#include "ListsOfGoodRuns.h"

// Selections to set:
Int_t cut_fVertexContrib = 2;
Double_t cut_fVertexZ = 10.;
Double_t cut_fY = 0.8;
Double_t cut_fEta = 0.8;
Double_t fN_DSCB = 10.;
// Options that will be set by the choice of iAnalysis in InitAnalysis:
vector<Int_t> runList_18q;
vector<Int_t> runList_18r;
Int_t nRuns_18q; 
Int_t nRuns_18r; 
Bool_t isPass3;
Bool_t isPIDCalibrated; // if NSigmas in MC data were shifted to zeros
Bool_t isNParInDSCBFixed; // if the values of the tail parameters "N" in DSCB are fixed to N_DSCB
Bool_t areBinYieldsUniform; // see the macro BinsThroughMassFit.C
// Array containing pT bin boundaries (will be created in SetPtBinning.h):
Double_t *ptBoundaries = NULL;
Int_t nPtBins;
// Array containing pT bin boundaries for pT fit (will be created in PtFit_SubtractBackground.h):
Double_t *ptBoundaries_PtFit = NULL;
Double_t *tBoundaries_PtFit = NULL;
Int_t nPtBins_PtFit;

// variables for both pass1 and pass3:
Int_t fRunNumber;
TString *fTriggerName = NULL;
Bool_t fTriggerInputsMC[11];
Double_t fTrk1SigIfMu, fTrk1SigIfEl, fTrk2SigIfMu, fTrk2SigIfEl;
Double_t fPt, fM, fY, fPhi;
Double_t fPt1, fPt2, fEta1, fEta2, fPhi1, fPhi2, fQ1, fQ2;
Double_t fZNA_energy, fZNC_energy;
Double_t fZNA_time[4], fZNC_time[4];
Int_t fV0A_dec, fV0C_dec, fADA_dec, fADC_dec;
Bool_t fMatchingSPD;
Double_t fPtGen, fYGen, fMGen, fPhiGen;
// only for pass1:
Double_t fV0A_time, fV0C_time, fADA_time, fADC_time;
// only for pass3:
Double_t fTrk1dEdx, fTrk2dEdx, fVertexZ;
Int_t fVertexContrib;
// only for pass3 & psi(2s) datasets:
Double_t fPtGen_Psi2s;

void SetReducedRunList(Bool_t pass3)
{
    if(!pass3){
        nRuns_18q = sizeof(DPG_LHC18q_pass1_reduced) / sizeof(DPG_LHC18q_pass1_reduced[0]); // 123 runs
        for(Int_t i = 0; i < nRuns_18q; i++) runList_18q.push_back(DPG_LHC18q_pass1_reduced[i]);
        nRuns_18r = sizeof(DPG_LHC18r_pass1_reduced) / sizeof(DPG_LHC18r_pass1_reduced[0]); // 96 runs
        for(Int_t i = 0; i < nRuns_18r; i++) runList_18r.push_back(DPG_LHC18r_pass1_reduced[i]);
    } else {
        nRuns_18q = sizeof(DPG_LHC18q_pass3_reduced) / sizeof(DPG_LHC18q_pass3_reduced[0]); // 122 runs
        for(Int_t i = 0; i < nRuns_18q; i++) runList_18q.push_back(DPG_LHC18q_pass3_reduced[i]);
        nRuns_18r = sizeof(DPG_LHC18r_pass3_reduced) / sizeof(DPG_LHC18r_pass3_reduced[0]); // 96 runs
        for(Int_t i = 0; i < nRuns_18r; i++) runList_18r.push_back(DPG_LHC18r_pass3_reduced[i]);
    }
    Printf("Number of runs in LHC18q run list: %i (%i)", nRuns_18q, (Int_t)runList_18q.size());
    Printf("Number of runs in LHC18r run list: %i (%i)", nRuns_18r, (Int_t)runList_18r.size());

    return;
}

void ConnectTreeVariables(TTree *t)
{
    // Set branch addresses
    // Basic things:
    t->SetBranchAddress("fRunNumber", &fRunNumber);
    t->SetBranchAddress("fTriggerName", &fTriggerName);
    // PID, sigmas:
    t->SetBranchAddress("fTrk1SigIfMu", &fTrk1SigIfMu);
    t->SetBranchAddress("fTrk1SigIfEl", &fTrk1SigIfEl);
    t->SetBranchAddress("fTrk2SigIfMu", &fTrk2SigIfMu);
    t->SetBranchAddress("fTrk2SigIfEl", &fTrk2SigIfEl);
    // Kinematics:
    t->SetBranchAddress("fPt", &fPt);
    t->SetBranchAddress("fPhi", &fPhi);
    t->SetBranchAddress("fY", &fY);
    t->SetBranchAddress("fM", &fM);
    // Two tracks:
    t->SetBranchAddress("fPt1", &fPt1);
    t->SetBranchAddress("fPt2", &fPt2);
    t->SetBranchAddress("fEta1", &fEta1);
    t->SetBranchAddress("fEta2", &fEta2);
    t->SetBranchAddress("fPhi1", &fPhi1);
    t->SetBranchAddress("fPhi2", &fPhi2);
    t->SetBranchAddress("fQ1", &fQ1);
    t->SetBranchAddress("fQ2", &fQ2);
    // ZDC:
    t->SetBranchAddress("fZNA_energy", &fZNA_energy);
    t->SetBranchAddress("fZNC_energy", &fZNC_energy);
    t->SetBranchAddress("fZNA_time", &fZNA_time);
    t->SetBranchAddress("fZNC_time", &fZNC_time);
    // V0:
    t->SetBranchAddress("fV0A_dec", &fV0A_dec);
    t->SetBranchAddress("fV0C_dec", &fV0C_dec);
    // AD:
    t->SetBranchAddress("fADA_dec", &fADA_dec);
    t->SetBranchAddress("fADC_dec", &fADC_dec);
    // Matching SPD clusters with FOhits:
    t->SetBranchAddress("fMatchingSPD", &fMatchingSPD);
    // if pass3
    if(isPass3){
        t->SetBranchAddress("fVertexZ", &fVertexZ);
        t->SetBranchAddress("fVertexContrib", &fVertexContrib);
        t->SetBranchAddress("fTrk1dEdx", &fTrk1dEdx);
        t->SetBranchAddress("fTrk2dEdx", &fTrk2dEdx);
    // if not
    } else {
        t->SetBranchAddress("fV0A_time", &fV0A_time);
        t->SetBranchAddress("fV0C_time", &fV0C_time);
        t->SetBranchAddress("fADA_time", &fADA_time);
        t->SetBranchAddress("fADC_time", &fADC_time);
    }

    Printf("Variables from %s connected.", t->GetName());
    return;
}

void ConnectTreeVariablesMCRec(TTree *t, Bool_t isPsi2sDataset = kFALSE)
{
    // Set branch addresses
    // Basic things:
    t->SetBranchAddress("fRunNumber", &fRunNumber);
    t->SetBranchAddress("fTriggerInputsMC", &fTriggerInputsMC);
    // PID, sigmas:
    t->SetBranchAddress("fTrk1SigIfMu", &fTrk1SigIfMu);
    t->SetBranchAddress("fTrk1SigIfEl", &fTrk1SigIfEl);
    t->SetBranchAddress("fTrk2SigIfMu", &fTrk2SigIfMu);
    t->SetBranchAddress("fTrk2SigIfEl", &fTrk2SigIfEl);
    // Kinematics:
    t->SetBranchAddress("fPt", &fPt);
    t->SetBranchAddress("fPhi", &fPhi);
    t->SetBranchAddress("fY", &fY);
    t->SetBranchAddress("fM", &fM);
    // Two tracks:
    t->SetBranchAddress("fPt1", &fPt1);
    t->SetBranchAddress("fPt2", &fPt2);
    t->SetBranchAddress("fEta1", &fEta1);
    t->SetBranchAddress("fEta2", &fEta2);
    t->SetBranchAddress("fPhi1", &fPhi1);
    t->SetBranchAddress("fPhi2", &fPhi2);
    t->SetBranchAddress("fQ1", &fQ1);
    t->SetBranchAddress("fQ2", &fQ2);
    // ZDC:
    t->SetBranchAddress("fZNA_energy", &fZNA_energy);
    t->SetBranchAddress("fZNC_energy", &fZNC_energy);
    t->SetBranchAddress("fZNA_time", &fZNA_time);
    t->SetBranchAddress("fZNC_time", &fZNC_time);
    // V0:
    t->SetBranchAddress("fV0A_dec", &fV0A_dec);
    t->SetBranchAddress("fV0C_dec", &fV0C_dec);
    // AD:
    t->SetBranchAddress("fADA_dec", &fADA_dec);
    t->SetBranchAddress("fADC_dec", &fADC_dec);
    // Matching SPD clusters with FOhits:
    t->SetBranchAddress("fMatchingSPD", &fMatchingSPD);
    // MC kinematics on generator level
    t->SetBranchAddress("fPtGen", &fPtGen);
    t->SetBranchAddress("fPhiGen", &fPhiGen);
    t->SetBranchAddress("fYGen", &fYGen);
    t->SetBranchAddress("fMGen", &fMGen);
    // if pass3
    if(isPass3){
        t->SetBranchAddress("fVertexZ", &fVertexZ);
        t->SetBranchAddress("fVertexContrib", &fVertexContrib);
        t->SetBranchAddress("fTrk1dEdx", &fTrk1dEdx);
        t->SetBranchAddress("fTrk2dEdx", &fTrk2dEdx);
        if(isPsi2sDataset){
            t->SetBranchAddress("fPtGen_Psi2s", &fPtGen_Psi2s);
        }
    // if not
    } else {
        t->SetBranchAddress("fV0A_time", &fV0A_time);
        t->SetBranchAddress("fV0C_time", &fV0C_time);
        t->SetBranchAddress("fADA_time", &fADA_time);
        t->SetBranchAddress("fADC_time", &fADC_time);
    }

    Printf("Variables from %s connected.", t->GetName());
    return;
}

void ConnectTreeVariablesMCGen(TTree *t)
{
    // Set branch addresses
    // Basic things:
    t->SetBranchAddress("fRunNumber", &fRunNumber);
    // MC kinematics on generator level
    t->SetBranchAddress("fPtGen", &fPtGen);
    t->SetBranchAddress("fPhiGen", &fPhiGen);
    t->SetBranchAddress("fYGen", &fYGen);
    t->SetBranchAddress("fMGen", &fMGen);

    Printf("Variables from %s connected.", t->GetName());
    return;
}

Bool_t RunNumberInListOfGoodRuns()
{
    // Run number in the GoodHadronPID lists published by DPG
    Bool_t GoodRunNumber = kFALSE;
    if(std::count(runList_18q.begin(), runList_18q.end(), fRunNumber) > 0) GoodRunNumber = kTRUE;
    if(std::count(runList_18r.begin(), runList_18r.end(), fRunNumber) > 0) GoodRunNumber = kTRUE;
    if(!GoodRunNumber){
        //Printf("Wrong run number: %i.", fRunNumber);
        return kFALSE;
    } else {
        return kTRUE;
    }
}

Bool_t EventPassed(Int_t iMassCut, Int_t iPtCut)
{
    // Run number in the GoodHadronPID lists published by DPG
    if(!RunNumberInListOfGoodRuns()) return kFALSE;

    // if pass1
    if(!isPass3){
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) At least two tracks associated with the vertex
        // 2) Distance from the IP lower than 15 cm
        // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        // 4) Central UPC trigger CCUP31:
        // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
        // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD

    // if pass3
    } else {
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        // 2) Central UPC trigger CCUP31:
        // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
        // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD

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

    // 12) Invariant mass cut
    Bool_t bMassCut = kFALSE;
    switch(iMassCut){
        case -1: // no inv mass cut
            bMassCut = kTRUE;
            break;
        case 0: // m between 2.2 and 4.5 GeV/c^2
            if(fM > 2.2 && fM < 4.5) bMassCut = kTRUE;
            break; 
        case 1: // m between 3.0 and 3.2 GeV/c^2
            if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
            break;
        case 2: // m between 1.5 and 7.0 GeV/c^2 (syst uncertainties in inv mass fit)
            if(fM > 1.5 && fM < 7.0) bMassCut = kTRUE;
            break;
    }
    if(!bMassCut) return kFALSE;

    // 13) Transverse momentum cut
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case -1: // no pt cut
            bPtCut = kTRUE;
            break;
        case 0: // 'inc': incoherent-enriched sample
            if(fPt > 0.20) bPtCut = kTRUE;
            break;
        case 1: // 'coh': coherent-enriched sample (~ Roman)
            if(fPt < 0.11) bPtCut = kTRUE;
            break;
        case 2: // 'all': total sample (pT < 2.0 GeV/c)
            if(fPt < 2.00) bPtCut = kTRUE;
            break;
        case 3: // 'allbins': sample with pT from 0.2 to 1 GeV/c 
            if(fPt > 0.20 && fPt < 1.00) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCRec(Int_t iMassCut, Int_t iPtCut, Int_t iPtBin = -1)
{
    // Run number in the GoodHadronPID lists published by DPG
    if(!RunNumberInListOfGoodRuns()) return kFALSE;

    // if pass1
    if(!isPass3){
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) At least two tracks associated with the vertex
        // 2) Distance from the IP lower than 15 cm
        // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

    // if pass3
    } else {
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

        // 2) At least two tracks associated with the vertex
        if(fVertexContrib < cut_fVertexContrib) return kFALSE;

        // 3) Distance from the IP lower than cut_fVertexZ
        if(fVertexZ > cut_fVertexZ) return kFALSE;
    }

    // 4) Central UPC trigger CCUP31:
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

    // 12) Invariant mass cut
    Bool_t bMassCut = kFALSE;
    switch(iMassCut){
        case -1: // no inv mass cut
            bMassCut = kTRUE;
            break;
        case 0: // m between 2.2 and 4.5 GeV/c^2
            if(fM > 2.2 && fM < 4.5) bMassCut = kTRUE;
            break; 
        case 1: // m between 3.0 and 3.2 GeV/c^2
            if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
            break;
    }
    if(!bMassCut) return kFALSE;

    // 13) Transverse momentum cut
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case -1: // no pT cut
            bPtCut = kTRUE;
            break;
        case 0: // 'inc': incoherent-enriched sample
            if(fPt > 0.20) bPtCut = kTRUE;
            break;
        case 1: // 'coh': coherent-enriched sample (~ Roman)
            if(fPt < 0.11) bPtCut = kTRUE;
            break;
        case 2: // 'all': total sample (pT < 2.0 GeV/c)
            if(fPt < 2.00) bPtCut = kTRUE;
            break;
        case 3: // 'allbins': sample with pT from 0.2 to 1 GeV/c 
            if(fPt > 0.20 && fPt < 1.00) bPtCut = kTRUE;
            break;
        case 4: // pT bins (4 or 5)
            if(fPt > ptBoundaries[iPtBin-1] && fPt <= ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCGen(Int_t iPtCut = -1, Int_t iPtBin = -1)
{
    // 1) Dilepton rapidity |y| < cut_fY
    if(!(abs(fYGen) < cut_fY)) return kFALSE;

    // 2) Transverse momentum cut (default: none)
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case -1: // no pT cut
            bPtCut = kTRUE;
            break;
        case 0: // no pT cut
            bPtCut = kTRUE;
            break;
        case 2: // no pT cut
            bPtCut = kTRUE;
            break;
        case 3: // sample with pT from 0.2 to 1 GeV/c 
            if(fPtGen > 0.20 && fPtGen < 1.00) bPtCut = kTRUE;
            break;
        case 4: // pT bins (4 or 5)
            if(fPtGen > ptBoundaries[iPtBin-1] && fPtGen <= ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}