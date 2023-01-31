// AxE_PtDep.cxx
// David Grund, Mar 28, 2022
// To investigate the pT dependence of AxE

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TColor.h"
#include "TPad.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "AxE_Utilities.h"

Int_t nBins = 0;
Float_t edges[25] = {0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.68, 0.76, 0.84, 0.92, 1.00, 1.10, 1.20, 1.30, 1.40};

const Int_t nCuts = 13;
Bool_t cuts[nCuts] = {
    0,  // 0) pt cut rec (kFALSE) or gen (kTRUE)
    1,  // 1) !0VBA (no signal in the V0A) && !0VBC (no signal in the V0C)
    1,  // 2) !0UBA (no signal in the ADA) && !0UBC (no signal in the ADC)
    1,  // 3) 0STG (SPD topological)
    1,  // 4) 0OMU (TOF two hits topology)
    1,  // 5) AD offline veto
    1,  // 6) V0 offline veto
    1,  // 7) SPD cluster matches FOhits
    1,  // 8) rapidity
    1,  // 9) pseudorapidity
    1,  // 10) opposite charges
    1,  // 11) muons only
    1   // 12) inv mass
};
// for CalculateRatiosOfNRec:
Bool_t CutsToBePlotted[nCuts-1] = {
    1,  // 1) !0VBA && !0VBC 
    1,  // 2) !0UBA && !0UBC 
    1,  // 3) 0STG
    1,  // 4) 0OMU
    1,  // 5) AD offline veto
    1,  // 6) V0 offline veto
    1,  // 7) SPD cluster matches FOhits
    1,  // 8) rapidity
    1,  // 9) pseudorapidity
    1,  // 10) opposite charges
    1,  // 11) muons only
    1   // 12) inv mass
};
TString CutsLabels[nCuts-1] = {
    "!0VBA && !0VBC",
    "!0UBA && !0UBC",
    "0STG",
    "0OMU",
    "AD offline",
    "V0 offline",
    "SPD matching",
    "|#it{y}| < 0.8",
    "|#eta_{1,2}| < 0.8",
    "#it{Q}_{1}#it{Q}_{2} < 0",
    "muons only",
    "2.2 < #it{m} < 4.5 GeV/#it{c}^{2}"
};

void CalculateAxEPtDep();
void FillHistNRec();
void FillHistNGen();
void CalculateRatiosOfNRec();
Bool_t EventPassedMCRec_AxEPtDep(Int_t iMassCut);
void SetPad(TPad* p);
TString ConvertCutsToString();
TString ConvertCutsToBePlottedToString();

void AxE_PtDep(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    //#########################################################################
    // AxE plots:
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_PtDep/fig/");

    // vs pT rec, all selections:
    cuts[0] = 0;
    for(Int_t i = 1; i < nCuts; i++) cuts[i] = 1;
    CalculateAxEPtDep();

    // vs pT rec, no other selection:
    for(Int_t i = 0; i < nCuts; i++) cuts[i] = 0;
    CalculateAxEPtDep();

    // vs pT rec, only 0STG:
    for(Int_t i = 0; i < nCuts; i++) cuts[i] = 0;
    cuts[3] = 1; // 0STG
    CalculateAxEPtDep();

    // vs pT rec, only 0OMU:
    for(Int_t i = 0; i < nCuts; i++) cuts[i] = 0;
    cuts[4] = 1; // 0OMU
    CalculateAxEPtDep();

    // vs pT rec, only SPD matching:
    for(Int_t i = 0; i < nCuts; i++) cuts[i] = 0;
    cuts[7] = 1; // SPD matching
    CalculateAxEPtDep();

    // vs pT rec, 0STG, 0OMU and SPD matching:
    for(Int_t i = 0; i < nCuts; i++) cuts[i] = 0;
    cuts[3] = 1; // 0STG
    cuts[4] = 1; // 0OMU
    cuts[7] = 1; // SPD matching
    CalculateAxEPtDep();

    //#########################################################################
    // ratios:
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_PtDep/ratios/");

    // ratios vs pt rec
    cuts[0] = 0;
    CalculateRatiosOfNRec();

    // ratios vs pt gen
    cuts[0] = 1;
    CalculateRatiosOfNRec();
    //#########################################################################

    return;
}

void CalculateAxEPtDep()
{
    // get number of bins
    nBins = sizeof(edges) / sizeof(edges[0]) - 1;
    Printf("%i pt bins defined.", nBins);
    // define the histograms
    hRec = new TH1F("hRec","N rec per bin",nBins,edges);
    hGen = new TH1F("hGen","N gen per bin",nBins,edges);
    // fill them
    FillHistNRec();
    FillHistNGen();
    // calculate AxE
    hAxE = (TH1F*)hRec->Clone("hAxE");
    hAxE->SetTitle("#it{N}_{rec}^{MC}/#it{N}_{gen}^{MC}");
    hAxE->Sumw2();
    hAxE->Divide(hGen);
    // prepare the figure
    TString sCutConfig = ConvertCutsToString();
    TCanvas* c = CreateCanvas("c" + sCutConfig);
    hAxE->SetBit(TH1::kNoStats);
    SetTH1<TH1F>(hAxE,kBlue,kFullCircle,0.8);
    TLegend *l = CreateLegendAxE();
    // axis title and range
    if(cuts[0] == 0) hAxE->GetXaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    if(cuts[0] == 1) hAxE->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    Bool_t onlyPtCut = kTRUE;
    for(Int_t i = 0; i < nCuts; i++) if(cuts[i] == kTRUE) onlyPtCut = kFALSE;
    if(onlyPtCut) hAxE->GetYaxis()->SetRangeUser(hAxE->GetMaximum()*0.75,hAxE->GetMaximum()*1.25);
    // save it
    c->cd();
    hAxE->Draw("E0");
    l->Draw();
    TString path = "Results/" + str_subfolder + "AxE_PtDep/fig/" + sCutConfig;
    c->Print((path + ".pdf").Data());

    // save the figures and print the results to txt file
    SaveToFile(path + ".txt",hAxE,4);

    // compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kFALSE;
    if(DebugErrors){
        Float_t ErrRoot = 0;
        Float_t ErrBayes = 0;    
        for(Int_t i = 1; i <= nPtBins; i++) {
            ErrRoot = hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(hRec->GetBinContent(i),hGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // cross-check: calculate the total value of AxE
    Float_t NRec_tot_val = 0;
    Float_t NGen_tot_val = 0;
    for(Int_t i = 1; i <= nPtBins; i++){
        NRec_tot_val += hRec->GetBinContent(i);
        NGen_tot_val += hGen->GetBinContent(i);
    }
    Float_t AxE_tot_val = NRec_tot_val / NGen_tot_val;
    Float_t AxE_tot_err = CalculateErrorBayes(NGen_tot_val, NGen_tot_val);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxE_tot_val*100., AxE_tot_err*100.);
    delete hRec;
    delete hGen;
    delete hAxE;
    return;
}

void FillHistNRec()
{
    // check if the corresponding text file already exists
    TString sCutConfig = ConvertCutsToString();
    TString file(("Results/" + str_subfolder + "AxE_PtDep/" + sCutConfig + ".txt").Data());

    ifstream ifs;
    ifs.open(file);
    if(!ifs.fail()) {
        // this configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // fill hRec with data from the text file
        Int_t bin;
        Float_t val, err;
        for(Int_t iBin = 0; iBin < nBins; iBin++) {
            ifs >> bin >> val >> err;
            hRec->SetBinContent(bin, val);
            hRec->SetBinError(bin, err);
        }
        ifs.close(); 
        return;
    } else {
        // this configuration needs to be calculated
        Printf("*** Calculating N_rec per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fRec) Printf("MC rec file loaded.");
        TTree *tRec = dynamic_cast<TTree*> (fRec->Get(str_in_MC_tree_rec.Data()));
        if(tRec) Printf("MC rec tree loaded.");
        ConnectTreeVariablesMCRec(tRec);

        // load the ratio to re-weight the spectra
        TH1F* hRatios = GetRatioHisto();
        TAxis* xAxis = hRatios->GetXaxis();
        // loop over tree entries
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) {
            tRec->GetEntry(iEntry);
            if(EventPassedMCRec_AxEPtDep(0)) {
                Int_t iBinGen = xAxis->FindBin(fPtGen);
                Float_t weight = hRatios->GetBinContent(iBinGen);
                if(cuts[0] == kFALSE) hRec->Fill(fPt, weight);
                if(cuts[0] == kTRUE)  hRec->Fill(fPtGen, weight);
            }
        }
        Printf("*** Finished! ***");
        SaveToFile(file,hRec);
        return;
    }
}

void FillHistNGen()
{
    // Check if the corresponding text file already exists
    TString file(("Results/" + str_subfolder + "AxE_PtDep/NGen.txt").Data());

    ifstream ifs;
    ifs.open(file);
    if(!ifs.fail()){
        // this configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // fill hGen with data from the text file
        Int_t bin;
        Float_t val, err;
        while(!ifs.eof()){
            ifs >> bin >> val >> err; // fist and second column
            hGen->SetBinContent(bin, val);
            hGen->SetBinError(bin, err);
        }
        ifs.close(); 
        return;
    } else {
        // this configuration needs to be calculated
        Printf("*** Calculating N_gen per bin for %s... ***", file.Data());

        TFile *fGen = TFile::Open((str_in_MC_fldr_gen + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fGen) Printf("MC gen file loaded.");
        TTree *tGen = dynamic_cast<TTree*> (fGen->Get(str_in_MC_tree_gen.Data()));
        if(tGen) Printf("MC gen tree loaded.");
        ConnectTreeVariablesMCGen(tGen);

        // load the ratio to re-weight the spectra
        TH1F* hRatios = GetRatioHisto();
        TAxis* xAxis = hRatios->GetXaxis();
        // loop over tree entries
        for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++) {
            tGen->GetEntry(iEntry);
            if(EventPassedMCGen()) {
                Int_t iBinGen = xAxis->FindBin(fPtGen);
                Float_t weight = hRatios->GetBinContent(iBinGen);
                hGen->Fill(fPtGen, weight);
            }
        }
        Printf("*** Finished! ***");
        SaveToFile(file,hGen);
        return;
    }
}

void CalculateRatiosOfNRec()
{
    gStyle->SetOptTitle(0);
    // Get number of bins
    nBins = sizeof(edges) / sizeof(edges[0]) - 1;
    Printf("%i pt bins defined.", nBins);
    // Define the histogram hRec
    hRec = new TH1F("hRec","N rec per bin",nBins,edges);

    // The first cut (pt_rec or pt_gen) must be selected manually
    // Turn off all the remaining cuts
    for(Int_t i = 1; i < nCuts; i++) cuts[i] = 0;
    // Create an array of histograms to store the results
    TH1F *hRecRatios[nCuts] = { NULL };
    // Calculate NRec per bin just for the pt cut
    FillHistNRec();
    // Store the results in the zeroth component of hRecRatios
    hRecRatios[0] = (TH1F*)hRec->Clone("hRecRatios_0");
    delete hRec;
    // Calculate NRec per bin for every other i-th cut (except the 1st one ofc)
    for(Int_t i = 1; i < nCuts; i++){
        cuts[i] = 1;
        hRec = new TH1F("hRec","N rec per bin",nBins,edges);
        FillHistNRec();
        TString hName("hRecRatios_%i", i);
        hRecRatios[i] = (TH1F*)hRec->Clone(hName.Data());
        delete hRec;
        cuts[i] = 0;
    }
    // Plot the results
    TH1F* hOne = new TH1F("hOne","ones",nBins,edges);
    for(Int_t i = 1; i <= nBins; i++){
        hOne->SetBinContent(i, 1.);
    }
    TCanvas *cRatios = new TCanvas("cRatios","cRatios",1000,600);
    TPad *PadL = new TPad("PadL", "PadL",0.00,0.0,0.72,1.0);
    TPad *PadR = new TPad("PadR", "PadR",0.72,0.0,1.00,1.0);
    PadL->Draw();
    PadR->Draw();
    SetPad(PadL);
    // Left pad
    PadL->cd();
    // Plot the histogram with ones
    // Marker and line
    hOne->SetMarkerStyle(21);
    hOne->SetMarkerColor(kBlack);
    hOne->SetMarkerSize(1.0);
    hOne->SetLineStyle(kDashed);
    hOne->SetLineColor(kBlack);
    hOne->SetLineWidth(2.0);
    // Vertical axis
    hOne->GetYaxis()->SetTitle("#it{N}_{rec}^{MC}[#it{p}_{T} sel. + add. sel.]/#it{N}_{rec}^{MC}[#it{p}_{T} sel.]");
    hOne->GetYaxis()->SetTitleSize(0.056);
    hOne->GetYaxis()->SetTitleOffset(1.1);
    hOne->GetYaxis()->SetLabelSize(0.056);
    hOne->GetYaxis()->SetDecimals(3);
    hOne->GetYaxis()->SetRangeUser(0.0,hOne->GetBinContent(1)*1.1);
    // Horizontal axis
    if(cuts[0] == 0) hOne->GetXaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    if(cuts[0] == 1) hOne->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hOne->GetXaxis()->SetTitleSize(0.056);
    hOne->GetXaxis()->SetTitleOffset(1.2);
    hOne->GetXaxis()->SetLabelSize(0.056);
    hOne->GetXaxis()->SetLabelOffset(0.015);
    hOne->GetXaxis()->SetDecimals(1);
    // Eventually draw it
    hOne->Draw("HIST");
    // Plot the rest of ratios
    Int_t colors[nCuts] = {2,   // red
                           634, // dark red
                           7,   // cyan
                           4,   // blue
                           8,   // green
                           209, // dark green
                           6,   // pink
                           222, // magenta
                           12,  // grey
                           92,  // gold
                           807, // orange
                           402, // dark yellow 
                           159};// 
    for(Int_t i = 1; i < nCuts; i++){
        hRecRatios[i]->Sumw2();
        hRecRatios[i]->Divide(hRecRatios[0]);
        if(CutsToBePlotted[i-1]){
            if(i % 2 == 1) hRecRatios[i]->SetMarkerStyle(22);
            else           hRecRatios[i]->SetMarkerStyle(23);
            hRecRatios[i]->SetMarkerColor(colors[i-1]);
            hRecRatios[i]->SetMarkerSize(1);
            //hRecRatios[i]->SetLineStyle(kDashed);
            hRecRatios[i]->SetLineColor(colors[i-1]);
            hRecRatios[i]->SetLineWidth(1);
            hRecRatios[i]->Draw("SAME E1 P");
        }
    }
    // Right pad: plot the legend
    PadR->cd();
    TLegend *l = new TLegend(0.0,0.1,1.0,0.98);
    //l->AddEntry(hist[0],"only S9");
    l->AddEntry(hOne,"only #it{p}_{T} sel. = 1.0","L");
    l->AddEntry((TObject*)0,Form("with additional:"),"");
    for(Int_t iHist = 1; iHist < nCuts; iHist++){
        l->AddEntry(hRecRatios[iHist],Form("%s", (CutsLabels[iHist-1]).Data()),"LP");
    }
    l->SetTextSize(0.085);
    l->SetBorderSize(0);
    l->SetFillStyle(0); 
    l->Draw();
    // Save the figures
    TString PlotConfiguration = ConvertCutsToBePlottedToString();
    TString Path("Results/" + str_subfolder + "AxE_PtDep/ratios/");
    cRatios->Print((Path + PlotConfiguration + ".pdf").Data());

    return;
}

Bool_t EventPassedMCRec_AxEPtDep(Int_t iMassCut = 0)
{
    // Run number in the GoodHadronPID lists published by DPG
    if(!RunNumberInListOfGoodRuns()) return kFALSE;

    // if pass1
    if(!isPass3){
        // selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) At least two tracks associated with the vertex
        // 2) Distance from the IP lower than 15 cm
        // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

    // if pass3
    } else {
        // selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

        // 2) At least two tracks associated with the vertex
        if(fVertexContrib < cut_fVertexContrib) return kFALSE;

        // 3) Distance from the IP lower than cut_fVertexZ
        if(fVertexZ > cut_fVertexZ) return kFALSE;
    }

    // 4) central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    if(cuts[1]){ 
        // !0VBA (no signal in the V0A)
        // and
        // !0VBC (no signal in the V0C)
        if(fTriggerInputsMC[0] || fTriggerInputsMC[1]) return kFALSE;
    }
    if(cuts[2]){ 
        // !0UBA (no signal in the ADA)
        // and
        // !0UBC (no signal in the ADC)
        if(fTriggerInputsMC[2] || fTriggerInputsMC[3]) return kFALSE;
    }
    if(cuts[3]){ 
        // 0STG (SPD topological)
        if(!fTriggerInputsMC[10]) return kFALSE;
    }
    if(cuts[4]){ 
        // 0OMU (TOF two hits topology)
        if(!fTriggerInputsMC[4]) return kFALSE;
    }

    // 5) AD offline veto (negligible effect on MC)
    if(cuts[5]){
        if(!(fADA_dec == 0 && fADC_dec == 0)) return kFALSE;
    }

    // 6) V0 offline veto (negligible effect on MC)
    if(cuts[6]){
        if(!(fV0A_dec == 0 && fV0C_dec == 0)) return kFALSE;
    }

    // 7) SPD cluster matches FOhits
    if(cuts[7]){
        if(!(fMatchingSPD == kTRUE)) return kFALSE;
    }

    // 8) muon pairs only
    if(cuts[8]){
        if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    }

    // 9) dilepton rapidity |y| < cut_fY
    if(cuts[9]){
        if(!(abs(fY) < cut_fY)) return kFALSE;
    }

    // 10) pseudorapidity of both tracks |eta| < cut_fEta
    if(cuts[10]){
        if(!(abs(fEta1) < cut_fEta && abs(fEta2) < cut_fEta)) return kFALSE;
    }

    // 11) tracks have opposite charges
    if(cuts[11]){
        if(!(fQ1 * fQ2 < 0)) return kFALSE;
    }

    // 12) invariant mass between 2.2 and 4.5 GeV/c^2
    if(cuts[12]){
        Bool_t bMassCut = kFALSE;
        switch(iMassCut){
            case -1: // No inv mass cut
                bMassCut = kTRUE;
                break;
            case 0:
                if(fM > 2.2 && fM < 4.5) bMassCut = kTRUE;
                break;
            case 1:
                if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
                break;
        }
        if(!bMassCut) return kFALSE;
    }

    // Event passed all the selections =>
    return kTRUE;
}

void SetPad(TPad* p)
{
    p->SetTopMargin(0.02);
    p->SetBottomMargin(0.15);
    p->SetRightMargin(0.03);
    p->SetLeftMargin(0.145);

    //gStyle->SetOptTitle(0);
    //gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);
    //gStyle->SetPaintTextFormat("4.2f");

    return;
}

TString ConvertCutsToString()
{
    TString s("cuts_");
    for(Int_t iCut = 0; iCut < nCuts; iCut++){
        if(cuts[iCut] == kTRUE) s.Append("1");
        else s.Append("0");
    }
    return s;
}

TString ConvertCutsToBePlottedToString()
{
    TString s("ratios_");
    if(cuts[0] == 0) s.Append("0");
    if(cuts[0] == 1) s.Append("1");
    for(Int_t iCut = 0; iCut < nCuts-1; iCut++){
        if(CutsToBePlotted[iCut] == kTRUE) s.Append("1");
        else s.Append("0");
    }
    return s;
}