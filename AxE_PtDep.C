// AxE_PtDep.C
// David Grund, Mar 28, 2022
// To investigate the pT dependence of AxE

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h" // gStyle
#include "TLegend.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

Int_t nBins = 0;
Double_t edges[25] = {0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.68, 0.76, 0.84, 0.92, 1.00, 1.15, 1.30, 1.45, 1.60};
TH1D *hNRec = NULL;
TH1D *hNGen = NULL;
TH1D* hAxE = NULL;

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
    1,  // 8) Rapidity
    1,  // 9) Pseudorapidity
    1,  // 10) Opposite charges
    1,  // 11) Muons only
    1   // 12) Inv mass
};

// For CalculateRatiosOfNRec:
Bool_t CutsToBePlotted[nCuts-1] = {
    1,  // 1) !0VBA && !0VBC 
    1,  // 2) !0UBA && !0UBC 
    1,  // 3) 0STG
    1,  // 4) 0OMU
    1,  // 5) AD offline veto
    1,  // 6) V0 offline veto
    1,  // 7) SPD cluster matches FOhits
    1,  // 8) Rapidity
    1,  // 9) Pseudorapidity
    1,  // 10) Opposite charges
    1,  // 11) Muons only
    1   // 12) Inv mass
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
Bool_t EventPassedMCRec_AxEPtDep(Int_t iMassCut, Int_t iPtBin);
Bool_t EventPassedMCGen_AxEPtDep(Int_t iPtBin);
void SetPad(TPad* p);
void SaveToFile(TH1D* hist, TString name);
TString ConvertCutsToString();
TString ConvertCutsToBePlottedToString();
Double_t CalculateErrorBayes(Double_t k, Double_t n);

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
    // Get number of bins
    nBins = sizeof(edges) / sizeof(edges[0]) - 1;
    Printf("%i pt bins defined.", nBins);
    // Define the histograms
    hNRec = new TH1D("hNRec","N rec per bin",nBins,edges);
    hNGen = new TH1D("hNGen","N gen per bin",nBins,edges);

    FillHistNRec();
    FillHistNGen();

    hAxE = (TH1D*)hNRec->Clone("hAxE");
    hAxE->SetTitle("AxE per bin");
    hAxE->Sumw2();
    hAxE->Divide(hNGen);

    // Draw the histogram:
    TCanvas *c = new TCanvas("c", "c", 900, 600);
    TPad *p = new TPad("p", "p",0.0,0.0,1.0,1.0);
    p->Draw();
    p->cd();
    SetPad(p);
    // Marker and line
    hAxE->SetMarkerStyle(21);
    hAxE->SetMarkerColor(kBlue);
    hAxE->SetMarkerSize(1.0);
    hAxE->SetLineColor(kBlue);
    hAxE->SetLineWidth(1.0);
    // Vertical axis
    hAxE->GetYaxis()->SetTitle("#it{N}_{rec}/#it{N}_{gen}");
    hAxE->GetYaxis()->SetTitleSize(0.056);
    hAxE->GetYaxis()->SetTitleOffset(1.3);
    hAxE->GetYaxis()->SetLabelSize(0.056);
    hAxE->GetYaxis()->SetDecimals(3);
    //hAxE->GetYaxis()->SetRangeUser(0.10,hAxE->GetBinContent(1)*1.5); // for pt cut only
    // Horizontal axis
    if(cuts[0] == 0) hAxE->GetXaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    if(cuts[0] == 1) hAxE->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hAxE->GetXaxis()->SetTitleSize(0.056);
    hAxE->GetXaxis()->SetTitleOffset(1.2);
    hAxE->GetXaxis()->SetLabelSize(0.056);
    hAxE->GetXaxis()->SetLabelOffset(0.015);
    hAxE->GetXaxis()->SetDecimals(1);
    // Eventually draw it
    hAxE->Draw("P E1");
    // Legend
    TLegend *l = new TLegend(0.52,0.77,0.85,0.97);
    l->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    l->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l->SetTextSize(0.056);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();
    // Legend 2
    TLegend *l2 = new TLegend(0.15,0.17,0.35,0.32);
    l2->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),""); 
    l2->AddEntry((TObject*)0,Form("2.2 < #it{m} < 4.5 GeV/#it{c}^{2}"),"");
    l2->SetTextSize(0.056);
    l2->SetBorderSize(0); // no border
    l2->SetFillStyle(0);  // legend is transparent
    l2->Draw();

    // Save the figures and print the results to txt file
    TString CutConfiguration = ConvertCutsToString();
    TString path(("Results/" + str_subfolder + "AxE_PtDep/fig/" + CutConfiguration).Data());
    c->Print((path + ".pdf").Data());
    c->Print((path + ".png").Data());
    ofstream outfile((path + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << "Bin \tPtLow \tPtUpp \tAxE [%%] \tAxE_err [%%] \n";
    for(Int_t i = 1; i <= nBins; i++){
        outfile << i << "\t" << edges[i-1] << "\t" << edges[i] << "\t" << hAxE->GetBinContent(i)*100 << "\t\t" << hAxE->GetBinError(i)*100 << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", (path + ".txt").Data());

    // Compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kFALSE;
    if(DebugErrors){
        Double_t ErrRoot = 0;
        Double_t ErrBayes = 0;    
        for(Int_t i = 1; i <= nPtBins; i++){
            ErrRoot = hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(hNRec->GetBinContent(i),hNGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // Cross-check: calculate the total value of AxE
    Double_t NRecTot = 0;
    Double_t NGenTot = 0;
    for(Int_t i = 1; i <= nPtBins; i++){
        NRecTot += hNRec->GetBinContent(i);
        NGenTot += hNGen->GetBinContent(i);
    }
    Double_t AxETot = NRecTot / NGenTot;
    Double_t AxETot_err = CalculateErrorBayes(NRecTot, NGenTot);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxETot*100, AxETot_err*100);

    return;
}

void FillHistNRec()
{
    // Check if the corresponding text file already exists
    TString CutConfiguration = ConvertCutsToString();
    TString file(("Results/" + str_subfolder + "AxE_PtDep/" + CutConfiguration + ".txt").Data());

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNRec with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNRec->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N rec per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fRec) Printf("MC rec file loaded.");

        TTree *tRec = dynamic_cast<TTree*> (fRec->Get(str_in_MC_tree_rec.Data()));
        if(tRec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCRec(tRec);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nBins; iPtBin++){
            Int_t NRec = 0;
            for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
                tRec->GetEntry(iEntry);
                if(EventPassedMCRec_AxEPtDep(0, iPtBin)) NRec++;
            }
            hNRec->SetBinContent(iPtBin, NRec);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hNRec, file);

        return;
    }
}

void FillHistNGen()
{
    // Check if the corresponding text file already exists
    TString file(("Results/" + str_subfolder + "AxE_PtDep/NGen.txt").Data());

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNGen with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNGen->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N gen per bin for %s... ***", file.Data());

        TFile *fGen = TFile::Open((str_in_MC_fldr_gen + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fGen) Printf("MC gen file loaded.");

        TTree *tGen = dynamic_cast<TTree*> (fGen->Get(str_in_MC_tree_gen.Data()));
        if(tGen) Printf("MC gen tree loaded.");
        
        ConnectTreeVariablesMCGen(tGen);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nBins; iPtBin++){
            Int_t NGen = 0;
            for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
                tGen->GetEntry(iEntry);
                if(EventPassedMCGen_AxEPtDep(iPtBin)) NGen++;
            }
            hNGen->SetBinContent(iPtBin, NGen);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hNGen, file);

        return;
    }
}

void CalculateRatiosOfNRec()
{
    // Get number of bins
    nBins = sizeof(edges) / sizeof(edges[0]) - 1;
    Printf("%i pt bins defined.", nBins);
    // Define the histogram hNRec
    hNRec = new TH1D("hNRec","N rec per bin",nBins,edges);

    // The first cut (pt_rec or pt_gen) must be selected manually
    // Turn off all the remaining cuts
    for(Int_t i = 1; i < nCuts; i++) cuts[i] = 0;
    // Create an array of histograms to store the results
    TH1D *hNRecRatios[nCuts] = { NULL };
    // Calculate NRec per bin just for the pt cut
    FillHistNRec();
    // Store the results in the zeroth component of hNRecRatios
    hNRecRatios[0] = (TH1D*)hNRec->Clone("hNRecRatios_0");
    // Calculate NRec per bin for every other i-th cut (except the 1st one ofc)
    for(Int_t i = 1; i < nCuts; i++){
        cuts[i] = 1;
        FillHistNRec();
        TString hName("hNRecRatios_%i", i);
        hNRecRatios[i] = (TH1D*)hNRec->Clone(hName.Data());
        cuts[i] = 0;
    }
    // Plot the results
    TH1D* hOne = new TH1D("hOne","ones",nBins,edges);
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
    hOne->GetYaxis()->SetTitle("#it{N}_{rec}[#it{p}_{T} sel. + add. sel.]/#it{N}_{rec}[#it{p}_{T} sel.]");
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
        hNRecRatios[i]->Sumw2();
        hNRecRatios[i]->Divide(hNRecRatios[0]);
        if(CutsToBePlotted[i-1]){
            if(i % 2 == 1) hNRecRatios[i]->SetMarkerStyle(22);
            else           hNRecRatios[i]->SetMarkerStyle(23);
            hNRecRatios[i]->SetMarkerColor(colors[i-1]);
            hNRecRatios[i]->SetMarkerSize(1);
            //hNRecRatios[i]->SetLineStyle(kDashed);
            hNRecRatios[i]->SetLineColor(colors[i-1]);
            hNRecRatios[i]->SetLineWidth(1);
            hNRecRatios[i]->Draw("SAME E1 P");
        }
    }
    // Right pad: plot the legend
    PadR->cd();
    TLegend *l = new TLegend(0.0,0.1,1.0,0.98);
    //l->AddEntry(hist[0],"only S9");
    l->AddEntry(hOne,"only #it{p}_{T} sel. = 1.0");
    l->AddEntry((TObject*)0,Form("with additional:"),"");
    for(Int_t iHist = 1; iHist < nCuts; iHist++){
        l->AddEntry(hNRecRatios[iHist],Form("%s", (CutsLabels[iHist-1]).Data()));
    }
    l->SetTextSize(0.085);
    l->SetBorderSize(0);
    l->SetFillStyle(0); 
    l->Draw();
    // Save the figures
    TString PlotConfiguration = ConvertCutsToBePlottedToString();
    TString Path("Results/" + str_subfolder + "AxE_PtDep/ratios/");
    cRatios->Print((Path + PlotConfiguration + ".pdf").Data());
    cRatios->Print((Path + PlotConfiguration + ".png").Data());

    return;
}

Bool_t EventPassedMCRec_AxEPtDep(Int_t iMassCut = 0, Int_t iPtBin = 0)
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

        // 3) Distance from the IP lower than 15 cm
        if(fVertexZ > cut_fVertexZ) return kFALSE;
    }

    // 4) Central UPC trigger CCUP31:
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

    // 8) Muon pairs only
    if(cuts[8]){
        if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    }

    // 9) Dilepton rapidity |y| < 0.8
    if(cuts[9]){
        if(!(abs(fY) < 0.8)) return kFALSE;
    }

    // 10) Pseudorapidity of both tracks |eta| < 0.8
    if(cuts[10]){
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    }

    // 11) Tracks have opposite charges
    if(cuts[11]){
        if(!(fQ1 * fQ2 < 0)) return kFALSE;
    }

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
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

    // 13) Transverse momentum cut
    if(cuts[0] == kFALSE){ // vs PtRec
        if(!(fPt > edges[iPtBin-1] && fPt <= edges[iPtBin])) return kFALSE;
    } else if(cuts[0] == kTRUE){ // vs PtGen
        if(!(fPtGen > edges[iPtBin-1] && fPtGen < edges[iPtBin])) return kFALSE;
    }

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCGen_AxEPtDep(Int_t iPtBin = 0)
{
    // 1) Dilepton rapidity |y| < 0.8
    if(!(abs(fYGen) < 0.8)) return kFALSE;

    // 2) Transverse momentum cut
    if(!(fPtGen > edges[iPtBin-1] && fPtGen < edges[iPtBin])) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

void SetPad(TPad* p)
{
    p->SetTopMargin(0.02);
    p->SetBottomMargin(0.15);
    p->SetRightMargin(0.03);
    p->SetLeftMargin(0.145);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    return;
}

void SaveToFile(TH1D* hist, TString name)
{
    ofstream outfile (name.Data());
    for(Int_t iBin = 1; iBin <= hist->GetNbinsX(); iBin++){
        outfile << iBin << "\t" << hist->GetBinContent(iBin) << "\n";
    }
    outfile.close();
    Printf("*** Saved to %s.***", name.Data());
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

Double_t CalculateErrorBayes(Double_t k, Double_t n)
{ 
    // k = NRec, n = NGen
    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}