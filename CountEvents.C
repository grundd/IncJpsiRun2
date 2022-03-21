// CountEvents.C
// David Grund, Mar 20, 2022
// To calculate the number of events passing the list of selections

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TList.h"
#include "TH1.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

Int_t counterData[18] = { 0 };

Bool_t CountEvents_EventPassed();

void CountEvents(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
    if(t_in) Printf("Input tree loaded.");

    ConnectTreeVariables(t_in);

    TList *l_in = dynamic_cast<TList*> (f_in->Get("AnalysisOutput/fOutputList"));
    if(l_in) Printf("Input list loaded.");

    TH1F *hCounterCuts = (TH1F*)l_in->FindObject("hCounterCuts");
    if(hCounterCuts) Printf("Histogram hCounterCuts loaded.");

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);
        CountEvents_EventPassed();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // Print the numbers:
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CountEvents/");
    TString name = "Results/" + str_subfolder + "CountEvents/cuts.txt";
    ofstream outfile (name.Data());
    outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
    outfile << "0) event non-empty:  \t" << hCounterCuts->GetBinContent(1) << "\n";
    // if pass1
    if(!isPass3){
        outfile << "1) vertex # contribs:\t" << hCounterCuts->GetBinContent(2) << "\n";
        outfile << "2) vertex Z distance:\t" << hCounterCuts->GetBinContent(3) << "\n";
        outfile << "3) two good tracks:  \t" << hCounterCuts->GetBinContent(4) << "\n";
        outfile << "4) CCUP31 trigger:   \t" << hCounterCuts->GetBinContent(5) << " (check: " << counterData[0] << ")\n";
    // if pass3
    } else {
        outfile << "1) two good tracks:  \t" << hCounterCuts->GetBinContent(2) << "\n";
        outfile << "2) CCUP31 trigger:   \t" << hCounterCuts->GetBinContent(3) << " (check: " << counterData[0] << ")\n";
        outfile << "3) vertex # contribs:\t" << counterData[1] << "\n";
        outfile << "4) vertex Z distance:\t" << counterData[2] << "\n";
    }
    outfile << "5) good run numbers: \t" << counterData[3] << "\n";
    outfile << "6a) ADA offline veto:\t" << counterData[4] << "\n";
    outfile << "6b) ADC offline veto:\t" << counterData[5] << "\n";
    outfile << "7a) V0A offline veto:\t" << counterData[6] << "\n";
    outfile << "7b) V0C offline veto:\t" << counterData[7] << "\n";
    outfile << "8) SPD match FOhits: \t" << counterData[8] << "\n";
    outfile << "9) muon pairs only:  \t" << counterData[9] << "\n";
    outfile << "10) dilept |y| < 0.8:\t" << counterData[10] << "\n";
    outfile << "11) trks |eta| < 0.8:\t" << counterData[11] << "\n";
    outfile << "12) opposite charges:\t" << counterData[12] << "\n";
    outfile << "13) mass 2.2 to 4.5: \t" << counterData[13] << "\n";
    outfile << "14) p_T 0.2 to 1.0:  \t" << counterData[14] << "\n";
    outfile << "15) mass 3.0 to 3.2: \t" << counterData[15] << "\n";

    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    // Print just the numbers (to be read by _CompareCountsPass1Pass3.C)
    name = "Results/" + str_subfolder + "CountEvents/cuts_numbersOnly.txt";
    outfile.open(name.Data());
    for(Int_t i = 0; i < 18; i++) outfile << counterData[i] << "\n";
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}

Bool_t CountEvents_EventPassed(){

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
        counterData[0]++;

    // if pass3
    } else {
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        // 2) Central UPC trigger CCUP31:
        // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
        // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
        counterData[0]++;

        // 3) At least two tracks associated with the vertex
        if(fVertexContrib < cut_fVertexContrib) return kFALSE;
        counterData[1]++;

        // 4) Distance from the IP lower than 15 cm
        if(fVertexZ > cut_fVertexZ) return kFALSE;
        counterData[2]++;
    }

    // 5) Run numbers from the DPG list
    if(!RunNumberInListOfGoodRuns()) return kFALSE;
    counterData[3]++;

    // 6a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counterData[4]++;

    // 6b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counterData[5]++;

    // 7a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counterData[6]++;

    // 7b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counterData[7]++;

    // 78) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counterData[8]++;

    // 9) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counterData[9]++;

    // 10) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < cut_fY)) return kFALSE;
    counterData[10]++;

    // 11) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < cut_fEta && abs(fEta2) < cut_fEta)) return kFALSE;
    counterData[11]++;

    // 12) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counterData[12]++;

    // 13) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counterData[13]++;

    // 14) Transverse momentum cut
    if(!(fPt > 0.20 && fPt < 1.00)) return kFALSE;
    counterData[14]++;

    // 15) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counterData[15]++;

    // Event passed all the selections =>
    return kTRUE;
}