// CountEvents_MC.c
// David Grund, Feb 26, 2022
// To calculate the number of MC events passing the list of selections

Int_t counterMC[19] = { 0 };

Bool_t CountEvents_MC_EventPassed();

void CountEvents_MC_main(){

    TFile *f_in = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_MC_tree_rec.Data()));
    if(t_in) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(t_in);

    TList *l_in = dynamic_cast<TList*> (f_in->Get("AnalysisOutput/fOutputList"));
    if(l_in) Printf("Input list loaded.");

    TH1F *hCounterCuts = (TH1F*)l_in->FindObject("hCounterCuts");
    if(hCounterCuts) Printf("Histogram hCounterCuts loaded.");

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);
        CountEvents_MC_EventPassed();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // Print the numbers:
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CountEvents_MC/");
    TString name = "Results/" + str_subfolder + "CountEvents_MC/cuts.txt";
    ofstream outfile (name.Data());
    outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
    outfile << "0) event non-empty:  \t" << hCounterCuts->GetBinContent(1) << "\n";
    // if pass1
    if(!isPass3){
        outfile << "1) vertex # contribs:\t" << hCounterCuts->GetBinContent(2) << "\n";
        outfile << "2) vertex Z distance:\t" << hCounterCuts->GetBinContent(3) << "\n";
        outfile << "3) two good tracks:  \t" << hCounterCuts->GetBinContent(4) << " (check: " << counterMC[0] << ")\n";
    // if pass3
    } else {
        outfile << "1) two good tracks:  \t" << hCounterCuts->GetBinContent(2) << " (check: " << counterMC[0] << ")\n";
        outfile << "2) vertex # contribs:\t" << counterMC[1] << "\n";
        outfile << "3) vertex Z distance:\t" << counterMC[2] << "\n";
    }
    outfile << "4) CCUP31 trigger:   \t" << counterMC[3] << "\n";
    outfile << "5) good run numbers: \t" << counterMC[4] << "\n";
    outfile << "6a) ADA offline veto:\t" << counterMC[5] << "\n";
    outfile << "6b) ADC offline veto:\t" << counterMC[6] << "\n";
    outfile << "7a) V0A offline veto:\t" << counterMC[7] << "\n";
    outfile << "7b) V0C offline veto:\t" << counterMC[8] << "\n";
    outfile << "8) SPD match FOhits: \t" << counterMC[9] << "\n";
    outfile << "9) muon pairs only:  \t" << counterMC[10] << "\n";
    outfile << "10) dilept |y| < 0.8:\t" << counterMC[11] << "\n";
    outfile << "11) trks |eta| < 0.8:\t" << counterMC[12] << "\n";
    outfile << "12) opposite charges:\t" << counterMC[13] << "\n";
    outfile << "13) mass 2.2 to 4.5: \t" << counterMC[14] << "\n";
    outfile << "14) p_T 0.2 to 1.0:  \t" << counterMC[15] << "\n";
    outfile << "15) mass 3.0 to 3.2: \t" << counterMC[16] << "\n";
    if(isZNcut){
        outfile << "16) ZNA < 10.5 neut: \t" << counterMC[17] << "\n";
        outfile << "17) ZNC < 10.5 neut: \t" << counterMC[18] << "\n";
    }

    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    // Print just the numbers (to be read by _CompareCountsPass1Pass3.C)
    name = "Results/" + str_subfolder + "CountEvents_MC/cuts_numbersOnly.txt";
    outfile.open(name.Data());
    for(Int_t i = 0; i < 19; i++) outfile << counterMC[i] << "\n";
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}

Bool_t CountEvents_MC_EventPassed(){

    // if pass1
    if(!isPass3){
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) At least two tracks associated with the vertex
        // 2) Distance from the IP lower than 15 cm
        // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        counterMC[0]++;

    // if pass3
    } else {
        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        counterMC[0]++;

        // 2) At least two tracks associated with the vertex
        if(fVertexContrib < cut_fVertexContrib) return kFALSE;
        counterMC[1]++;

        // 3) Distance from the IP lower than 15 cm
        if(fVertexZ > cut_fVertexZ) return kFALSE;
        counterMC[2]++;
    }

    // 4) Central UPC trigger CCUP31
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
    counterMC[3]++;

    // 5) Run numbers from the DPG list
    if(!RunNumberInListOfGoodRuns()) return kFALSE;
    counterMC[4]++;

    // 6a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counterMC[5]++;

    // 6b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counterMC[6]++;

    // 7a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counterMC[7]++;

    // 7b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counterMC[8]++;

    // 8) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counterMC[9]++;

    // 9) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counterMC[10]++;

    // 10) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < cut_fY)) return kFALSE;
    counterMC[11]++;

    // 11) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < cut_fEta && abs(fEta2) < cut_fEta)) return kFALSE;
    counterMC[12]++;

    // 12) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counterMC[13]++;

    // 13) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counterMC[14]++;

    // 14) Transverse momentum cut
    if(!(fPt > 0.20)) return kFALSE;
    counterMC[15]++;

    // 15) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counterMC[16]++;

    if(isZNcut){
        Bool_t fZNA_hit = kFALSE;
        Bool_t fZNC_hit = kFALSE;
        Double_t fZNA_n = fZNA_energy / 2510.;
        Double_t fZNC_n = fZNC_energy / 2510.;
        for(Int_t i = 0; i < 4; i++){
            // hit in ZNA
            if(TMath::Abs(fZNA_time[i]) < 2) fZNA_hit = kTRUE;
            // hit in ZNC
            if(TMath::Abs(fZNC_time[i]) < 2) fZNC_hit = kTRUE;
        }    
        // 16) If ZNA signal, then max 10.5 neutrons
        if(fZNA_hit && fZNA_n > cut_fZN_neutrons) return kFALSE;
        counterMC[17]++;

        // 17) If ZNC signal, then max 10.5 neutrons
        if(fZNC_hit && fZNC_n > cut_fZN_neutrons) return kFALSE;
        counterMC[18]++;
    }

    // Event passed all the selections =>
    return kTRUE;
}