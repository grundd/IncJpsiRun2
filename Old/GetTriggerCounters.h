// GetTriggerCounters.h
// David Grund, Feb 26, 2022
// To extract the counters of CCUP31-triggered events per run for both periods

void GetTriggerCounters_main(){

    vector<Int_t> Counts18q;
    vector<Int_t> Counts18r;

    Printf("Selected run list for 18q contains %i runs.", nRuns_18q);
    Printf("Selected run list for 18r contains %i runs.", nRuns_18r);

    TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
    if(f_in) Printf("Input data loaded.");  

    TList *l_in = dynamic_cast<TList*> (f_in->Get("AnalysisOutput/fOutputList"));
    if(l_in) Printf("Input list loaded.");

    TH1F *hCounterTrigger = (TH1F*)l_in->FindObject("hCounterTrigger");
    if(hCounterTrigger) Printf("Histogram hCounterTrigger loaded.");

    hCounterTrigger->Draw();

    Printf("Number of entries in the histogram: %.0f", hCounterTrigger->Integral());
    Double_t counter = 0;
    Int_t iCurrentRun = 0;

    for(Int_t iBin = 1; iBin <= hCounterTrigger->GetNbinsX(); iBin++){
        if(iCurrentRun < nRuns_18q){ // 2018q
            if(hCounterTrigger->GetBinCenter(iBin) == runList_18q[iCurrentRun]){
                Counts18q.push_back(hCounterTrigger->GetBinContent(iBin));
                iCurrentRun++;
                counter += hCounterTrigger->GetBinContent(iBin);
                // Debugging:
                Printf("LHC18q, iRun: %i, Run no.: %.0f, Counter: %.0f", iCurrentRun, hCounterTrigger->GetBinCenter(iBin), hCounterTrigger->GetBinContent(iBin));
            }
        } else { // 2018r
            if(hCounterTrigger->GetBinCenter(iBin) == runList_18r[iCurrentRun-nRuns_18q]){
                Counts18r.push_back(hCounterTrigger->GetBinContent(iBin));
                iCurrentRun++;
                counter += hCounterTrigger->GetBinContent(iBin);
                // Debugging:
                Printf("LHC18r, iRun: %i, Run no.: %.0f, Counter: %.0f", iCurrentRun, hCounterTrigger->GetBinCenter(iBin), hCounterTrigger->GetBinContent(iBin));                
            }
        }     
    }
    Printf("Total number of events found: %.0f", counter);

    // Print the results
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "GetTriggerCounters/");
    TString name = "Results/" + str_subfolder + "GetTriggerCounters/trigger_counters.txt";
    ofstream outfile (name.Data());

    outfile << Form("Counts18q (%i runs):\n", nRuns_18q);
    for(Int_t i = 0; i < nRuns_18q; i++){
        outfile << i+1 << "\t" << runList_18q[i] << "\t" << Counts18q[i] << "\n";
    }

    outfile << Form("\n\nCounts18r (%i runs):\n", nRuns_18r);
    for(Int_t i = 0; i < nRuns_18r; i++){
        outfile << i+1 << "\t" << runList_18r[i] << "\t" << Counts18r[i] << "\n";
    }

    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    // Print the results separately for both periods (to be read by the lumi macro)
    // LHC18q
    TString name_18q = "Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18q.txt";
    ofstream out_f_18q (name_18q.Data());
    for(Int_t i = 0; i < nRuns_18q; i++){
        out_f_18q << Counts18q[i] << "\n";
    }
    out_f_18q.close();
    Printf("*** Results printed to %s.***", name_18q.Data());

    // LHC18r
    TString name_18r = "Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18r.txt";
    ofstream out_f_18r (name_18r.Data());
    for(Int_t i = 0; i < nRuns_18r; i++){
        out_f_18r << Counts18r[i] << "\n";
    }
    out_f_18r.close();
    Printf("*** Results printed to %s.***", name_18r.Data());

    return;
}