// RunListCheck.h
// David Grund, Feb 27, 2021
// To check if the run numbers of analysed data match the official list

void RunListCheck_main(){

    Int_t nRuns_18qr = nRuns_18q + nRuns_18r;
    vector<Bool_t> RunNumbersFound_18q(nRuns_18q, kFALSE);
    vector<Bool_t> RunNumbersFound_18r(nRuns_18r, kFALSE);
    vector<Int_t> RunNumbersNotFound;
    
    TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
    if(t_in) Printf("Input tree loaded.");

    ConnectTreeVariables(t_in);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "RunListCheck/");
    ofstream outfile("Results/" + str_subfolder + "RunListCheck/status.txt");

    // Check if the current run number is on the list
    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;
    Int_t nRunNumbersNotFound = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);

        // Run number from the GoodHadronPID lists published by DPG
        // https://www.techiedelight.com/find-index-element-vector-cpp/
        Bool_t isRunIn18q = kFALSE;
        Bool_t isRunIn18r = kFALSE;
        std::vector<Int_t>::iterator itr_q = std::find(runList_18q.begin(), runList_18q.end(), fRunNumber);
        std::vector<Int_t>::iterator itr_r = std::find(runList_18r.begin(), runList_18r.end(), fRunNumber);
        if(itr_q != runList_18q.cend()) isRunIn18q = kTRUE;
        if(itr_r != runList_18r.cend()) isRunIn18r = kTRUE;

        if(isRunIn18q == kTRUE && isRunIn18r == kTRUE){
            outfile << Form("Error. Terminating...\n");
        }

        // if run number found
        if(isRunIn18q || isRunIn18r){
            if(isRunIn18q) RunNumbersFound_18q[std::distance(runList_18q.begin(), itr_q)] = kTRUE;
            if(isRunIn18r) RunNumbersFound_18r[std::distance(runList_18r.begin(), itr_r)] = kTRUE;
        }

        // if run number not found
        if(isRunIn18q == kFALSE && isRunIn18r == kFALSE){
            Bool_t isMissingRunAlreadyKnown = kFALSE;
            std::vector<Int_t>::iterator itr = std::find(RunNumbersNotFound.begin(), RunNumbersNotFound.end(), fRunNumber);
            if(itr != RunNumbersNotFound.cend()) isMissingRunAlreadyKnown = kTRUE;
            if(!isMissingRunAlreadyKnown){
                RunNumbersNotFound.push_back(fRunNumber);
                outfile << Form("Run number %i not in the list.\n", fRunNumber);
                nRunNumbersNotFound++;
            }
        }

        // number of events analyzed
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // check how many run numbers that are not in the list was found
    if(nRunNumbersNotFound == 0){
        outfile << Form("No run numbers outside the list of good runs found.\n");
    } else {
        outfile << Form("%i run numbers outside the list of good runs found:\n", nRunNumbersNotFound);
        for(Int_t i = 0; i < (Int_t)RunNumbersNotFound.size(); i++){
            outfile << RunNumbersNotFound[i] << endl;
        }
    }

    // check if all the run numbers from the list were found
    Bool_t allRunNumbersFound = kTRUE;
    for(Int_t i = 0; i < nRuns_18q; i++){
        if(RunNumbersFound_18q[i] == kFALSE){
            allRunNumbersFound = kFALSE;
            outfile << Form("Run number %i from the list LHC18q was not found in the tree.\n", runList_18q[i]);       
        } 
    }
    for(Int_t i = 0; i < nRuns_18r; i++){
        if(RunNumbersFound_18r[i] == kFALSE){
            allRunNumbersFound = kFALSE;
            outfile << Form("Run number %i from the list LHC18r was not found in the tree.\n", runList_18r[i]);   
        } 
    }
    if(allRunNumbersFound) outfile << Form("All run numbers appear in the tree.");

    outfile.close();

    return;
}