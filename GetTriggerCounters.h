// RunListCheck.c
// David Grund, Feb 26, 2022
// To extract the counters of CCUP31-triggered events per run for both periods

void GetTriggerCounters_main(){

    Int_t RunList18q[123] = {
        295585, 295586, 295588, 295589, 295610, 295611, 295612, 295615, 295666, 295667, 295668, 295673, 295675, 
        295676, 295712, 295714, 295717, 295718, 295719, 295721, 295723, 295725, 295754, 295755, 295758, 295759, 
        295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829, 295831, 
        295853, 295854, 295855, 295856, 295859, 295860, 295861, 295909, 295910, 295913, 295936, 295937, 295941, 
        295942, 296016, 296060, 296062, 296063, 296065, 296066, 296123, 296132, 296133, 296134, 296135, 296142, 
        296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 
        296246, 296247, 296269, 296270, 296273, 296279, 296280, 296303, 296304, 296309, 296312, 296377, 296378, 
        296379, 296380, 296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472, 296509, 
        296510, 296511, 296512, 296516, 296547, 296548, 296549, 296550, 296551, 296552, 296553, 296594, 296615, 
        296616, 296618, 296619, 296621, 296622, 296623
    };

    Int_t RunList18r[96] = {
        296690, 296691, 296693, 296694, 296749, 296750, 296781, 296784, 296785, 296786, 296787, 296790, 296793, 
        296794, 296799, 296835, 296836, 296838, 296839, 296848, 296849, 296850, 296851, 296852, 296890, 296894, 
        296899, 296900, 296903, 296930, 296931, 296932, 296934, 296935, 296938, 296941, 296966, 297029, 297031, 
        297035, 297085, 297117, 297118, 297119, 297123, 297124, 297128, 297129, 297132, 297133, 297193, 297194, 
        297195, 297196, 297218, 297219, 297221, 297222, 297278, 297310, 297311, 297317, 297332, 297333, 297335, 
        297336, 297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414, 297415, 297441, 
        297442, 297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540, 297541, 297542, 
        297544, 297558, 297588, 297590, 297595
    };

    Int_t Counts18q[123] = { 0 };
    Int_t Counts18r[96] = { 0 };

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
        if(iCurrentRun < 123){ // 2018q
            if(hCounterTrigger->GetBinCenter(iBin) == RunList18q[iCurrentRun]){
                Counts18q[iCurrentRun] = hCounterTrigger->GetBinContent(iBin);
                iCurrentRun++;
                counter += hCounterTrigger->GetBinContent(iBin);
                // Debugging:
                Printf("LHC18q, iRun: %i, Run no.: %.0f, Counter: %.0f", iCurrentRun, hCounterTrigger->GetBinCenter(iBin), hCounterTrigger->GetBinContent(iBin));
            }
        } else { // 2018r
            if(hCounterTrigger->GetBinCenter(iBin) == RunList18r[iCurrentRun-123]){
                Counts18r[iCurrentRun-123] = hCounterTrigger->GetBinContent(iBin);
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
    outfile << "Counts18q (123 runs):\n";
    Int_t RunsPerLine = 13;
    for(Int_t i = 0; i < 123; i++){
        if((i+1) % RunsPerLine == 0){
            outfile << Counts18q[i] << ", \n";
        } else if(i == 122){
            outfile << Counts18q[i];
        } else {
            outfile << Counts18q[i] << ", ";
        }
    }
    outfile << "\n\nCounts18r (96 runs):\n";
    for(Int_t i = 0; i < 96; i++){
        if((i+1) % RunsPerLine == 0){
            outfile << Counts18r[i] << ", \n";
        } else if(i == 95){
            outfile << Counts18r[i];
        } else {
            outfile << Counts18r[i] << ", ";
        }
    }
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    // Print the results separately for both periods (to be read by the lumi macro)
    // LHC18q
    TString name_18q = "Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18q.txt";
    ofstream out_f_18q (name_18q.Data());
    for(Int_t i = 0; i < 123; i++){
        out_f_18q << Counts18q[i] << "\n";
    }
    out_f_18q.close();
    Printf("*** Results printed to %s.***", name_18q.Data());

    // LHC18r
    TString name_18r = "Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18r.txt";
    ofstream out_f_18r (name_18r.Data());
    for(Int_t i = 0; i < 96; i++){
        out_f_18r << Counts18r[i] << "\n";
    }
    out_f_18r.close();
    Printf("*** Results printed to %s.***", name_18r.Data());

    return;
}