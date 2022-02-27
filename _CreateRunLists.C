// RunListCheck.C
// David Grund, Feb 27, 2022
// To check if the run numbers of analysed data match the official list

#include "TFile.h"
#include <fstream> // print output to txt file
#include <vector>

const Int_t nRuns_18q_pass1_GoodHadronPID = 126;
const Int_t nRuns_18r_pass1_GoodHadronPID = 90;
const Int_t nRuns_18r_pass1_TimeRangeCuts = 7;
const Int_t nRuns_18q_pass3_GoodHadronPID = 125;
const Int_t nRuns_18r_pass3_GoodHadronPID = 89;
const Int_t nRuns_18r_pass3_TimeRangeCuts = 8;
const Int_t nRuns_18q_ExcludedRuns = 3;
const Int_t nRuns_18r_ExcludedRuns = 1;

ofstream *outfile = NULL;

Int_t RunListDPG_18q_pass1_GoodHadronPID[nRuns_18q_pass1_GoodHadronPID] = {
    // https://twiki.cern.ch/twiki/pub/ALICE/AliDPGRunList18q/RunList_LHC18q_pass1_CentralBarrelTracking_hadronPID.txt
    296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551, 296550, 296549, 
    296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420, 296419, 
    296415, 296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304, 
    296303, 296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241, 296240, 
    296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142, 296135, 296134, 296133, 296132, 
    296123, 296074, 296066, 296065, 296063, 296062, 296060, 296016, 295942, 295941, 295937, 295936, 295913, 
    295910, 295909, 295861, 295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825, 
    295822, 295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 295759, 295758, 295755, 295754, 
    295725, 295723, 295721, 295719, 295718, 295717, 295714, 295712, 295676, 295675, 295673, 295668, 295667, 
    295666, 295615, 295612, 295611, 295610, 295589, 295588, 295586, 295585
};

Int_t RunListDPG_18r_pass1_GoodHadronPID[nRuns_18r_pass1_GoodHadronPID] = {
    // https://twiki.cern.ch/twiki/pub/ALICE/AliDPGRunList18r/RunList_LHC18r_pass1_CentralBarrelTracking_hadronPID_7.txt
    297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483, 297481, 297479, 
    297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380, 297379, 
    297372, 297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297311, 297310, 297278, 297222, 
    297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 
    297117, 297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934, 296932, 296931, 296930, 296903, 
    296900, 296899, 296894, 296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794, 
    296793, 296790, 296787, 296786, 296785, 296784, 296781, 296752, 296694, 296693, 296691, 296690
};

Int_t RunListDPG_18r_pass1_TimeRangeCuts[nRuns_18r_pass1_TimeRangeCuts] = {
    // https://twiki.cern.ch/twiki/pub/ALICE/AliDPGRunList18r/RunList_LHC18r_pass1_CentralBarrelTracking_hadronPID_7.txt
    // these can be used in pass1 only with TimeRangeCuts
    297219, 297194, 297029, 296890, 296849, 296750, 296749
};

Int_t RunListDPG_18q_pass3_GoodHadronPID[nRuns_18q_pass3_GoodHadronPID] = {
    // https://twiki.cern.ch/twiki/pub/ALICE/AliDPGRunList18q/RunList_LHC18q_pass3_CentralBarrelTracking_hadronPID.txt 
    // 296549 exluded w.r.t. to pass1
    296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551, 296550, 296548, 
    296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420, 296419, 296415, 
    296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304, 296303, 
    296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241, 296240, 296198, 
    296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142, 296135, 296134, 296133, 296132, 296123, 
    296074, 296066, 296065, 296063, 296062, 296060, 296016, 295942, 295941, 295937, 295936, 295913, 295910, 
    295909, 295861, 295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825, 295822, 
    295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 295759, 295758, 295755, 295754, 295725, 
    295723, 295721, 295719, 295718, 295717, 295714, 295712, 295676, 295675, 295673, 295668, 295667, 295666, 
    295615, 295612, 295611, 295610, 295589, 295588, 295586, 295585
};

Int_t RunListDPG_18r_pass3_GoodHadronPID[nRuns_18r_pass3_GoodHadronPID] = {
    // https://twiki.cern.ch/twiki/pub/ALICE/AliDPGRunList18r/Runlist_LHC18r_pass3_CentralBarrelTracking_hadronPID.txt
    // 297481 exluded w.r.t. to pass1 (in pass3 moved to the list requiring TimeRangeCuts)
    297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483, 297479, 297452, 
    297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380, 297379, 297372, 
    297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297311, 297310, 297278, 297222, 297221, 
    297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 297117, 
    297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 
    296899, 296894, 296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794, 296793, 
    296790, 296787, 296786, 296785, 296784, 296781, 296752, 296694, 296693, 296691, 296690

};

Int_t RunListDPG_18r_pass3_TimeRangeCuts[nRuns_18r_pass3_TimeRangeCuts] = {
    // https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGRunList18r
    296749, 296750, 296849, 296890, 297029, 297194, 297219, 297481
};

Int_t RunList_18q_ExcludedRuns[nRuns_18q_ExcludedRuns] = {
    296074, // LHC18q, no CCUP31 events
    296375, // LHC18q, AD set to collect cosmics
    296376  // LHC18q, AD set to collect cosmics
};

Int_t RunList_18r_ExcludedRuns[nRuns_18r_ExcludedRuns] = {
    296752  // LHC18r, no CCUP31 events
};

void CreateAndSortRunLists();
void PrintRunList(ofstream &ofs, vector<Int_t>runList);

void _CreateRunLists(){

    CreateAndSortRunLists();

    // note: for pass1, there are none events in LHC18r with run no. 297219 because all jobs for this run failed in a LEGO train

    return;
}

void CreateAndSortRunLists(){
    // create total run lists for both periods (separately)
    // sort the run numbers in an increasing order

    // check if the lengths of arrays match the anticipated numbers of runs
    Printf("Checking the lenghts of arrays with run numbers:");

    if(sizeof(RunListDPG_18q_pass1_GoodHadronPID) / sizeof(RunListDPG_18q_pass1_GoodHadronPID[0]) == nRuns_18q_pass1_GoodHadronPID) Printf("LHC18q_pass1_GoodHadronPID: OK.");
    else Printf("LHC18q_pass1_GoodHadronPID: Dimensions don't match!.");

    if(sizeof(RunListDPG_18r_pass1_GoodHadronPID) / sizeof(RunListDPG_18r_pass1_GoodHadronPID[0]) == nRuns_18r_pass1_GoodHadronPID) Printf("LHC18r_pass1_GoodHadronPID: OK.");
    else Printf("LHC18r_pass1_GoodHadronPID: Dimensions don't match!.");

    if(sizeof(RunListDPG_18r_pass1_TimeRangeCuts) / sizeof(RunListDPG_18r_pass1_TimeRangeCuts[0]) == nRuns_18r_pass1_TimeRangeCuts) Printf("LHC18r_pass1_TimeRangeCuts: OK.");
    else Printf("LHC18r_pass1_TimeRangeCuts: Dimensions don't match!.");

    if(sizeof(RunListDPG_18q_pass3_GoodHadronPID) / sizeof(RunListDPG_18q_pass3_GoodHadronPID[0]) == nRuns_18q_pass3_GoodHadronPID) Printf("LHC18q_pass3_GoodHadronPID: OK.");
    else Printf("LHC18q_pass3_GoodHadronPID: Dimensions don't match!.");

    if(sizeof(RunListDPG_18r_pass3_GoodHadronPID) / sizeof(RunListDPG_18r_pass3_GoodHadronPID[0]) == nRuns_18r_pass3_GoodHadronPID) Printf("LHC18r_pass3_GoodHadronPID: OK.");
    else Printf("LHC18r_pass3_GoodHadronPID: Dimensions don't match!.");

    if(sizeof(RunListDPG_18r_pass3_TimeRangeCuts) / sizeof(RunListDPG_18r_pass3_TimeRangeCuts[0]) == nRuns_18r_pass3_TimeRangeCuts) Printf("LHC18r_pass3_TimeRangeCuts: OK.");
    else Printf("LHC18r_pass3_TimeRangeCuts: Dimensions don't match!.");

    // create new arrays 
    // https://stackoverflow.com/questions/28625465/c-creating-an-array-with-a-size-entered-by-the-user
    vector<Int_t> TotalRunList_18q_pass1;
    vector<Int_t> TotalRunList_18r_pass1;
    vector<Int_t> TotalRunList_18q_pass3;
    vector<Int_t> TotalRunList_18r_pass3;

    // fill the new arrays
    // LHC18q pass1
    for(Int_t i = 0; i < nRuns_18q_pass1_GoodHadronPID; i++){
        TotalRunList_18q_pass1.push_back(RunListDPG_18q_pass1_GoodHadronPID[i]);
    } 
    // LHC18r pass1
    for(Int_t i = 0; i < nRuns_18r_pass1_GoodHadronPID; i++){
        TotalRunList_18r_pass1.push_back(RunListDPG_18r_pass1_GoodHadronPID[i]);
    } 
    for(Int_t i = nRuns_18r_pass1_GoodHadronPID; i < nRuns_18r_pass1_GoodHadronPID + nRuns_18r_pass1_TimeRangeCuts; i++){
        TotalRunList_18r_pass1.push_back(RunListDPG_18r_pass1_TimeRangeCuts[i-nRuns_18r_pass1_GoodHadronPID]);
    }
    // LHC18q pass3
    for(Int_t i = 0; i < nRuns_18q_pass3_GoodHadronPID; i++){
        TotalRunList_18q_pass3.push_back(RunListDPG_18q_pass3_GoodHadronPID[i]);
    } 
    // LHC18r pass3
    for(Int_t i = 0; i < nRuns_18r_pass3_GoodHadronPID; i++){
        TotalRunList_18r_pass3.push_back(RunListDPG_18r_pass3_GoodHadronPID[i]);
    } 
    for(Int_t i = nRuns_18r_pass3_GoodHadronPID; i < nRuns_18r_pass3_GoodHadronPID + nRuns_18r_pass3_TimeRangeCuts; i++){
        TotalRunList_18r_pass3.push_back(RunListDPG_18r_pass3_TimeRangeCuts[i-nRuns_18r_pass3_GoodHadronPID]);
    }

    // sort the new arrays
    // sort the new arrays
    std::sort(TotalRunList_18q_pass1.begin(), TotalRunList_18q_pass1.end());
    std::sort(TotalRunList_18r_pass1.begin(), TotalRunList_18r_pass1.end());
    std::sort(TotalRunList_18q_pass3.begin(), TotalRunList_18q_pass3.end());
    std::sort(TotalRunList_18r_pass3.begin(), TotalRunList_18r_pass3.end());

    // create new arrays if some runs are excluded
    vector<Int_t> TotalRunList_18q_pass1_exc;
    vector<Int_t> TotalRunList_18r_pass1_exc;
    vector<Int_t> TotalRunList_18q_pass3_exc;
    vector<Int_t> TotalRunList_18r_pass3_exc;

    // fill the new arrays
    // LHC18q pass1
    for(Int_t i = 0; i < (Int_t)TotalRunList_18q_pass1.size(); i++){
        Bool_t RunExcluded = std::find(std::begin(RunList_18q_ExcludedRuns), std::end(RunList_18q_ExcludedRuns), TotalRunList_18q_pass1[i]) != std::end(RunList_18q_ExcludedRuns);
        if(!RunExcluded) TotalRunList_18q_pass1_exc.push_back(TotalRunList_18q_pass1[i]);
    } 
    // LHC18r pass1
    for(Int_t i = 0; i < (Int_t)TotalRunList_18r_pass1.size(); i++){
        Bool_t RunExcluded = std::find(std::begin(RunList_18r_ExcludedRuns), std::end(RunList_18r_ExcludedRuns), TotalRunList_18r_pass1[i]) != std::end(RunList_18r_ExcludedRuns);
        if(!RunExcluded) TotalRunList_18r_pass1_exc.push_back(TotalRunList_18r_pass1[i]);
    }
    // LHC18q pass3
    for(Int_t i = 0; i < (Int_t)TotalRunList_18q_pass3.size(); i++){
        Bool_t RunExcluded = std::find(std::begin(RunList_18q_ExcludedRuns), std::end(RunList_18q_ExcludedRuns), TotalRunList_18q_pass3[i]) != std::end(RunList_18q_ExcludedRuns);
        if(!RunExcluded) TotalRunList_18q_pass3_exc.push_back(TotalRunList_18q_pass3[i]);
    } 
    // LHC18r pass3
    for(Int_t i = 0; i < (Int_t)TotalRunList_18r_pass3.size(); i++){
        Bool_t RunExcluded = std::find(std::begin(RunList_18r_ExcludedRuns), std::end(RunList_18r_ExcludedRuns), TotalRunList_18r_pass3[i]) != std::end(RunList_18r_ExcludedRuns);
        if(!RunExcluded) TotalRunList_18r_pass3_exc.push_back(TotalRunList_18r_pass3[i]);
    }

    // sort the new arrays
    std::sort(TotalRunList_18q_pass1_exc.begin(), TotalRunList_18q_pass1_exc.end());
    std::sort(TotalRunList_18r_pass1_exc.begin(), TotalRunList_18r_pass1_exc.end());
    std::sort(TotalRunList_18q_pass3_exc.begin(), TotalRunList_18q_pass3_exc.end());
    std::sort(TotalRunList_18r_pass3_exc.begin(), TotalRunList_18r_pass3_exc.end());

    // *************************************************************************
    // pass1: print the results:
    TString name_pass1 = "Results/_CreateRunLists/RunLists_pass1.txt";
    ofstream outfile_pass1 (name_pass1.Data());

    // LHC18q
    outfile_pass1 << "LHC18q_pass1 (" << TotalRunList_18q_pass1.size() << " runs):\n";   
    PrintRunList(outfile_pass1, TotalRunList_18q_pass1);

    // LHC18r
    outfile_pass1 << "LHC18r_pass1 (" << TotalRunList_18r_pass1.size() << " runs):\n"; 
    PrintRunList(outfile_pass1, TotalRunList_18r_pass1);

    // LHC18q, excluded runs
    outfile_pass1 << "LHC18q_pass1, excluded runs (" << TotalRunList_18q_pass1_exc.size() << " runs):\n";   
    PrintRunList(outfile_pass1, TotalRunList_18q_pass1_exc);

    // LHC18r, excluded runs
    outfile_pass1 << "LHC18r_pass1, excluded runs (" << TotalRunList_18r_pass1_exc.size() << " runs):\n"; 
    PrintRunList(outfile_pass1, TotalRunList_18r_pass1_exc);

    outfile_pass1.close();
    Printf("*** Results printed to %s.***", name_pass1.Data());

    // *************************************************************************
    // pass3: print the results:
    TString name_pass3 = "Results/_CreateRunLists/RunLists_pass3.txt";
    ofstream outfile_pass3 (name_pass3.Data());

    // LHC18q
    outfile_pass3 << "LHC18q_pass3 (" << TotalRunList_18q_pass3.size() << " runs):\n";   
    PrintRunList(outfile_pass3, TotalRunList_18q_pass3);

    // LHC18r
    outfile_pass3 << "LHC18r_pass3 (" << TotalRunList_18r_pass3.size() << " runs):\n"; 
    PrintRunList(outfile_pass3, TotalRunList_18r_pass3);

    // LHC18q, excluded runs
    outfile_pass3 << "LHC18q_pass3, excluded runs (" << TotalRunList_18q_pass3_exc.size() << " runs):\n";   
    PrintRunList(outfile_pass3, TotalRunList_18q_pass3_exc);

    // LHC18r, excluded runs
    outfile_pass3 << "LHC18r_pass3, excluded runs (" << TotalRunList_18r_pass3_exc.size() << " runs):\n"; 
    PrintRunList(outfile_pass3, TotalRunList_18r_pass3_exc);

    outfile_pass3.close();
    Printf("*** Results printed to %s.***", name_pass3.Data());

    return;
}

void PrintRunList(ofstream &ofs, vector<Int_t>runList){

    Int_t RunsPerLine = 12;

    for(Int_t i = 0; i < (Int_t)runList.size(); i++){
        if(i == (Int_t)runList.size() - 1){
            ofs << runList[i] << "\n\n";
        } else if((i+1) % RunsPerLine == 0){
            ofs << runList[i] << ", \n";
        } else {
            ofs << runList[i] << ", ";
        }
    }

    return;
}