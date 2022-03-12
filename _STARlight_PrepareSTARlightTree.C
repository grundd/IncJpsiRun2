// _STARlight_PrepareSTARlightTree.C
// David Grund, Mar 07, 2022

// root headers
#include "TSystem.h"
// my headers
#include "_STARlight_Utilities.h"

void PrepareTrees_ChooseDataset(Int_t iDS);

void _STARlight_PrepareSTARlightTree()
{
    // R_A = 7.350 vs. 6.624 fm:
    // IncJ
    PrepareTrees_ChooseDataset(1);
    PrepareTrees_ChooseDataset(2);
    // CohP
    PrepareTrees_ChooseDataset(3);
    PrepareTrees_ChooseDataset(4);
    // IncP
    PrepareTrees_ChooseDataset(5);
    PrepareTrees_ChooseDataset(6);

    return;
}

void PrepareTrees_ChooseDataset(Int_t iDS){

    Int_t nGenEv = 0;
    TString folder_in = "";
    TString folder_out = "";

    if(iDS == 1){
        // Mar 07, 2022
        // IncJ, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncJ_6.624/";
        folder_out = "Trees/STARlight/IncJ_6.624/";
    } else if(iDS == 2){
        // Mar 07, 2022
        // IncJ, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncJ_7.350/";
        folder_out = "Trees/STARlight/IncJ_7.350/";
    } else if(iDS == 3){
        // Mar 12, 2022
        // CohP, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/CohP_6.624/";
        folder_out = "Trees/STARlight/CohP_6.624/";
    } else if(iDS == 4){
        // Mar 12, 2022
        // CohP, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/CohP_7.350/";
        folder_out = "Trees/STARlight/CohP_7.350/";
    } else if(iDS == 5){
        // Mar 12, 2022
        // IncP, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncP_6.624/";
        folder_out = "Trees/STARlight/IncP_6.624/";
    } else if(iDS == 6){
        // Mar 12, 2022
        // IncP, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncP_7.350/";
        folder_out = "Trees/STARlight/IncP_7.350/";
    } 

    gSystem->Exec("mkdir -p " + folder_out);

    ConvertStarlightAsciiToTree(nGenEv, folder_in, folder_out);

    return;
}