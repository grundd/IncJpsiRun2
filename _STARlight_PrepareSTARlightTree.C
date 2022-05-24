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
    // CohJ
    PrepareTrees_ChooseDataset(10);
    PrepareTrees_ChooseDataset(11);
    // IncJ
    PrepareTrees_ChooseDataset(20);
    PrepareTrees_ChooseDataset(21);
    // CohP
    PrepareTrees_ChooseDataset(30);
    PrepareTrees_ChooseDataset(31);
    // IncP
    PrepareTrees_ChooseDataset(40);
    PrepareTrees_ChooseDataset(41);

    // IncJ: |t| vs. pT of the J_psi
    PrepareTrees_ChooseDataset(68);

    return;
}

void PrepareTrees_ChooseDataset(Int_t iDS){

    Int_t nGenEv = 0;
    TString folder_in = "";
    TString folder_out = "";


    if(iDS == 10){
        // Mar 13, 2022
        // CohJ, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/CohJ_6.624/";
        folder_out = "Trees/STARlight/CohJ_6.624/";
    } else if(iDS == 11){
        // Mar 13, 2022
        // CohJ, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/CohJ_7.350/";
        folder_out = "Trees/STARlight/CohJ_7.350/";
    } else if(iDS == 20){
        // Mar 07, 2022
        // IncJ, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncJ_6.624/";
        folder_out = "Trees/STARlight/IncJ_6.624/";
    } else if(iDS == 21){
        // Mar 07, 2022
        // IncJ, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncJ_7.350/";
        folder_out = "Trees/STARlight/IncJ_7.350/";
    } else if(iDS == 30){
        // Mar 12, 2022
        // CohP, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/CohP_6.624/";
        folder_out = "Trees/STARlight/CohP_6.624/";
    } else if(iDS == 31){
        // Mar 12, 2022
        // CohP, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/CohP_7.350/";
        folder_out = "Trees/STARlight/CohP_7.350/";
    } else if(iDS == 40){
        // Mar 12, 2022
        // IncP, R_A set to 6.624 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncP_6.624/";
        folder_out = "Trees/STARlight/IncP_6.624/";
    } else if(iDS == 41){
        // Mar 12, 2022
        // IncP, R_A set to 7.350 fm in src/nucleus.cpp 
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncP_7.350/";
        folder_out = "Trees/STARlight/IncP_7.350/";
    } else if(iDS == 68){
        // May 24, 2022
        // IncJ, R_A = 6.624 fm, PtGammaVMPom.txt produced
        nGenEv = 6000000;
        folder_in = "STARlight_src/installation/IncJ_tVsPt/";
        folder_out = "Trees/STARlight/IncJ_tVsPt/";
        PrepareTreesPtGammaVMPom(nGenEv, folder_in, folder_out);
    } 

    gSystem->Exec("mkdir -p " + folder_out);

    ConvertStarlightAsciiToTree(nGenEv, folder_in, folder_out);

    return;
}