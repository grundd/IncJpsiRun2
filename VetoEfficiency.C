// VetoEfficiency.C
// David Grund, June 18, 2022

// root headers
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

// ******** options to set: **********
const Int_t nBinsN = 5;
const Int_t nBinsPt = 5;
Double_t fBkgM_low = 1.8; // GeV
Double_t fBkgM_upp = 2.8; // GeV
// arrays:
Double_t fNumberOfN[nBinsN+1] = {0.0, 1.5, 5.5, 10.5, 20.5, 50.5};
TString  sNumberOfN[nBinsN] = {"0-1", "2-5", "6-10", "11-20", "21-50"};
Double_t fVetoIneff_A[nBinsN] = {0.085, 0.157, 0.265, 0.413, 0.579};
Double_t fVetoIneff_C[nBinsN] = {0.146, 0.333, 0.425, 0.542, 0.844};
// ***********************************
// tree variables:
Bool_t fZNA_hit, fZNC_hit;
Double_t fZNA_n, fZNC_n;
// **********************************************************************
// numbers of background events
// total numbers in neutron classes
Int_t nBkg_0n0n(0), nBkg_Xn0n(0), nBkg_0nXn(0), nBkg_XnXn(0);
// per neutron class and pT bin
Int_t nBkg_0n0n_BinsPt[nBinsPt] = { 0 }; // 0n0n class, index = pT bin
Int_t nBkg_Xn0n_BinsPt[nBinsPt] = { 0 }; // Xn0n class, index = pT bin
Int_t nBkg_0nXn_BinsPt[nBinsPt] = { 0 }; // 0nXn class, index = pT bin
Int_t nBkg_XnXn_BinsPt[nBinsPt] = { 0 }; // XnXn class, index = pT bin
// per neutron class and neutron bin
Int_t nBkg_Xn0n_BinsN[nBinsN] = { 0 };   // Xn0n class, index = neutron bin (A)
Int_t nBkg_0nXn_BinsN[nBinsN] = { 0 };   // 0nXn class, index = neutron bin (C)
Int_t nBkg_XnXn_BinsN[nBinsN][nBinsN] = { 0 }; // XnXn class, first index = neutron bin (A), second index = neutron bin (C)
// per neutron class, pT bin and neutron bin
Int_t nBkg_Xn0n_BinsPtN[nBinsPt][nBinsN] = { 0 }; // Xn0n class, first index = pT bin, second index = neutron bin (A)
Int_t nBkg_0nXn_BinsPtN[nBinsPt][nBinsN] = { 0 }; // 0nXn class, first index = pT bin, second index = neutron bin (C)
Int_t nBkg_XnXn_BinsPtN[nBinsPt][nBinsN][nBinsN] = { 0 }; // XnXn class, first index = pT bin, second index = neutron bin (A), third index = neutron bin (C)
// **********************************************************************
// fraction of background events in each channel
// total fractions per neutron classes
Double_t fBkg_0n0n(0.), fBkg_Xn0n(0.), fBkg_0nXn(0.), fBkg_XnXn(0.);
// fractions per neutron class and pT bin
Double_t fBkg_0n0n_BinsPt[nBinsPt] = { 0. }; // 0n0n class, index = pT bin
Double_t fBkg_Xn0n_BinsPt[nBinsPt] = { 0. }; // Xn0n class, index = pT bin
Double_t fBkg_0nXn_BinsPt[nBinsPt] = { 0. }; // 0nXn class, index = pT bin
Double_t fBkg_XnXn_BinsPt[nBinsPt] = { 0. }; // XnXn class, index = pT bin
// fractions per neutron class and neutron bin
Double_t fBkg_Xn0n_BinsN[nBinsN] = { 0. };   // Xn0n class, index = neutron bin (A)
Double_t fBkg_0nXn_BinsN[nBinsN] = { 0. };   // 0nXn class, index = neutron bin (C)
Double_t fBkg_XnXn_BinsN[nBinsN][nBinsN] = { 0. }; // XnXn class, first index = neutron bin (A), second index = neutron bin (C)
// fractions per neutron class, pT bin and neutron bin
Double_t fBkg_Xn0n_BinsPtN[nBinsPt][nBinsN] = { 0. }; // Xn0n class, first index = pT bin, second index = neutron bin (A)
Double_t fBkg_0nXn_BinsPtN[nBinsPt][nBinsN] = { 0. }; // 0nXn class, first index = pT bin, second index = neutron bin (C)
Double_t fBkg_XnXn_BinsPtN[nBinsPt][nBinsN][nBinsN] = { 0. }; // XnXn class, first index = pT bin, second index = neutron bin (A), third index = neutron bin (C)
// **********************************************************************

void VetoEfficiency_ClassifyEventsToChannels(Int_t opt);
void VetoEfficiency_PrepareTree();

void VetoEfficiency(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VetoEfficiency/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/");

    // prepare the tree containing information about mass, pT and ZN signal
    VetoEfficiency_PrepareTree();

    return;
}

void VetoEfficiency_ClassifyEventsToChannels(Int_t opt)
// opt == 0 => background (mass range: fBkgM_low to fBkgM_upp; pT bins)
// opt == 1 =>
{


    return;
}

void VetoEfficiency_PrepareTree()
{
    TString name = "Trees/" + str_subfolder + "VetoEfficiency/tNeutrons.root";

    TFile *file = TFile::Open(name.Data(),"read");
    if(file)
    {
        Printf("Tree already created.");
        return;
    } 
    else 
    {
        Printf("Tree will be created.");

        // data
        TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
        if(f_in) Printf("Input data loaded.");

        TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
        if(t_in) Printf("Input tree loaded.");

        ConnectTreeVariables(t_in);

        // Create new data tree 
        file = new TFile(name.Data(),"RECREATE");

        TTree *t_out = new TTree("tNeutrons","tNeutrons");
        t_out->Branch("fPt", &fPt, "fPt/D");
        t_out->Branch("fM", &fM, "fM/D");
        t_out->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
        t_out->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
        t_out->Branch("fZNA_hit", &fZNA_hit, "fZNA_hit/O");
        t_out->Branch("fZNC_hit", &fZNC_hit, "fZNC_hit/O");
        t_out->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
        t_out->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
        t_out->Branch("fZNA_n", &fZNA_n, "fZNA_n/D");
        t_out->Branch("fZNC_n", &fZNC_n, "fZNC_n/D");

        Printf("%lli entries found in the tree.", t_in->GetEntries());
        Int_t nEntriesAnalysed = 0;
    
        for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
        {
            t_in->GetEntry(iEntry);
            // no mass cut, pT in 0.2 to 1.0 GeV/c, then mass between 1.6 GeV and 3.2 GeV
            if(EventPassed(-1, 3) && fM > 1.6 && fM < 3.2)
            {
                fZNA_hit = kFALSE;
                fZNC_hit = kFALSE;
                for(Int_t i = 0; i < 4; i++)
                {
                    // hit in ZNA
                    if(TMath::Abs(fZNA_time[i]) < 2) fZNA_hit = kTRUE;
                    // hit in ZNC
                    if(TMath::Abs(fZNC_time[i]) < 2) fZNC_hit = kTRUE;
                }
                fZNA_n = fZNA_energy / 2510.;
                fZNC_n = fZNC_energy / 2510.;
                t_out->Fill();
            }

            if((iEntry+1) % 100000 == 0)
            {
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        Printf("Tree %s filled with %lli entries.", t_out->GetName(), t_out->GetEntries());

        file->Write("",TObject::kWriteDelete);

        return;
    }
}