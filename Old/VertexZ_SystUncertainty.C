// VertexZ_SystUncertainty.C
// David Grund, Jun 07, 2022
// we suppose that values of the parameters fD and fC do not change significantly

// cpp headers
#include <fstream>
// root headers
#include "TFile.h"
#include "TString.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "InvMassFit_MC_Utilities.h"
#include "InvMassFit_Utilities.h"
#include "AxE_PtBins_Utilities.h"

void NewCutZ_InvMassFits(Double_t fCutZ);
void NewCutZ_AxE_PtBins(Double_t fCutZ);
void NewCutZ_PrepareTree_InvMassFit_MC(Double_t fCutZ);
void NewCutZ_PrepareTree_InvMassFit(Double_t fCutZ);

void VertexZ_SystUncertainty(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VertexZ_SystUncertainty/");

    // calculate new yields of J/psi in the same pT bins
    NewCutZ_InvMassFits(12.5);

    // calculate new AxE in the pT bins
    NewCutZ_AxE_PtBins(12.5);

    return;
}

void NewCutZ_InvMassFits(Double_t fCutZ)
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit_MC/", fCutZ));
    gSystem->Exec("mkdir -p Results/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit/", fCutZ));
    // prepare the tree for new inv mass fits
    NewCutZ_PrepareTree_InvMassFit_MC(fCutZ);
    NewCutZ_PrepareTree_InvMassFit(fCutZ);
    // go over pT bins and do new fits
    for(Int_t i = 0; i < nPtBins; i++)
    {
        // mc
        TString str_mc_new = "Results/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit_MC/bin%i", fCutZ, i+1);
        //InvMassFit_MC_DoFit(4+i, str_mc_new, kTRUE, fCutZ);
        // data
        // CB tail parameters: load MC values
        Double_t fAlpha_L;
        Double_t fAlpha_R;
        Double_t fN_L;
        Double_t fN_R;
        // open the file and load the values
        char name[20];
        Double_t values[4];
        Double_t errors[4];

        Bool_t KeepOldMCValues = kTRUE;
        TString str_mc_old = "Results/" + str_subfolder + Form("InvMassFit_MC/%ibins/bin%i", nPtBins, i+1);
        ifstream ifs;
        if(KeepOldMCValues) ifs.open((str_mc_old + ".txt").Data());
        else                ifs.open((str_mc_new + ".txt").Data());
        if(ifs.fail()){
            Printf("\n");
            Printf("*** Warning! ***");
            Printf("*** MC values for tail parameters not found. Terminating... *** \n");
            return;
        } else {
            Int_t i_line = 0;
            while(!ifs.eof()){
                ifs >> name >> values[i_line] >> errors[i_line];
                i_line++;
            }
        }
        ifs.close();
        // set the variables to loaded values
        fAlpha_L = values[0];
        fAlpha_R = values[1];
        fN_L = values[2];
        fN_R = values[3];
        // perform the fit
        TString str_out2 = "Results/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit/bin%i", fCutZ, i+1);
        InvMassFit_DoFit(4+i, 2.2, 4.5, fAlpha_L, fAlpha_R, fN_L, fN_R, str_out2, kTRUE, fCutZ);
    }

    return;
}

void NewCutZ_AxE_PtBins(Double_t fCutZ)
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_AxE_PtBins/", fCutZ));
    AxE_PtBins_Calculate(fCutZ);

    return;
}

void NewCutZ_PrepareTree_InvMassFit_MC(Double_t fCutZ)
{
    TString name = "Trees/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit_MC.root", fCutZ);

    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("Tree for InvMassFit_MC already created.");
        return;

    } else { 

        Printf("Tree for InvMassFit_MC will be created.");

        // kIncohJpsiToMu
        TFile *f_in = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(f_in) Printf("Input data loaded.");

        TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_MC_tree_rec.Data()));
        if(t_in) Printf("Input tree loaded.");

        ConnectTreeVariablesMCRec(t_in);

        // Create new data tree with applied cuts
        file = new TFile(name.Data(),"RECREATE");

        TTree *tIncEnrSample = new TTree("tIncEnrSample", "tIncEnrSample");
        tIncEnrSample->Branch("fPt", &fPt, "fPt/D");
        tIncEnrSample->Branch("fM", &fM, "fM/D");
        tIncEnrSample->Branch("fY", &fY, "fY/D");

        Printf("%lli entries found in the tree.", t_in->GetEntries());
        Int_t nEntriesAnalysed = 0;

        // save the original value of cut_fVertexZ
        Printf("Original cut on vertex Z: %.1f", cut_fVertexZ);
        Double_t fCutZ_orig = cut_fVertexZ;
        // set the new value of cut_fVertexZ
        cut_fVertexZ = fCutZ;
        Printf("New cut on vertex Z: %.1f", cut_fVertexZ);

        for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
        {
            t_in->GetEntry(iEntry);
            // inv mass cut: 2.2 < m < 4.5, pT cut: all
            if(EventPassed(0, 2)) tIncEnrSample->Fill();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // set back the original value of cut_fVertexZ
        cut_fVertexZ = fCutZ_orig;
        Printf("Restoring the original cut on vertex Z: %.1f", cut_fVertexZ);        

        file->Write("",TObject::kWriteDelete);

        return;
    }
}

void NewCutZ_PrepareTree_InvMassFit(Double_t fCutZ)
{
    TString name = "Trees/" + str_subfolder + Form("VertexZ_SystUncertainty/Zcut%.1f_InvMassFit.root", fCutZ);

    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("Tree for InvMassFit already created.");
        return;

    } else { 

        Printf("Tree for InvMassFit will be created.");

        // data
        TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
        if(f_in) Printf("Input data loaded.");

        TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
        if(t_in) Printf("Input tree loaded.");

        ConnectTreeVariables(t_in);  

        // Create new data tree with applied cuts
        file = new TFile(name.Data(),"RECREATE");

        TTree *tIncEnrSample = new TTree("tIncEnrSample", "tIncEnrSample");
        tIncEnrSample->Branch("fPt", &fPt, "fPt/D");
        tIncEnrSample->Branch("fM", &fM, "fM/D");
        tIncEnrSample->Branch("fY", &fY, "fY/D");

        Printf("%lli entries found in the tree.", t_in->GetEntries());
        Int_t nEntriesAnalysed = 0;

        // save the original value of cut_fVertexZ
        Printf("Original cut on vertex Z: %.1f", cut_fVertexZ);
        Double_t fCutZ_orig = cut_fVertexZ;
        // set the new value of cut_fVertexZ
        cut_fVertexZ = fCutZ;
        Printf("New cut on vertex Z: %.1f", cut_fVertexZ);

        for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
        {
            t_in->GetEntry(iEntry);
            // inv mass cut: 2.2 < m < 4.5, pT cut: inc (pT > 0.2)
            if(EventPassed(0, 0)) tIncEnrSample->Fill();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // set back the original value of cut_fVertexZ
        cut_fVertexZ = fCutZ_orig;
        Printf("Restoring the original cut on vertex Z: %.1f", cut_fVertexZ);        

        file->Write("",TObject::kWriteDelete);

        return;
    }
}