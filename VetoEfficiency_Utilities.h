// VetoEfficiency_Utilities.h
// David Grund, June 26, 2022

// cpp headers
#include <stdio.h> // printf
#include <iostream> // cout, cin
// root headers
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TRandom3.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"

// tree variables:
Bool_t fZNA_hit, fZNC_hit;
Double_t fZNA_n, fZNC_n;

// neutron bins:
const Int_t nBinsN = 5; 
Double_t fNumberOfN[nBinsN+1] = {0.0, 1.5, 5.5, 10.5, 20.5, 50.5};
TString  sNumberOfN[nBinsN+1] = {"none","0.0,1.5","1.5,5.5","5.5,10.5","10.5,20.5","20.5,50.5"};
Double_t nEv_SelAD_A[nBinsN] = {4, 11, 13, 19, 55};
Double_t nEv_SelAD_C[nBinsN] = {6, 23, 17, 26, 76};
Double_t nEv_Sel_A[nBinsN] = {47, 70, 49, 46, 95};
Double_t nEv_Sel_C[nBinsN] = {41, 69, 40, 48, 90};
Double_t fInf_A_val[nBinsN+1] = { 0.0 };
Double_t fInf_A_err[nBinsN+1] = { 0.0 };
Double_t fInf_C_val[nBinsN+1] = { 0.0 };
Double_t fInf_C_err[nBinsN+1] = { 0.0 };
Double_t fEff_A_val[nBinsN+1] = { 0 };
Double_t fEff_A_err[nBinsN+1] = { 0 };
Double_t fEff_C_val[nBinsN+1] = { 0 };
Double_t fEff_C_err[nBinsN+1] = { 0 };
Double_t fInf_val[nBinsN+1] = { 0 };
Double_t fInf_err[nBinsN+1] = { 0 };
Double_t fEff_val[nBinsN+1] = { 0 };
Double_t fEff_err[nBinsN+1] = { 0 };
Double_t SampledEff_A[nBinsN+1] = { 0 };
Double_t SampledEff_C[nBinsN+1] = { 0 };
Double_t SampledEff[nBinsN+1] = { 0 };

TH1D *hSampledEffPartial_A[nBinsN] = { NULL }; 
TH1D *hSampledEffPartial_C[nBinsN] = { NULL }; 
TH1D *hSampledEffPartial[nBinsN] = { NULL }; 
TH1D *hSampledEffTotal = new TH1D("hSampledEffTotal","hSampledEffTotal",100,0.,1.);

Double_t CalculateErrorBinomial(Double_t k, Double_t n)
{
    Double_t var = k * (n - k) / n / n / n;
    Double_t err = TMath::Sqrt(var);

    return err;
}

void VetoEff_CalcEfficiencies()
{
    // calculate
    // we leave first bins to be zero
    fInf_A_val[0] = 0.;
    fInf_A_err[0] = 0.;
    fInf_C_val[0] = 0.;
    fInf_C_err[0] = 0.;
    fEff_A_val[0] = 1.;
    fEff_A_err[0] = 0.;
    fEff_C_val[0] = 1.;
    fEff_C_err[0] = 0.;
    fInf_val[0] = 0.;
    fInf_err[0] = 0.;
    fEff_val[0] = 1.;
    fEff_err[0] = 0.;
    for(Int_t iBinN = 0; iBinN < nBinsN; iBinN++)
    {
        fInf_A_val[iBinN+1] = nEv_SelAD_A[iBinN] / nEv_Sel_A[iBinN];
        fInf_A_err[iBinN+1] = CalculateErrorBinomial(nEv_SelAD_A[iBinN], nEv_Sel_A[iBinN]); 
        fInf_C_val[iBinN+1] = nEv_SelAD_C[iBinN] / nEv_Sel_C[iBinN];
        fInf_C_err[iBinN+1] = CalculateErrorBinomial(nEv_SelAD_C[iBinN], nEv_Sel_C[iBinN]); 
        fEff_A_val[iBinN+1] = (nEv_Sel_A[iBinN] - nEv_SelAD_A[iBinN]) / nEv_Sel_A[iBinN];
        fEff_A_err[iBinN+1] = fInf_A_err[iBinN+1];
        fEff_C_val[iBinN+1] = (nEv_Sel_C[iBinN] - nEv_SelAD_C[iBinN]) / nEv_Sel_C[iBinN];
        fEff_C_err[iBinN+1] = fInf_C_err[iBinN+1];
        fInf_val[iBinN+1] = (nEv_SelAD_A[iBinN] + nEv_SelAD_C[iBinN]) / (nEv_Sel_A[iBinN] + nEv_Sel_C[iBinN]);
        fInf_err[iBinN+1] = CalculateErrorBinomial(nEv_SelAD_A[iBinN] + nEv_SelAD_C[iBinN], nEv_Sel_A[iBinN] + nEv_Sel_C[iBinN]); 
        fEff_val[iBinN+1] = (nEv_Sel_A[iBinN] + nEv_Sel_C[iBinN] - nEv_SelAD_A[iBinN] - nEv_SelAD_C[iBinN]) / (nEv_Sel_A[iBinN] + nEv_Sel_C[iBinN]);
        fEff_err[iBinN+1] = fInf_err[iBinN+1];
    }
    // print the results
    ofstream outfile;
    TString str_out = "Results/" + str_subfolder + "VetoEfficiency/efficiencies.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(3);
    for(Int_t iBinN = 0; iBinN < nBinsN+1; iBinN++)
    {
        outfile << fInf_A_val[iBinN] << "\t" << fInf_A_err[iBinN] << "\t" << fInf_C_val[iBinN] << "\t" << fInf_C_err[iBinN] << "\t"
                << fEff_A_val[iBinN] << "\t" << fEff_A_err[iBinN] << "\t" << fEff_C_val[iBinN] << "\t" << fEff_C_err[iBinN] << "\t"
                << fInf_val[iBinN] << "\t" << fInf_err[iBinN] << "\t" << fEff_val[iBinN] << "\t" << fEff_err[iBinN] << "\n";
    }
    outfile.close();
    Printf("*** Efficiencies printed to %s. ***", str_out.Data());
    return;
}

void VetoEff_SetEfficiencies(Bool_t sample)
{
    TRandom3 *ran = new TRandom3();
    ran->SetSeed(0);
    if(!sample){
        for(Int_t iBinN = 0; iBinN < nBinsN+1; iBinN++){
            SampledEff_A[iBinN] = fEff_A_val[iBinN];
            SampledEff_C[iBinN] = fEff_C_val[iBinN];
            SampledEff[iBinN] = fEff_val[iBinN];
        }
    } else {
        SampledEff_A[0] = 1.;
        SampledEff_C[0] = 1.;
        SampledEff[0] = 1.;
        for(Int_t iBinN = 1; iBinN < nBinsN+1; iBinN++){
            // sample the values until both positive
            Double_t val_A = ran->Gaus(fEff_A_val[iBinN],fEff_A_err[iBinN]);
            Double_t val_C = ran->Gaus(fEff_C_val[iBinN],fEff_C_err[iBinN]);
            Double_t val = ran->Gaus(fEff_val[iBinN],fEff_err[iBinN]);
            while(val_A < 0.0 || val_C < 0.0 || val_A > 1.0 || val_C > 1.0){
                Printf("Sampled values < 0 or > 1. Retrying.");
                val_A = ran->Gaus(fEff_A_val[iBinN],fEff_A_err[iBinN]);
                val_C = ran->Gaus(fEff_C_val[iBinN],fEff_C_err[iBinN]);
                val = ran->Gaus(fEff_val[iBinN],fEff_err[iBinN]);
            }
            SampledEff_A[iBinN] = val_A;
            SampledEff_C[iBinN] = val_C;
            SampledEff[iBinN] = val;
            hSampledEffPartial_A[iBinN-1]->Fill(val_A);
            hSampledEffPartial_C[iBinN-1]->Fill(val_C);
            hSampledEffPartial[iBinN-1]->Fill(val);
        }
    }
    ofstream outfile;
    TString str_out = "Results/" + str_subfolder + "VetoEfficiency/efficiencies_sampled.txt";
    outfile.open(str_out.Data());
    outfile << std::fixed << std::setprecision(3);
    for(Int_t iBinN = 0; iBinN < nBinsN+1; iBinN++) outfile << SampledEff_A[iBinN] << "\t" << SampledEff_C[iBinN] << "\t" << SampledEff[iBinN] << "\n";
    outfile.close();
    Printf("*** Efficiencies printed to %s. ***", str_out.Data());   
    return;
}

void ConnectTreeVariables_tNeutrons(TTree *t)
{
    // Set branch addresses
    t->SetBranchAddress("fPt", &fPt);
    t->SetBranchAddress("fM", &fM);
    t->SetBranchAddress("fZNA_time", &fZNA_time);
    t->SetBranchAddress("fZNC_time", &fZNC_time);
    t->SetBranchAddress("fZNA_hit", &fZNA_hit);
    t->SetBranchAddress("fZNC_hit", &fZNC_hit);
    t->SetBranchAddress("fZNA_energy", &fZNA_energy);
    t->SetBranchAddress("fZNC_energy", &fZNC_energy);
    t->SetBranchAddress("fZNA_n", &fZNA_n);
    t->SetBranchAddress("fZNC_n", &fZNC_n);

    Printf("Variables from %s connected.", t->GetName());

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

Double_t VetoEffiency_LoadBkg(Int_t iPt)
// iPt == 0 => allbins (0.2 to 1.0 GeV/c)
//     == i => bin i (1, 2, 3, 4, 5)
{
    TString str_bkg = "";
    if(iPt == 0) str_bkg = "Results/" + str_subfolder + "InvMassFit/allbins/allbins_bkg.txt";
    if(iPt >= 1) str_bkg = "Results/" + str_subfolder + "InvMassFit/" + Form("%ibins/bin%i_bkg.txt",nPtBins,iPt);
    Double_t fBkg_val(-1.), fBkg_err(-1.);
    ifstream ifs;
    ifs.open(str_bkg.Data());
    if(!ifs.fail())
    {
        // Read data from the file
        ifs >> fBkg_val >> fBkg_err;
        Printf("Loaded number of bkg events: %.1f.", fBkg_val);
    } else {
        Printf("Cannot open the file %s.", str_bkg.Data());
    }
    ifs.close();
    return fBkg_val;
}

class NeutronMatrix
{
    // first index (rows) = neutron bin in A
    // second index (columns) = neutron bin in C
    // first bin = no neutrons, then five neutron bins
    public:
        NeutronMatrix();
        ~NeutronMatrix();
        void     AddEvent(Int_t iBinA, Int_t iBinC) {fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] + 1;}
        Double_t GetBinContent(Int_t iBinA, Int_t iBinC) {return fEv_neutron_bins[iBinA][iBinC];}
        Double_t CountEvents_tot();
        Double_t CountEvents_0n0n() {return fEv_neutron_bins[0][0];}
        Double_t CountEvents_Xn0n();
        Double_t CountEvents_0nXn();
        Double_t CountEvents_XnXn();
        void     Multiply(Double_t x);
        void     SubtractMatrix(NeutronMatrix *nm);
        void     ApplyEfficiencies_AC();
        void     ApplyEfficiencies_combined1();
        void     ApplyEfficiencies_combined2();
        void     Plot(TString path);
        void     PrintToConsole();
        void     PrintToFile(TString name, Int_t precision = 0);
        void     LoadFromFile(TString name);
    private:
        Double_t fEv_neutron_bins[nBinsN+1][nBinsN+1];
};

NeutronMatrix::NeutronMatrix()
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) fEv_neutron_bins[iBinA][iBinC] = 0;
    }
    Printf("Neutron matrix created.");
}

NeutronMatrix::~NeutronMatrix(){}

Double_t NeutronMatrix::CountEvents_tot()
{
    Double_t sum = 0;
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) sum += fEv_neutron_bins[iBinA][iBinC];
    }
    return sum;
}

Double_t NeutronMatrix::CountEvents_Xn0n()
{
    Double_t sum = 0;
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++) sum += fEv_neutron_bins[iBinA][0];
    return sum;
}

Double_t NeutronMatrix::CountEvents_0nXn()
{
    Double_t sum = 0;
    for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) sum += fEv_neutron_bins[0][iBinC];
    return sum;
}

Double_t NeutronMatrix::CountEvents_XnXn()
{
    Double_t sum = 0;
    for(Int_t iBinA = 1; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 1; iBinC < nBinsN+1; iBinC++) sum += fEv_neutron_bins[iBinA][iBinC];
    }
    return sum;
}

void NeutronMatrix::Multiply(Double_t x)
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] * x;
    }
    return;
}

void NeutronMatrix::SubtractMatrix(NeutronMatrix *nm)
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) 
            if(fEv_neutron_bins[iBinA][iBinC] < nm->GetBinContent(iBinA,iBinC)) fEv_neutron_bins[iBinA][iBinC] = 0;
            else fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] - nm->GetBinContent(iBinA,iBinC);
    }
    return;
}

void NeutronMatrix::ApplyEfficiencies_AC()
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) 
            fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] / SampledEff_A[iBinA] / SampledEff_C[iBinC];
    }
    return;
}

void NeutronMatrix::ApplyEfficiencies_combined1()
{
    // 0n0n => we no correct at all
    // 0nXn
    for(Int_t iBinC = 1; iBinC < nBinsN+1; iBinC++){
        fEv_neutron_bins[0][iBinC] = fEv_neutron_bins[0][iBinC] / SampledEff[iBinC];
    }
    // Xn + whatever
    for(Int_t iBinA = 1; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) 
            fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] / SampledEff[iBinA];
    }
    return;
}

void NeutronMatrix::ApplyEfficiencies_combined2()
{
    // 0n0n => we no correct at all
    // Xn0n
    for(Int_t iBinA = 1; iBinA < nBinsN+1; iBinA++){
        fEv_neutron_bins[iBinA][0] = fEv_neutron_bins[iBinA][0] / SampledEff[iBinA];
    }
    // whatever + Xn
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 1; iBinC < nBinsN+1; iBinC++) 
            fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] / SampledEff[iBinC];
    }
    return;
}

void NeutronMatrix::Plot(TString path)
{
    TH2D *h = new TH2D("h","h",nBinsN+1,0.,6.,nBinsN+1,0.,6.);
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) h->SetBinContent(iBinA+1, iBinC+1, fEv_neutron_bins[iBinA][iBinC]);
    }
    TCanvas *c = new TCanvas("c","c",1050,900);
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.12);
    c->SetLeftMargin(0.16);
    c->cd();
    h->SetStats(0);
    h->SetTitle(0);
    gStyle->SetPaintTextFormat(".1f");
    h->GetXaxis()->SetTitle("A side, # of neutrons");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitle("C side, # of neutrons");
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.65);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetZaxis()->SetLabelSize(0.05);
    h->SetMarkerSize(2.5);
    h->Draw("colz,text");
    for(Int_t i = 0; i < nBinsN+1; i++)
    {
        h->GetXaxis()->SetBinLabel(i+1,sNumberOfN[i]);
        h->GetYaxis()->SetBinLabel(i+1,sNumberOfN[i]);
    }
    h->Draw("COLZ,TEXT");
    c->Print(path.Data());
    delete h;
    delete c;
    return;
}

void NeutronMatrix::PrintToConsole()
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) cout << fEv_neutron_bins[iBinA][iBinC] << "\t";
        cout << endl;
    }
    return;
}

void NeutronMatrix::PrintToFile(TString name, Int_t precision = 0)
{
    ofstream outfile;
    outfile.open(name.Data());
    outfile << std::fixed << std::setprecision(precision);
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) outfile << fEv_neutron_bins[iBinA][iBinC] << "\t";
        outfile << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", name.Data());
    return;
}

void NeutronMatrix::LoadFromFile(TString name)
{
    ifstream ifs;
    ifs.open(name.Data());
    if(!ifs.fail())
    {
        // Read data from the file
        for(Int_t iRow = 0; iRow < nBinsN+1; iRow++){
            for(Int_t iCol = 0; iCol < nBinsN+1; iCol++){
                ifs >> fEv_neutron_bins[iRow][iCol];
            }
        }
        Printf("Neutron matrix loaded from %s.", name.Data());
    } else {
        Printf("Cannot open the file %s.", name.Data());
    }
    ifs.close();
    return;
}