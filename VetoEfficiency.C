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
// histograms in the number of neutrons
const Int_t nBins = 200;
Double_t n_low = 0.; // number of neutrons
Double_t n_upp = 50.;
// arrays:
Double_t fNumberOfN[nBinsN+1] = {0.0, 1.5, 5.5, 10.5, 20.5, 50.5};
TString  sNumberOfN[nBinsN] = {"0-1", "2-5", "6-10", "11-20", "21-50"};
Double_t fVetoIneff_A[nBinsN] = {0.085, 0.157, 0.265, 0.413, 0.579};
Double_t fVetoIneff_C[nBinsN] = {0.146, 0.333, 0.425, 0.542, 0.844};
// ***********************************
// tree variables:
Bool_t fZNA_hit, fZNC_hit;
Double_t fZNA_n, fZNC_n;
// normalization when calculating fractions:
Double_t NORM;

void VetoEfficiency_ClassifyEventsToChannels(Int_t mass_range, Bool_t normalized);
void VetoEfficiency_PrepareTree();
TCanvas* PlotNeutronDistribution(const char* name, TH1 *hZNA, TH1 *hZNC, Double_t fPtMin, Double_t fPtMax, Double_t fMMin, Double_t fMMax);
void ConnectTreeVariables_tNeutrons(TTree *t);
void normalize(Double_t &val){ val = val / NORM; return; }

void VetoEfficiency(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VetoEfficiency/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/");

    // prepare the tree containing information about mass, pT and ZN signal
    VetoEfficiency_PrepareTree();

    // classify events into classes: background (mass from 1.8 to 2.8 GeV)
    VetoEfficiency_ClassifyEventsToChannels(0, kFALSE);
    VetoEfficiency_ClassifyEventsToChannels(0, kTRUE);

    // classify events into classes: signal+bkg (mass from 3.0 to 3.2 GeV)
    VetoEfficiency_ClassifyEventsToChannels(1, kFALSE);
    VetoEfficiency_ClassifyEventsToChannels(1, kTRUE);

    return;
}

void VetoEfficiency_ClassifyEventsToChannels(Int_t mass_range, Bool_t normalized)
// mass_range == 0 => background (mass range: fBkgM_low to fBkgM_upp; only 1 pT bin)
//            == 1 => signal+bkg (mass range: 3.0 to 3.2 GeV; only 1 pT bin)
{
    // **********************************************************************
    // numbers or fractions of events 
    // total number of events
    Double_t fEv_tot(0); 
    // per ZN channel summing over the corresponding bin types
    Double_t fEv_ch_Pt[4] = { 0 };
    Double_t fEv_ch_N[4] = { 0 };
    Double_t fEv_ch_PtN[4] = { 0 };
    // per pT bin
    Double_t fEv_BinsPt[nBinsPt] = { 0 };
    // in neutron classes
    Double_t fEv_0n0n(0), fEv_Xn0n(0), fEv_0nXn(0), fEv_XnXn(0);
    // per neutron class and pT bin
    Double_t fEv_0n0n_BinsPt[nBinsPt] = { 0 }; // 0n0n class, index = pT bin
    Double_t fEv_Xn0n_BinsPt[nBinsPt] = { 0 }; // Xn0n class, index = pT bin
    Double_t fEv_0nXn_BinsPt[nBinsPt] = { 0 }; // 0nXn class, index = pT bin
    Double_t fEv_XnXn_BinsPt[nBinsPt] = { 0 }; // XnXn class, index = pT bin
    // per neutron class and neutron bin
    Double_t fEv_Xn0n_BinsN[nBinsN] = { 0 };   // Xn0n class, index = neutron bin (A)
    Double_t fEv_0nXn_BinsN[nBinsN] = { 0 };   // 0nXn class, index = neutron bin (C)
    Double_t fEv_XnXn_BinsN[nBinsN][nBinsN] = { 0 }; // XnXn class, first index = neutron bin (A), second index = neutron bin (C)
    // per neutron class, pT bin and neutron bin
    Double_t fEv_Xn0n_BinsPtN[nBinsPt][nBinsN] = { 0 }; // Xn0n class, first index = pT bin, second index = neutron bin (A)
    Double_t fEv_0nXn_BinsPtN[nBinsPt][nBinsN] = { 0 }; // 0nXn class, first index = pT bin, second index = neutron bin (C)
    Double_t fEv_XnXn_BinsPtN[nBinsPt][nBinsN][nBinsN] = { 0 }; // XnXn class, first index = pT bin, second index = neutron bin (A), third index = neutron bin (C)
    // **********************************************************************

    Double_t m_low(0.), m_upp(0.);
    if(mass_range == 0){
        m_low = fBkgM_low;
        m_upp = fBkgM_upp;
    } else if(mass_range == 1){
        m_low = 3.0;
        m_upp = 3.2;
    }
    TString str_mass_subfolder = Form("mass_%.2fto%.2f/", m_low, m_upp);
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/" + str_mass_subfolder);

    TString name = "Trees/" + str_subfolder + "VetoEfficiency/tNeutrons.root";
    TFile *f_in = new TFile(name.Data(),"read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("tNeutrons"));
    if(t_in) Printf("Input tree %s loaded.", t_in->GetName());

    ConnectTreeVariables_tNeutrons(t_in);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    gROOT->cd();
    TH1D *hZNA[6] = { NULL }; // 0 = pT from 0.2 to 1.0 GeV/c, then pT bins
    TH1D *hZNC[6] = { NULL };
    TCanvas *c[6] = { NULL };

    for(Int_t i = 0; i < 6; i++){
        hZNA[i] = new TH1D(Form("hZNA%i",i),Form("hZNA%i",i),nBins,n_low,n_upp);
        hZNC[i] = new TH1D(Form("hZNC%i",i),Form("hZNC%i",i),nBins,n_low,n_upp);
    }

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
    {
        t_in->GetEntry(iEntry);

        if(!(fM > m_low && fM < m_upp)) continue;

        if(fZNA_hit == kTRUE && fZNA_n > 50.5){ Printf("Ev %i: More than 50 neutrons on A side (%.2f). Skipping...", iEntry, fZNA_n); continue; }
        if(fZNC_hit == kTRUE && fZNC_n > 50.5){ Printf("Ev %i: More than 50 neutrons on C side (%.2f). Skipping...", iEntry, fZNC_n); continue; }
        // find index of the neutron bin
        Int_t iBinN_A(0), iBinN_C(0);
        if(fZNA_hit) while(fZNA_n > fNumberOfN[iBinN_A+1]) iBinN_A++;
        if(fZNC_hit) while(fZNC_n > fNumberOfN[iBinN_C+1]) iBinN_C++;
        // find index of the pT bin
        Int_t iBinPt(0);
        while(fPt > ptBoundaries[iBinPt+1]) iBinPt++;
        // 0n0n class
        if(fZNA_hit == kFALSE && fZNC_hit == kFALSE)
        {
            fEv_0n0n++;
            fEv_0n0n_BinsPt[iBinPt]++;
        }
        // Xn0n class
        if(fZNA_hit == kTRUE && fZNC_hit == kFALSE)
        {
            fEv_Xn0n++;
            fEv_Xn0n_BinsPt[iBinPt]++;
            fEv_Xn0n_BinsN[iBinN_A]++;
            fEv_Xn0n_BinsPtN[iBinPt][iBinN_A]++;
        }
        // 0nXn class
        if(fZNA_hit == kFALSE && fZNC_hit == kTRUE)
        {
            fEv_0nXn++;
            fEv_0nXn_BinsPt[iBinPt]++;
            fEv_0nXn_BinsN[iBinN_C]++;
            fEv_0nXn_BinsPtN[iBinPt][iBinN_C]++;
        }
        // XnXn class
        if(fZNA_hit == kTRUE && fZNC_hit == kTRUE)
        {
            fEv_XnXn++;
            fEv_XnXn_BinsPt[iBinPt]++;
            fEv_XnXn_BinsN[iBinN_A][iBinN_C]++;
            fEv_XnXn_BinsPtN[iBinPt][iBinN_A][iBinN_C]++;        
        }
        // fill the histograms
        if(fZNA_hit) hZNA[0]->Fill(fZNA_n); 
        if(fZNC_hit) hZNC[0]->Fill(fZNC_n); 
        if(fZNA_hit) hZNA[iBinPt+1]->Fill(fZNA_n); 
        if(fZNC_hit) hZNC[iBinPt+1]->Fill(fZNC_n); 
    }

    f_in->Close();    

    // ##########################################################################################################
    // sum over classes in pT bins
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++) fEv_BinsPt[iBinPt] = fEv_0n0n_BinsPt[iBinPt] + fEv_Xn0n_BinsPt[iBinPt] + fEv_0nXn_BinsPt[iBinPt] + fEv_XnXn_BinsPt[iBinPt];
    // total sum 
    fEv_tot = fEv_0n0n + fEv_Xn0n + fEv_0nXn + fEv_XnXn;
    // total sum per ZN channel over pT bins
    for(Int_t i = 0; i < nPtBins; i++){
        fEv_ch_Pt[0] += fEv_0n0n_BinsPt[i];
        fEv_ch_Pt[1] += fEv_Xn0n_BinsPt[i];
        fEv_ch_Pt[2] += fEv_0nXn_BinsPt[i];
        fEv_ch_Pt[3] += fEv_XnXn_BinsPt[i];
    }
    // total sum per ZN channel over neutron bins
    fEv_ch_N[0] = fEv_0n0n;
    for(Int_t i = 0; i < 5; i++){
        fEv_ch_N[1] += fEv_Xn0n_BinsN[i];
        fEv_ch_N[2] += fEv_0nXn_BinsN[i];
        for(Int_t j = 0; j < 5; j++) fEv_ch_N[3] += fEv_XnXn_BinsN[i][j];
    }
    // total sum per ZN channel over pT bins and neutron bins
    fEv_ch_PtN[0] = fEv_0n0n;
    for(Int_t i = 0; i < 5; i++){
        for(Int_t j = 0; j < 5; j++){
            fEv_ch_PtN[1] += fEv_Xn0n_BinsPtN[i][j];
            fEv_ch_PtN[2] += fEv_0nXn_BinsPtN[i][j];
            for(Int_t k = 0; k < 5; k++) fEv_ch_PtN[3] += fEv_XnXn_BinsPtN[i][j][k];
        }
    }

    // ##########################################################################################################
    // calculate the fractions: normalize all numbers by the total number of events
    if(!normalized) NORM = 1.0;
    else            NORM = fEv_tot;
    normalize(fEv_tot); normalize(fEv_0n0n); normalize(fEv_Xn0n); normalize(fEv_0nXn); normalize(fEv_XnXn);
    for(Int_t i = 0; i < 4; i++){ 
        normalize(fEv_ch_Pt[i]); 
        normalize(fEv_ch_N[i]); 
        normalize(fEv_ch_PtN[i]);
    }
    for(Int_t i = 0; i < nPtBins; i++){
        normalize(fEv_BinsPt[i]);
        normalize(fEv_0n0n_BinsPt[i]);
        normalize(fEv_Xn0n_BinsPt[i]);
        normalize(fEv_0nXn_BinsPt[i]);
        normalize(fEv_XnXn_BinsPt[i]);
        for(Int_t j = 0; j < nBinsN; j++){
            normalize(fEv_Xn0n_BinsPtN[i][j]);
            normalize(fEv_0nXn_BinsPtN[i][j]);
            for(Int_t k = 0; k < nBinsN; k++) normalize(fEv_XnXn_BinsPtN[i][j][k]);
        }
    }
    for(Int_t i = 0; i < nBinsN; i++){
        normalize(fEv_Xn0n_BinsN[i]);
        normalize(fEv_0nXn_BinsN[i]);
        for(Int_t j = 0; j < nBinsN; j++) normalize(fEv_XnXn_BinsN[i][j]);
    }
    // ##########################################################################################################
    // plot neutron distribution in allbins
    c[0] = PlotNeutronDistribution("c0",hZNA[0],hZNC[0],0.2,1.0,m_low,m_upp);
    c[0]->Draw();
    TString str_out = "Results/" + str_subfolder + "VetoEfficiency/" + str_mass_subfolder + "ZN_n_allbins";
    c[0]->Print((str_out + ".pdf").Data());
    c[0]->Print((str_out + ".png").Data());
    // plots neutron distribution in bins
    for(Int_t i = 1; i < nPtBins+1; i++){
        c[i] = PlotNeutronDistribution(Form("c%i",i),hZNA[i],hZNC[i],ptBoundaries[i-1],ptBoundaries[i],m_low,m_upp);
        c[i]->Draw();
        str_out = Form("Results/%sVetoEfficiency/%sZN_n_bin%i", str_subfolder.Data(), str_mass_subfolder.Data(), i);
        c[i]->Print((str_out + ".pdf").Data());
        c[i]->Print((str_out + ".png").Data());
    }
    // ##########################################################################################################
    // print the numbers
    Int_t precision(0);
    if(normalized) precision = 3;
    // numbers of events per neutron class and pT bin
    ofstream outfile;
    if(!normalized) str_out = Form("Results/%sVetoEfficiency/%snEv_binsPt.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    else            str_out = Form("Results/%sVetoEfficiency/%sfEv_binsPt.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    outfile.open(str_out.Data());
    outfile << "pT_low\tpT_upp\t0n0n\tXn0n\t0nXn\tXnXn\ttotal\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++)
    {
        outfile << std::fixed << std::setprecision(3)
                << ptBoundaries[iBinPt] << "\t" << ptBoundaries[iBinPt+1] << "\t"
                << std::fixed << std::setprecision(precision)
                << fEv_0n0n_BinsPt[iBinPt] << "\t" << fEv_Xn0n_BinsPt[iBinPt] << "\t" << fEv_0nXn_BinsPt[iBinPt] << "\t" << fEv_XnXn_BinsPt[iBinPt] << "\t" << fEv_BinsPt[iBinPt] << "\n";
    }
    outfile << "sum\t\t" << fEv_ch_Pt[0] << "\t" << fEv_ch_Pt[1] << "\t" << fEv_ch_Pt[2] << "\t" << fEv_ch_Pt[3] << "\t" << fEv_tot; 
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data());

    // numbers of events per neutron class and neutron bin
    if(!normalized) str_out = Form("Results/%sVetoEfficiency/%snEv_binsN.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    else            str_out = Form("Results/%sVetoEfficiency/%sfEv_binsN.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    outfile.open(str_out.Data());
    // 0n0n
    outfile << std::fixed << std::setprecision(precision)
            << "0n0n:\n"
            << "total: " << fEv_ch_N[0] << "\n\n";
    // Xn0n
    outfile << "Xn0n:\n" 
            << "total: " << fEv_ch_N[1] << "\n"
            << "0-1\t2-5\t6-10\t11-20\t21-50\n"
            << fEv_Xn0n_BinsN[0] << "\t" 
            << fEv_Xn0n_BinsN[1] << "\t" 
            << fEv_Xn0n_BinsN[2] << "\t" 
            << fEv_Xn0n_BinsN[3] << "\t" 
            << fEv_Xn0n_BinsN[4] << "\n\n";
    // 0nXn
    outfile << "0nXn:\n" 
            << "total: " << fEv_ch_N[2] << "\n"
            << "0-1\t2-5\t6-10\t11-20\t21-50\n"
            << fEv_0nXn_BinsN[0] << "\t" 
            << fEv_0nXn_BinsN[1] << "\t" 
            << fEv_0nXn_BinsN[2] << "\t" 
            << fEv_0nXn_BinsN[3] << "\t" 
            << fEv_0nXn_BinsN[4] << "\n\n";
    // XnXn
    outfile << "XnXn:\n"
            << "total: " << fEv_ch_N[3] << "\n"
            << "row: fZNA_n\ncol: fZNC_n\n"
            << "\t0-1\t2-5\t6-10\t11-20\t21-50\n";
    for(Int_t iBinN = 0; iBinN < 5; iBinN++){
        outfile << sNumberOfN[iBinN] << "\t" << fEv_XnXn_BinsN[0][iBinN] << "\t" 
                                             << fEv_XnXn_BinsN[1][iBinN] << "\t" 
                                             << fEv_XnXn_BinsN[2][iBinN] << "\t" 
                                             << fEv_XnXn_BinsN[3][iBinN] << "\t" 
                                             << fEv_XnXn_BinsN[4][iBinN] << "\n";
    }        
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data()); 

    // numbers of events per neutron class, pT bin and neutron bin 
    if(!normalized) str_out = Form("Results/%sVetoEfficiency/%snEv_binsPtN.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    else            str_out = Form("Results/%sVetoEfficiency/%sfEv_binsPtN.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    outfile.open(str_out.Data());
    // 0n0n
    outfile << std::fixed << std::setprecision(precision)
            << "Xn0n:\n"
            << "total: " << fEv_ch_PtN[0] << "\n" 
            << "pT_low\tpT_upp\tfEv\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << std::fixed << std::setprecision(3)
                << ptBoundaries[iBinPt] << "\t" << ptBoundaries[iBinPt+1] << "\t"
                << std::fixed << std::setprecision(precision)
                << fEv_0n0n_BinsPt[iBinPt] << "\n";
    }
    outfile << "\n";
    // Xn0n
    outfile << "Xn0n:\n"
            << "total: " << fEv_ch_PtN[1] << "\n" 
            << "pT_low\tpT_upp\t0-1\t2-5\t6-10\t11-20\t21-50\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << std::fixed << std::setprecision(3)
                << ptBoundaries[iBinPt] << "\t" << ptBoundaries[iBinPt+1] << "\t"
                << std::fixed << std::setprecision(precision)
                << fEv_Xn0n_BinsPtN[iBinPt][0] << "\t"
                << fEv_Xn0n_BinsPtN[iBinPt][1] << "\t"
                << fEv_Xn0n_BinsPtN[iBinPt][2] << "\t"
                << fEv_Xn0n_BinsPtN[iBinPt][3] << "\t"
                << fEv_Xn0n_BinsPtN[iBinPt][4] << "\n";
    }
    outfile << "\n";
    // 0nXn
    outfile << "0nXn:\n"
            << "total: " << fEv_ch_PtN[2] << "\n" 
            << "pT_low\tpT_upp\t0-1\t2-5\t6-10\t11-20\t21-50\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << std::fixed << std::setprecision(3)
                << ptBoundaries[iBinPt] << "\t" << ptBoundaries[iBinPt+1] << "\t"
                << std::fixed << std::setprecision(precision)
                << fEv_0nXn_BinsPtN[iBinPt][0] << "\t"
                << fEv_0nXn_BinsPtN[iBinPt][1] << "\t"
                << fEv_0nXn_BinsPtN[iBinPt][2] << "\t"
                << fEv_0nXn_BinsPtN[iBinPt][3] << "\t"
                << fEv_0nXn_BinsPtN[iBinPt][4] << "\n";
    }
    outfile << "\n";
    // XnXn
    outfile << "XnXn:\n"
            << "total: " << fEv_ch_PtN[3] << "\n"
            << "row: fZNA_n\ncol: fZNC_n\n";
    // over pT bins
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << Form("pT_low = %.3f, pT_upp = %.3f\n", ptBoundaries[iBinPt], ptBoundaries[iBinPt+1])
                << "\t0-1\t2-5\t6-10\t11-20\t21-50\n";
        for(Int_t i = 0; i < 5; i++){
            outfile << sNumberOfN[i] << "\t";
            for(Int_t j = 0; j < 5; j++) outfile << fEv_XnXn_BinsPtN[iBinPt][i][j] << "\t";
            outfile << "\n";
        }
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", str_out.Data()); 

    // ##########################################################################################################
    
    for(Int_t i = 0; i < 6; i++) {delete c[i]; delete hZNA[i]; delete hZNC[i];}

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

TCanvas* PlotNeutronDistribution(const char* name, TH1 *hZNA, TH1 *hZNC, Double_t fPtMin, Double_t fPtMax, Double_t fMMin, Double_t fMMax)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas(name,name,1500,500);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.07);
    c->SetRightMargin(0.015);
    //c->SetLogy();
    // X-axis
    hZNA->GetXaxis()->SetTitle("# of neutrons (ZN energy/2510 GeV)");
    hZNA->GetXaxis()->SetTitleSize(0.06);
    hZNA->GetXaxis()->SetLabelSize(0.06);
    // Y-axis
    hZNA->GetYaxis()->SetTitle("Counts per 0.25");
    hZNA->GetYaxis()->SetTitleSize(0.06);
    hZNA->GetYaxis()->SetLabelSize(0.06);
    hZNA->GetYaxis()->SetTitleOffset(0.53);
    // Style hist ZNA
    hZNA->SetLineColor(kRed);
    hZNA->SetLineWidth(1);
    hZNA->SetMarkerStyle(21);
    hZNA->SetMarkerColor(kRed);
    hZNA->SetMarkerSize(0.5);
    // Draw
    hZNA->Draw("HIST");
    if(hZNC){
        // Style hist ZNC
        hZNC->SetLineColor(kBlue);
        hZNC->SetLineWidth(1);
        hZNC->SetMarkerStyle(21);
        hZNC->SetMarkerColor(kBlue);
        hZNC->SetMarkerSize(0.5);
        // Draw
        hZNC->Draw("SAME HIST");
    }
    // legend1
    TLegend *l1 = new TLegend(0.60,0.72,0.99,0.96);
    l1->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    l1->AddEntry((TObject*)0,Form("#it{m}_{#mu#mu} #in (%.2f,%.2f) GeV/#it{c}^{2}",fMMin,fMMax),"");
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}",fPtMin,fPtMax),"");
    l1->SetTextSize(0.06);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();
    // legend2
    TLegend *l2 = new TLegend(0.70,0.56,0.92,0.70);
    l2->AddEntry(hZNA,Form("ZNA (total: %.0f events)", hZNA->Integral()),"L");
    l2->AddEntry(hZNC,Form("ZNC (total: %.0f events)", hZNC->Integral()),"L");
    l2->SetTextSize(0.06);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    return c;
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