// _VetoInefficiencies.C
// David Grund, Mar 20, 2022
// use pass3 data and defined pT boundaries (to be found in Results/pass3_4bins_calibPID/)
// then calculate a weighted veto inefficiency for each neutron class and neutron bin number
// use the values from Guillermo's presentation in Documents/VetoIncoherent.pdf

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h" 

Double_t *fPtBins = NULL;
Double_t fPtBins_4[5] = {0.200, 0.282, 0.381, 0.571, 1.000};
//Double_t fPtBins_5[6] = {0.200, 1.000};

Double_t fNumberOfN[6] = {0.0, 1.5, 5.5, 10.5, 20.5, 50.5};
TString sNumberOfN[5] = {"0-1", "2-5", "6-10", "11-20", "21-50"};
Double_t fVetoIneff_A[5] = {0.085, 0.157, 0.265, 0.413, 0.579};
Double_t fVetoIneff_C[5] = {0.146, 0.333, 0.425, 0.542, 0.844};
Double_t fVetoEff_A[5] = { 0 };
Double_t fVetoEff_C[5] = { 0 };

Bool_t fZNA_hit, fZNC_hit;
Double_t fZNA_n, fZNC_n;

const Int_t nBins = 100;
Double_t n_low = 0.0; // number of neutrons
Double_t n_upp = 20.0;

// total number of events
Int_t nEv_tot(0); 
// cross-check: total number of events per ZN channel summing over the corresponding bin types
Int_t nEv_ch_Pt[4] = { 0 };
Int_t nEv_ch_N[4] = { 0 };
Int_t nEv_ch_PtN[4] = { 0 };
// total numbers per pT bin
Int_t nEv_BinsPt[5] = { 0 };
// total numbers per neutron class
Int_t nEv_0n0n(0), nEv_Xn0n(0), nEv_0nXn(0), nEv_XnXn(0);
// numbers per neutron class and pT bin
Int_t nEv_0n0n_BinsPt[5] = { 0 }; // 0n0n class, index = pT bin
Int_t nEv_Xn0n_BinsPt[5] = { 0 }; // Xn0n class, index = pT bin
Int_t nEv_0nXn_BinsPt[5] = { 0 }; // 0nXn class, index = pT bin
Int_t nEv_XnXn_BinsPt[5] = { 0 }; // XnXn class, index = pT bin
// numbers per neutron class and neutron bin
Int_t nEv_Xn0n_BinsN[5] = { 0 };    // Xn0n class, index = neutron bin
Int_t nEv_0nXn_BinsN[5] = { 0 };    // 0nXn class, index = neutron bin
Int_t nEv_XnXn_BinsN[5][5] = { 0 }; // XnXn class, first index = neutron bin (A), second index = neutron bin (C)
// numbers per neutron class, pT bin and neutron bin
Int_t nEv_Xn0n_BinsPtN[5][5] = { 0 };   // Xn0n class, first index = pT bin, second index = neutron bin (A)
Int_t nEv_0nXn_BinsPtN[5][5] = { 0 };   // 0nXn class, first index = pT bin, second index = neutron bin (C)
Int_t nEv_XnXn_BinsPtN[5][5][5] = { 0 };// XnXn class, first index = pT bin, second index = neutron bin (A), third index = neutron bin (C)
// the total veto inefficiency
Double_t fVetoIneff_direct = 0.;
Double_t fVetoIneff_weight = 0.;
// veto inefficiencies per pT bins
Double_t fVetoIneff_BinsPt[5] = { 0 };
// the total veto efficiency
Double_t fVetoEff_direct = 0.;
Double_t fVetoEff_weight = 0.;
// veto efficiencies per pT bins
Double_t fVetoEff_BinsPt[5] = { 0 };

// veto inefficiency
void VetoIneff_Calculate(Int_t nPtBins); // 4 or 5
void VetoIneff_PrintNumbers();
void VetoIneff_PrepareTree();
// distribution of neutrons in same-sign events
void SameSignEv_MakePlots();
void SameSignEv_PrepareTree();
// opposite-sign events, low pT
void LowPt_MakePlots();
void LowPt_PrepareTree();
// support function
void ConnectTreeVariables_tNeutrons(TTree *t);
TCanvas* PlotNeutronDistribution(const char* name, TH1 *hZNA, TH1 *hZNC, Double_t fPtMin, Double_t fPtMax);

void _VetoInefficiencies()
{
    gSystem->Exec("mkdir -p Trees/_VetoInefficiencies/");
    gSystem->Exec("mkdir -p Results/_VetoInefficiencies/");

    if(kTRUE)
    {
        //VetoIneff_PrepareTree();

        VetoIneff_PrintNumbers();

        VetoIneff_Calculate(4);

        //VetoIneff_Calculate(5);

        //SameSignEv_PrepareTree();

        SameSignEv_MakePlots();

        //LowPt_PrepareTree();

        LowPt_MakePlots();
    }

    return;
}

void VetoIneff_Calculate(Int_t nPtBins)
{
    if(nPtBins == 4)      fPtBins = &fPtBins_4[0];
    //else if(nPtBins == 5) fPtBins = &fPtBins_5[0];

    TFile *f_in = new TFile("Trees/_VetoInefficiencies/tNeutrons.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("tNeutrons"));
    if(t_in) Printf("Input tree %s loaded.", t_in->GetName());

    ConnectTreeVariables_tNeutrons(t_in);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    gROOT->cd();
    TH1D *hZNA[6] = { NULL }; // 0 = pT from 0.2 to 1.0 GeV/c
    TH1D *hZNC[6] = { NULL };
    TCanvas *c[6] = { NULL };

    for(Int_t i = 0; i < 6; i++){
        hZNA[i] = new TH1D(Form("hZNA%i",i),Form("hZNA%i",i),nBins,n_low,n_upp);
        hZNC[i] = new TH1D(Form("hZNC%i",i),Form("hZNC%i",i),nBins,n_low,n_upp);
    }

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
    {
        t_in->GetEntry(iEntry);

        if(fZNA_hit == kTRUE && fZNA_n > 50.5){ Printf("Ev %i: More than 50 neutrons on A side (%.2f). Skipping...", iEntry, fZNA_n); continue; }
        if(fZNC_hit == kTRUE && fZNC_n > 50.5){ Printf("Ev %i: More than 50 neutrons on C side (%.2f). Skipping...", iEntry, fZNC_n); continue; }
        // find index of the neutron bin
        Int_t iBinN_A(0), iBinN_C(0);
        if(fZNA_hit) while(fZNA_n > fNumberOfN[iBinN_A+1]) iBinN_A++;
        if(fZNC_hit) while(fZNC_n > fNumberOfN[iBinN_C+1]) iBinN_C++;
        // find index of the pT bin
        Int_t iBinPt(0);
        while(fPt > fPtBins[iBinPt+1]) iBinPt++;
        // 0n0n class
        if(fZNA_hit == kFALSE && fZNC_hit == kFALSE)
        {
            nEv_0n0n++;
            nEv_0n0n_BinsPt[iBinPt]++;
        }
        // Xn0n class
        if(fZNA_hit == kTRUE && fZNC_hit == kFALSE)
        {
            nEv_Xn0n++;
            nEv_Xn0n_BinsPt[iBinPt]++;
            nEv_Xn0n_BinsN[iBinN_A]++;
            nEv_Xn0n_BinsPtN[iBinPt][iBinN_A]++;
        }
        // 0nXn class
        if(fZNA_hit == kFALSE && fZNC_hit == kTRUE)
        {
            nEv_0nXn++;
            nEv_0nXn_BinsPt[iBinPt]++;
            nEv_0nXn_BinsN[iBinN_C]++;
            nEv_0nXn_BinsPtN[iBinPt][iBinN_C]++;
        }
        // XnXn class
        if(fZNA_hit == kTRUE && fZNC_hit == kTRUE)
        {
            nEv_XnXn++;
            nEv_XnXn_BinsPt[iBinPt]++;
            nEv_XnXn_BinsN[iBinN_A][iBinN_C]++;
            nEv_XnXn_BinsPtN[iBinPt][iBinN_A][iBinN_C]++;          
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
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++) nEv_BinsPt[iBinPt] = nEv_0n0n_BinsPt[iBinPt] + nEv_Xn0n_BinsPt[iBinPt] + nEv_0nXn_BinsPt[iBinPt] + nEv_XnXn_BinsPt[iBinPt];
    // total sum 
    nEv_tot = nEv_0n0n + nEv_Xn0n + nEv_0nXn + nEv_XnXn;
    // total sum per ZN channel over pT bins
    for(Int_t i = 0; i < nPtBins; i++){
        nEv_ch_Pt[0] += nEv_0n0n_BinsPt[i];
        nEv_ch_Pt[1] += nEv_Xn0n_BinsPt[i];
        nEv_ch_Pt[2] += nEv_0nXn_BinsPt[i];
        nEv_ch_Pt[3] += nEv_XnXn_BinsPt[i];
    }
    // total sum per ZN channel over neutron bins
    nEv_ch_N[0] = nEv_0n0n;
    for(Int_t i = 0; i < 5; i++){
        nEv_ch_N[1] += nEv_Xn0n_BinsN[i];
        nEv_ch_N[2] += nEv_0nXn_BinsN[i];
        for(Int_t j = 0; j < 5; j++) nEv_ch_N[3] += nEv_XnXn_BinsN[i][j];
    }
    // total sum per ZN channel over pT bins and neutron bins
    nEv_ch_PtN[0] = nEv_0n0n;
    for(Int_t i = 0; i < 5; i++){
        for(Int_t j = 0; j < 5; j++){
            nEv_ch_PtN[1] += nEv_Xn0n_BinsPtN[i][j];
            nEv_ch_PtN[2] += nEv_0nXn_BinsPtN[i][j];
            for(Int_t k = 0; k < 5; k++) nEv_ch_PtN[3] += nEv_XnXn_BinsPtN[i][j][k];
        }
    }

    // ##########################################################################################################
    // print the results
    gSystem->Exec(Form("mkdir -p Results/_VetoInefficiencies/%ibins/", nPtBins));
    // plot in allbins
    c[0] = PlotNeutronDistribution("c0",hZNA[0],hZNC[0],0.2,1.0);
    c[0]->Draw();
    TString str = "Results/_VetoInefficiencies/ZN_n_allbins";
    c[0]->Print((str + ".pdf").Data());
    c[0]->Print((str + ".png").Data());
    // plots in bins
    for(Int_t i = 1; i < nPtBins+1; i++){
        c[i] = PlotNeutronDistribution(Form("c%i",i),hZNA[i],hZNC[i],fPtBins[i-1],fPtBins[i]);
        c[i]->Draw();
        TString str = Form("Results/_VetoInefficiencies/%ibins/ZN_n_bin%i", nPtBins, i);
        c[i]->Print((str + ".pdf").Data());
        c[i]->Print((str + ".png").Data());
    }
    // numbers of events per neutron class and pT bin
    ofstream outfile;
    outfile.open(Form("Results/_VetoInefficiencies/%ibins/nEvents_binsPt.txt", nPtBins));
    outfile << "pT_low\tpT_upp\t0n0n\tXn0n\t0nXn\tXnXn\ttotal\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++)
    {
        outfile << std::fixed << std::setprecision(3)
                << fPtBins[iBinPt] << "\t" << fPtBins[iBinPt+1] << "\t"
                << nEv_0n0n_BinsPt[iBinPt] << "\t" << nEv_Xn0n_BinsPt[iBinPt] << "\t" << nEv_0nXn_BinsPt[iBinPt] << "\t" << nEv_XnXn_BinsPt[iBinPt] << "\t" << nEv_BinsPt[iBinPt] << "\n";
    }
    outfile << "sum\t\t" << nEv_ch_Pt[0] << "\t" << nEv_ch_Pt[1] << "\t" << nEv_ch_Pt[2] << "\t" << nEv_ch_Pt[3] << "\t" << nEv_tot; 
    outfile.close();
    Printf("*** Results printed to Results/_VetoInefficiencies/%ibins/nEvents_binsPt.txt.***", nPtBins);

    // numbers of events per neutron class and neutron bin
    outfile.open("Results/_VetoInefficiencies/nEvents_binsN.txt");
    // 0n0n
    outfile << "0n0n:\n"
            << "total: " << nEv_ch_N[0] << "\n\n";
    // Xn0n
    outfile << "Xn0n:\n" 
            << "total: " << nEv_ch_N[2] << "\n"
            << "0-1\t2-5\t6-10\t11-20\t21-50\n"
            << nEv_Xn0n_BinsN[0] << "\t" 
            << nEv_Xn0n_BinsN[1] << "\t" 
            << nEv_Xn0n_BinsN[2] << "\t" 
            << nEv_Xn0n_BinsN[3] << "\t" 
            << nEv_Xn0n_BinsN[4] << "\n\n";
    // 0nXn
    outfile << "0nXn:\n" 
            << "total: " << nEv_ch_N[1] << "\n"
            << "0-1\t2-5\t6-10\t11-20\t21-50\n"
            << nEv_0nXn_BinsN[0] << "\t" 
            << nEv_0nXn_BinsN[1] << "\t" 
            << nEv_0nXn_BinsN[2] << "\t" 
            << nEv_0nXn_BinsN[3] << "\t" 
            << nEv_0nXn_BinsN[4] << "\n\n";
    // XnXn
    outfile << "XnXn:\n"
            << "total: " << nEv_ch_N[3] << "\n"
            << "row: fZNA_n\ncol: fZNC_n\n"
            << "\t0-1\t2-5\t6-10\t11-20\t21-50\n";
    for(Int_t iBinN = 0; iBinN < 5; iBinN++){
        outfile << sNumberOfN[iBinN] << "\t" << nEv_XnXn_BinsN[0][iBinN] << "\t" 
                                             << nEv_XnXn_BinsN[1][iBinN] << "\t" 
                                             << nEv_XnXn_BinsN[2][iBinN] << "\t" 
                                             << nEv_XnXn_BinsN[3][iBinN] << "\t" 
                                             << nEv_XnXn_BinsN[4][iBinN] << "\n";
    }        
    outfile.close();
    Printf("*** Results printed to Results/_VetoInefficiencies/nEvents_binsN.txt.***"); 

    // numbers of events per neutron class, pT bin and neutron bin 
    outfile.open(Form("Results/_VetoInefficiencies/%ibins/nEvents_binsPtN.txt", nPtBins));
    // 0n0n
    outfile << "Xn0n:\n"
            << "total: " << nEv_ch_PtN[0] << "\n" 
            << "pT_low\tpT_upp\tnEv\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << std::fixed << std::setprecision(3)
                << fPtBins[iBinPt] << "\t" << fPtBins[iBinPt+1] << "\t"
                << nEv_0n0n_BinsPt[iBinPt] << "\n";
    }
    outfile << "\n";
    // Xn0n
    outfile << "Xn0n:\n"
            << "total: " << nEv_ch_PtN[1] << "\n" 
            << "pT_low\tpT_upp\t0-1\t2-5\t6-10\t11-20\t21-50\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << std::fixed << std::setprecision(3)
                << fPtBins[iBinPt] << "\t" << fPtBins[iBinPt+1] << "\t"
                << nEv_Xn0n_BinsPtN[iBinPt][0] << "\t"
                << nEv_Xn0n_BinsPtN[iBinPt][1] << "\t"
                << nEv_Xn0n_BinsPtN[iBinPt][2] << "\t"
                << nEv_Xn0n_BinsPtN[iBinPt][3] << "\t"
                << nEv_Xn0n_BinsPtN[iBinPt][4] << "\n";
    }
    outfile << "\n";
    // 0nXn
    outfile << "0nXn:\n"
            << "total: " << nEv_ch_PtN[2] << "\n" 
            << "pT_low\tpT_upp\t0-1\t2-5\t6-10\t11-20\t21-50\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << std::fixed << std::setprecision(3)
                << fPtBins[iBinPt] << "\t" << fPtBins[iBinPt+1] << "\t"
                << nEv_0nXn_BinsPtN[iBinPt][0] << "\t"
                << nEv_0nXn_BinsPtN[iBinPt][1] << "\t"
                << nEv_0nXn_BinsPtN[iBinPt][2] << "\t"
                << nEv_0nXn_BinsPtN[iBinPt][3] << "\t"
                << nEv_0nXn_BinsPtN[iBinPt][4] << "\n";
    }
    outfile << "\n";
    // XnXn
    outfile << "XnXn:\n"
            << "total: " << nEv_ch_PtN[3] << "\n"
            << "row: fZNA_n\ncol: fZNC_n\n";
    // over pT bins
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++){
        outfile << Form("pT_low = %.3f, pT_upp = %.3f\n", fPtBins[iBinPt], fPtBins[iBinPt+1])
                << "\t0-1\t2-5\t6-10\t11-20\t21-50\n";
        for(Int_t i = 0; i < 5; i++){
            outfile << sNumberOfN[i] << "\t";
            for(Int_t j = 0; j < 5; j++) outfile << nEv_XnXn_BinsPtN[iBinPt][i][j] << "\t";
            outfile << "\n";
        }
    }

    outfile.close();
    Printf("*** Results printed to Results/_VetoInefficiencies/%ibins/nEvents_binsPtN.txt.***", nPtBins);  

    // ##########################################################################################################
    // calculate veto inefficiencies from Guillermo's numbers
    // in pT bins
    // and the total number as a weighted average over pT bins
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++)
    {
        Double_t Ineff_Xn0n(0), Ineff_0nXn(0), Ineff_XnXn(0);
        for(Int_t iBinN1 = 0; iBinN1 < 5; iBinN1++){
            Ineff_Xn0n += fVetoIneff_A[iBinN1] * nEv_Xn0n_BinsPtN[iBinPt][iBinN1];
            Ineff_0nXn += fVetoIneff_C[iBinN1] * nEv_0nXn_BinsPtN[iBinPt][iBinN1];
            for(Int_t iBinN2 = 0; iBinN2 < 5; iBinN2++){
                Ineff_XnXn += fVetoIneff_A[iBinN1] * fVetoIneff_C[iBinN2] * nEv_XnXn_BinsPtN[iBinPt][iBinN1][iBinN2];
            }
        }
        fVetoIneff_BinsPt[iBinPt] = (Ineff_Xn0n + Ineff_0nXn + Ineff_XnXn) / nEv_BinsPt[iBinPt];
        fVetoIneff_weight += fVetoIneff_BinsPt[iBinPt] * nEv_BinsPt[iBinPt];
    }
    fVetoIneff_weight = fVetoIneff_weight / (Double_t)nEv_tot;
    // the total number: direct calculation
    Double_t Ineff_Xn0n(0), Ineff_0nXn(0), Ineff_XnXn(0);
    for(Int_t iBinN1 = 0; iBinN1 < 5; iBinN1++){
        Ineff_Xn0n += fVetoIneff_A[iBinN1] * nEv_Xn0n_BinsN[iBinN1];
        Ineff_0nXn += fVetoIneff_C[iBinN1] * nEv_0nXn_BinsN[iBinN1];
        for(Int_t iBinN2 = 0; iBinN2 < 5; iBinN2++){
            Ineff_XnXn += fVetoIneff_A[iBinN1] * fVetoIneff_C[iBinN2] * nEv_XnXn_BinsN[iBinN1][iBinN2];
        }
    }
    fVetoIneff_direct = (Ineff_Xn0n + Ineff_0nXn + Ineff_XnXn) / nEv_tot;
    // ##########################################################################################################
    // print the results
    outfile.open(Form("Results/_VetoInefficiencies/%ibins/VetoIneff.txt", nPtBins));
    outfile << "total value:\n"
            << Form("directly calculated: %.4f = %.2f%%\n", fVetoIneff_direct, fVetoIneff_direct*100)
            << Form("weighted avg over pT bins: %.4f = %.2f%%\n", fVetoIneff_weight, fVetoIneff_weight*100)
            << "in pT bins:\n"
            << "pT_low\tpT_upp\tineff\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++)
    {
        outfile << std::fixed << std::setprecision(3)
                << fPtBins[iBinPt] << "\t" 
                << fPtBins[iBinPt+1] << "\t"
                << fVetoIneff_BinsPt[iBinPt] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to Results/_VetoInefficiencies/%ibins/VetoIneff.txt.***", nPtBins); 
    // ##########################################################################################################
    // calculate veto efficiencies from Guillermo's numbers
    // in pT bins
    // and the total number as a weighted average over pT bins
    for(Int_t i = 0; i < 5; i++)
    {
        fVetoEff_A[i] = 1 - fVetoIneff_A[i];
        fVetoEff_C[i] = 1 - fVetoIneff_C[i];
    }

    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++)
    {
        Double_t Eff_0n0n(0), Eff_Xn0n(0), Eff_0nXn(0), Eff_XnXn(0);
        for(Int_t iBinN1 = 0; iBinN1 < 5; iBinN1++){
            Eff_Xn0n += fVetoEff_A[iBinN1] * nEv_Xn0n_BinsPtN[iBinPt][iBinN1];
            Eff_0nXn += fVetoEff_C[iBinN1] * nEv_0nXn_BinsPtN[iBinPt][iBinN1];
            for(Int_t iBinN2 = 0; iBinN2 < 5; iBinN2++){
                Eff_XnXn += fVetoEff_A[iBinN1] * fVetoEff_C[iBinN2] * nEv_XnXn_BinsPtN[iBinPt][iBinN1][iBinN2];
            }
        }
        Eff_0n0n = 1.0 * nEv_0n0n_BinsPt[iBinPt];
        fVetoEff_BinsPt[iBinPt] = (Eff_0n0n + Eff_Xn0n + Eff_0nXn + Eff_XnXn) / nEv_BinsPt[iBinPt];
        fVetoEff_weight += fVetoEff_BinsPt[iBinPt] * nEv_BinsPt[iBinPt];
    }
    fVetoEff_weight = fVetoEff_weight / (Double_t)nEv_tot;
    // the total number: direct calculation
    Double_t Eff_0n0n(0), Eff_Xn0n(0), Eff_0nXn(0), Eff_XnXn(0);
    for(Int_t iBinN1 = 0; iBinN1 < 5; iBinN1++){
        Eff_Xn0n += fVetoEff_A[iBinN1] * nEv_Xn0n_BinsN[iBinN1];
        Eff_0nXn += fVetoEff_C[iBinN1] * nEv_0nXn_BinsN[iBinN1];
        for(Int_t iBinN2 = 0; iBinN2 < 5; iBinN2++){
            Eff_XnXn += fVetoEff_A[iBinN1] * fVetoEff_C[iBinN2] * nEv_XnXn_BinsN[iBinN1][iBinN2];
        }
    }
    Eff_0n0n = 1.0 * nEv_0n0n;
    fVetoEff_direct = (Eff_0n0n + Eff_Xn0n + Eff_0nXn + Eff_XnXn) / nEv_tot;
    // ##########################################################################################################
    // print the results
    outfile.open(Form("Results/_VetoInefficiencies/%ibins/VetoEff.txt", nPtBins));
    outfile << "total value:\n"
            << Form("directly calculated: %.4f = %.2f%%\n", fVetoEff_direct, fVetoEff_direct*100)
            << Form("weighted avg over pT bins: %.4f = %.2f%%\n", fVetoEff_weight, fVetoEff_weight*100)
            << "in pT bins:\n"
            << "pT_low\tpT_upp\teff\n";
    for(Int_t iBinPt = 0; iBinPt < nPtBins; iBinPt++)
    {
        outfile << std::fixed << std::setprecision(3)
                << fPtBins[iBinPt] << "\t" 
                << fPtBins[iBinPt+1] << "\t"
                << fVetoEff_BinsPt[iBinPt] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to Results/_VetoInefficiencies/%ibins/VetoEff.txt.***", nPtBins); 
    // ##########################################################################################################

    for(Int_t i = 0; i < 6; i++) delete c[i];

    return;
}

void VetoIneff_PrintNumbers()
{
    TFile *f_in = new TFile("Trees/_VetoInefficiencies/tNeutrons.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("tNeutrons"));
    if(t_in) Printf("Input tree %s loaded.", t_in->GetName());

    ConnectTreeVariables_tNeutrons(t_in);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    // Go over the entries and print the information to the text file
    ofstream outfile("Results/_VetoInefficiencies/tNeutrons_eventInfo.txt");
    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
    {
        t_in->GetEntry(iEntry);
        outfile << std::fixed << std::setprecision(3)
                << fPt << "\t" << fM << "\t" << fZNA_hit << "\t" << fZNC_hit << "\t" 
                << std::fixed << std::setprecision(2)
                << fZNA_n << "\t" << fZNC_n << "\n";
    }
    outfile.close();
    Printf("*** Results printed to Results/_VetoInefficiencies/tNeutrons_eventInfo.txt.***");

    f_in->Close();

    return;
}

void VetoIneff_PrepareTree()
{
    isPass3 = kTRUE;

    SetReducedRunList(kTRUE);

    TFile *f_in = new TFile("Trees/AnalysisData_pass3/AnalysisResults.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("AnalysisOutput/fTreeJpsi"));
    if(t_in) Printf("Input tree %s loaded.", t_in->GetName());

    ConnectTreeVariables(t_in);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    gROOT->cd();
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

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
    {
        t_in->GetEntry(iEntry);
        // m in 3.0 to 3.2 GeV/c^2, pT in 0.2 to 1.0 GeV/c
        if(EventPassed(1, 3))
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

    // create new file
    TFile *f_out = new TFile("Trees/_VetoInefficiencies/tNeutrons.root","RECREATE");
    // open the file
    f_out->cd();
    // write the list and tree to this directory
    t_out->Write("tNeutrons",TObject::kSingleKey);
    // list the contents of the file
    f_out->ls();
    // close the file
    f_out->Close();

    f_in->Close();

    return;
}

void SameSignEv_PrepareTree()
{
    isPass3 = kTRUE;

    SetReducedRunList(kTRUE);

    TFile *f_in = new TFile("Trees/AnalysisData_pass3/AnalysisResults.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("AnalysisOutput/fTreeJpsi"));
    if(t_in) Printf("Input tree %s loaded.", t_in->GetName());

    ConnectTreeVariables(t_in);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    gROOT->cd();
    TTree *tSameSign = new TTree("tSameSign","tSameSign");
    tSameSign->Branch("fPt", &fPt, "fPt/D");
    tSameSign->Branch("fM", &fM, "fM/D");
    tSameSign->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    tSameSign->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    tSameSign->Branch("fZNA_hit", &fZNA_hit, "fZNA_hit/O");
    tSameSign->Branch("fZNC_hit", &fZNC_hit, "fZNC_hit/O");
    tSameSign->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    tSameSign->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    tSameSign->Branch("fZNA_n", &fZNA_n, "fZNA_n/D");
    tSameSign->Branch("fZNC_n", &fZNC_n, "fZNC_n/D");

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
    {
        t_in->GetEntry(iEntry);

        if((iEntry+1) % 100000 == 0)
        {
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

        // Run number in the GoodHadronPID lists published by DPG
        if(!RunNumberInListOfGoodRuns()) continue;

        // Selections applied on the GRID:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        // 2) Central UPC trigger CCUP31:
        // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
        // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD

        // 3) At least two tracks associated with the vertex
        if(fVertexContrib < cut_fVertexContrib) continue;

        // 4) Distance from the IP lower than 15 cm
        if(fVertexZ > cut_fVertexZ) continue;
    
        // 5a) ADA offline veto (no effect on MC)
        if(!(fADA_dec == 0)) continue;

        // 5b) ADC offline veto (no effect on MC)
        if(!(fADC_dec == 0)) continue;

        // 6a) V0A offline veto (no effect on MC)
        if(!(fV0A_dec == 0)) continue;

        // 6b) V0C offline veto (no effect on MC)
        if(!(fV0C_dec == 0)) continue;

        // 7) SPD cluster matches FOhits
        if(!(fMatchingSPD == kTRUE)) continue;

        // 8) Muon pairs only
        if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) continue;

        // 9) Dilepton rapidity |y| < 0.8
        if(!(abs(fY) < 0.8)) continue;

        // 10) Pseudorapidity of both tracks |eta| < 0.8
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;

        // 11) Tracks have same charges (!)
        if(fQ1 * fQ2 < 0) continue;

        // 12) Invariant mass cut
        if(!(fM > 3.0 && fM < 3.2)) continue;

        // prepate ZDC info
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

        // 13) Transverse momentum cut
        if(fPt < 1.00) tSameSign->Fill();
    }

    Printf("Tree %s filled with %lli entries.", tSameSign->GetName(), tSameSign->GetEntries());

    // create new file
    TFile *f_out = new TFile("Trees/_VetoInefficiencies/tSameSign.root","RECREATE");
    // open the file
    f_out->cd();
    // write the list and tree to this directory
    tSameSign->Write("tSameSign",TObject::kSingleKey);
    // list the contents of the file
    f_out->ls();
    // close the file
    f_out->Close();

    f_in->Close();

    return;
}

void SameSignEv_MakePlots()
{
    TFile *f_in = new TFile("Trees/_VetoInefficiencies/tSameSign.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *tSameSign = dynamic_cast<TTree*> (f_in->Get("tSameSign"));
    if(tSameSign) Printf("Input tree %s loaded.", tSameSign->GetName());

    ConnectTreeVariables_tNeutrons(tSameSign);

    TH1D *hZNA[2] = { NULL };
    TH1D *hZNC[2] = { NULL };
    TCanvas *c[2] = { NULL };
    TString name[2] = {"coh","inc"};
    Double_t fPtBoundaries[3] = {0.0, 0.2, 1.0};

    for(Int_t i = 0; i < 2; i++){
        hZNA[i] = new TH1D(("hZNA_" + name[i]).Data(),("hZNA_" + name[i]).Data(),nBins,n_low,n_upp);
        hZNC[i] = new TH1D(("hZNC_" + name[i]).Data(),("hZNC_" + name[i]).Data(),nBins,n_low,n_upp);
    }

    Printf("Tree %s filled with %lli entries.", tSameSign->GetName(), tSameSign->GetEntries());

    for(Int_t iEntry = 0; iEntry < tSameSign->GetEntries(); iEntry++){
        tSameSign->GetEntry(iEntry);
        for(Int_t i = 0; i < 2; i++){
            if(fPt > fPtBoundaries[i] && fPt <= fPtBoundaries[i+1]){
                if(fZNA_hit) {hZNA[i]->Fill(fZNA_n); Printf("Filling %.3f to hZNA_%s.", fZNA_n, name[i].Data());}
                if(fZNC_hit) {hZNC[i]->Fill(fZNC_n); Printf("Filling %.3f to hZNC_%s.", fZNC_n, name[i].Data());}
            }
        }
    }

    for(Int_t i = 0; i < 2; i++){
        c[i] = PlotNeutronDistribution(("c_" + name[i]).Data(),hZNA[i],hZNC[i],fPtBoundaries[i],fPtBoundaries[i+1]);
        c[i]->Draw();
        TString str = "Results/_VetoInefficiencies/SameSign_" + name[i];
        c[i]->Print((str + ".pdf").Data());
        c[i]->Print((str + ".png").Data());
    }

    for(Int_t i = 0; i < 2; i++) delete c[i];

    return;
}

void LowPt_PrepareTree()
{
    isPass3 = kTRUE;

    SetReducedRunList(kTRUE);

    TFile *f_in = new TFile("Trees/AnalysisData_pass3/AnalysisResults.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("AnalysisOutput/fTreeJpsi"));
    if(t_in) Printf("Input tree %s loaded.", t_in->GetName());

    ConnectTreeVariables(t_in);

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    gROOT->cd();
    TTree *t_out = new TTree("tLowPt","tLowPt");
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

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++)
    {
        t_in->GetEntry(iEntry);
        // m in 3.0 to 3.2 GeV/c^2, no pT cut
        if(EventPassed(1, -1) && fPt < 0.1)
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

    // create new file
    TFile *f_out = new TFile("Trees/_VetoInefficiencies/tLowPt.root","RECREATE");
    // open the file
    f_out->cd();
    // write the list and tree to this directory
    t_out->Write("tLowPt",TObject::kSingleKey);
    // list the contents of the file
    f_out->ls();
    // close the file
    f_out->Close();

    f_in->Close();

    return;
}

void LowPt_MakePlots()
{
    TFile *f_in = new TFile("Trees/_VetoInefficiencies/tLowPt.root","read");
    if(f_in) Printf("Input file %s loaded.", f_in->GetName());

    TTree *tLowPt = dynamic_cast<TTree*> (f_in->Get("tLowPt"));
    if(tLowPt) Printf("Input tree %s loaded.", tLowPt->GetName());

    ConnectTreeVariables_tNeutrons(tLowPt);

    TH1D *hZNA[3] = { NULL };
    TH1D *hZNC[3] = { NULL };
    TCanvas *c[3] = { NULL };
    Double_t fPtBoundaries[4] = {0.00, 0.02, 0.05, 0.10};

    for(Int_t i = 0; i < 3; i++){
        hZNA[i] = new TH1D(Form("hZNA%i",i+1),Form("hZNA%i",i+1),nBins,n_low,n_upp);
        hZNC[i] = new TH1D(Form("hZNC%i",i+1),Form("hZNC%i",i+1),nBins,n_low,n_upp);
    }

    Printf("Tree %s filled with %lli entries.", tLowPt->GetName(), tLowPt->GetEntries());

    for(Int_t iEntry = 0; iEntry < tLowPt->GetEntries(); iEntry++){
        tLowPt->GetEntry(iEntry);
        for(Int_t i = 0; i < 3; i++){
            if(fPt > fPtBoundaries[i] && fPt <= fPtBoundaries[i+1]){
                if(fZNA_hit) {hZNA[i]->Fill(fZNA_n);}
                if(fZNC_hit) {hZNC[i]->Fill(fZNC_n);}
            }
        }
    }

    for(Int_t i = 0; i < 3; i++){
        c[i] = PlotNeutronDistribution(Form("c%i",i+1),hZNA[i],hZNC[i],fPtBoundaries[i],fPtBoundaries[i+1]);
        c[i]->Draw();
        TString str = Form("Results/_VetoInefficiencies/LowPt_%i", i+1);
        c[i]->Print((str + ".pdf").Data());
        c[i]->Print((str + ".png").Data());
    }

    for(Int_t i = 0; i < 3; i++) delete c[i];

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

TCanvas* PlotNeutronDistribution(const char* name, TH1 *hZNA, TH1 *hZNC, Double_t fPtMin, Double_t fPtMax)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas(name,name,900,600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.11);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.04);
    //c->SetLogy();
    // X-axis
    hZNA->GetXaxis()->SetTitle("# of neutrons (ZN energy/2510 GeV)");
    hZNA->GetXaxis()->SetTitleSize(0.05);
    hZNA->GetXaxis()->SetLabelSize(0.05);
    // Y-axis
    hZNA->GetYaxis()->SetTitle("Counts per 0.2");
    hZNA->GetYaxis()->SetTitleSize(0.05);
    hZNA->GetYaxis()->SetLabelSize(0.05);
    hZNA->GetYaxis()->SetTitleOffset(1.0);
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
    TLegend *l1 = new TLegend(0.37,0.75,0.99,0.96);
    l1->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    l1->AddEntry((TObject*)0,"#it{m}_{#mu#mu} #in (3.0,3.2) GeV/#it{c}^{2}","");
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}",fPtMin,fPtMax),"");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();
    // legend2
    TLegend *l2 = new TLegend(0.53,0.60,0.95,0.72);
    l2->AddEntry(hZNA,Form("ZNA (total: %.0f events)", hZNA->Integral()),"L");
    l2->AddEntry(hZNC,Form("ZNC (total: %.0f events)", hZNC->Integral()),"L");
    l2->SetTextSize(0.05);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    return c;
}