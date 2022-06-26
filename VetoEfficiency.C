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
#include "VetoEfficiency_Utilities.h"

// ******** options to set: **********
const Int_t nBinsPt = 5;
Double_t fBkgM_low = 1.8; // GeV
Double_t fBkgM_upp = 2.8; // GeV
// histograms in the number of neutrons
const Int_t nBins = 200;
Double_t n_low = 0.; // number of neutrons
Double_t n_upp = 50.;
// ***********************************
// tree variables:
Bool_t fZNA_hit, fZNC_hit;
Double_t fZNA_n, fZNC_n;
// normalization when calculating fractions:
Double_t NORM;

void VetoEfficiency_ClassifyEventsToChannels(Int_t mass_range, Bool_t normalized);
void VetoEfficiency_SubtractBackground();
void VetoEfficiency_Calculate();
void VetoEfficiency_PrepareTree();
TCanvas* PlotNeutronDistribution(const char* name, TH1 *hZNA, TH1 *hZNC, Double_t fPtMin, Double_t fPtMax, Double_t fMMin, Double_t fMMax);
void ConnectTreeVariables_tNeutrons(TTree *t);

void VetoEfficiency(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VetoEfficiency/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/");

    // prepare the tree containing information about mass, pT and ZN signal
    VetoEfficiency_PrepareTree();

    /*
    // classify events into classes: background (mass from 1.8 to 2.8 GeV)
    VetoEfficiency_ClassifyEventsToChannels(0, kFALSE);
    VetoEfficiency_ClassifyEventsToChannels(0, kTRUE);

    // classify events into classes: signal+bkg (mass from 3.0 to 3.2 GeV)
    VetoEfficiency_ClassifyEventsToChannels(1, kFALSE);
    VetoEfficiency_ClassifyEventsToChannels(1, kTRUE);
    */

    VetoEfficiency_SubtractBackground();

    VetoEfficiency_Calculate();

    return;
}

void VetoEfficiency_ClassifyEventsToChannels(Int_t mass_range, Bool_t normalized)
// mass_range == 0 => background (mass range: fBkgM_low to fBkgM_upp; only 1 pT bin)
//            == 1 => signal+bkg (mass range: 3.0 to 3.2 GeV; only 1 pT bin)
{
    NeutronMatrix *nEv = new NeutronMatrix();
    NeutronMatrix *nEv_PtBins[nBinsPt] = { NULL };
    for(Int_t i = 0; i < nBinsPt; i++) nEv_PtBins[i] = new NeutronMatrix();

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
        if(fZNA_hit){
            iBinN_A = 1;
            while(fZNA_n > fNumberOfN[iBinN_A]) iBinN_A++;
        } 
        if(fZNC_hit){
            iBinN_C = 1;
            while(fZNC_n > fNumberOfN[iBinN_C]) iBinN_C++;
        } 
        // find index of the pT bin
        Int_t iBinPt(0);
        while(fPt > ptBoundaries[iBinPt+1]) iBinPt++;

        nEv->AddEvent(iBinN_A,iBinN_C);
        nEv_PtBins[iBinPt]->AddEvent(iBinN_A,iBinN_C);

        // fill the histograms
        if(fZNA_hit) hZNA[0]->Fill(fZNA_n); 
        if(fZNC_hit) hZNC[0]->Fill(fZNC_n); 
        if(fZNA_hit) hZNA[iBinPt+1]->Fill(fZNA_n); 
        if(fZNC_hit) hZNC[iBinPt+1]->Fill(fZNC_n); 
    }

    f_in->Close();    

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
    if(normalized){
        // normalize by the total number of events
        nEv->Multiply(1/nEv->CountEvents_tot());
        for(Int_t i = 0; i < nBinsPt; i++) nEv_PtBins[i]->Multiply(1/nEv_PtBins[i]->CountEvents_tot());
        precision = 4;
    } 
    // total pT range
    if(!normalized) str_out = Form("Results/%sVetoEfficiency/%snEv_PtAll.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    else            str_out = Form("Results/%sVetoEfficiency/%snormalized_PtAll.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    nEv->PrintToFile(str_out, precision);
    // in pT bins
    for(Int_t i = 0; i < nBinsPt; i++){
        if(!normalized) str_out = Form("Results/%sVetoEfficiency/%snEv_PtBin%i.txt", str_subfolder.Data(), str_mass_subfolder.Data(), i+1);
        else            str_out = Form("Results/%sVetoEfficiency/%snormalized_PtBin%i.txt", str_subfolder.Data(), str_mass_subfolder.Data(), i+1);
        nEv_PtBins[i]->PrintToFile(str_out, precision);
    }

    // ##########################################################################################################
    
    for(Int_t i = 0; i < 6; i++) {delete c[i]; delete hZNA[i]; delete hZNC[i];}

    return;
}

void VetoEfficiency_SubtractBackground()
{
    // in full pT range
    // load fractions of bkg events
    NeutronMatrix *nEv_bkg = new NeutronMatrix();
    nEv_bkg->LoadFromFile("Results/" + str_subfolder + "VetoEfficiency/mass_1.80to2.80/normalized_PtAll.txt");
    Double_t nBkg = 205.; // from the invariant mass fit in allbins
    nEv_bkg->Multiply(nBkg);
    // first load all events (sig + bkg)
    NeutronMatrix *nEv_sig = new NeutronMatrix();
    nEv_sig->LoadFromFile("Results/" + str_subfolder + "VetoEfficiency/mass_3.00to3.20/nEv_PtAll.txt");
    // subtract background
    nEv_sig->SubtractMatrix(nEv_bkg);
    nEv_sig->PrintToConsole();
    Printf("Remaining number of events: %.2f", nEv_sig->CountEvents_tot());
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/");
    nEv_bkg->PrintToFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_bkg.txt",1);
    nEv_sig->PrintToFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_sig.txt",1);

    // in pT bins
    Double_t nBkg_PtBins[nBinsPt] = {26.,25.,40.,56.,57.}; // from the invariant mass fit in 5bins
    for(Int_t i = 1; i < nBinsPt+1; i++){
        NeutronMatrix *nEv_bkg = new NeutronMatrix();
        nEv_bkg->LoadFromFile(Form("Results/" + str_subfolder + "VetoEfficiency/mass_1.80to2.80/normalized_PtBin%i.txt", i));
        nEv_bkg->Multiply(nBkg_PtBins[i-1]);
        // first load all events (sig + bkg)
        NeutronMatrix *nEv_sig = new NeutronMatrix();
        nEv_sig->LoadFromFile(Form("Results/" + str_subfolder + "VetoEfficiency/mass_3.00to3.20/nEv_PtBin%i.txt", i));
        // subtract background
        nEv_sig->SubtractMatrix(nEv_bkg);
        nEv_sig->PrintToConsole();
        Printf("Remaining number of events: %.2f", nEv_sig->CountEvents_tot());
        gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/");
        nEv_bkg->PrintToFile(Form("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_PtBin%i_bkg.txt", i),1);
        nEv_sig->PrintToFile(Form("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_PtBin%i_sig.txt", i),1);    
    }

    return;
}

void VetoEfficiency_Calculate()
{
    // calculate efficiencies from inefficiencies
    for(Int_t i = 0; i < nBinsN+1; i++){
        fVetoEff_A[i] = 1. - fVetoIneff_A[i];
        fVetoEff_C[i] = 1. - fVetoIneff_C[i];
    }
    // in full pT range
    NeutronMatrix *nEv_sig = new NeutronMatrix();
    nEv_sig->LoadFromFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_sig.txt");
    Double_t nEv_uncorr = nEv_sig->CountEvents_tot();
    nEv_sig->ApplyEfficiencies();
    nEv_sig->PrintToFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_corrected.txt",1);
    Double_t nEv_corrected = nEv_sig->CountEvents_tot();
    Double_t fEff_total = nEv_uncorr / nEv_corrected;
    Printf("Total pile-up efficiency: %.3f", fEff_total);
    // in pT bins
    for(Int_t i = 1; i < nBinsPt+1; i++){
        NeutronMatrix *nEv_sig = new NeutronMatrix();
        nEv_sig->LoadFromFile(Form("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_PtBin%i_sig.txt", i));
        Double_t nEv_uncorr = nEv_sig->CountEvents_tot();
        nEv_sig->ApplyEfficiencies();
        nEv_sig->PrintToFile(Form("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_PtBin%i_corrected.txt",i),1);
        Double_t nEv_corrected = nEv_sig->CountEvents_tot();
        Double_t fEff_total = nEv_uncorr / nEv_corrected;
        Printf("Total pile-up efficiency in bin %i: %.3f", i, fEff_total);
    }

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