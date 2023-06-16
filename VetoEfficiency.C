// VetoEfficiency.C
// David Grund, June 18, 2022

// root headers
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
// my headers
#include "VetoEfficiency_Utilities.h"

// ******** options to set: **********
const Int_t nBinsPt = 5;
Double_t fBkgM_low = 1.8; // GeV
Double_t fBkgM_upp = 2.8; // GeV
// histograms in the number of neutrons
const Int_t nBins = 200;
Double_t n_low = 0.; // number of neutrons
Double_t n_upp = 50.;

void VetoEff_ClassifyEvents(Int_t mass_range, Bool_t normalized);
void VetoEff_SubtractBkg();
Double_t VetoEff_Calculate(Int_t iEff, Bool_t SystUncr = kFALSE);
// iEff == 0 => weight both the A and C side (XnXn)
//      == 1 => weight as 0nXn and XnYn(whatever) (combined1)
//      == 2 => weight as Xn0n and (whatever)YnXn (combined2)
void VetoEff_SystUncertainty();
TCanvas* PlotNeutronDistribution(const char* name, TH1 *hZNA, TH1 *hZNC, Double_t fPtMin, Double_t fPtMax, Double_t fMMin, Double_t fMMax);
TCanvas* Plot2DNeutronDistribution(const char* name, TH2 *hZN, Double_t fPtMin, Double_t fPtMax, Double_t fMMin, Double_t fMMax);

void VetoEfficiency(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "VetoEfficiency/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/");

    // prepare the tree containing information about mass, pT and ZN signal
    VetoEfficiency_PrepareTree();

    if(kTRUE)
    {
        // classify events into classes: background (mass from 1.8 to 2.8 GeV)
        VetoEff_ClassifyEvents(0, kFALSE);
        VetoEff_ClassifyEvents(0, kTRUE);
        // classify events into classes: signal+bkg (mass from 3.0 to 3.2 GeV)
        VetoEff_ClassifyEvents(1, kFALSE);
        VetoEff_ClassifyEvents(1, kTRUE);
    }

    // calculate partial (in)efficiencies in neutron bins
    VetoEff_CalcEfficiencies();
    // subtract bkg from signal in 3.0 to 3.2 GeV in mass
    VetoEff_SubtractBkg();
    // calculate the total efficiency
    VetoEff_Calculate(0,kFALSE);
    VetoEff_Calculate(1,kFALSE);
    VetoEff_Calculate(2,kFALSE);
    // calculate systematic uncertainties
    if(kTRUE) VetoEff_SystUncertainty();

    return;
}

void VetoEff_ClassifyEvents(Int_t mass_range, Bool_t normalized)
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
    // 0 = pT from 0.2 to 1.0 GeV/c, then pT bins
    TH2D *hZN[6] = { NULL };
    TH2D *hZN_hits[6] = { NULL };
    TH1D *hZNA[6] = { NULL }; 
    TH1D *hZNC[6] = { NULL };
    TCanvas *c2d[6] = { NULL };
    TCanvas *c2d_hits[6] = { NULL };
    TCanvas *c[6] = { NULL };

    for(Int_t i = 0; i < 6; i++){
        hZN[i] = new TH2D(Form("hZN%i",i),Form("hZN%i",i),48,-4.,8.,48,-4.,8.);
        hZN_hits[i] = new TH2D(Form("hZNhits%i",i),Form("hZN%i",i),48,-4.,8.,48,-4.,8.);
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
        if(fZNA_hit || fZNC_hit) hZN_hits[0]->Fill(fZNA_n*2.510,fZNC_n*2.510);
        hZN[0]->Fill(fZNA_n*2.510,fZNC_n*2.510);
        if(fZNA_hit) hZNA[0]->Fill(fZNA_n); 
        if(fZNC_hit) hZNC[0]->Fill(fZNC_n); 
        if(fZNA_hit || fZNC_hit) hZN_hits[iBinPt+1]->Fill(fZNA_n*2.510,fZNC_n*2.510);
        hZN[iBinPt+1]->Fill(fZNA_n*2.510,fZNC_n*2.510);
        if(fZNA_hit) hZNA[iBinPt+1]->Fill(fZNA_n); 
        if(fZNC_hit) hZNC[iBinPt+1]->Fill(fZNC_n); 
    }

    f_in->Close();    

    // ##########################################################################################################
    // plot neutron distribution in allbins
    // 1d
    c[0] = PlotNeutronDistribution("c0",hZNA[0],hZNC[0],0.2,1.0,m_low,m_upp);
    c[0]->Draw();
    TString str_out = "Results/" + str_subfolder + "VetoEfficiency/" + str_mass_subfolder + "ZN_n_all";
    c[0]->Print((str_out + ".pdf").Data());
    // 2d
    // at least one ZN hit
    c2d_hits[0] = Plot2DNeutronDistribution("c2d0",hZN_hits[0],0.2,1.0,m_low,m_upp);
    c2d_hits[0]->Draw();
    str_out = "Results/" + str_subfolder + "VetoEfficiency/" + str_mass_subfolder + "ZN_2dHits_n_all";
    c2d_hits[0]->Print((str_out + ".pdf").Data());
    // everything
    c2d[0] = Plot2DNeutronDistribution("c2d0",hZN[0],0.2,1.0,m_low,m_upp);
    c2d[0]->Draw();
    str_out = "Results/" + str_subfolder + "VetoEfficiency/" + str_mass_subfolder + "ZN_2d_n_all";
    c2d[0]->Print((str_out + ".pdf").Data());
    // plots neutron distribution in bins
    for(Int_t i = 1; i < nPtBins+1; i++){
        // 1d
        c[i] = PlotNeutronDistribution(Form("c%i",i),hZNA[i],hZNC[i],ptBoundaries[i-1],ptBoundaries[i],m_low,m_upp);
        c[i]->Draw();
        str_out = Form("Results/%sVetoEfficiency/%sZN_n_bin%i", str_subfolder.Data(), str_mass_subfolder.Data(), i);
        c[i]->Print((str_out + ".pdf").Data());
        // 2d
        // at least one ZN hit
        c2d_hits[i] = Plot2DNeutronDistribution(Form("c2d%i",i),hZN_hits[i],ptBoundaries[i-1],ptBoundaries[i],m_low,m_upp);
        c2d_hits[i]->Draw();
        str_out = Form("Results/%sVetoEfficiency/%sZN_2dHits_n_bin%i", str_subfolder.Data(), str_mass_subfolder.Data(), i);
        c2d_hits[i]->Print((str_out + ".pdf").Data()); 
        // everything
        c2d[i] = Plot2DNeutronDistribution(Form("c2d%i",i),hZN[i],ptBoundaries[i-1],ptBoundaries[i],m_low,m_upp);
        c2d[i]->Draw();
        str_out = Form("Results/%sVetoEfficiency/%sZN_2d_n_bin%i", str_subfolder.Data(), str_mass_subfolder.Data(), i);
        c2d[i]->Print((str_out + ".pdf").Data());
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
    if(!normalized) str_out = Form("Results/%sVetoEfficiency/%snEv_all.txt", str_subfolder.Data(), str_mass_subfolder.Data());
    else            str_out = Form("Results/%sVetoEfficiency/%snormalized_all.txt", str_subfolder.Data(), str_mass_subfolder.Data());
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

void VetoEff_SubtractBkg()
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/");
    // in full pT range
    // load fractions of bkg events
    NeutronMatrix *nEv_bkg = new NeutronMatrix();
    nEv_bkg->LoadFromFile("Results/" + str_subfolder + "VetoEfficiency/mass_1.80to2.80/normalized_all.txt");
    Double_t nBkg = VetoEffiency_LoadBkg(0); // from the invariant mass fit in allbins
    nEv_bkg->Multiply(nBkg);
    // first load all events (sig + bkg)
    NeutronMatrix *nEv_sig = new NeutronMatrix();
    nEv_sig->LoadFromFile("Results/" + str_subfolder + "VetoEfficiency/mass_3.00to3.20/nEv_all.txt");
    nEv_sig->Plot("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_all.pdf");
    // subtract background
    nEv_sig->SubtractMatrix(nEv_bkg);
    nEv_sig->PrintToConsole();
    Printf("Remaining number of events: %.2f", nEv_sig->CountEvents_tot());
    nEv_bkg->PrintToFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_bkg.txt",1);
    nEv_bkg->Plot("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_bkg.pdf");
    nEv_sig->PrintToFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_sig.txt",1);
    nEv_sig->Plot("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_sig.pdf");
    // in pT bins
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/PtBins/");
    NeutronMatrix *nEv_bkg_bins[5] = { NULL };
    NeutronMatrix *nEv_sig_bins[5] = { NULL };
    for(Int_t i = 0; i < nPtBins; i++)
    {
        // load fractions of bkg events
        nEv_bkg_bins[i] = new NeutronMatrix();
        nEv_bkg_bins[i]->LoadFromFile("Results/" + str_subfolder + Form("VetoEfficiency/mass_1.80to2.80/normalized_PtBin%i.txt",i+1));
        Double_t nBkgPtBin = VetoEffiency_LoadBkg(i+1); // from the invariant mass fit in allbins
        nEv_bkg_bins[i]->Multiply(nBkgPtBin);
        // first load all events (sig + bkg)
        nEv_sig_bins[i] = new NeutronMatrix();
        nEv_sig_bins[i]->LoadFromFile("Results/" + str_subfolder + Form("VetoEfficiency/mass_3.00to3.20/nEv_PtBin%i.txt",i+1));
        nEv_sig_bins[i]->Plot("Results/" + str_subfolder + Form("VetoEfficiency/bkg_subtracted/PtBins/nEv_bin%i.pdf",i+1));
        // subtract background
        nEv_sig_bins[i]->SubtractMatrix(nEv_bkg_bins[i]);
        nEv_sig_bins[i]->PrintToConsole();
        Printf("Remaining number of events: %.2f", nEv_sig_bins[i]->CountEvents_tot());
        //nEv_bkg_bins[i]->PrintToFile("Results/" + str_subfolder + Form("VetoEfficiency/bkg_subtracted/PtBins/nEv_bkg_bin%i.txt",i+1),1);
        nEv_bkg_bins[i]->Plot("Results/" + str_subfolder + Form("VetoEfficiency/bkg_subtracted/PtBins/nEv_bkg_bin%i.pdf",i+1));
        //nEv_sig_bins[i]->PrintToFile("Results/" + str_subfolder + Form("VetoEfficiency/bkg_subtracted/PtBins/nEv_sig_bin%i.txt",i+1),1);
        nEv_sig_bins[i]->Plot("Results/" + str_subfolder + Form("VetoEfficiency/bkg_subtracted/PtBins/nEv_sig_bin%i.pdf",i+1));
    }

    return;
}

Double_t VetoEff_Calculate(Int_t iEff, Bool_t SystUncr)
{
    VetoEff_SetEfficiencies(SystUncr);
    // in full pT range
    NeutronMatrix *nEv_sig = new NeutronMatrix();
    nEv_sig->LoadFromFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_sig.txt");
    Double_t nEv_uncorr = nEv_sig->CountEvents_tot();
    TString name = "";
    switch(iEff)
    {
        case 0:
            nEv_sig->ApplyEfficiencies_AC();
            name = "XnXn";
            break;
        case 1:
            nEv_sig->ApplyEfficiencies_combined1();
            name = "XnYn";
            break;
        case 2:
            nEv_sig->ApplyEfficiencies_combined2();
            name = "YnXn";
            break;
    }
    // calculate the veto eff
    Double_t nEv_corr = nEv_sig->CountEvents_tot();
    Double_t fEff_total = nEv_uncorr / nEv_corr;
    Printf("nEv uncorr: %.1f corr: %.1f", nEv_uncorr, nEv_corr);
    Printf("Total pile-up efficiency: %.3f", fEff_total);
    // if not the syst uncr calculation
    if(!SystUncr){
        // save the matrix containing corrected nEv
        nEv_sig->PrintToFile("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_corr_" + name + ".txt",1);
        nEv_sig->Plot("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/nEv_corr_" + name + ".pdf");
        // print the result to a text file
        ofstream outfile;
        outfile.open("Results/" + str_subfolder + "VetoEfficiency/bkg_subtracted/VetoEff_" + name + ".txt");
        outfile << std::fixed << std::setprecision(3);
        outfile << fEff_total;
        outfile.close();
    }

    return fEff_total;
}

void VetoEff_SystUncertainty()
{
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "VetoEfficiency/SystUncertainty/");

    for(Int_t i = 0; i < nBinsN; i++){
        hSampledEffPartial_A[i] = new TH1D(Form("hSampledEffPartial_A%i",i+1),Form("hSampledEffPartial_A%i",i+1),100,0.,1.);
        hSampledEffPartial_C[i] = new TH1D(Form("hSampledEffPartial_C%i",i+1),Form("hSampledEffPartial_C%i",i+1),100,0.,1.);
        hSampledEffPartial[i] = new TH1D(Form("hSampledEffPartial%i",i+1),Form("hSampledEffPartial%i",i+1),100,0.,1.);
    } 

    for(Int_t i = 0; i < 1e4; i++)
    {
        Double_t fEff = VetoEff_Calculate(1,kTRUE);
        hSampledEffTotal->Fill(fEff);
    }

    TCanvas *cA[5] = { NULL };
    TCanvas *cC[5] = { NULL };
    TCanvas *c[5] = { NULL };
    for(Int_t i = 0; i < nBinsN; i++){
        cA[i] = new TCanvas(Form("cA%i", i+1), Form("cA%i", i+1), 900, 600);
        cC[i] = new TCanvas(Form("cC%i", i+1), Form("cC%i", i+1), 900, 600);
        c[i] = new TCanvas(Form("c%i", i+1), Form("c%i", i+1), 900, 600);
        cA[i]->cd();
        hSampledEffPartial_A[i]->Draw();
        cC[i]->cd();
        hSampledEffPartial_C[i]->Draw();
        c[i]->cd();
        hSampledEffPartial[i]->Draw();
        // print the canvases
        cA[i]->Print("Results/" + str_subfolder + Form("VetoEfficiency/SystUncertainty/hSampledEff_A%i.pdf",i+1));
        cC[i]->Print("Results/" + str_subfolder + Form("VetoEfficiency/SystUncertainty/hSampledEff_C%i.pdf",i+1));
        c[i]->Print("Results/" + str_subfolder + Form("VetoEfficiency/SystUncertainty/hSampledEff%i.pdf",i+1));
    }
    // fit the gaussian peak
    TF1 *fGauss = new TF1("fGauss", "gaus", 0.0, 1.0);
    hSampledEffTotal->Fit(fGauss);
    hSampledEffTotal->SetStats(0);
    fGauss->SetLineWidth(3);
    fGauss->SetLineColor(kRed);
    // plot the results
    TCanvas *cTotal = new TCanvas("cTotal","cTotal",900,800);
    cTotal->cd();
    cTotal->SetTopMargin(0.06);
    cTotal->SetBottomMargin(0.12);
    cTotal->SetRightMargin(0.045);
    cTotal->SetLeftMargin(0.14);

    hSampledEffTotal->SetTitle(";#varepsilon^{veto}_{diss} (-);Counts");
    hSampledEffTotal->SetLineWidth(3);
    hSampledEffTotal->SetLineColor(kBlue);
    // Vertical axis
    hSampledEffTotal->GetYaxis()->SetTitleSize(0.05);
    hSampledEffTotal->GetYaxis()->SetTitleOffset(1.3);
    hSampledEffTotal->GetYaxis()->SetLabelSize(0.05);
    hSampledEffTotal->GetYaxis()->SetMaxDigits(3);
    hSampledEffTotal->GetYaxis()->SetDecimals(1);
    // Horizontal axis
    hSampledEffTotal->GetXaxis()->SetTitleSize(0.05);
    hSampledEffTotal->GetXaxis()->SetTitleOffset(1.1);
    hSampledEffTotal->GetXaxis()->SetLabelOffset(0.01);
    hSampledEffTotal->GetXaxis()->SetLabelSize(0.05);
    hSampledEffTotal->GetXaxis()->SetRangeUser(0.5,0.8);
    hSampledEffTotal->GetXaxis()->SetDecimals(2);
    hSampledEffTotal->Draw();
    fGauss->Draw("SAME");

    TLegend *l = new TLegend(0.51,0.74,0.95,0.92);
    l->AddEntry((TObject*)0,"gaussian fit:","");
    l->AddEntry((TObject*)0,Form("#mu = (%.2f #pm %.2f)%%", fGauss->GetParameter(1) * 1e2, fGauss->GetParError(1) * 1e2),"");
    l->AddEntry((TObject*)0,Form("#sigma = (%.2f #pm %.2f)%%", fGauss->GetParameter(2) * 1e2, fGauss->GetParError(2) * 1e2),"");
    l->SetTextSize(0.045);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    cTotal->Print("Results/" + str_subfolder + "VetoEfficiency/SystUncertainty/hSampledEffTotal.pdf");

    TLegend *ltw = new TLegend(0.22,0.83,0.38,0.92);
    ltw->AddEntry((TObject*)0,"#bf{This work}","");
    ltw->SetMargin(0.);
    ltw->SetTextSize(0.05);
    ltw->SetBorderSize(0);
    ltw->SetFillStyle(0);
    ltw->Draw();

    cTotal->Print("Results/" + str_subfolder + "_rozprava/eff_veto_syst.pdf");

    return;
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

TCanvas* Plot2DNeutronDistribution(const char* name, TH2 *hZN, Double_t fPtMin, Double_t fPtMax, Double_t fMMin, Double_t fMMax)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas(name,name,700,700);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.11);
    c->SetLeftMargin(0.09);
    c->SetRightMargin(0.11);
    c->SetGrid();
    // X-axis
    hZN->GetXaxis()->SetTitle("ZNA energy (TeV)");
    hZN->GetXaxis()->SetTitleSize(0.05);
    hZN->GetXaxis()->SetLabelSize(0.05);
    // Y-axis
    hZN->GetYaxis()->SetTitle("ZNC energy (TeV)");
    hZN->GetYaxis()->SetTitleSize(0.05);
    hZN->GetYaxis()->SetLabelSize(0.05);
    hZN->GetYaxis()->SetTitleOffset(0.7);
    // Z-axis
    hZN->GetZaxis()->SetLabelSize(0.042);
    hZN->GetZaxis()->SetDecimals(1);
    // Draw
    hZN->Draw("COLZ");
    // Legend 1
    TLegend *l1 = new TLegend(0.12,0.93,0.6,1.0);
    l1->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();
    // Legend 2
    TLegend *l2 = new TLegend(0.40,0.15,0.90,0.25);
    l2->AddEntry((TObject*)0,Form("#it{m}_{#mu#mu} #in (%.2f,%.2f) GeV/#it{c}^{2}",fMMin,fMMax),"");
    l2->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}",fPtMin,fPtMax),"");
    l2->SetTextSize(0.042);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->SetMargin(0.);
    l2->Draw();
    return c;
}