// AxE_PtBins_Utilities.h
// David Grund, Jun 07, 2022

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"

TH1F* hRec = NULL; 
TH1F* hGen = NULL; 
TH1F* hAxE = NULL;
Float_t NRec_tot_val, NRec_tot_err;
Float_t NGen_tot_val, NGen_tot_err;

void SaveToFile(TString path, TH1F* h, Int_t prec = 0, Float_t N_tot_val = -1, Float_t N_tot_err = -1, Float_t fact = 1.) 
{
    ofstream of(path.Data());
    of << std::fixed << std::setprecision(prec);
    if(N_tot_val > 0) {
        of << "0\t" << fact * N_tot_val 
           << "\t" << fact * N_tot_err 
           << "\n";
    }
    for(Int_t iBin = 1; iBin <= h->GetNbinsX(); iBin++) {
        of << iBin << "\t" << fact * h->GetBinContent(iBin) 
           << "\t" << fact * h->GetBinError(iBin)
           << "\n";
    }
    of.close();
    Printf("*** File saved in %s.***", path.Data());
}

Float_t CalculateErrorBayes(Float_t k, Float_t n){ // k = NRec, n = NGen

    Float_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Float_t err = TMath::Sqrt(var);

    return err;
}

template <typename TH>
void SetTH1(TH* h, Color_t c, Marker_t m, Float_t mSize = 0.8, TString xTitle = "", Int_t lStyle = 1)
{
    gStyle->SetOptStat(0);
    h->SetBit(TH::kNoTitle);
    // line style
    h->SetLineColor(c);
    h->SetLineWidth(2);
    h->SetLineStyle(lStyle);
    // marker style
    h->SetMarkerColor(c);
    h->SetMarkerSize(mSize);
    h->SetMarkerStyle(m);
    // x-axis
    h->GetXaxis()->SetTitle(xTitle.Data());
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetDecimals(1);
    // y-axis
    h->GetYaxis()->SetTitle(h->GetTitle());
    h->GetYaxis()->SetTitleOffset(1.42);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetMaxDigits(3);
}

TCanvas* CreateCanvas(TString name, bool log = kFALSE)
{
    TCanvas* c = new TCanvas("c","c",700,600);
    c->SetLeftMargin(0.135);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.12);
    if(log) c->SetLogy();
    return c;
}

TLegend* CreateLegendAxE()
{
    TLegend *l = new TLegend(0.50,0.69,0.90,0.94);
    l->AddEntry((TObject*)0,Form("ALICE Simulation"),"");
    l->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),"");
    l->AddEntry((TObject*)0,Form("2.2 < #it{m} < 4.5 GeV/#it{c}^{2}"),"");
    l->SetTextSize(0.045);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    return l;
}

void PlotHistos(TString folder, TString name, TString opt, Bool_t log, Float_t mSize, TH1F* h1, TH1F* h2 = NULL, TLegend* l = NULL)
{
    TCanvas* c = CreateCanvas(name,log);
    SetTH1<TH1F>(h1,kBlue,kFullCircle,mSize,"#it{p}_{T} (GeV/#it{c})");
    if(log) {
        Float_t minimum = h1->GetMinimum();
        if(minimum < 0.1) minimum = 0.1;
        h1->GetYaxis()->SetRangeUser(0.5*minimum,1.5*h1->GetMaximum());
    }
    if(h2) {
        h1->SetBit(TH1::kNoStats);
        h2->SetBit(TH1::kNoStats);
        SetTH1<TH1F>(h2,kRed,kFullCircle,mSize);
    }
    // draw it
    h1->Draw(opt.Data());
    if(h2) h2->Draw(Form("%s SAME",opt.Data()));
    // legend
    if(l) l->Draw();
    c->Print("Results/" + str_subfolder + folder + name + ".pdf");
    delete c;
    return;
}

TH1F* GetRatioHisto()
{
    TFile* f = TFile::Open("Results/" + str_subfolder + "AxE_Dissociative/incJpsi/ratios.root");
    if(f) Printf("File %s loaded.", f->GetName()); 
    TList *l = (TList*) f->Get("HistList");
    if(l) Printf("List %s loaded.", l->GetName()); 
    TH1F* h = (TH1F*)l->FindObject("hRatios");
    if(h) Printf("Histo %s loaded.", h->GetName());
    return h;
}

void AxE_PtBins_FillHistNRec(Bool_t reWeight, Float_t fCutZ)
{
    // check if the corresponding text file already exists
    TString sOut;
    if(fCutZ == cut_fVertexZ) sOut = "Results/" + str_subfolder + "AxE_PtBins/";
    else                      sOut = "Results/" + str_subfolder + Form("VertexZ_SystUncertainties/Zcut%.1f_AxE_PtBins/", fCutZ);
    if(reWeight) sOut += "reweighted_";
    sOut.Append(Form("NRec_%ibins.txt", nPtBins));

    ifstream ifs;
    ifs.open(sOut);
    if(!ifs.fail())
    {
        // this configuration has already been calculated
        Printf("*** The file %s already exists. ***", sOut.Data());
        // fill hRec with data from the file
        Int_t i_bin; Float_t NRec_val, NRec_err;
        // fiducial
        ifs >> i_bin >> NRec_tot_val >> NRec_tot_err;
        // pT bins
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++) {
            ifs >> i_bin >> NRec_val >> NRec_err;
            hRec->SetBinContent(iBin, NRec_val);
        }
        ifs.close(); 
        return;
    } 
    else 
    {
        // this configuration is yet to be calculated
        Printf("*** Calculating N rec per bin for %s... ***", sOut.Data());
        TFile *fRec = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fRec) Printf("MC rec file loaded.");
        TTree *tRec = dynamic_cast<TTree*> (fRec->Get(str_in_MC_tree_rec.Data()));
        if(tRec) Printf("MC rec tree loaded.");
        ConnectTreeVariablesMCRec(tRec);

        // |> *********** for VertexZ_SystUncertainties.C ***********
        // save the original value of cut_fVertexZ
        Printf("Original cut on vertex Z: %.1f", cut_fVertexZ);
        Float_t fCutZ_orig = cut_fVertexZ;
        if(fCutZ != cut_fVertexZ) {
            // set the new value of cut_fVertexZ
            cut_fVertexZ = fCutZ;
            Printf("New cut on vertex Z: %.1f", cut_fVertexZ);
        }
        // <| *****************************************************

        // load the ratio to re-weight the spectra
        TH1F* hRatios = GetRatioHisto();
        TAxis* xAxis = hRatios->GetXaxis();
        // go over tree entries and calculate NRec in the total range and in bins
        Float_t NRec_tot(0.);
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) {
            tRec->GetEntry(iEntry);
            Int_t iBinGen = xAxis->FindBin(fPtGen);
            Float_t weight = 1.0; 
            if(reWeight) weight = hRatios->GetBinContent(iBinGen);
            // m between 2.2 and 4.5 GeV/c^2
            // & pT from 0.2 to 1.0 GeV/c
            if(EventPassedMCRec(0, 3)) {
                NRec_tot += weight;
                hRec->Fill(fPt, weight);
            }
        }
        Printf("*** Finished! ***");

        // |> *********** for VertexZ_SystUncertainties.C ***********
        if(cut_fVertexZ != fCutZ_orig)
        {
            // set back the original value of cut_fVertexZ
            cut_fVertexZ = fCutZ_orig;
            Printf("Restoring the original cut on vertex Z: %.1f", cut_fVertexZ);  
        }
        // <| *****************************************************
        
        NRec_tot_val = NRec_tot;
        NRec_tot_err = TMath::Sqrt(NRec_tot);
        SaveToFile(sOut,hRec,0,NRec_tot_val,NRec_tot_err);
        return;
    }
}

void AxE_PtBins_FillHistNGen(Bool_t reWeight)
{
    // check if the corresponding text file already exists
    TString sOut = "Results/" + str_subfolder + "AxE_PtBins/";
    if(reWeight) sOut += "reweighted_";
    sOut.Append(Form("NGen_%ibins.txt", nPtBins));

    ifstream ifs;
    ifs.open(sOut);
    if(!ifs.fail())
    {
        // this configuration has already been calculated
        Printf("*** The file %s already exists. ***", sOut.Data());
        // fill hGen with data from the text file
        Int_t i_bin; Float_t NGen_val, NGen_err;
        // fiducial
        ifs >> i_bin >> NGen_tot_val >> NGen_tot_err;
        // pT bins
        for(Int_t iBin = 1; iBin <= nPtBins; iBin++) {
            ifs >> i_bin >> NGen_val >> NGen_err;
            hGen->SetBinContent(iBin, NGen_val);
        }
        ifs.close(); 
        return;
    } 
    else 
    {
        // this configuration is yet to be calculated
        Printf("*** Calculating N gen per bin for %s... ***", sOut.Data());
        TFile *fGen = TFile::Open((str_in_MC_fldr_gen + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fGen) Printf("MC gen file loaded.");
        TTree *tGen = dynamic_cast<TTree*> (fGen->Get(str_in_MC_tree_gen.Data()));
        if(tGen) Printf("MC gen tree loaded.");
        ConnectTreeVariablesMCGen(tGen);

        // load the ratio to re-weight the spectra
        TH1F* hRatios = GetRatioHisto();
        TAxis* xAxis = hRatios->GetXaxis();
        // go over tree entries and calculate NRec in the total range and in bins
        Float_t NGen_tot(0.);
        for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++) {
            tGen->GetEntry(iEntry);
            Int_t iBinGen = xAxis->FindBin(fPtGen);
            Float_t weight = 1.0; 
            if(reWeight) weight = hRatios->GetBinContent(iBinGen);
            // m between 2.2 and 4.5 GeV/c^2
            // & pT from 0.2 to 1.0 GeV/c
            if(EventPassedMCGen(3)) {
                NGen_tot += weight;
                hGen->Fill(fPtGen, weight);
            }
        }
        Printf("*** Finished! ***");
        
        NGen_tot_val = NGen_tot;
        NGen_tot_err = TMath::Sqrt(NGen_tot);
        SaveToFile(sOut,hGen,0,NGen_tot_val,NGen_tot_err);
        return;
    }
    return;
}

void AxE_PtBins_Calculate(Bool_t reWeight, Float_t fCutZ)
{
    hRec = new TH1F("hRec","N_{rec} per bin",nPtBins,ptBoundaries);
    hGen = new TH1F("hGen","N_{gen} per bin",nPtBins,ptBoundaries);

    AxE_PtBins_FillHistNRec(reWeight,fCutZ);
    AxE_PtBins_FillHistNGen(reWeight);

    hAxE = (TH1F*)hRec->Clone("hAxE");
    hAxE->SetTitle("#it{N}_{rec}^{MC}/#it{N}_{gen}^{MC}");
    hAxE->Sumw2();
    hAxE->Divide(hGen);
    hAxE->SetBit(TH1::kNoTitle);
    hAxE->SetBit(TH1::kNoStats);
    hAxE->GetYaxis()->SetRangeUser(0.015,0.035);
    TLegend *l = CreateLegendAxE();
    // save the figures and print the results to txt file
    TString folder(""), name("");
    if(fCutZ == cut_fVertexZ) folder = "AxE_PtBins/";
    else                      folder = Form("VertexZ_SystUncertainties/Zcut%.1f_AxE_PtBins/",fCutZ);
    if(reWeight) name += "reweighted_";
    name += Form("AxE_%ibins",nPtBins);
    PlotHistos(folder,name,"E0",kFALSE,0.8,hAxE,NULL,l);

    // compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kTRUE;
    if(DebugErrors){
        Float_t ErrRoot = 0;
        Float_t ErrBayes = 0;
        for(Int_t i = 1; i <= nPtBins; i++) {
            ErrRoot = hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(hRec->GetBinContent(i),hGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // calculate the total value of AxE
    Float_t AxE_tot_val = NRec_tot_val / NGen_tot_val;
    Float_t AxE_tot_err = CalculateErrorBayes(NRec_tot_val, NGen_tot_val);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxE_tot_val*100., AxE_tot_err*100.);

    // print the results to a text file
    SaveToFile("Results/" + str_subfolder + folder + name + ".txt",hAxE,3,AxE_tot_val,AxE_tot_err,100.);
    delete hRec;
    delete hGen;
    delete hAxE;
    return;
}