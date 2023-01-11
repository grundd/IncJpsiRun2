// AxE_Dissociative.C
// David Grund, Jan 11, 2023

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

void TH1_SetStyle(TH1F* h, Color_t c, Int_t style = 1)
{
    h->SetLineColor(c);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
}

void PlotHistos(TH1F* h1, TString name, TString opt, TString yTitle, Bool_t log = kFALSE)
{
    TCanvas c("c","c",700,600);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.03);
    c.SetTopMargin(0.085);
    c.SetBottomMargin(0.12);
    if(log) c.SetLogy();
    TH1_SetStyle(h1,kBlue);
    // x-axis
    h1->GetXaxis()->SetTitle("#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})");
    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetDecimals(1);
    // y-axis
    if(log) {
        Float_t minimum = h1->GetMinimum();
        if(minimum < 0.1) minimum = 0.1;
        h1->GetYaxis()->SetRangeUser(0.5*minimum,1.5*h1->GetMaximum());
    }
    h1->GetYaxis()->SetTitle(yTitle.Data());
    h1->GetYaxis()->SetTitleOffset(1.2);
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetMaxDigits(3);
    // draw it
    h1->Draw(Form("%s",opt.Data()));
    c.Print(("Results/" + str_subfolder + "AxE_Dissociative/" + name + ".pdf").Data());
    return;
}

void ReweightIncPtShape()
{
    TFile *fGen = TFile::Open("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root","read");
    if(fGen) Printf("MC gen file loaded.");
    TTree *tGen = dynamic_cast<TTree*> (fGen->Get("AnalysisOutput/fTreeJpsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    ConnectTreeVariablesMCGen(tGen);
    TH1F* hGenOld = new TH1F("hGenOld","#it{N}_{gen}^{old}",260,0.0,1.2);
    TH1F* hGenNew = new TH1F("hGenNew","#it{N}_{gen}^{new}",260,0.0,1.2);

    // go over generated events
    Float_t N_gen_tot = tGen->GetEntries();
    Float_t N_gen_rap = 0;
    for(Int_t iEntry = 0; iEntry < N_gen_tot; iEntry++) 
    {
        tGen->GetEntry(iEntry);
        if(TMath::Abs(fYGen) < 1.0) { 
            N_gen_rap++;
            hGenOld->Fill(fPtGen);
        }
    }
    // now take 70% of hGenOld
    for(Int_t i = 0; i < 0.7*N_gen_rap; i++) hGenNew->Fill(hGenOld->GetRandom());
    // and add 30% from the H1 parametrization
    TF1 *fDissH1 = new TF1("fDissH1","x*pow((1 + x*x*[0]/[1]),-[1])",0.0,1.2);
    fDissH1->SetParameter(0,1.79);
    fDissH1->SetParameter(1,3.58);
    for(Int_t i = 0; i < 0.3*N_gen_rap; i++) hGenNew->Fill(fDissH1->GetRandom());
    // calculate the ratios
    TH1F* hRatios = (TH1F*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("#it{N}_{gen}^{new}/#it{N}_{gen}^{old}");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);
    TAxis *xAxis = hRatios->GetXaxis();
    // plot everything
    PlotHistos(hGenOld,"hGenOld","E0","#it{N}_{gen}^{old}");
    PlotHistos(hGenNew,"hGenNew","E0","#it{N}_{gen}^{new}");
    PlotHistos(hRatios,"hRatios","E0","#it{N}_{gen}^{new}/#it{N}_{gen}^{old}");

    TFile *fRec = TFile::Open("Trees/AnalysisDataMC_pass3/PIDCalibrated/AnalysisResults_MC_kIncohJpsiToMu.root","read");
    if(fRec) Printf("MC rec file loaded.");
    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJpsi"));
    if(tRec) Printf("MC rec tree loaded.");
    ConnectTreeVariablesMCRec(tRec);

    TH1F* hRecOld = new TH1F("hRecOld","#it{N}_{rec}^{old}",260,0.0,1.2);
    TH1F* hRecNew = new TH1F("hRecNew","#it{N}_{rec}^{new}",260,0.0,1.2);
    // go over reconstructed events and apply ratios
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) 
    {
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(1, 2))
        {
            // fill the original rec histogram
            hRecOld->Fill(fPt);
            // fill the new one using the ratios as weights
            // find index of the bin to which the current fPtGen corresponds
            Int_t iBinGen = xAxis->FindBin(fPtGen);
            // find index of the bin to which the current fPt corresponds
            Int_t iBinRec = xAxis->FindBin(fPt);
            // scale the J/psi entry by the ratio with the index iBinPsi2s
            Float_t ratio = hRatios->GetBinContent(iBinGen);
            // add the entry to hRecNew
            hRecNew->SetBinContent(iBinRec,hRecNew->GetBinContent(iBinRec)+ratio);
        }
    }
    PlotHistos(hRecOld,"hRecOld","E0","#it{N}_{rec}^{old}");
    PlotHistos(hRecNew,"hRecNew","E0","#it{N}_{rec}^{new}");

    return;
}

void AxE_Dissociative(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_Dissociative/");

    ReweightIncPtShape();

    return;
}