// IntegratedLuminosity.C
// David Grund, Mar 20, 2022
// To calculate the integrated luminosity of the analyzed samples using the trending files

// cpp headers
#include <fstream>
// root headers
#include "TFile.h"
#include "TSystem.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
// aliroot headers
#include "AliTriggerClass.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"

// Arrays containing the lists of good run numbers for LHC18qr
// see ListsOfGoodRuns.h
// particular lists selected in AnalysisConfig.h

// Arrays containing counts of fired triggers per run
// LHC18q:
// In increasing order (from 295585 to 296623)
vector<Int_t> Counts18q;
// LHC18r:
// In increasing order (from 296690 to 297595)
vector<Int_t> Counts18r;
// will be loaded from the folder Results/GetTriggerCounters/...

Int_t nRunsInList = 0;
vector<Int_t> RunList;
vector<Int_t> CountsList;

TString PeriodName[2] = {"18q", "18r"};

void CalculateLumi(Int_t period);
void SumLumi(Int_t period);
void SetLumiHisto(TH1D* h, Int_t period, Color_t color);

void IntegratedLuminosity(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "Lumi/");

    // Load trigger counters
    ifstream ifs;
    ifs.open(("Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18q.txt").Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < 123; i++){
            Int_t iCount;
            ifs >> iCount;
            Counts18q.push_back(iCount);
        }
        ifs.close(); 
    }
    ifs.open(("Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18r.txt").Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < 96; i++){
            Int_t iCount;
            ifs >> iCount;
            Counts18r.push_back(iCount);
        }
        ifs.close(); 
    }

    // LHC18q
    CalculateLumi(0);
    // LHC18r
    CalculateLumi(1);

    return;
}

void CalculateLumi(Int_t period)
{
    // Choose the period
    if(period == 0){ 
        // LHC18q
        RunList = runList_18q;
        CountsList = Counts18q;
        nRunsInList = nRuns_18q;
    } else if(period == 1){ 
        // LHC18r
        RunList = runList_18r;
        CountsList = Counts18r;
        nRunsInList = nRuns_18r;
    }

    const Int_t nRunsMax = 130;

    TString ClassName1 = "CCUP31-B-NOPF-CENTNOTRD"; // for run number < 295881
    TString ClassName2 = "CCUP31-B-SPD2-CENTNOTRD"; // for run number >= 295881

    // Load the trending file and tree
    TFile *fTrendFile = TFile::Open("Trees/Lumi/trending_merged_PbPb_2018.root", "read");
    if(fTrendFile) Printf("Input file %s loaded.", fTrendFile->GetName());

    TTree *fTree = dynamic_cast<TTree*> (fTrendFile->Get("trending"));
    if(fTree) Printf("Input tree %s loaded.", fTree->GetName());

    // Connect the branch addresses
    TObjArray* classes = new TObjArray();
    Double_t  lumi_seen[nRunsMax] = {0};
    Double_t  class_lumi[nRunsMax] = {0};
    Double_t  class_ds[nRunsMax] = {0};
    ULong64_t class_l2a[nRunsMax] = {0};
    Int_t run;
    Double_t mu = 0;
    fTree->SetBranchAddress("mu",&mu);
    fTree->SetBranchAddress("run",&run);
    fTree->SetBranchAddress("classes",&classes);
    fTree->SetBranchAddress("lumi_seen",&lumi_seen);
    fTree->SetBranchAddress("class_lumi",&class_lumi);
    fTree->SetBranchAddress("class_ds",&class_ds);
    fTree->SetBranchAddress("class_l2a",&class_l2a);
    fTree->BuildIndex("run");

    TH1D* hLumi = new TH1D("hLumi","",nRunsInList,0,nRunsInList); // Recorded luminosity per run
    TH1D* hLumiS = new TH1D("hLumiS","",nRunsInList,0,nRunsInList); // Seen luminosity per run (in my analysis)
    TH1D* hScale = new TH1D("hScale","",nRunsInList,0,nRunsInList); // Scale between seen and recorded lumi (<= 1.0) per run
    TH1D* hCCUP31ds = new TH1D("hCCUP31ds","CCUP31 downscaling",nRunsInList,0,nRunsInList); // Downscale of the trigger class per run

    // Calculate seen luminosity for a specified trigger class
    Int_t iBadScale = 0;

    for (Int_t i = 0; i < nRunsInList; i++){
        Int_t iRun = RunList[i];
        char* sRun = Form("%i",iRun); // Convert run number from int to char/string
        //Printf("%s %i %i",sRun, i, iRun);
        fTree->GetEntryWithIndex(iRun);

        // Check trigger class name
        AliTriggerClass* cl;
        if(iRun < 295881){
        cl = (AliTriggerClass*) classes->FindObject(ClassName1.Data());
        }
        if(iRun >= 295881){
        cl = (AliTriggerClass*) classes->FindObject(ClassName2.Data());
        }    
        if (!cl) continue;

        Int_t iClass = classes->IndexOf(cl);
        //Printf("%i %i %s",iRun, iClass, cl->GetName());
        Double_t l2a = (Double_t) class_l2a[iClass];
        //Printf("%.llu",class_l2a[iClass]);
        //Printf("%.10f",class_lumi[iClass]);
        Double_t scale = CountsList[i]/l2a;

        if(scale > 1.0){
            Printf("In run %i the scale is %.2f", iRun, scale);
            iBadScale++;
        }

        hScale->Fill(sRun,scale);
        hLumiS->Fill(sRun,scale*class_lumi[iClass]);
        hLumi->Fill(sRun,class_lumi[iClass]);
        hCCUP31ds->Fill(sRun,class_ds[iClass]);
    }

    // Write the results to the output root file:
    TFile* fOutputLumiHisto = new TFile("Results/" + str_subfolder + "Lumi/LumiHisto.root","recreate");
    hLumi->Write();
    hLumiS->Write();
    hScale->Write();
    hCCUP31ds->Write();
    fOutputLumiHisto->Close();

    if(iBadScale == 0){
        Printf("No scale factors above 1.00.");
    } else if(iBadScale > 0){
        Printf("%i suspicious runs.", iBadScale);
    }

    // Calculate the integrated luminosity (seen and official)
    SumLumi(period);
}

void SumLumi(Int_t period)
{
    // Here the total integrated luminosity is calculated
    // Open input file and read lumi per run
    TFile *fInputLumiHisto = TFile::Open("Results/" + str_subfolder + "Lumi/LumiHisto.root", "read");
    TH1D *hLumi = dynamic_cast<TH1D*> (fInputLumiHisto->Get("hLumi"));
    TH1D *hLumiS = dynamic_cast<TH1D*> (fInputLumiHisto->Get("hLumiS"));
    Double_t int_lumi_ana(0.), int_lumi_rec(0.);
    for(Int_t i(0); i < hLumiS->GetNbinsX(); i++){
        int_lumi_ana += hLumiS->GetBinContent(i+1);
        int_lumi_rec += hLumi->GetBinContent(i+1);
    }

    // Print the values to a text file
    TString str_outfile = "Results/" + str_subfolder + "Lumi/lumi_" + PeriodName[period] + ".txt";
    ofstream outfile(str_outfile.Data());
    outfile << int_lumi_ana << "\n" << int_lumi_rec;
    outfile.close();
    Printf("*** Results printed to %s.***", str_outfile.Data());

    // Also, we want to sum the luminosity corresponding to ClassName1 (CCUP31-B-NOPF-CENTNOTRD), i.e., 
    // to the situation without the past-future protection
    Double_t int_lumi_ana_NOPF(0.), int_lumi_rec_NOPF(0.);
    for(Int_t i(0); i < nRunsInList; i++){
        if(RunList[i] < 295881){
            int_lumi_ana_NOPF += hLumiS->GetBinContent(hLumiS->GetXaxis()->FindBin(Form("%i",RunList[i])));
            int_lumi_rec_NOPF += hLumi->GetBinContent(hLumi->GetXaxis()->FindBin(Form("%i",RunList[i])));
        }
    }
    TString str_outfile_NOPF = "Results/" + str_subfolder + "Lumi/lumi_NOPF_" + PeriodName[period] + ".txt";
    ofstream outfile_NOPF(str_outfile_NOPF.Data());
    outfile_NOPF << "Integrated luminosity corresponding to CCUP31-B-NOPF-CENTNOTRD (run number < 295881):\n";
    outfile_NOPF << Form("analyzed: %.3f (%.1f percent of the total lumi for this period) \n",int_lumi_ana_NOPF, int_lumi_ana_NOPF / int_lumi_ana * 100.);
    outfile_NOPF << Form("recorded: %.3f (%.1f percent of the total lumi for this period) \n",int_lumi_rec_NOPF, int_lumi_rec_NOPF / int_lumi_rec * 100.);
    outfile_NOPF.close();
    Printf("*** Results printed to %s.***", str_outfile_NOPF.Data());

    SetLumiHisto(hLumi, period, kBlue); // recorded lumi
    SetLumiHisto(hLumiS, period, kRed); // seen lumi

    TCanvas *cLumi = new TCanvas("cLumi","cLumi",1300,400);
    // Set margins
    cLumi->SetTopMargin(0.03);
    cLumi->SetRightMargin(0.01);
    //cLumi->SetLeftMargin(0.07);
    cLumi->SetBottomMargin(0.15);
    // Set histograms
    hLumi->GetYaxis()->SetTitle("L_{int} [#mub^{-1}]");
    hLumi->GetXaxis()->SetDecimals(1);
    hLumi->Draw();
    hLumiS->Draw("sameP0");
    // Legend
    TLegend* legLumi = new TLegend(0.295,0.74,0.52,0.96);
    legLumi->SetFillColor(kWhite);
    if(period == 0){
        legLumi->SetHeader("CCUP31 trigger class, LHC18q","l");
        cLumi->SetLeftMargin(0.07);
    } else if(period == 1){
        legLumi->SetHeader("CCUP31 trigger class, LHC18r","l");
        cLumi->SetLeftMargin(0.05);
    }
    legLumi->AddEntry(hLumi,Form("Total lumi rec.: %.3f #mub^{-1}",int_lumi_rec),"l");
    legLumi->AddEntry(hLumiS,Form("Total lumi ana.: %.3f #mub^{-1}",int_lumi_ana),"l");
    //legLumi->AddEntry((TObject*)0,Form("#bf{This thesis}"),"");
    legLumi->SetTextSize(0.055);
    legLumi->SetBorderSize(0);
    legLumi->SetFillStyle(0);
    legLumi->Draw();

    /*    TLegend *ltw = new TLegend(0.75,0.88,0.85,0.94);
    ltw->AddEntry((TObject*)0,"#bf{This work}","");
    ltw->SetMargin(0.);
    ltw->SetTextSize(0.055);
    ltw->SetBorderSize(0);
    ltw->SetFillStyle(0);
    ltw->Draw();
    */

    cLumi->SaveAs("Results/" + str_subfolder + "Lumi/lumi_ccup31_" + PeriodName[period] + ".pdf");
    cLumi->SaveAs("Results/" + str_subfolder + "Lumi/lumi_ccup31_" + PeriodName[period] + ".png");
}

void SetLumiHisto(TH1D* h, Int_t period, Color_t color)
{
    // A function to set the properties of the final histograms
    //gStyle->SetperiodStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    h->SetTitleFont(43);
    h->SetTitleSize(25);
    h->GetYaxis()->SetTitleFont(43);
    h->GetXaxis()->SetLabelFont(43);
    h->GetYaxis()->SetLabelFont(43);
    h->GetYaxis()->SetTitleSize(28);
    h->GetXaxis()->SetLabelSize(13);
    h->GetYaxis()->SetLabelSize(28);
    h->GetYaxis()->SetTickLength(0.01);
    if(period == 0){
        h->GetYaxis()->SetTitleOffset(1.50);
    } else if(period == 1){
        h->GetYaxis()->SetTitleOffset(1.00);
    }
    h->GetYaxis()->SetDecimals(1);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
    h->SetMarkerSize(0.7);
    h->SetMarkerStyle(kFullCross);
    //h->Labelsperiodion("v");
    h->SetMinimum(0);
    h->SetLineWidth(2);
    h->Sumw2(kFALSE);
}