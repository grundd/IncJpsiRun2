// _FiducialCrossSec.C
// David Grund, Jan 05, 2023

// cpp headers
#include <fstream>
// root headers
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TSystem.h"
#include "TStyle.h"

Int_t iCount(0);

void TH1_SetStyle(TH1F* h, Color_t c, Int_t style = 1)
{
    h->SetLineColor(c);
    h->SetLineWidth(3);
    h->SetLineStyle(style);
}

void TF1_SetStyle(TF1 *f, Color_t c, Int_t style = 1)
{
    f->SetLineColor(c);
    f->SetLineWidth(3);
    f->SetLineStyle(style);
}

TF1* PtShape(TString name, Float_t par0, Float_t par1)
{
    TF1* f = new TF1(name.Data(),"[0] * abs(x) * exp(-[1] * pow(x,2))",0.,2.); // x = pT
    f->SetParameter(0,par0);
    f->SetParameter(1,par1);
    return f;
}

TF1* Pt2Shape(TString name, Float_t par0, Float_t par1)
{
    TF1* f = new TF1(name.Data(),"[0] * exp(-[1] * x)",0.,4.); // x = pT^2
    f->SetParameter(0,par0);
    f->SetParameter(1,par1);
    return f;
}

void PlotHistos(TString name, TString opt, Bool_t log, TH1F* h1, TF1* f1 = NULL)
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
    h1->GetYaxis()->SetTitle("value");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetMaxDigits(3);
    // draw it
    h1->Draw(Form("%s",opt.Data()));
    if(f1) {
        TF1_SetStyle(f1,kGreen+1,9);
        f1->Draw("SAME"); 
        TLegend* l = new TLegend(0.50,0.78,0.80,0.92);
        l->AddEntry(f1,"fit: #it{N} exp(#minus#it{b}#it{p}_{T}^{2})","L");
        l->AddEntry((TObject*)0,Form("#it{N} = %.0f #pm %.0f",f1->GetParameter(0),f1->GetParError(0)),"");
        l->AddEntry((TObject*)0,Form("#it{b} = %.3f #pm %.3f",f1->GetParameter(1),f1->GetParError(1)),"");
        l->SetTextSize(0.03);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetMargin(0.25);
        l->Draw();
    }
    c.Print(("Results/_FiducialCrossSec/" + name + ".pdf").Data());
    return;
}

void CalculateAxE(TString name, TH1F* hRec, TH1F* hGen)
{
    TH1F* hAxE = (TH1F*)hRec->Clone("hAxE");
    hAxE->SetTitle("Acc #times eff: #it{N}_{rec}/#it{N}_{gen}");
    hAxE->Sumw2();
    hAxE->Divide(hGen);
    PlotHistos(name,"E0",kFALSE,hAxE);
    return;
}

void PlotFilledHist(TString subfolder, TTree* t, TH1F* h)
{
    t->Draw(Form("fPt>>%s",h->GetName()));
    PlotHistos(Form("%s%s",subfolder.Data(),h->GetName()),"E0",kFALSE,h);
    return;
}

void GenerateEvents(TString subfolder, Float_t par_rec[], Float_t par_gen[], Float_t N_data, Bool_t useH1 = kTRUE)
{
    TF1* fSTARlight = Pt2Shape("fSTARlight",1.,1.);
    TString s = "Results/_FiducialCrossSec/" + subfolder + "generatedEvs.root";
    TFile* f = new TFile(s.Data(),"RECREATE");
    Float_t fPt;
    // rec
    TTree* tRec = new TTree("tRec","");
    tRec->Branch("fPt", &fPt, "fPt/F"); 
    fSTARlight->SetParameter(1,par_rec[1]);
    for(Int_t i = 0; i < par_rec[0]; i++) { 
        fPt = 0;
        while(fPt < 0.04 || fPt > 1.00) fPt = fSTARlight->GetRandom(); 
        tRec->Fill();
    }
    // gen
    TTree* tGen = new TTree("tGen","");
    tGen->Branch("fPt", &fPt, "fPt/F"); 
    fSTARlight->SetParameter(1,par_gen[1]);
    for(Int_t i = 0; i < par_gen[0]; i++) { 
        fPt = 0;
        while(fPt < 0.04 || fPt > 1.00) fPt = fSTARlight->GetRandom(); 
        tGen->Fill(); 
    }
    // data
    Float_t N_data_SL = N_data;
    Float_t N_data_H1 = 0;
    if(useH1) {
        N_data_SL = N_data * 0.7;
        N_data_H1 = N_data * 0.3;
    }
    TTree* tData = new TTree("tData","");
    tData->Branch("fPt", &fPt, "fPt/F"); 
    // first mimic inc from STARlight
    // according to the pT from our analysis, this is ~ 281/(281+121) ~ 70% of events
    fSTARlight->SetParameter(1,par_rec[1]);
    for(Int_t i = 0; i < N_data_SL; i++) { 
        fPt = 0;
        while(fPt < 0.04 || fPt > 1.00) fPt = fSTARlight->GetRandom(); 
        tData->Fill(); 
    }
    // now mimic inc dissociative using the H1 shape
    TF1 *fDissH1 = new TF1("fDissH1","x*pow((1 + x*x*[0]/[1]),-[1])",0.,4.);
    fDissH1->SetParameter(0,1.79);
    fDissH1->SetParameter(1,3.58);
    for(Int_t i = 0; i < N_data_H1; i++) { 
        fPt = 0;
        while(fPt < 0.04 || fPt > 1.00) fPt = fDissH1->GetRandom(); 
        tData->Fill();
    }
    // plot the histograms
    TH1F* hRnd_rec = new TH1F("hRnd_rec","simulated #it{N}_{rec} vs #it{p}_{T}^{2}",100,0.04,1.);
    PlotFilledHist(subfolder,tRec,hRnd_rec);
    TH1F* hRnd_gen = new TH1F("hRnd_gen","simulated #it{N}_{gen} vs #it{p}_{T}^{2}",100,0.04,1.);
    PlotFilledHist(subfolder,tGen,hRnd_gen);
    TH1F* hRnd_data = new TH1F("hRnd_data","simulated pseudo-data vs #it{p}_{T}^{2}",100,0.04,1.);
    PlotFilledHist(subfolder,tData,hRnd_data);
    // AxE
    CalculateAxE(subfolder + "hRnd_AxE",hRnd_rec,hRnd_gen);
    f->Write("",TObject::kWriteDelete);
    f->Close();
    return;
}

void TryBinning(TString subfolder, Int_t nBins, Float_t &fiducialDir, Float_t &fiducialInt, Float_t* fBins = NULL)
{
    TString s = "Results/_FiducialCrossSec/" + subfolder + "generatedEvs.root";
    TFile* f = new TFile(s.Data(),"read");
    TTree *tRec = dynamic_cast<TTree*> (f->Get("tRec"));
    TTree *tGen = dynamic_cast<TTree*> (f->Get("tGen"));
    TTree *tData = dynamic_cast<TTree*> (f->Get("tData"));

    TH1F* hRec = NULL;
    TH1F* hGen = NULL;
    TH1F* hData = NULL;
    TH1F* hCS = NULL;
    if(fBins == NULL) {
        cout << "Using uniform binning.\n";
        hRec = new TH1F(Form("%02i_hRec",iCount),"",nBins,0.04,1.00);
        hGen = new TH1F(Form("%02i_hGen",iCount),"",nBins,0.04,1.00);
        hData = new TH1F(Form("%02i_hData",iCount),"",nBins,0.04,1.00);
        hCS = new TH1F(Form("%02i_hCS",iCount),"",nBins,0.04,1.00);
    } else {
        cout << "Using variable bins sizes.\n";
        hRec = new TH1F(Form("%02i_hRec",iCount),"",nBins,fBins);
        hGen = new TH1F(Form("%02i_hGen",iCount),"",nBins,fBins);
        hData = new TH1F(Form("%02i_hData",iCount),"",nBins,fBins);
        hCS = new TH1F(Form("%02i_hCS",iCount),"",nBins,fBins);
    }
    hRec->SetTitle("#it{p}_{T}^{2} dist of simulated #it{N}_{rec}");
    hGen->SetTitle("#it{p}_{T}^{2} dist of simulated #it{N}_{gen}");
    hData->SetTitle("#it{p}_{T}^{2} dist of simulated pseudo-data");
    hCS->SetTitle("#it{p}_{T}^{2} dist of simulated d#sigma_{#gammaPb}/d|#it{t}|");
    PlotFilledHist(subfolder,tRec,hRec);
    PlotFilledHist(subfolder,tGen,hGen);
    PlotFilledHist(subfolder,tData,hData);

    // AxE
    CalculateAxE(Form("%s%02i_hAxE",subfolder.Data(),iCount),hRec,hGen);
    // cross section in |t| bins
    Float_t num_fact = 4.44e-4;
    for(Int_t i = 1; i <= nBins; i++) {
        hCS->SetBinContent(i, num_fact 
            * hData->GetBinContent(i) // N_inc
            * hGen->GetBinContent(i) // N_gen
            / hRec->GetBinContent(i) // N_rec
            / (hCS->GetBinLowEdge(i+1) - hCS->GetBinLowEdge(i))); // Delta|t|
        hCS->SetBinError(i, hCS->GetBinContent(i)
            * TMath::Sqrt(hData->GetBinContent(i))
            / hData->GetBinContent(i));
    }
    PlotHistos(Form("%s%02i_hCS",subfolder.Data(),iCount),"E0",kTRUE,hCS);
    // cross section (per the full interval of |t|) - direct calculation
    Float_t CS = num_fact * hData->Integral() * hGen->Integral() / hRec->Integral() / (hCS->GetBinLowEdge(nBins+1) - hCS->GetBinLowEdge(1));
    // fiducial cross section
    // directly
    fiducialDir = CS * (hCS->GetBinLowEdge(nBins+1) - hCS->GetBinLowEdge(1));
    // integral
    fiducialInt = 0;
    for(Int_t i = 1; i <= nBins; i++) {
        fiducialInt += hCS->GetBinContent(i) * (hCS->GetBinLowEdge(i+1) - hCS->GetBinLowEdge(i));
    }
    // compare:
    ofstream outfile(Form("Results/_FiducialCrossSec/%s%02i_log.txt",subfolder.Data(),iCount));
    outfile << "differential cross section dsigma_gPb/d|t| (" << nBins << " bins):\n"
            << "range \tN_rec \tN_gen \tN_data \tcross sec. (mub)\n"
            << "full\t" 
            << std::fixed << std::setprecision(0)
            << hRec->Integral() << "\t"
            << hGen->Integral() << "\t"
            << hData->Integral() << "\t"
            << std::fixed << std::setprecision(2)
            << CS << "\n";
    for(Int_t i = 1; i <= nBins; i++) {
        outfile << "bin " << i << "\t"
                << std::fixed << std::setprecision(0)
                << hRec->GetBinContent(i) << "\t"
                << hGen->GetBinContent(i) << "\t"
                << hData->GetBinContent(i) << "\t"
                << std::fixed << std::setprecision(2)
                << hCS->GetBinContent(i) << "\n";
    }            
    outfile << "\nfiducial cross section sigma_gPb:\n"
            << "directly: " << fiducialDir << " mub\n"
            << "integral: " << fiducialInt << " mub\n";
    outfile.close();
    iCount++; // counter of tries
    return;
}

void TryConfiguration(Bool_t useH1, Float_t b_rec, Float_t b_gen)
{
    TString subfolder = "";
    if(useH1) subfolder += "H1_";
    else      subfolder += "noH1_";
    subfolder += Form("bRec%.1f_",b_rec);
    subfolder += Form("bGen%.1f/",b_gen);
    gSystem->Exec("mkdir -p Results/_FiducialCrossSec/" + subfolder);
    
    // generate new events according to the fits
    // pars: N, b
    Float_t par_rec[2] = {81000, b_rec};
    Float_t par_gen[2] = {3254000, b_gen};
    Float_t N_data = 405;
    GenerateEvents(subfolder,par_rec,par_gen,N_data,useH1);
    
    // try various binnings
    Float_t fiducialDir[11] = { 0 };
    Float_t fiducialInt[11] = { 0 };
    // binning from the original analysis
    Float_t bins[6] = {0.04, 0.08, 0.15, 0.26, 0.48, 1.00};
    TryBinning(subfolder,5,fiducialDir[0],fiducialInt[0],bins);
    // new binnings: 1 to 10 bins in total
    for(Int_t i = 1; i <= 10; i++) TryBinning(subfolder,i,fiducialDir[i],fiducialInt[i]);

    ofstream outfile("Results/_FiducialCrossSec/" + subfolder + "summary.txt");
    outfile << "binning\tdirect \tintegral\n" << std::fixed << std::setprecision(3)
            << "orig.\t" << fiducialDir[0] << "\t" << fiducialInt[0] << "\n";
    for(Int_t i = 1; i <= 10; i++) {
        outfile << i << " bins\t" << fiducialDir[i] << "\t" << fiducialInt[i] << "\n";
    }
    outfile.close();
    return;
}

void _FiducialCrossSec()
{
    // arguments: useH1, b_rec, b_gen
    // b_rec and b_gen from the fits
    TryConfiguration(kTRUE,5.3,4.0);
    TryConfiguration(kFALSE,5.3,4.0);
    // b_rec same as b_gen
    TryConfiguration(kTRUE,4.0,4.0);
    TryConfiguration(kFALSE,4.0,4.0);
    return;
}