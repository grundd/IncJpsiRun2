// AxE_Dissociative.C
// David Grund, Jan 11, 2023

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
// roofit
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
// my headers
#include "AnalysisManager.h"
#include "AnalysisConfig.h"
#include "SetPtBinning.h"
#include "AxE_PtBins_Utilities.h"

using namespace RooFit;

const Int_t nBins(240);
Float_t fPtLow(0.0);
Float_t fPtUpp(1.2);
Float_t perc_gen_H1 = 57.;

void TH1_SetStyle(TH1F* h, Color_t c, Int_t style = 1)
{
    h->SetLineColor(c);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
}

void PlotHistos(TString opt, Bool_t log, TH1F* h1, TH1F* h2 = NULL)
{
    TCanvas c("c","c",700,600);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.03);
    c.SetTopMargin(0.085);
    c.SetBottomMargin(0.12);
    if(log) c.SetLogy();
    TH1_SetStyle(h1,kBlue);
    if(h2) TH1_SetStyle(h2,kRed);
    // x-axis
    h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetLabelOffset(0.01);
    h1->GetXaxis()->SetDecimals(1);
    // y-axis
    if(log) {
        Float_t minimum = h1->GetMinimum();
        if(minimum < 0.1) minimum = 0.1;
        h1->GetYaxis()->SetRangeUser(0.5*minimum,1.5*h1->GetMaximum());
    }
    h1->GetYaxis()->SetTitle(h1->GetTitle());
    h1->GetYaxis()->SetTitleOffset(1.2);
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetLabelOffset(0.01);
    h1->GetYaxis()->SetMaxDigits(3);
    // draw it
    h1->Draw(Form("%s",opt.Data()));
    if(h2) h2->Draw(Form("%s SAME",opt.Data()));
    // legend
    // (...)
    c.Print(("Results/" + str_subfolder + "AxE_Dissociative/" + h1->GetName() + ".pdf").Data());
    return;
}

void ReweightIncPtShape()
{
    TFile *fGen = TFile::Open("Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root","read");
    if(fGen) Printf("MC gen file loaded.");
    TTree *tGen = dynamic_cast<TTree*> (fGen->Get("AnalysisOutput/fTreeJpsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    ConnectTreeVariablesMCGen(tGen);
    TH1F* hGenOld = new TH1F("hGenOld","#it{N}_{gen}^{old}",nBins,fPtLow,fPtUpp);
    TH1F* hGenNew = new TH1F("hGenNew","#it{N}_{gen}^{new}",nBins,fPtLow,fPtUpp);

    // go over generated events
    Float_t NGen_tot = tGen->GetEntries();
    Float_t NGen_rap = 0;
    for(Int_t iEntry = 0; iEntry < NGen_tot; iEntry++) 
    {
        tGen->GetEntry(iEntry);
        if(TMath::Abs(fYGen) < 1.0) { 
            NGen_rap++;
            hGenOld->Fill(fPtGen);
        }
    }
    // now take 70% of hGenOld
    Float_t howMany_SL = NGen_rap * (100. - perc_gen_H1) / 100.;
    for(Int_t i = 0; i < howMany_SL; i++) hGenNew->Fill(hGenOld->GetRandom());
    // and add 30% from the H1 parametrization
    TF1 *fH1 = new TF1("fH1","x*pow((1 + x*x*[0]/[1]),-[1])",fPtLow,fPtUpp);
    fH1->SetParameter(0,1.79);
    fH1->SetParameter(1,3.58);
    Float_t howMany_H1 = NGen_rap * (perc_gen_H1) / 100.;
    for(Int_t i = 0; i < howMany_H1; i++) hGenNew->Fill(fH1->GetRandom());
    // calculate the ratios
    TH1F* hRatios = (TH1F*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("#it{N}_{gen}^{new}/#it{N}_{gen}^{old}");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);
    TAxis *xAxis = hRatios->GetXaxis();
    // plot everything
    PlotHistos("E0",kFALSE,hGenOld,hGenNew);
    PlotHistos("E0",kFALSE,hRatios);

    TFile *fRec = TFile::Open("Trees/AnalysisDataMC_pass3/PIDCalibrated/AnalysisResults_MC_kIncohJpsiToMu.root","read");
    if(fRec) Printf("MC rec file loaded.");
    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJpsi"));
    if(tRec) Printf("MC rec tree loaded.");
    ConnectTreeVariablesMCRec(tRec);

    TH1F* hRecOld = new TH1F("hRecOld","#it{N}_{rec}^{old}",nBins,fPtLow,fPtUpp);
    TH1F* hRecNew = new TH1F("hRecNew","#it{N}_{rec}^{new}",nBins,fPtLow,fPtUpp);
    TH1F* hRecOld_fit = new TH1F("hRecOld_fit","#it{N}_{rec}^{old}",60,fPtLow,fPtUpp);
    TH1F* hRecNew_fit = new TH1F("hRecNew_fit","#it{N}_{rec}^{new}",60,fPtLow,fPtUpp);
    // go over reconstructed events and apply ratios
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) 
    {
        tRec->GetEntry(iEntry);
        // 3.0 < m < 3.2 GeV/c^2, pT < 2 GeV/c
        if(EventPassedMCRec(1, 2))
        {
            // fill the original rec histogram
            hRecOld->Fill(fPt);
            hRecOld_fit->Fill(fPt);
            // fill the new one using the ratios as weights
            Int_t iBinGen = xAxis->FindBin(fPtGen);
            Float_t weight = hRatios->GetBinContent(iBinGen);
            hRecNew->Fill(fPt,weight);
            hRecNew_fit->Fill(fPt,weight);
        }
    }
    PlotHistos("E0",kFALSE,hRecOld,hRecNew);
    PlotHistos("E0",kFALSE,hRecOld_fit,hRecNew_fit);

    // fit hNRecNew to get the fraction of events coming from H1
    RooRealVar vPt("vPt", "vPt", fPtLow, fPtUpp);
    RooDataHist dhData("dhData","dhData",vPt,hRecNew_fit);
    RooRealVar vb_SL("vb_SL","vb_SL",4.,1.,10.);
    RooRealVar vb_H1("vb_H1","vb_H1",1.79,1.,10.);
    RooRealVar vn_H1("vn_H1","vn_H1",3.58,1.,10.);
    vb_H1.setConstant(kTRUE);
    vn_H1.setConstant(kTRUE);
    RooGenericPdf pdfSL("pdfSL","","vPt*exp(-vb_SL*pow(vPt,2))",RooArgSet(vPt,vb_SL));
    RooGenericPdf pdfH1("pdfH1","","vPt*pow((1 + pow(vPt,2)*vb_H1/vn_H1),-vn_H1)",RooArgSet(vPt,vb_H1,vn_H1));
    RooDataHist dhSL("dhSL","dhSL",vPt,hRecOld_fit);
    RooHistPdf  hpdfSL("hpdfSL","hpdfSL",vPt,dhSL,0);
    Float_t NRec_tot = hRecNew->Integral();
    RooRealVar vN_SL("vN_SL","vN_SL",NRec_tot*0.7,NRec_tot*0.1,NRec_tot*1.);
    RooRealVar vN_H1("vN_H1","vN_H1",NRec_tot*0.3,NRec_tot*0.1,NRec_tot*1.);
    RooAddPdf* CombinedPDF = NULL;
    // fit with hpdfSL instead of pdfSL?
    Bool_t useHPdf = kTRUE;
    if(useHPdf) CombinedPDF = new RooAddPdf("CombinedPDF","", RooArgList(hpdfSL,pdfH1),RooArgList(vN_SL,vN_H1));
    else        CombinedPDF = new RooAddPdf("CombinedPDF","", RooArgList(pdfSL,pdfH1),RooArgList(vN_SL,vN_H1));

    RooFitResult* fResFit = CombinedPDF->fitTo(dhData,SumW2Error(kFALSE),Extended(kTRUE),Save());
    TCanvas c("c","c",700,600);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.03);
    c.SetTopMargin(0.085);
    c.SetBottomMargin(0.12);
    RooPlot* fPlot = vPt.frame(Title("fit of the #it{p}_{T} distribution of #it{N}_{rec}^{new}")); 
    dhData.plotOn(fPlot,Name("dhData"),MarkerStyle(kFullCircle),MarkerSize(0.8),MarkerColor(kBlack),LineColor(kBlack),LineWidth(2));
    if(useHPdf) CombinedPDF->plotOn(fPlot,Name("hpdfSL"),Components(hpdfSL),Range(0.001,fPtUpp),LineColor(kRed),LineWidth(3),LineStyle(2));
    else        CombinedPDF->plotOn(fPlot,Name("pdfSL"),Components(pdfSL),Range(0.001,fPtUpp),LineColor(kRed),LineWidth(3),LineStyle(2));
    CombinedPDF->plotOn(fPlot,Name("pdfH1"),Components(pdfH1),Range(0.001,fPtUpp),LineColor(kViolet),LineWidth(3),LineStyle(2));
    CombinedPDF->plotOn(fPlot,Name("CombinedPDF"),Range(0.001,fPtUpp),LineColor(kBlue),LineWidth(3),LineStyle(9));
    Float_t chi2 = fPlot->chiSquare("CombinedPDF","dhData",fResFit->floatParsFinal().getSize());
    // x-axis
    fPlot->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fPlot->GetXaxis()->SetTitleOffset(1.1);
    fPlot->GetXaxis()->SetTitleSize(0.045);
    fPlot->GetXaxis()->SetLabelSize(0.045);
    fPlot->GetXaxis()->SetLabelOffset(0.01);
    fPlot->GetXaxis()->SetDecimals(1);
    // y-axis
    fPlot->GetYaxis()->SetRangeUser(0.0,hRecNew_fit->GetMaximum()*1.05);
    fPlot->GetYaxis()->SetTitle("counts per 15 MeV/#it{c}");
    fPlot->GetYaxis()->SetTitleOffset(1.2);
    fPlot->GetYaxis()->SetTitleSize(0.045);
    fPlot->GetYaxis()->SetLabelSize(0.045);
    fPlot->GetYaxis()->SetLabelOffset(0.01);
    fPlot->GetYaxis()->SetMaxDigits(3);
    fPlot->Draw();
    // legend
    vPt.setRange("rPt_all",fPtLow,fPtUpp);
    RooAbsReal* arN_SL = NULL;
    if(useHPdf) arN_SL = hpdfSL.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
    else        arN_SL = pdfSL.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
    RooAbsReal* arN_H1 = pdfH1.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
    Float_t N_SL = arN_SL->getVal()*vN_SL.getVal();
    Float_t N_H1 = arN_H1->getVal()*vN_H1.getVal();
    Float_t perc_SL = 100. * N_SL / NRec_tot;
    Float_t perc_H1 = 100. * N_H1 / NRec_tot;
    Int_t nRows(4);
    TLegend l(0.55,0.90-nRows*0.05,0.90,0.90);
    l.AddEntry((TObject*)0,Form("#chi^{2}: %.3f",chi2),"");
    l.AddEntry((TObject*)0,Form("total #it{N}_{rec} = %.0f",NRec_tot),"");
    l.AddEntry((TObject*)0,Form("#it{N}_{rec}^{SL} = %.0f (%.1f%%)",N_SL,perc_SL),"");
    l.AddEntry((TObject*)0,Form("#it{N}_{rec}^{H1} = %.0f (%.1f%%)",N_H1,perc_H1),"");
    l.SetTextSize(0.045);
    l.SetBorderSize(0);
    l.SetFillStyle(0);
    l.Draw();
    if(useHPdf) c.Print("Results/" + str_subfolder + "AxE_Dissociative/fit_pdf.pdf");
    else        c.Print("Results/" + str_subfolder + "AxE_Dissociative/fit_hPdf.pdf");

    // AxE in pT bins
    TH1F* hNGen = new TH1F("hNGen","#it{N}_{gen}",nPtBins,ptBoundaries);
    TH1F* hNRec = new TH1F("hNRec","#it{N}_{rec}",nPtBins,ptBoundaries);
    // go over generated events
    Float_t NGen_all = 0;
    for(Int_t iEntry = 0; iEntry < NGen_tot; iEntry++) 
    {
        tGen->GetEntry(iEntry);
        // |y| < 0.8, 0.2 < pT < 1.0 GeV/c
        if(EventPassedMCGen(3)) 
        {
            Int_t iBinGen = xAxis->FindBin(fPtGen);
            Float_t weight = hRatios->GetBinContent(iBinGen);
            NGen_all += weight;
            hNGen->Fill(fPtGen,weight);
        }
    }
    PlotHistos("E0",kFALSE,hNGen);
    AxE_PtBins_SaveToFile(NGen_all,hNGen,"Results/" + str_subfolder + Form("AxE_Dissociative/NGen_%ibins.txt",nPtBins));
    // go over reconstructed events
    Float_t NRec_all = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++) 
    {
        tRec->GetEntry(iEntry);
        // 2.2 < m < 4.5 GeV/c^2, 0.2 < pT 1.0 GeV/c
        if(EventPassedMCRec(0, 3))
        {
            Int_t iBinGen = xAxis->FindBin(fPtGen);
            Float_t weight = hRatios->GetBinContent(iBinGen);
            NRec_all += weight;
            hNRec->Fill(fPt,weight);
        }
    }
    PlotHistos("E0",kFALSE,hNRec);
    AxE_PtBins_SaveToFile(NRec_all,hNRec,"Results/" + str_subfolder + Form("AxE_Dissociative/NRec_%ibins.txt",nPtBins));
    // calculate AxE
    TH1F* hAxE = (TH1F*)hNRec->Clone("hAxE");
    hAxE->SetTitle("(Acc#times#varepsilon)_{MC} = #it{N}_{rec}/#it{N}_{gen}");
    hAxE->Sumw2();
    hAxE->Divide(hNGen);
    PlotHistos("E0",kFALSE,hAxE);

    return;
}

void AxE_Dissociative(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_Dissociative/");

    ReweightIncPtShape();

    return;
}