// AxE_Dissociative.cxx
// David Grund, Jan 11, 2023

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
#include "SetPtBinning_PtFit.h"
#include "AxE_Utilities.h"

using namespace RooFit;

const Int_t nBins(280);
Float_t fPtLow(0.0);
Float_t fPtUpp(1.4);
Float_t perc_gen_H1;

void ReweightIncPtShape(TString process = "IncJ")
{    
    TString prefix = "";
    TString sGen = "Trees/AnalysisDataMC_pass3/";
    TString sRec = "Trees/AnalysisDataMC_pass3/PIDCalibrated/";
    TString sOut = "Results/" + str_subfolder + "AxE_Dissociative/";
    if(process == "IncJ") {
        perc_gen_H1 = 61.2;
        prefix = "incJpsi/";
        sGen += "AnalysisResults_MC_kIncohJpsiToMu.root";
        sRec += "AnalysisResults_MC_kIncohJpsiToMu.root";
    } else if(process == "IncP") {
        perc_gen_H1 = 49.8;
        prefix = "incPsi2s/";
        sGen += "AnalysisResults_MC_kIncohPsi2sToMuPi.root";
        sRec += "AnalysisResults_MC_kIncohPsi2sToMuPi.root";
    } else return;

    // check if the ratio file already exists
    sOut += prefix;
    gSystem->Exec("mkdir -p " + sOut);
    sOut += "ratios.root";
    Bool_t fileExists = !gSystem->AccessPathName(sOut.Data());
    if(fileExists) {
        cout << "Re-weighting already performed. Skipping... \n";
        return;
    } else {
        cout << "Re-weighting will be performed: \n";

        TFile *fGen = TFile::Open(sGen.Data(),"read");
        if(fGen) Printf("MC gen file loaded.");
        TTree *tGen = dynamic_cast<TTree*> (fGen->Get("AnalysisOutput/fTreeJpsiMCGen"));
        if(tGen) Printf("MC gen tree loaded.");
        ConnectTreeVariablesMCGen(tGen);
        TH1F* hGenOld = new TH1F("hGenOld","#it{N}^{gen}_{MC}",nBins,fPtLow,fPtUpp);
        TH1F* hGenNew = new TH1F("hGenNew","#it{N}^{gen}_{MC}",nBins,fPtLow,fPtUpp);
        TH1F* hH1_fit = new TH1F("hH1_fit","#it{N}^{rec}_{H1}",60,fPtLow,fPtUpp);

        // go over generated events
        Float_t NGen_tot = tGen->GetEntries();
        for(Int_t iEntry = 0; iEntry < NGen_tot; iEntry++) 
        {
            tGen->GetEntry(iEntry);
            if(TMath::Abs(fYGen) < 1.0) hGenOld->Fill(fPtGen);
        }
        // now take 70% of hGenOld
        Float_t NGen_old = hGenOld->Integral();
        cout << "NGen old events: " << NGen_old << "\n";
        Float_t howMany_SL = NGen_old * (100. - perc_gen_H1) / 100.;
        cout << "NGen SL events: " << howMany_SL << "\n";
        for(Int_t i = 0; i < howMany_SL; i++) hGenNew->Fill(hGenOld->GetRandom());
        // and add 30% from the H1 parametrization
        TF1 *fH1 = new TF1("fH1","x*pow((1 + x*x*[0]/[1]),-[1])",fPtLow,fPtUpp);
        fH1->SetParameter(0,1.79);
        fH1->SetParameter(1,3.58);
        Float_t howMany_H1 = NGen_old * (perc_gen_H1) / 100.;
        cout << "NGen H1 events: " << howMany_H1 << "\n";
        for(Int_t i = 0; i < howMany_H1; i++) {
            Float_t fPtH1 = fH1->GetRandom();
            hGenNew->Fill(fPtH1);
            hH1_fit->Fill(fPtH1);
        }
    
        // calculate the ratios
        TH1F* hRatios = (TH1F*)hGenNew->Clone("hRatios");
        hRatios->SetTitle("#it{R} = (#it{N}^{gen}_{MC})_{new}/(#it{N}^{gen}_{MC})_{old}");
        hRatios->Sumw2();
        hRatios->Divide(hGenOld);
        TAxis *xAxis = hRatios->GetXaxis();
        // plot everything
        int nRowsLeg = 2;
        bool tw = false;
        if(tw) nRowsLeg++;
        TLegend lGen(0.47,0.95-nRowsLeg*0.05,0.95,0.95);
        lGen.AddEntry(hGenOld,Form("original dist.: %.0f ev.",hGenOld->Integral()),"L"); // (#it{N}^{gen}_{MC})_{old}
        lGen.AddEntry(hGenNew,Form("re-weighted dist.: %.0f ev.",hGenNew->Integral()),"L"); // (#it{N}^{gen}_{MC})_{new}
        if(tw) lGen.AddEntry((TObject*)0,"#bf{This work}","");
        lGen.SetMargin(0.175);
        lGen.SetTextSize(0.040);
        lGen.SetBorderSize(0);
        lGen.SetFillStyle(0);
        lGen.Draw();

        PlotHistos("AxE_Dissociative/",prefix+"hGen_oldNew","E0",kFALSE,0.,hGenOld,hGenNew,&lGen);
        PlotHistos("AxE_Dissociative/",prefix+"ratios","E0",kFALSE,0.,hRatios);
        // save the ratio file
        TFile* fRat = new TFile(sOut.Data(),"RECREATE");
        TList *lRat = new TList();
        lRat->Add(hRatios);
        lRat->Write("HistList", TObject::kSingleKey);
        lRat->ls();
        fRat->ls();
        fRat->Close();

        TFile *fRec = TFile::Open(sRec.Data(),"read");
        if(fRec) Printf("MC rec file loaded.");
        TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJpsi"));
        if(tRec) Printf("MC rec tree loaded.");
        ConnectTreeVariablesMCRec(tRec);
        TH1F* hRecOld = new TH1F("hRecOld","#it{N}^{rec}_{MC}",nBins,fPtLow,fPtUpp);
        TH1F* hRecNew = new TH1F("hRecNew","#it{N}^{rec}_{MC}",nBins,fPtLow,fPtUpp);
        TH1F* hRecOld_fit = new TH1F("hRecOld_fit","#it{N}^{rec}_{MC}",60,fPtLow,fPtUpp);
        TH1F* hRecNew_fit = new TH1F("hRecNew_fit","#it{N}^{rec}_{MC}",60,fPtLow,fPtUpp);
        TH1D* hRec_ptFit = new TH1D("hRec_ptFit","",nPtBins_PtFit,ptBoundaries_PtFit);

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
                hRec_ptFit->Fill(fPt,weight);
            }
        }
        TLegend lRec(0.47,0.95-nRowsLeg*0.05,0.95,0.95);
        lRec.AddEntry(hRecOld,Form("original dist.: %.0f ev.",hRecOld->Integral()),"L"); // (#it{N}^{rec}_{MC})_{old}
        lRec.AddEntry(hRecNew,Form("re-weighted dist.: %.0f ev.",hRecNew->Integral()),"L"); // (#it{N}^{rec}_{MC})_{new}
        if(tw) lRec.AddEntry((TObject*)0,"#bf{This work}","");
        lRec.SetMargin(0.175);
        lRec.SetTextSize(0.040);
        lRec.SetBorderSize(0);
        lRec.SetFillStyle(0);
        lRec.Draw();

        PlotHistos("AxE_Dissociative/",prefix+"hRec_oldNew","E0",kFALSE,0.,hRecOld,hRecNew,&lRec);
        PlotHistos("AxE_Dissociative/",prefix+"hRec_oldNew_fit","E0",kFALSE,0.,hRecOld_fit,hRecNew_fit);
        // rozprava:
        PlotHistos("_rozprava/","hGen_oldNew","E0",kFALSE,0.,hGenOld,hGenNew,&lGen);
        PlotHistos("_rozprava/","hRec_oldNew","E0",kFALSE,0.,hRecOld,hRecNew,&lRec);

        // save the new template for the pT fit
        if(process == "IncJ") 
        {
            TFile* f = new TFile("Results/" + str_subfolder + "AxE_Dissociative/" + prefix + "incTemplate.root","RECREATE");
            TList *l = new TList();
            l->Add(hRec_ptFit);
            l->Write("HistList", TObject::kSingleKey);
            l->ls();
            f->ls();
            f->Close();
        }

        // fit hNRecNew to get the fraction of events coming from H1
        RooRealVar vPt("vPt","",fPtLow,fPtUpp);
        RooDataHist dhData("dhData","",vPt,hRecNew_fit);
        RooRealVar vb_SL("vb_SL","",4.,1.,10.);
        RooRealVar vb_H1("vb_H1","",1.79,1.,10.);
        RooRealVar vn_H1("vn_H1","",3.58,1.,10.);
        vb_H1.setConstant(kTRUE);
        vn_H1.setConstant(kTRUE);
        RooGenericPdf pdfSL("pdfSL","","vPt*exp(-vb_SL*pow(vPt,2))",RooArgSet(vPt,vb_SL));
        RooGenericPdf pdfH1("pdfH1","","vPt*pow((1 + pow(vPt,2)*vb_H1/vn_H1),-vn_H1)",RooArgSet(vPt,vb_H1,vn_H1));
        RooDataHist dhSL("dhSL","",vPt,hRecOld_fit);
        RooHistPdf  hpdfSL("hpdfSL","",vPt,dhSL,0);
        RooDataHist dhH1("dhH1","",vPt,hH1_fit);
        RooHistPdf  hpdfH1("hpdfH1","",vPt,dhH1,0);
        Float_t NRec_tot = hRecNew->Integral();
        RooRealVar vN_SL("vN_SL","",NRec_tot*0.7,NRec_tot*0.1,NRec_tot*1.);
        RooRealVar vN_H1("vN_H1","",NRec_tot*0.3,NRec_tot*0.1,NRec_tot*1.);
        RooAddPdf* CombinedPDF = NULL;
        // fit with hpdfSL instead of pdfSL?
        Bool_t useHistPdfs = kTRUE;
        if(useHistPdfs) CombinedPDF = new RooAddPdf("CombinedPDF","", RooArgList(hpdfSL,hpdfH1),RooArgList(vN_SL,vN_H1));
        else            CombinedPDF = new RooAddPdf("CombinedPDF","", RooArgList(pdfSL,pdfH1),RooArgList(vN_SL,vN_H1));

        RooFitResult* fResFit = CombinedPDF->fitTo(dhData,SumW2Error(kFALSE),Extended(kTRUE),Save());
        TCanvas c("c","c",700,600);
        c.SetLeftMargin(0.12);
        c.SetRightMargin(0.03);
        c.SetTopMargin(0.085);
        c.SetBottomMargin(0.12);
        RooPlot* fPlot = vPt.frame(Title("fit of the #it{p}_{T} distribution of #it{N}_{rec}^{new}")); 
        dhData.plotOn(fPlot,Name("dhData"),MarkerStyle(kFullCircle),MarkerSize(0.8),MarkerColor(kBlack),LineColor(kBlack),LineWidth(2));
        if(useHistPdfs) {
            CombinedPDF->plotOn(fPlot,Name("hpdfSL"),Components(hpdfSL),Range(0.001,fPtUpp),LineColor(kRed),LineWidth(3),LineStyle(2));
            CombinedPDF->plotOn(fPlot,Name("hpdfH1"),Components(hpdfH1),Range(0.001,fPtUpp),LineColor(kViolet),LineWidth(3),LineStyle(2));
        } else {
            CombinedPDF->plotOn(fPlot,Name("pdfSL"),Components(pdfSL),Range(0.001,fPtUpp),LineColor(kRed),LineWidth(3),LineStyle(2));
            CombinedPDF->plotOn(fPlot,Name("pdfH1"),Components(pdfH1),Range(0.001,fPtUpp),LineColor(kViolet),LineWidth(3),LineStyle(2));
        }
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
        fPlot->Draw("][");
        // legend
        vPt.setRange("rPt_all",fPtLow,fPtUpp);
        RooAbsReal *arN_SL(NULL), *arN_H1(NULL);
        if(useHistPdfs) {
            arN_SL = hpdfSL.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
            arN_H1 = hpdfH1.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
        } else {
            arN_SL = pdfSL.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
            arN_H1 = pdfH1.createIntegral(vPt,NormSet(vPt),Range("rPt_all"));
        }        
        Float_t N_SL = arN_SL->getVal()*vN_SL.getVal();
        Float_t N_H1 = arN_H1->getVal()*vN_H1.getVal();
        Float_t perc_SL = 100. * N_SL / NRec_tot;
        Float_t perc_H1 = 100. * N_H1 / NRec_tot;
        Int_t nRows(4);
        TLegend l(0.55,0.90-nRows*0.05,0.90,0.90);
        l.AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2),"");
        l.AddEntry("CombinedPDF",Form("total #it{N}_{rec}^{new} = %.0f",NRec_tot),"L");
        if(useHistPdfs) {
            l.AddEntry("hpdfSL",Form("#it{N}_{rec}^{SL} = %.0f (%.1f%%)",N_SL,perc_SL),"L");
            l.AddEntry("hpdfH1",Form("#it{N}_{rec}^{H1} = %.0f (%.1f%%)",N_H1,perc_H1),"L");
        } else {
            l.AddEntry("pdfSL",Form("#it{N}_{rec}^{SL} = %.0f (%.1f%%)",N_SL,perc_SL),"L");
            l.AddEntry("pdfH1",Form("#it{N}_{rec}^{H1} = %.0f (%.1f%%)",N_H1,perc_H1),"L");
        }
        l.SetTextSize(0.045);
        l.SetBorderSize(0);
        l.SetFillStyle(0);
        l.Draw();
        if(useHistPdfs) c.Print("Results/" + str_subfolder + "AxE_Dissociative/" + prefix + "fit_histPdfs.pdf");
        else            c.Print("Results/" + str_subfolder + "AxE_Dissociative/" + prefix + "fit_pdfs.pdf");
        return;
    }
}

void AxE_Dissociative(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();
    SetPtBinning_PtFit();

    ReweightIncPtShape("IncJ");
    ReweightIncPtShape("IncP");
    return;
}