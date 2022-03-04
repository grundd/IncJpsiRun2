// PtFit_PrepareMCTemplates.h
// David Grund, Mar 04, 2022

TString NamesPDFs[6] = {"hCohJ","hIncJ","hCohP","hIncP","hBkgr","hDiss"};

//Double_t fPtGenerated;
//Bool_t bStopWeight = kFALSE;
//Double_t StopWeight = 0.2;

void PtFit_PreparePDFs();
void PtFit_FillHistogramsMC(Int_t iMC, TH1D *hist);

void PtFit_PrepareMCTemplates_main(){

    // prepare histograms that will be normalized to PDFs
    // variable sized bins introduced in PtFit_SubtractBackground.h used
    // official STARlight data used for hCohJ, hIncJ, hCohP, hIncP and hBkgr
    // H1 parametrization used for hDiss
    PtFit_PreparePDFs();

    return;
}

void PtFit_PreparePDFs(){

    TString name = "Trees/" + str_subfolder + "PtFit/MCTemplates.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("PDFs from MC data already created.");
        return;

    } else {   

        // ***************************************************************
        // Go over MC data
        // Define output histograms with predefined binning to create PDFs
        TList *l = new TList();
        TH1D *HistPDFs[6] = { NULL };
        HistPDFs[0] = new TH1D(NamesPDFs[0].Data(), NamesPDFs[0].Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[1] = new TH1D(NamesPDFs[1].Data(), NamesPDFs[1].Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[2] = new TH1D(NamesPDFs[2].Data(), NamesPDFs[2].Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[3] = new TH1D(NamesPDFs[3].Data(), NamesPDFs[3].Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[4] = new TH1D(NamesPDFs[4].Data(), NamesPDFs[4].Data(), nPtBins_PtFit, ptBoundaries_PtFit);

        for(Int_t i = 0; i < 5; i++){
            PtFit_FillHistogramsMC(i, HistPDFs[i]);
            l->Add(HistPDFs[i]);
        }
        // ***************************************************************
        // Create the dissociative PDF
        HistPDFs[5] = new TH1D(NamesPDFs[5].Data(), NamesPDFs[5].Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        TF1 *fDissH1 = new TF1("fDissH1","x*pow((1 + x*x*[0]/[1]),-[1])",fPtCutLow_PtFit,fPtCutUpp_PtFit);
        Double_t b_pd = 1.79;
        Double_t n_pd = 3.58;
        fDissH1->SetParameter(0,b_pd);
        fDissH1->SetParameter(1,n_pd);
        Int_t i = 0;
        while(i < 1e6){
            HistPDFs[5]->Fill(fDissH1->GetRandom());
            i++;
        }
        //HistPDFs[5]->Draw();
        l->Add(HistPDFs[5]);

        // ***************************************************************
        // Save results to the output file
        // Create the output file
        TFile *f = new TFile(name.Data(),"RECREATE");
        l->Write("HistList", TObject::kSingleKey);
        f->ls();

        return;
    }
}

void PtFit_FillHistogramsMC(Int_t iMC, TH1D *hist){

    // Load the data
    TFile *file = NULL;
    switch(iMC){
        case 0:
            file = TFile::Open((str_in_MC_fldr + "AnalysisResults_MC_kCohJpsiToMu.root").Data(), "read");
            break;
        case 1:
            file = TFile::Open((str_in_MC_fldr + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
            break;
        case 2:
            file = TFile::Open((str_in_MC_fldr + "AnalysisResults_MC_kCohPsi2sToMuPi.root").Data(), "read");
            break;
        case 3:
            file = TFile::Open((str_in_MC_fldr + "AnalysisResults_MC_kIncohPsi2sToMuPi.root").Data(), "read");
            break;
        case 4:
            file = TFile::Open((str_in_MC_fldr + "AnalysisResults_MC_kTwoGammaToMuMedium.root").Data(), "read");
            break;
    }
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tRec = dynamic_cast<TTree*> (file->Get(str_in_MC_tree.Data()));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    Printf("Tree %s has %lli entries.", tRec->GetName(), tRec->GetEntries());

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);

        // m between 3.0 and 3.2 GeV/c^2, pT cut: all
        if(EventPassedMCRec(1, 2)){
            nEvPassed++;
            hist->Fill(fPt);
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("%i events passed the selections.", nEvPassed);

    return;
}

/*
void PreparePDFs_modRA_CohJ(){

    TString str_name = "";
    if(bStopWeight) str_name = Form("%sPDFs_MC_modRA_Binning%i_StopWeight.root", OutputPDFs.Data(), BinningOpt);
    else            str_name = Form("%sPDFs_MC_modRA_Binning%i.root", OutputPDFs.Data(), BinningOpt);

    TFile *file = TFile::Open(str_name.Data(),"read");
    if(file){
        Printf("MC PDFs for modRA with this binning already created.");
        return;
    }

    TH1D *hRec = new TH1D("hRec","hRec",nBins,ptEdges);
    // Fill the histogram with reconstructed events
    PtFit_FillHistogramsMC(0, hRec);

    // ***************************************************************
    // Go over MC data
    // Define output histograms with predefined binning to create PDFs
    TList *l = new TList();
    TString str_modRA[14] = {""};
    str_modRA[0] = "hCohJ_modRA_7.53";
    Double_t R_A = 6.60;
    for(Int_t i = 0; i < 13; i++){
        str_modRA[i+1] = Form("hCohJ_modRA_%.2f", R_A);
        R_A += 0.10;
    }
    TH1D *hCohJ_modRA[14] = { NULL };
    TH1D *hGenOld[14] = { NULL };
    TH1D *hGenNew[14] = { NULL };
    TH1D *hRatios[14] = { NULL };

    for(Int_t i = 0; i < 14; i++){

        if(i == 0) Printf("Now calculating the histogram with R_A = 7.53 fm.");
        else Printf("Now calculating the histogram with R_A = %.2f fm.", 6.6+(i-1)*0.1);

        hCohJ_modRA[i] = (TH1D*)hRec->Clone(str_modRA[i].Data());
        hCohJ_modRA[i]->SetTitle(str_modRA[i].Data());
        
        // Correct the shape => calculate the ratios
        hGenOld[i] = new TH1D(Form("hGenOld%i",i), Form("hGenOld%i",i), nBins, ptEdges);
        hGenNew[i] = new TH1D(Form("hGenNew%i",i), Form("hGenNew%i",i), nBins, ptEdges); 

        TString str_fGenOld = "";
        TString str_fGenNew = "";
        if(i == 0){
            str_fGenOld = "Trees/STARlight/CoherentShape/tGenOld_6000k.root";
            str_fGenNew = "Trees/STARlight/CoherentShape/tGenNew_6000k_7.530.root";
        } else {
            str_fGenOld = Form("Trees/STARlight/OptimalRA/tGenOld_6000k.root");
            str_fGenNew = Form("Trees/STARlight/OptimalRA/tGenNew_6000k_%.3f.root", 6.6+(i-1)*0.1);
        }

        // Open fGenOld file and get the tree
        TFile *fGenOld = TFile::Open(str_fGenOld.Data(),"read");
        if(!fGenOld){
            Printf("File fGenOld not found! Terminating...");
            return;
        }

        TTree *tGenOld = (TTree*)fGenOld->Get("tGenOld");
        if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName());
        tGenOld->SetBranchAddress("fPtGen", &fPtGenerated);

        // Open fGenNew file and get the tree
        TFile *fGenNew = TFile::Open(str_fGenNew.Data(),"read");
        if(!fGenNew){
            Printf("File fGenNew not found! Terminating...");
            return;
        }

        TTree *tGenNew = (TTree*)fGenNew->Get("tGenNew");
        if(tGenNew) Printf("Tree %s loaded.", tGenNew->GetName());
        tGenNew->SetBranchAddress("fPtGen", &fPtGenerated);

        // Loop over tGenOld entries
        for(Int_t iEntry = 0; iEntry < tGenOld->GetEntries(); iEntry++){
            tGenOld->GetEntry(iEntry);
            hGenOld[i]->Fill(fPtGenerated);
        }
        // Loop over tGenNew entries
        for(Int_t iEntry = 0; iEntry < tGenNew->GetEntries(); iEntry++){
            tGenNew->GetEntry(iEntry);
            hGenNew[i]->Fill(fPtGenerated);
        }
        // Calculate the ratios
        hRatios[i] = (TH1D*)hGenNew[i]->Clone(Form("hRatios%i",i));
        hRatios[i]->SetTitle(Form("hRatios%i",i));
        hRatios[i]->Sumw2();
        hRatios[i]->Divide(hGenOld[i]);

        // Stop weighting at: Double_t StopWeight    
        for(Int_t iBin = 1; iBin <= nBins; iBin++){
            if(hRatios[i]->GetBinCenter(iBin) > StopWeight) hRatios[i]->SetBinContent(iBin, 1.0);
        }    

        // Print the results to console
        Printf("\n");
        Printf("*****");
        Printf("Output:");
        Printf("pT_low\tpT_upp\tnEvOld\tnEvNew\tRatio");
        for(Int_t iBin = 1; iBin <= nBins; iBin++){
            Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
                hRatios[i]->GetBinLowEdge(iBin), hRatios[i]->GetBinLowEdge(iBin+1), 
                hGenOld[i]->GetBinContent(iBin), hGenNew[i]->GetBinContent(iBin), hRatios[i]->GetBinContent(iBin));
        }
        Printf("*****");
        Printf("\n");  

        // Correct the shape of reconstructed events by the ratios
        hCohJ_modRA[i]->Multiply(hRatios[i]);

        l->Add(hCohJ_modRA[i]);

    }
    // ***************************************************************
    // Save results to the output file
    // Create the output file
    TFile *f = new TFile(str_name.Data(),"RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f->ls();

    return;
}

void PreparePDFs_modRA_all(){

    return;
}

void ConnectTreeVariablesPt(TTree *t){

    t->SetBranchAddress("fPt", &fPt);

    Printf("Variables from %s connected.", t->GetName());

    return;
}

void SetCanvas(TCanvas *c, Bool_t isLogScale){

    if(isLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.02);

    return;
}

void DrawCorrelationMatrix(TCanvas *cCM, RooFitResult* ResFit){

    // Set margins
    cCM->SetTopMargin(0.03);
    cCM->SetBottomMargin(0.11);
    cCM->SetRightMargin(0.17);
    cCM->SetLeftMargin(0.15);
    // Get 2D corr hist
    TH2* hCorr = ResFit->correlationHist();
    // Set X axis
    hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
    hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
    hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{inc}");
    // Set Y axis
    hCorr->GetYaxis()->SetBinLabel(1,"#it{N}_{inc}");
    hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{diss}");
    hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{coh}");
    // Set corr hist and draw it
    hCorr->SetMarkerSize(3.6);
    hCorr->GetXaxis()->SetLabelSize(0.13);
    hCorr->GetYaxis()->SetLabelSize(0.13);
    hCorr->GetZaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}
*/