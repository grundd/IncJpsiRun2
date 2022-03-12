// PtFit_PrepareMCTemplates.h
// David Grund, Mar 04, 2022

TString NamesPDFs[6] = {"hCohJ","hIncJ","hCohP","hIncP","hBkgr","hDiss"};
Double_t fPtStopWeigh = 0.2; // GeV/c

Double_t fPtGenerated_PtFit;

void PtFit_FillHistogramsMC(Int_t iMC, TH1D *hist);
void PtFit_PreparePDFs();
void PtFit_PreparePDFs_CohJmodRA(Bool_t bStopWeigh);
void PtFit_SetCanvas(TCanvas *c, Bool_t isLogScale);

void PtFit_PrepareMCTemplates_main()
{
    // prepare histograms that will be normalized to PDFs
    // variable sized bins introduced in PtFit_SubtractBackground.h used
    // official STARlight data used for hCohJ, hIncJ, hCohP, hIncP and hBkgr
    // H1 parametrization used for hDiss
    PtFit_PreparePDFs();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_NoBkg/CohJmodRA_ratios/");
    // prepare PDFs for CohJ from STARlight data generated with modified values of R_A
    // weigh hRec over the whole pT range:
    PtFit_PreparePDFs_CohJmodRA(kFALSE);
    // stop weighing at fPtStopWeigh = 0.2 GeV/c
    PtFit_PreparePDFs_CohJmodRA(kTRUE);

    return;
}

void PtFit_PreparePDFs()
{
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
        f->Close();

        return;
    }
}

void PtFit_FillHistogramsMC(Int_t iMC, TH1D *hist)
{
    // Load the data
    TFile *file = NULL;
    switch(iMC){
        case 0:
            file = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kCohJpsiToMu.root").Data(), "read");
            break;
        case 1:
            file = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
            break;
        case 2:
            file = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kCohPsi2sToMuPi.root").Data(), "read");
            break;
        case 3:
            file = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohPsi2sToMuPi.root").Data(), "read");
            break;
        case 4:
            file = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kTwoGammaToMuMedium.root").Data(), "read");
            break;
    }
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tRec = dynamic_cast<TTree*> (file->Get(str_in_MC_tree_rec.Data()));
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

void PtFit_PreparePDFs_CohJmodRA(Bool_t bStopWeigh)
{
    TString name;
    if(!bStopWeigh) name = "Trees/" + str_subfolder + "PtFit/MCTemplates_CohJmodRA.root";
    else            name = "Trees/" + str_subfolder + "PtFit/MCTemplates_CohJmodRA_StopWeigh.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("PDFs for CohJ with modified R_A already created.");
        return;

    } else { 

        // Fill the histogram with reconstructed events of CohJ
        TH1D *hRec = new TH1D("hRec","hRec",nPtBins_PtFit,ptBoundaries_PtFit);
        PtFit_FillHistogramsMC(0, hRec);

        // Define output histograms with predefined binning to create PDFs
        TList *l = new TList();
        Double_t RA[14] = { 0 };
        TString str_modRA[14] = { "" };
        RA[0] = 7.53;
        for(Int_t i = 0; i < 13; i++) RA[i+1] = 6.60 + i*0.10;
        if(!bStopWeigh) for(Int_t i = 0; i < 14; i++) str_modRA[i] = Form("hCohJmodRA_%.2f", RA[i]);
        else            for(Int_t i = 0; i < 14; i++) str_modRA[i] = Form("hCohJmodRA_StopWeigh_%.2f", RA[i]);
        
        TH1D *hCohJmodRA[14] = { NULL };
        TH1D *hGenOld[14] = { NULL };
        TH1D *hGenNew[14] = { NULL };
        TH1D *hRatios[14] = { NULL };

        // Go over MC data
        // add the input file for 0 !!!!, then start from 0
        for(Int_t i = 1; i < 14; i++){
            Printf("Now calculating hRec for R_A = %.2f fm.", RA[i]);

            hCohJmodRA[i] = (TH1D*)hRec->Clone(str_modRA[i].Data());
            hCohJmodRA[i]->SetTitle(str_modRA[i].Data());

            // Correct the shape => calculate the ratios
            hGenOld[i] = new TH1D(Form("hGenOld%i",i), Form("hGenOld%i",i), nPtBins_PtFit, ptBoundaries_PtFit);
            hGenNew[i] = new TH1D(Form("hGenNew%i",i), Form("hGenNew%i",i), nPtBins_PtFit, ptBoundaries_PtFit); 

            // Open fGenOld file and get the tree
            TFile *fGenOld = TFile::Open("Trees/STARlight/tGen_CohJ_RA_6.624.root","read");
            if(!fGenOld){
                Printf("File fGenOld not found! Terminating...");
                return;
            }

            TTree *tGenOld = (TTree*)fGenOld->Get("tGenOld");
            if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName());
            tGenOld->SetBranchAddress("fPtGen", &fPtGenerated_PtFit);

            // Open fGenNew file and get the tree
            TFile *fGenNew = TFile::Open(Form("Trees/STARlight/tGen_CohJ_RA_%.3f.root", RA[i]),"read");
            if(!fGenNew){
                Printf("File fGenNew not found! Terminating...");
                return;
            }

            TTree *tGenNew = (TTree*)fGenNew->Get("tGenNew");
            if(tGenNew) Printf("Tree %s loaded.", tGenNew->GetName());
            tGenNew->SetBranchAddress("fPtGen", &fPtGenerated_PtFit);

            // Loop over tGenOld entries
            for(Int_t iEntry = 0; iEntry < tGenOld->GetEntries(); iEntry++){
                tGenOld->GetEntry(iEntry);
                hGenOld[i]->Fill(fPtGenerated_PtFit);
            }
            // Loop over tGenNew entries
            for(Int_t iEntry = 0; iEntry < tGenNew->GetEntries(); iEntry++){
                tGenNew->GetEntry(iEntry);
                hGenNew[i]->Fill(fPtGenerated_PtFit);
            }
            // Calculate the ratios
            hRatios[i] = (TH1D*)hGenNew[i]->Clone(Form("hRatios%i",i));
            hRatios[i]->SetTitle(Form("hRatios%i",i));
            hRatios[i]->Sumw2();
            hRatios[i]->Divide(hGenOld[i]);

            // If we want to stop weighing at fPtStopWeigh, set all ratios above this value to 1.0
            if(bStopWeigh){
                for(Int_t iBin = 1; iBin <= nPtBins_PtFit; iBin++){
                    if(hRatios[i]->GetBinCenter(iBin) > fPtStopWeigh) hRatios[i]->SetBinContent(iBin, 1.0);
                }  

            }

            // Print the results to the text file
            TString name_out = "";
            if(!bStopWeigh) name_out = "Results/" + str_subfolder + Form("PtFit_NoBkg/CohJmodRA_ratios/RA_%.3f.txt", RA[i]);
            else            name_out = "Results/" + str_subfolder + Form("PtFit_NoBkg/CohJmodRA_ratios/RA_StopWeigh_%.3f.txt", RA[i]);
            ofstream outfile(name_out.Data());
            outfile << "pT_low\tpT_upp\tnEvOld\tnEvNew\tratio\n";
            for(Int_t iBin = 1; iBin <= nPtBins_PtFit; iBin++){
                outfile << Form("%.3f\t%.3f\t%.0f\t%.0f\t%.3f\n",
                                hRatios[i]->GetBinLowEdge(iBin), hRatios[i]->GetBinLowEdge(iBin+1), 
                                hGenOld[i]->GetBinContent(iBin), hGenNew[i]->GetBinContent(iBin), hRatios[i]->GetBinContent(iBin));
            }
            outfile.close();

            // Correct the shape of reconstructed events by the ratios
            hCohJmodRA[i]->Multiply(hRatios[i]);
            // Add the histogram to the list
            l->Add(hCohJmodRA[i]);
        }
        // Save results to the output file
        // Create the output file
        TFile *f = new TFile(name.Data(),"RECREATE");
        l->Write("HistList", TObject::kSingleKey);
        f->ls();
        f->Close();

        return;
    }   
}

void PtFit_PreparePDFs_AllmodRA(Bool_t bStopWeigh)
{
    TString name;
    if(!bStopWeigh) name = "Trees/" + str_subfolder + "PtFit/MCTemplates_AllmodRA.root";
    else            name = "Trees/" + str_subfolder + "PtFit/MCTemplates_AllmodRA_StopWeigh.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("PDFs with R_A = 7.350 already created.");
        return;

    } else { 
    
        

    }

    return;
}

void PtFit_SetCanvas(TCanvas *c, Bool_t isLogScale)
{
    if(isLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.02);

    return;
}

/*
void PreparePDFs_modRA_all(){

    return;
}

void ConnectTreeVariablesPt(TTree *t){

    t->SetBranchAddress("fPt", &fPt);

    Printf("Variables from %s connected.", t->GetName());

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