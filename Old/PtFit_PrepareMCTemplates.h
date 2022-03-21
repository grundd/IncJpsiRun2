// PtFit_PrepareMCTemplates.h
// David Grund, Mar 04, 2022

TString NamesPDFs[6] = {"CohJ","IncJ","CohP","IncP","Bkgr","Diss"};
Double_t fPtStopWeigh[4] = {0.17, 1.2, 0.5, 1.2}; // GeV/c

Double_t fPtGenerated_PtFit;

void PtFit_FillHistogramsMC(Int_t iMC, TH1D *hist);
void PtFit_PreparePDFs();
void PtFit_PreparePDFs_modRA_CohJ(Bool_t bStopWeigh);
void PtFit_PreparePDFs_modRA_all(Bool_t bStopWeigh);
void PtFit_SetCanvas(TCanvas *c, Bool_t isLogScale);
TTree* PtFit_GetTreeMCRecPsi2s(Int_t iMC);

void PtFit_PrepareMCTemplates_main()
{
    // prepare histograms that will be normalized to PDFs
    // variable sized bins introduced in PtFit_SubtractBackground.h used
    // official STARlight data used for hCohJ, hIncJ, hCohP, hIncP and hBkgr
    // H1 parametrization used for hDiss
    PtFit_PreparePDFs();

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_NoBkg/modRA_CohJ_ratios/");
    // prepare PDFs for CohJ from STARlight data generated with modified values of R_A
    // weigh hRec over the whole pT range:
    PtFit_PreparePDFs_modRA_CohJ(kFALSE);
    // stop weighing at fPtStopWeigh[0] = 0.2 GeV/c
    PtFit_PreparePDFs_modRA_CohJ(kTRUE);

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_NoBkg/modRA_all_ratios/");
    // prepare PDFs for CohJ, IncJ, CohP and IncP from STARlight data generated with R_A = 7.350 fm (optimal value)
    // weigh hRec over the whole pT range:
    PtFit_PreparePDFs_modRA_all(kFALSE);
    // stop weighing at fPtStopWeigh[iMC]
    PtFit_PreparePDFs_modRA_all(kTRUE);

    return;
}

// #############################################################################################

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
        HistPDFs[0] = new TH1D(("h" + NamesPDFs[0]).Data(), ("h" + NamesPDFs[0]).Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[1] = new TH1D(("h" + NamesPDFs[1]).Data(), ("h" + NamesPDFs[1]).Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[2] = new TH1D(("h" + NamesPDFs[2]).Data(), ("h" + NamesPDFs[2]).Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[3] = new TH1D(("h" + NamesPDFs[3]).Data(), ("h" + NamesPDFs[3]).Data(), nPtBins_PtFit, ptBoundaries_PtFit);
        HistPDFs[4] = new TH1D(("h" + NamesPDFs[4]).Data(), ("h" + NamesPDFs[4]).Data(), nPtBins_PtFit, ptBoundaries_PtFit);

        for(Int_t i = 0; i < 5; i++){
            PtFit_FillHistogramsMC(i, HistPDFs[i]);
            l->Add(HistPDFs[i]);
        }
        // ***************************************************************
        // Create the dissociative PDF
        HistPDFs[5] = new TH1D(("h" + NamesPDFs[5]).Data(), ("h" + NamesPDFs[5]).Data(), nPtBins_PtFit, ptBoundaries_PtFit);
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

// #############################################################################################

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

// #############################################################################################

void PtFit_PreparePDFs_modRA_CohJ(Bool_t bStopWeigh)
{
    TString name;
    if(!bStopWeigh) name = "Trees/" + str_subfolder + "PtFit/MCTemplates_modRA_CohJ.root";
    else            name = "Trees/" + str_subfolder + "PtFit/MCTemplates_modRA_CohJ_stopWeigh.root";
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
        if(!bStopWeigh) for(Int_t i = 0; i < 14; i++) str_modRA[i] = Form("hCohJ_modRA_%.2f", RA[i]);
        else            for(Int_t i = 0; i < 14; i++) str_modRA[i] = Form("hCohJ_modRA_%.2f_stopWeigh", RA[i]);
        
        TH1D *hCohJ_modRA[14] = { NULL };
        TH1D *hGenOld[14] = { NULL };
        TH1D *hGenNew[14] = { NULL };
        TH1D *hRatios[14] = { NULL };

        // Go over MC data
        // add the input file for 0 !!!!, then start from 0
        for(Int_t i = 1; i < 14; i++){
            Printf("Now calculating hRec for R_A = %.2f fm.", RA[i]);

            hCohJ_modRA[i] = (TH1D*)hRec->Clone(str_modRA[i].Data());
            hCohJ_modRA[i]->SetTitle(str_modRA[i].Data());

            // Correct the shape => calculate the ratios
            hGenOld[i] = new TH1D(Form("hGenOld%i",i), Form("hGenOld%i",i), nPtBins_PtFit, ptBoundaries_PtFit);
            hGenNew[i] = new TH1D(Form("hGenNew%i",i), Form("hGenNew%i",i), nPtBins_PtFit, ptBoundaries_PtFit); 

            // Open fGenOld file and get the tree
            TFile *fGenOld = TFile::Open("Trees/STARlight/tGen_CohJ_RA_6.624.root","read");
            if(!fGenOld){
                Printf("File fGenOld not found! Terminating...");
                return;
            }

            TTree *tGenOld = (TTree*)fGenOld->Get("tGen");
            if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName());
            tGenOld->SetBranchAddress("fPtGen", &fPtGenerated_PtFit);

            // Open fGenNew file and get the tree
            TFile *fGenNew = TFile::Open(Form("Trees/STARlight/tGen_CohJ_RA_%.3f.root", RA[i]),"read");
            if(!fGenNew){
                Printf("File fGenNew not found! Terminating...");
                return;
            }

            TTree *tGenNew = (TTree*)fGenNew->Get("tGen");
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

            // If we want to stop weighing at fPtStopWeigh[0] (for CohJ), set all ratios above this value to 1.0
            if(bStopWeigh){
                for(Int_t iBin = 1; iBin <= nPtBins_PtFit; iBin++){
                    if(hRatios[i]->GetBinCenter(iBin) > fPtStopWeigh[0]) hRatios[i]->SetBinContent(iBin, 1.0);
                }  
            }

            // Print the results to the text file
            TString name_out = "";
            if(!bStopWeigh) name_out = "Results/" + str_subfolder + Form("PtFit_NoBkg/modRA_CohJ_ratios/RA_%.3f.txt", RA[i]);
            else            name_out = "Results/" + str_subfolder + Form("PtFit_NoBkg/modRA_CohJ_ratios/RA_%.3f_stopWeigh.txt", RA[i]);
            ofstream outfile(name_out.Data());
            outfile << "pT_low\tpT_upp\tnEvOld\tnEvNew\tratio\n";
            for(Int_t iBin = 1; iBin <= nPtBins_PtFit; iBin++){
                outfile << Form("%.3f\t%.3f\t%.0f\t%.0f\t%.3f\n",
                                hRatios[i]->GetBinLowEdge(iBin), hRatios[i]->GetBinLowEdge(iBin+1), 
                                hGenOld[i]->GetBinContent(iBin), hGenNew[i]->GetBinContent(iBin), hRatios[i]->GetBinContent(iBin));
            }
            outfile.close();

            // Correct the shape of reconstructed events by the ratios
            hCohJ_modRA[i]->Multiply(hRatios[i]);
            // Add the histogram to the list
            Printf("Adding %s to the list.", hCohJ_modRA[i]->GetName());
            l->Add(hCohJ_modRA[i]);
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

// #############################################################################################

void PtFit_PreparePDFs_modRA_all(Bool_t bStopWeigh)
{
    // Here the PDFs with R_A = 7.350 fm (optimal value) are created
    TString name;
    if(!bStopWeigh) name = "Trees/" + str_subfolder + "PtFit/MCTemplates_modRA_all.root";
    else            name = "Trees/" + str_subfolder + "PtFit/MCTemplates_modRA_all_stopWeigh.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("PDFs with R_A = 7.350 already created.");
        return;

    } else { 
    
        TString str_modRA[4] = {"hCohJ_modRA_7.350",
                                "hIncJ_modRA_7.350",
                                "hCohP_modRA_7.350",
                                "hIncP_modRA_7.350"};
        if(bStopWeigh) for(Int_t iMC = 0; iMC < 4; iMC++) str_modRA[iMC] += "_stopWeigh";

        TList *l = new TList();
        TH1D *hRec[2] = { NULL };
        TH1D *h_modRA[4] = { NULL };
        TH1D *hGenOld[4] = { NULL };
        TH1D *hGenNew[4] = { NULL };
        TH1D *hRatios[4] = { NULL };
     
        for(Int_t iMC = 0; iMC < 4; iMC++)
        {
            Printf("Now calculating h_modRA for %s.", str_modRA[iMC].Data());

            // Correct the shape => calculate the ratios
            // Load hGenOld, hGenNew
            hGenOld[iMC] = new TH1D(Form("hGenOld%i",iMC), Form("hGenOld%i",iMC), nPtBins_PtFit, ptBoundaries_PtFit);
            hGenNew[iMC] = new TH1D(Form("hGenNew%i",iMC), Form("hGenNew%i",iMC), nPtBins_PtFit, ptBoundaries_PtFit); 

            // Open fGenOld file and get the tree
            TString str_in_old = Form("Trees/STARlight/tGen_%s_RA_6.624.root", NamesPDFs[iMC].Data());
            TFile *fGenOld = TFile::Open(str_in_old.Data(),"read");
            if(!fGenOld){
                Printf("File fGenOld not found! Terminating...");
                return;
            }

            TTree *tGenOld = (TTree*)fGenOld->Get("tGen");
            if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName());
            tGenOld->SetBranchAddress("fPtGen", &fPtGenerated_PtFit);

            // Open fGenNew file and get the tree
            TString str_in_new = Form("Trees/STARlight/tGen_%s_RA_7.350.root", NamesPDFs[iMC].Data());
            TFile *fGenNew = TFile::Open(str_in_new.Data(),"read");
            if(!fGenNew){
                Printf("File fGenNew not found! Terminating...");
                return;
            }

            TTree *tGenNew = (TTree*)fGenNew->Get("tGen");
            if(tGenNew) Printf("Tree %s loaded.", tGenNew->GetName());
            tGenNew->SetBranchAddress("fPtGen", &fPtGenerated_PtFit);

            // Loop over tGenOld entries
            for(Int_t iEntry = 0; iEntry < tGenOld->GetEntries(); iEntry++){
                tGenOld->GetEntry(iEntry);
                hGenOld[iMC]->Fill(fPtGenerated_PtFit);
            }
            // Loop over tGenNew entries
            for(Int_t iEntry = 0; iEntry < tGenNew->GetEntries(); iEntry++){
                tGenNew->GetEntry(iEntry);
                hGenNew[iMC]->Fill(fPtGenerated_PtFit);
            }
            // Calculate the ratios
            hRatios[iMC] = (TH1D*)hGenNew[iMC]->Clone(Form("hRatios%i",iMC));
            hRatios[iMC]->SetTitle(Form("hRatios%i",iMC));
            hRatios[iMC]->Sumw2();
            hRatios[iMC]->Divide(hGenOld[iMC]);

            // CohJ and IncJ (iMC == 0,1)
            if(iMC == 0 || iMC == 1)
            {
                // Define the histogram with reconstructed events and fill it
                hRec[iMC] = new TH1D("hRec","hRec",nPtBins_PtFit,ptBoundaries_PtFit);
                PtFit_FillHistogramsMC(iMC, hRec[iMC]);

                h_modRA[iMC] = (TH1D*)hRec[iMC]->Clone(str_modRA[iMC].Data());
                h_modRA[iMC]->SetTitle(str_modRA[iMC].Data());   

                // If we want to stop weighing at fPtStopWeigh[iMC], set all ratios above this value to 1.0
                if(bStopWeigh){
                    for(Int_t iBin = 1; iBin <= nPtBins_PtFit; iBin++){
                        if(hRatios[iMC]->GetBinCenter(iBin) > fPtStopWeigh[iMC]) hRatios[iMC]->SetBinContent(iBin, 1.0);
                    }  
                }

                // Correct the shape of reconstructed events by the ratios
                h_modRA[iMC]->Multiply(hRatios[iMC]);
            }
            // CohP and IncP (iMC == 2,3)
            if(iMC == 2 || iMC == 3)
            {
                if(!isPass3){
                    Printf("This option is not supported. Skipping..."); 
                    continue;
                }
                // detine the output histogram with rec. events
                h_modRA[iMC] = new TH1D(str_modRA[iMC].Data(),str_modRA[iMC].Data(),nPtBins_PtFit,ptBoundaries_PtFit);
                // define the paths to the file and tRec
                TString str_f_in = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_";
                // choose MC dataset
                if(iMC == 2) str_f_in += "kCohPsi2sToMuPi_2.root";
                if(iMC == 3) str_f_in += "kIncohPsi2sToMuPi_2.root";
                // open the input file
                TFile *fRec = TFile::Open(str_f_in.Data(), "read");
                if(fRec) Printf("File %s loaded.", fRec->GetName());
                // get fOutputList
                TList *l = (TList*) fRec->Get("AnalysisOutput/fOutputListcharged");
                if(l) Printf("List %s loaded.", l->GetName()); 
                // get the MCRec tree
                TTree *tRec = (TTree*)l->FindObject("fTreeJpsi");
                if(tRec) Printf("Tree %s loaded.", tRec->GetName());
                // connect tree varibles
                ConnectTreeVariablesMCRec(tRec, kTRUE);

                TAxis *xAxis = hRatios[iMC]->GetXaxis();
                // run over reconstructed feed-down events
                Int_t iBinP = 0;
                Int_t iBinJ = 0;
                Double_t fJpsi = 0;
                // counters
                Int_t nEntriesAnalysed = 0;
                for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++)
                {
                    tRec->GetEntry(iEntry);
                    // m between 3.0 and 3.2 GeV/c^2, pT cut: all
                    if(EventPassedMCRec(1, 2))
                    {
                        // find index of the bin to which the current fPtGen_Psi2s corresponds
                        iBinP = xAxis->FindBin(fPtGen_Psi2s);
                        // find index of the bin to which the current fPt corresponds
                        iBinJ = xAxis->FindBin(fPt);
                        // scale the J/psi entry by the ratio with the index iBinPsi2s
                        fJpsi = hRatios[iMC]->GetBinContent(iBinP);
                        // add the entry to h_modRA[iMC]
                        h_modRA[iMC]->SetBinContent(iBinJ,h_modRA[iMC]->GetBinContent(iBinJ)+fJpsi);
                        // fill the histogram hRecOld
                        h_modRA[iMC]->Fill(fPt);
                    }
                    if((iEntry+1) % 100000 == 0)
                    {
                        nEntriesAnalysed += 100000;
                        Printf("%i entries analysed.", nEntriesAnalysed);
                    }
                } 
            }
            // Add the histogram to the list
            Printf("Adding %s to the list.", h_modRA[iMC]->GetName());
            l->Add(h_modRA[iMC]);

            // Print the results to the text file
            TString name_out = "Results/" + str_subfolder + Form("PtFit_NoBkg/modRA_all_ratios/%s_RA_7.350", NamesPDFs[iMC].Data());
            if(!bStopWeigh) name_out += ".txt";
            else            name_out += "_stopWeigh.txt"; 
            ofstream outfile(name_out.Data());
            outfile << "pT_low\tpT_upp\tnEvOld\tnEvNew\tratio\n";
            for(Int_t iBin = 1; iBin <= nPtBins_PtFit; iBin++){
                outfile << Form("%.3f\t%.3f\t%.0f\t%.0f\t%.3f\n",
                                hRatios[iMC]->GetBinLowEdge(iBin), hRatios[iMC]->GetBinLowEdge(iBin+1), 
                                hGenOld[iMC]->GetBinContent(iBin), hGenNew[iMC]->GetBinContent(iBin), hRatios[iMC]->GetBinContent(iBin));
            }
            outfile.close();
        }
        // Save results to the output file
        // Create the output file
        TFile *f = new TFile(name.Data(),"RECREATE");
        l->Write("HistList", TObject::kSingleKey);
        l->ls();
        f->ls();
        f->Close();

        return;
    }
}

// #############################################################################################