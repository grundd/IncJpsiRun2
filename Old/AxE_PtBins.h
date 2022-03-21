// AxE_PtBins.h
// David Grund, Feb 27, 2022
// To calculate the acceptance x efficiency from MC data in defined pt bins

TH1D* AxE_PtBins_hNRec = NULL; 
TH1D* AxE_PtBins_hNGen = NULL; 
TH1D* AxE_PtBins_hAxE = NULL;

void AxE_PtBins_FillHistNRec();
void AxE_PtBins_FillHistNGen();
void AxE_PtBins_SaveToFile(TH1D* hist, TString name);
Double_t CalculateErrorBayes(Double_t k, Double_t n);

void AxE_PtBins_main(){

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "AxE_PtBins/");

    AxE_PtBins_hNRec = new TH1D("AxE_PtBins_hNRec","N_rec per bin",nPtBins,ptBoundaries);
    AxE_PtBins_hNGen = new TH1D("AxE_PtBins_hNGen","N_gen per bin",nPtBins,ptBoundaries);

    AxE_PtBins_FillHistNRec();
    AxE_PtBins_FillHistNGen();

    AxE_PtBins_hAxE = (TH1D*)AxE_PtBins_hNRec->Clone("AxE_PtBins_hAxE");
    AxE_PtBins_hAxE->SetTitle("AxE per bin");
    AxE_PtBins_hAxE->Sumw2();
    AxE_PtBins_hAxE->Divide(AxE_PtBins_hNGen);

    // Draw the histogram:
    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.145);
    // gStyle
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");
    // Marker and line
    //AxE_PtBins_hAxE->SetMarkerStyle(21);
    //AxE_PtBins_hAxE->SetMarkerColor(kBlue);
    //AxE_PtBins_hAxE->SetMarkerSize(1.0);
    AxE_PtBins_hAxE->SetLineColor(kBlue);
    AxE_PtBins_hAxE->SetLineWidth(1.0);
    // Vertical axis
    AxE_PtBins_hAxE->GetYaxis()->SetTitle("#it{N}_{rec}/#it{N}_{gen}");
    AxE_PtBins_hAxE->GetYaxis()->SetTitleSize(0.056);
    AxE_PtBins_hAxE->GetYaxis()->SetTitleOffset(1.3);
    AxE_PtBins_hAxE->GetYaxis()->SetLabelSize(0.056);
    AxE_PtBins_hAxE->GetYaxis()->SetDecimals(3);
    AxE_PtBins_hAxE->GetYaxis()->SetRangeUser(0.0,AxE_PtBins_hAxE->GetBinContent(1)*1.1);
    // Horizontal axis
    AxE_PtBins_hAxE->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    AxE_PtBins_hAxE->GetXaxis()->SetTitleSize(0.056);
    AxE_PtBins_hAxE->GetXaxis()->SetTitleOffset(1.2);
    AxE_PtBins_hAxE->GetXaxis()->SetLabelSize(0.056);
    AxE_PtBins_hAxE->GetXaxis()->SetLabelOffset(0.015);
    AxE_PtBins_hAxE->GetXaxis()->SetDecimals(1);
    // Eventually draw it
    AxE_PtBins_hAxE->Draw("P E1");
    // Legend
    TLegend *l = new TLegend(0.52,0.77,0.85,0.97);
    l->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    l->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l->SetTextSize(0.056);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();
    // Legend 2
    TLegend *l2 = new TLegend(0.15,0.17,0.35,0.32);
    l2->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),""); 
    l2->AddEntry((TObject*)0,Form("2.2 < #it{m} < 4.5 GeV/#it{c}^{2}"),"");
    l2->SetTextSize(0.056);
    l2->SetBorderSize(0); // no border
    l2->SetFillStyle(0);  // legend is transparent
    l2->Draw();

    // Save the figures and print the results to txt file
    TString str = "Results/" + str_subfolder + Form("AxE_PtBins/AxE_%ibins", nPtBins);
    c->Print((str + ".pdf").Data());
    c->Print((str + ".png").Data());
    ofstream outfile((str + ".txt").Data());
    outfile << std::fixed << std::setprecision(5);
    //outfile << "Bin \tAxE [%%] \tAxE_err [%%] \n";
    for(Int_t i = 1; i <= nPtBins; i++){
        outfile << i << "\t" << AxE_PtBins_hAxE->GetBinContent(i)*100 << "\t\t" << AxE_PtBins_hAxE->GetBinError(i)*100 << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", (str + ".txt").Data());

    // Compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kFALSE;
    if(DebugErrors){
        Double_t ErrRoot = 0;
        Double_t ErrBayes = 0;    
        for(Int_t i = 1; i <= nPtBins; i++){
            ErrRoot = AxE_PtBins_hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(AxE_PtBins_hNRec->GetBinContent(i),AxE_PtBins_hNGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // Cross-check: calculate the total value of AxE
    Double_t NRecTot = 0;
    Double_t NGenTot = 0;
    for(Int_t i = 1; i <= nPtBins; i++){
        NRecTot += AxE_PtBins_hNRec->GetBinContent(i);
        NGenTot += AxE_PtBins_hNGen->GetBinContent(i);
    }
    Double_t AxETot = NRecTot / NGenTot;
    Double_t AxETot_err = CalculateErrorBayes(NRecTot, NGenTot);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxETot*100, AxETot_err*100);

    return;
}

void AxE_PtBins_FillHistNRec(){

    // Check if the corresponding text file already exists
    TString file("Results/" + str_subfolder + "AxE_PtBins/");
    file.Append(Form("NRec_%ibins.txt", nPtBins));

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill AxE_PtBins_hNRec with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            AxE_PtBins_hNRec->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N rec per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open((str_in_MC_fldr_rec + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fRec) Printf("MC rec file loaded.");

        TTree *tRec = dynamic_cast<TTree*> (fRec->Get(str_in_MC_tree_rec.Data()));
        if(tRec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCRec(tRec);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nPtBins; iPtBin++){
            Int_t NRec = 0;
            for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
                tRec->GetEntry(iEntry);
                // m between 2.2 and 4.5 GeV/c^2, pT in a given bin 
                if(EventPassedMCRec(0, 4, iPtBin)) NRec++;
            }
            AxE_PtBins_hNRec->SetBinContent(iPtBin, NRec);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        AxE_PtBins_SaveToFile(AxE_PtBins_hNRec, file);

        return;
    }
}

void AxE_PtBins_FillHistNGen(){
    // Check if the corresponding text file already exists
    TString file("Results/" + str_subfolder + "AxE_PtBins/");
    file.Append(Form("NGen_%ibins.txt", nPtBins));

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill AxE_PtBins_hNGen with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            AxE_PtBins_hNGen->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N gen per bin for %s... ***", file.Data());

        TFile *fGen = TFile::Open((str_in_MC_fldr_gen + "AnalysisResults_MC_kIncohJpsiToMu.root").Data(), "read");
        if(fGen) Printf("MC gen file loaded.");

        TTree *tGen = dynamic_cast<TTree*> (fGen->Get(str_in_MC_tree_gen.Data()));
        if(tGen) Printf("MC gen tree loaded.");
        
        ConnectTreeVariablesMCGen(tGen);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nPtBins; iPtBin++){
            Int_t NGen = 0;
            for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
                tGen->GetEntry(iEntry);
                // pT in a given bin 
                if(EventPassedMCGen(4, iPtBin)) NGen++;
            }
            AxE_PtBins_hNGen->SetBinContent(iPtBin, NGen);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        AxE_PtBins_SaveToFile(AxE_PtBins_hNGen, file);

        return;
    }
    return;
}

void AxE_PtBins_SaveToFile(TH1D* hist, TString name){
    ofstream outfile (name.Data());
    for(Int_t iBin = 1; iBin <= hist->GetNbinsX(); iBin++){
        outfile << iBin << "\t" << hist->GetBinContent(iBin) << "\n";
    }
    outfile.close();
    Printf("*** File saved in %s.***", name.Data());
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}