// PtFit_SubtractBkg.h
// David Grund, Feb 28, 2022

Double_t fPtCutLow_PtFit = 0.0;
Double_t fPtCutUpp_PtFit = 2.0;
RooBinning fPtBins_PtFit(fPtCutLow_PtFit, fPtCutUpp_PtFit);

vector<Double_t> BinsBoundaries_PtFit;
Double_t N_Jpsi_val, N_Jpsi_err;
Double_t N_Bkgr_val, N_Bkgr_err;

void PtFit_SetPtBinning();
void PtFit_PrepareData();
void PtFit_SubtractBkg();
void PtFit_DoInvMassFit(Double_t fPtCutLow, Double_t fPtCutUpp, Int_t iBin);

void PtFit_SubtractBkg_main(){

    gSystem->Exec("mkdir -p Trees/" + str_subfolder + "PtFit/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PtFit_SubtractBkg/InvMassFitsInBins/");

    // set pT binning
    PtFit_SetPtBinning();
    // prepare data tree before subtracting background
    PtFit_PrepareData();
    // subtract background in each of these bins
    PtFit_SubtractBkg();

    return;
}

void PtFit_SetPtBinning(){

    vector<Double_t> BinsWidths;
    vector<Double_t> BinsUpTo; // value of pT up to which we want the bins of a given width
    Int_t nBinsTypes = 0;
    // Define binning
    // Variable bin widths optimal for subtracting the background through inv mass fits
    nBinsTypes = 4;
    BinsWidths = {0.010, 0.020, 0.080, 0.400}; // GeV 
    BinsUpTo = {0.00, 0.20, 0.40, 1.20, 2.00}; // GeV 

    // Create a vector containing the calculated boundaries
    
    BinsBoundaries_PtFit.push_back(0.);
    Double_t fPtNow = 0.0;
    for(Int_t iBinType = 0; iBinType < nBinsTypes; iBinType++){
        Double_t nBinsThisType = (BinsUpTo[iBinType+1] - BinsUpTo[iBinType]) / BinsWidths[iBinType];
        Printf("Bin type %i: width %.3f GeV, from %.3f GeV to %.3f GeV, nBins = %.1f", 
            iBinType+1, BinsWidths[iBinType], BinsUpTo[iBinType], BinsUpTo[iBinType+1], nBinsThisType);
        Int_t iBin = 0;
        while(iBin < nBinsThisType){
            fPtNow += BinsWidths[iBinType];
            BinsBoundaries_PtFit.push_back(fPtNow);
            iBin++;
        }
    }
    // set the number of bins
    nPtBins_PtFit = BinsBoundaries_PtFit.size() - 1;
    // set the values of boundaries
    ptBoundaries_PtFit = &BinsBoundaries_PtFit[0];
    // set RooBinning
    for(Int_t i = 0; i < nPtBins_PtFit - 1; i++) fPtBins_PtFit.addBoundary(ptBoundaries_PtFit[i + 1]);
    // to print a detailed information
    if(kFALSE){
        Printf("%i bins created with the following boundaries:", nPtBins_PtFit);
        for(Int_t i = 0; i < nPtBins_PtFit + 1; i++) Printf("%.3f", ptBoundaries_PtFit[i]);
        Printf("RooBinning:");
        for(Int_t i = 0; i < nPtBins_PtFit; i++) Printf("%.3f", fPtBins_PtFit.binLow(i));
        Printf("%.3f", fPtBins_PtFit.binHigh(nPtBins_PtFit - 1));
    }

    return;
}

void PtFit_PrepareData()
{
    TString name = "Trees/" + str_subfolder + "PtFit/PtFit.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("Data tree already created.");
        return;

    } else {    

        Printf("Data tree will be created.");

        TFile *f_in = TFile::Open((str_in_DT_fldr + "AnalysisResults.root").Data(), "read");
        if(f_in) Printf("Input data loaded.");

        TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_in_DT_tree.Data()));
        if(t_in) Printf("Input tree loaded.");

        ConnectTreeVariables(t_in);

        // Create new data tree with applied cuts
        file = new TFile(name.Data(),"RECREATE");

        TTree *tPtFit = new TTree("tPtFit", "tPtFit");
        tPtFit->Branch("fPt", &fPt, "fPt/D");
        tPtFit->Branch("fM", &fM, "fM/D");
        tPtFit->Branch("fY", &fY, "fY/D");

        Printf("%lli entries found in the tree.", t_in->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
            t_in->GetEntry(iEntry);
            // m between 2.2 and 4.5 GeV/c^2, pT cut: all
            if(EventPassed(0,2)) tPtFit->Fill();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        file->Write("",TObject::kWriteDelete);

        return;
    }
}

void PtFit_SubtractBkg()
{
    TString name = "Trees/" + str_subfolder + "PtFit/SignalWithBkgSubtracted.root";
    TFile *file = TFile::Open(name.Data(),"read");
    if(file){
        Printf("Background already subtracted.");
        return;

    } else {
        
        Printf("Background will be subtracted.");

        TList *l = new TList();
        TH1D *hNSigPerBin = new TH1D("hNSigPerBin", "hNSigPerBin", nPtBins_PtFit, ptBoundaries_PtFit);
        TH1D *hNBkgPerBin = new TH1D("hNBkgPerBin", "hNBkgPerBin", nPtBins_PtFit, ptBoundaries_PtFit);
        l->Add(hNSigPerBin);
        l->Add(hNBkgPerBin);

        Double_t fPt = fPtCutLow_PtFit;
        Int_t iBin = 1;
        while(iBin <= nPtBins_PtFit){
            PtFit_DoInvMassFit(ptBoundaries_PtFit[iBin-1], ptBoundaries_PtFit[iBin], iBin);
            hNSigPerBin->SetBinContent(iBin, N_Jpsi_val);
            hNSigPerBin->SetBinError  (iBin, N_Jpsi_err);
            hNBkgPerBin->SetBinContent(iBin, N_Bkgr_val);
            hNBkgPerBin->SetBinError  (iBin, N_Bkgr_err);
            iBin++;
        }

        TCanvas *cHist_sig = new TCanvas("cHist_sig", "cHist_sig", 900, 600);
        cHist_sig->SetLogy();    
        hNSigPerBin->GetYaxis()->SetTitle("Counts per bin");
        hNSigPerBin->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
        hNSigPerBin->GetXaxis()->SetDecimals(1);
        // Draw the histogram
        //hNSigPerBin->Scale(1., "width");
        hNSigPerBin->Draw("E0");

        TCanvas *cHist_bkg = new TCanvas("cHist_bkg", "cHist_bkg", 900, 600);
        cHist_bkg->SetLogy();    
        hNBkgPerBin->GetYaxis()->SetTitle("Counts per bin");
        hNBkgPerBin->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
        hNBkgPerBin->GetXaxis()->SetDecimals(1);
        // Draw the histogram
        //hNSigPerBin->Scale(1., "width");
        hNBkgPerBin->Draw("E0");

        // Save results to the output file
        // Create the output file
        file = new TFile(name.Data(),"RECREATE");
        l->Write("HistList", TObject::kSingleKey);
        file->ls();

        TString name = "Results/" + str_subfolder + "PtFit_SubtractBkg/";
        cHist_sig->Print((name + "hNSigPerBins.pdf").Data());
        cHist_sig->Print((name + "hNSigPerBins.png").Data());
        cHist_bkg->Print((name + "hNBkgPerBins.pdf").Data());
        cHist_bkg->Print((name + "hNBkgPerBins.png").Data());

        delete cHist_sig;
        delete cHist_bkg;

        return;
    }
}

void PtFit_DoInvMassFit(Double_t fPtCutLow, Double_t fPtCutUpp, Int_t iBin)
{
    // Cuts:
    char fStrReduce[120];
    Double_t fYCut      = 0.80;
    Double_t fMCutLow   = 2.2;
    Double_t fMCutUpp   = 4.5;

    sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);

    // Binning:
    Int_t nBins = 115; // so that each bin between 2.2 and 4.5 GeV is 20 MeV wide
    RooBinning binM(nBins,fMCutLow,fMCutUpp);
    Double_t BinSizeDouble = (fMCutUpp - fMCutLow) * 1000 / nBins; // in MeV
    BinSizeDouble = BinSizeDouble + 0.5;
    // https://stackoverflow.com/questions/9695329/c-how-to-round-a-double-to-an-int
    Int_t BinSize = (Int_t)BinSizeDouble;

    Printf("\n");
    Printf("*** Bin size (double): %.3f ***", BinSizeDouble);
    Printf("*** Bin size (int): %i ***\n", BinSize);
  
    // Roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0,10.);
    RooRealVar fY("fY","fY",-0.8,0.8);

    // Get the data trees
    TFile *f_in = new TFile("Trees/" + str_subfolder + "PtFit/PtFit.root"); 
    TTree *t_in = NULL;
    f_in->GetObject("tPtFit",t_in);
        
    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*t_in));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Crystal Ball parameters from MC (to be fixed)
    Double_t fAlpha_L;
    Double_t fAlpha_R;
    Double_t fN_L;
    Double_t fN_R;

    char name[20];
    Double_t values[4];
    Double_t errors[4];

    TString pathMC;
    if(fPtCutLow < 0.20) pathMC = "Results/" + str_subfolder + "InvMassFit_MC/coh/coh.txt";
    else                 pathMC = "Results/" + str_subfolder + "InvMassFit_MC/inc/inc.txt";

    ifstream f_txt_in;
    f_txt_in.open(pathMC.Data());
    if(f_txt_in.fail()){
        Printf("\n");
        Printf("*** Warning! ***");
        Printf("*** MC values for tail parameters not found. Terminating... *** \n");
        return;
    } else {
        Int_t i_line = 0;
        while(!f_txt_in.eof()){
            f_txt_in >> name >> values[i_line] >> errors[i_line];
            i_line++;
        }
        f_txt_in.close();
    }

    fAlpha_L = values[0];
    fAlpha_R = values[1];
    fN_L = values[2];
    fN_R = values[3];

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,20.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,20.);
    alpha_L.setConstant(kTRUE);
    alpha_R.setConstant(kTRUE);
    n_L.setConstant(kTRUE);
    n_R.setConstant(kTRUE);

    // Crystal Ball for J/Psi
    RooRealVar mass_Jpsi("mass_Jpsi","J/psi mass",3.097,3.00,3.20);
    Double_t sigma;
    if(fPtCutLow < 0.20) sigma = 0.020;
    else                 sigma = 0.021;
    RooRealVar sigma_Jpsi("sigma_Jpsi","J/psi resolution",sigma,0.01,0.1);
    sigma_Jpsi.setConstant(kTRUE);
    RooGenericPdf mean_R("mean_R","J/psi mass","mass_Jpsi",RooArgSet(mass_Jpsi));
    RooGenericPdf sigma_R("sigma_R","J/psi resolution","sigma_Jpsi",RooArgSet(sigma_Jpsi));
    RooRealVar N_Jpsi("N_Jpsi","number of J/psi events",0.4*nEvents,0,nEvents);

    // Background
    RooRealVar lambda("lambda","background exp",-1.2,-10.,0.);
    RooRealVar N_bkg("N_bkg","number of background events",0.6*nEvents,0,nEvents);

    // Functions for fitting
    // J/psi:
    RooCBShape CB_left("CB_left","CB_left",fM,mass_Jpsi,sigma_Jpsi,alpha_L,n_L);
    RooCBShape CB_right("CB_right","CB_right",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_left,CB_right),RooArgList(frac));
    // Background:
    RooGenericPdf BkgPdf("BkgPdf","exp(fM*lambda)",RooArgSet(fM,lambda));

    // Create Model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    Double_t N_Jpsi_out[2];
    Double_t N_Bkgr_out[2];
    fM.setRange("JpsiMassRange",3.0,3.2);
    // Calculate the number of J/psi events
    RooAbsReal *intDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_Jpsi_out[0] = intDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_out[1] = intDSCB->getVal()*N_Jpsi.getError();
    N_Jpsi_val = N_Jpsi_out[0];
    N_Jpsi_err = N_Jpsi_out[1];
    // Calculate the number of Bkgr events
    RooAbsReal *intBkgr = BkgPdf.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_Bkgr_out[0] = intBkgr->getVal()*N_bkg.getVal();
    N_Bkgr_out[1] = intBkgr->getVal()*N_bkg.getError();
    N_Bkgr_val = N_Bkgr_out[0];
    N_Bkgr_err = N_Bkgr_out[1];

    // ##########################################################
    // Plot the results
    gStyle->SetOptTitle(0); // suppress title
    gStyle->SetOptStat(0);  // the type of information printed in the histogram statistics box
                            // 0 = no information

    // Draw histogram with fit results
    TCanvas *cHist = new TCanvas("cHist","cHist",800,600);
    cHist->SetTopMargin(0.055);
    cHist->SetBottomMargin(0.12);
    cHist->SetRightMargin(0.03);
    cHist->SetLeftMargin(0.11);

    RooPlot* fFrameM = fM.frame(Title("Mass fit")); 
    fDataSet->plotOn(fFrameM,Name("fDataSet"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DoubleSidedCB"),Components(DoubleSidedCB),LineColor(kBlack),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("BkgPdf"),Components(BkgPdf),LineColor(kRed),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DSCBAndBkgPdf"),LineColor(215),LineWidth(3));
    // Vertical axis
    fFrameM->GetYaxis()->SetTitle(Form("Counts per %i MeV/#it{c}^{2}", BinSize));
    fFrameM->GetYaxis()->SetTitleSize(0.05);
    fFrameM->GetYaxis()->SetTitleOffset(1.1);
    fFrameM->GetYaxis()->SetLabelSize(0.05);
    fFrameM->GetYaxis()->SetLabelOffset(0.01);
    fFrameM->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    fFrameM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    fFrameM->GetXaxis()->SetTitleSize(0.05);
    fFrameM->GetXaxis()->SetLabelSize(0.05);
    fFrameM->GetXaxis()->SetDecimals(1);
    fFrameM->Draw();

    // Get chi2 
    Double_t chi2 = fFrameM->chiSquare("DSCBAndBkgPdf","fDataSet",fResFit->floatParsFinal().getSize()); // last argument = number of parameters
    Printf("********************");
    Printf("chi2/NDF = %.3f", chi2);
    Printf("NDF = %i", fResFit->floatParsFinal().getSize());
    Printf("chi2/NDF = %.3f/%i", chi2*fResFit->floatParsFinal().getSize(), fResFit->floatParsFinal().getSize());
    Printf("********************");   

    // -------------------------------------------------------------------------------- 
    // Legend1
    TLegend *l1 = new TLegend(0.09,0.76,0.3,0.935);
    //l1->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < %.1f", fYCut),"");
    // Print the pt cut
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    l1->SetTextSize(0.040);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    TLegend *lTitle = new TLegend(0.325,0.88,0.95,0.935);
    lTitle->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    lTitle->SetTextSize(0.05);
    lTitle->SetBorderSize(0);
    lTitle->SetFillStyle(0);
    lTitle->Draw();

    // Legend2
    TLegend *l2 = new TLegend(0.52,0.29,0.95,0.87);
    l2->SetMargin(0.14);
    l2->AddEntry("DSCBAndBkgPdf","sum","L");
    l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l2->AddEntry("DoubleSidedCB","J/#psi signal","L");
    l2->AddEntry((TObject*)0,Form("#it{N}_{J/#psi} = %.0f #pm %.0f",N_Jpsi_out[0],N_Jpsi_out[1]),"");
    l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
    l2->AddEntry((TObject*)0,Form("#sigma = %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.2f", alpha_L.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{R} = %.2f", (-1)*(alpha_R.getVal())),"");
    l2->AddEntry("BkgPdf","background","L");
    l2->AddEntry((TObject*)0,Form("#lambda = %.2f #pm %.2f GeV^{-1}#it{c}^{2}",lambda.getVal(), lambda.getError()),"");
    l2->AddEntry((TObject*)0,"with #it{m}_{#mu#mu} #in (3.0,3.2) GeV/#it{c}^{2}:","");
    l2->AddEntry((TObject*)0,Form("#it{N}_{bkg} = %.0f #pm %.0f",N_Bkgr_out[0],N_Bkgr_out[1]),"");
    l2->SetTextSize(0.040); // was 0.042
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    TLegend *l3 = new TLegend(0.74,0.48,0.85,0.58);
    l3->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f", n_L.getVal()),"");
    l3->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f", n_R.getVal()),"");
    l3->SetTextSize(0.040); // was 0.042
    l3->SetBorderSize(0);
    l3->SetFillStyle(0);
    l3->Draw();
    // Prepare path
    TString path = Form("Results/%s/PtFit_SubtractBkg/InvMassFitsInBins/bin%i", str_subfolder.Data(), iBin);
    // Print the plots
    cHist->Print((path + ".pdf").Data()); 
    cHist->Print((path + ".png").Data()); 

    delete cHist;

    return;
}