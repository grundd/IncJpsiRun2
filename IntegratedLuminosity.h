// IntegratedLuminosity.h
// David Grund, Feb 26, 2022
// To calculate the integrated luminosity of the analyzed samples using the trending files

// root headers
#include "TObjArray.h"
// aliroot headers
#include "AliTriggerClass.h"

// Lists of good run numbers
// LHC18q, 123 runs
Int_t RunList18q[123] = {
    295585, 295586, 295588, 295589, 295610, 295611, 295612, 295615, 295666, 295667, 295668, 295673, 295675, 
    295676, 295712, 295714, 295717, 295718, 295719, 295721, 295723, 295725, 295754, 295755, 295758, 295759, 
    295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829, 295831, 
    295853, 295854, 295855, 295856, 295859, 295860, 295861, 295909, 295910, 295913, 295936, 295937, 295941, 
    295942, 296016, 296060, 296062, 296063, 296065, 296066, 296123, 296132, 296133, 296134, 296135, 296142, 
    296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 
    296246, 296247, 296269, 296270, 296273, 296279, 296280, 296303, 296304, 296309, 296312, 296377, 296378, 
    296379, 296380, 296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472, 296509, 
    296510, 296511, 296512, 296516, 296547, 296548, 296549, 296550, 296551, 296552, 296553, 296594, 296615, 
    296616, 296618, 296619, 296621, 296622, 296623
};
// LHC18r, 96 runs
Int_t RunList18r[96] = {
    296690, 296691, 296693, 296694, 296749, 296750, 296781, 296784, 296785, 296786, 296787, 296790, 296793, 
    296794, 296799, 296835, 296836, 296838, 296839, 296848, 296849, 296850, 296851, 296852, 296890, 296894, 
    296899, 296900, 296903, 296930, 296931, 296932, 296934, 296935, 296938, 296941, 296966, 297029, 297031, 
    297035, 297085, 297117, 297118, 297119, 297123, 297124, 297128, 297129, 297132, 297133, 297193, 297194, 
    297195, 297196, 297218, 297219, 297221, 297222, 297278, 297310, 297311, 297317, 297332, 297333, 297335, 
    297336, 297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414, 297415, 297441, 
    297442, 297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540, 297541, 297542, 
    297544, 297558, 297588, 297590, 297595
};

// Arrays with counters of fired triggers per run
// LHC18q:
// In increasing order (from 295585 to 296623)
Int_t Counts18q[123] = { 0 };
// LHC18r:
// In increasing order (from 296690 to 297595)
Int_t Counts18r[96] = { 0 };

Int_t nRunsInList = 0;
Int_t *RunList = NULL;
Int_t *CountsList = NULL;

void CalculateLumi(Int_t period);
void SumLumi(Int_t period);
void SetLumiHisto(TH1D* h, Int_t period, Color_t color);

void IntegratedLuminosity_main(){

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "Lumi/");

    // Load trigger counters
    ifstream ifs;
    ifs.open(("Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18q.txt").Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < 123; i++){
            Int_t iCount;
            ifs >> iCount;
            Counts18q[i] = iCount;
        }
        ifs.close(); 
    }
    ifs.open(("Results/" + str_subfolder + "GetTriggerCounters/trigger_counters_LHC18r.txt").Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < 96; i++){
            Int_t iCount;
            ifs >> iCount;
            Counts18r[i] = iCount;
        }
        ifs.close(); 
    }

    // LHC18q
    CalculateLumi(0);
    // LHC18r
    CalculateLumi(1);

    return;
}

void CalculateLumi(Int_t period){

    // Choose the period
    if(period == 0){ 
        // LHC18q
        RunList = RunList18q;
        CountsList = Counts18q;
        nRunsInList = 123;
    } else if(period == 1){ 
        // LHC18r
        RunList = RunList18r;
        CountsList = Counts18r;
        nRunsInList = 96;
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

void SumLumi(Int_t period){
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

    // Also, we want to sum the luminosity corresponding to ClassName1 (CCUP31-B-NOPF-CENTNOTRD), i.e., 
    // to the situation without the past-future protection
    Double_t int_lumi_ana_NOPF(0.), int_lumi_rec_NOPF(0.);
    for(Int_t i(0); i < nRunsInList; i++){
        if(RunList[i] < 295881){
            int_lumi_ana_NOPF += hLumiS->GetBinContent(hLumiS->GetXaxis()->FindBin(Form("%i",RunList[i])));
            int_lumi_rec_NOPF += hLumi->GetBinContent(hLumi->GetXaxis()->FindBin(Form("%i",RunList[i])));
        }
    }
    TString str_NOPF_file = "";
    if(period == 0) str_NOPF_file = "Results/" + str_subfolder + "Lumi/lumi_NOPF_18q.txt";
    if(period == 1) str_NOPF_file = "Results/" + str_subfolder + "Lumi/lumi_NOPF_18r.txt";
    ofstream NOPF_file(str_NOPF_file.Data());
    NOPF_file << "Integrated luminosity corresponding to CCUP31-B-NOPF-CENTNOTRD (run number < 295881):\n";
    NOPF_file << Form("analyzed: %.3f (%.1f percent of the total lumi for this period) \n",int_lumi_ana_NOPF, int_lumi_ana_NOPF / int_lumi_ana * 100.);
    NOPF_file << Form("recorded: %.3f (%.1f percent of the total lumi for this period) \n",int_lumi_rec_NOPF, int_lumi_rec_NOPF / int_lumi_rec * 100.);
    NOPF_file.close();
    Printf("*** Results printed to %s.***", str_NOPF_file.Data());


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
    if(period == 0){ 
        // LHC18q
        cLumi->SaveAs("Results/" + str_subfolder + "Lumi/lumi_ccup31_18q.pdf");
        cLumi->SaveAs("Results/" + str_subfolder + "Lumi/lumi_ccup31_18q.png");
    } else if(period == 1){
        // LHC18r
        cLumi->SaveAs("Results/" + str_subfolder + "Lumi/lumi_ccup31_18r.pdf");
        cLumi->SaveAs("Results/" + str_subfolder + "Lumi/lumi_ccup31_18r.png");   
    }
}

void SetLumiHisto(TH1D* h, Int_t period, Color_t color){
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
        h->GetYaxis()->SetTitleOffset(0.45);
    } else if(period == 1){
        h->GetYaxis()->SetTitleOffset(0.3);
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