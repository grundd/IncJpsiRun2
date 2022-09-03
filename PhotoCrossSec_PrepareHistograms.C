// PhotoCrossSec_PrepareHistograms.C
// David Grund, July 14, 2022

// root headers
#include "TGraph.h"
// my headers
#include "PhotoCrossSec_Utilities.h"
#include "SetPtBinning.h"

Double_t *tBoundaries = NULL;
Double_t tBoundaries_4bins[5] = { 0 };
Double_t tBoundaries_5bins[6] = { 0 };

TH1D *h_models[9] = { NULL };
TH1D *h_models_binned[9] = { NULL };
// order: SL, CCK_hs, CCK_n, MS_fl, MS_nf, GSZ_tot_max, GSZ_tot_min, GSZ_el_max, GSZ_el_min
TString path;

void PrepareModelHistograms(TString str_name, Int_t n_data, Double_t *abs_t_val, Double_t *sig_val);
void FindAverageAbsTPerBin(TH1D *h, TString str);
void PrintIntegralValue(Double_t integral, TString str);
void PrintHistogramContent(TH1D *h, TString str);
TH1D* LoadHistograms(TString str);

void PhotoCrossSec_PrepareHistograms(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    SetPtBinning();

    if(nPtBins == 4) tBoundaries = tBoundaries_4bins;
    if(nPtBins == 5) tBoundaries = tBoundaries_5bins;
    for(Int_t i = 0; i <= nPtBins; i++) tBoundaries[i] = ptBoundaries[i] * ptBoundaries[i];

    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/continuous/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/integrals/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/histogram_plots/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/print_histograms/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/avg_t/");

    // create and store the histograms
    ReadInput_data();
    Double_t integral_data = IntegrateData(1.0); // in micro barns
    integral_data = integral_data / 1e3; // in mili barns
    PrintIntegralValue(integral_data,str_data);

    ReadInput_SL(); 
    PrepareModelHistograms(str_models[0],n_SL,abs_t_SL,sig_SL);

    ReadInput_CCK();
    PrepareModelHistograms(str_models[1],n_CCK,abs_t_CCK,sig_CCK_inc_hs);
    PrepareModelHistograms(str_models[2], n_CCK,abs_t_CCK,sig_CCK_inc_n);

    ReadInput_MS(); 
    PrepareModelHistograms(str_models[3],n_MS,abs_t_MS,sig_MS_fluct);
    PrepareModelHistograms(str_models[4],n_MS,abs_t_MS,sig_MS_noflu);

    ReadInput_GSZ();
    for(Int_t i = 0; i < n_GSZ; i++){
        // GSZ uses nb instead of mb
        sig_GSZ_tot_min[i] = sig_GSZ_tot_min[i] / 1e6;
        sig_GSZ_tot_max[i] = sig_GSZ_tot_max[i] / 1e6;
        sig_GSZ_el_min[i] = sig_GSZ_el_min[i] / 1e6;
        sig_GSZ_el_max[i] = sig_GSZ_el_max[i] / 1e6;
    }
    PrepareModelHistograms(str_models[5],n_GSZ,abs_t_GSZ,sig_GSZ_tot_max);
    PrepareModelHistograms(str_models[6],n_GSZ,abs_t_GSZ,sig_GSZ_tot_min);
    PrepareModelHistograms(str_models[7],n_GSZ,abs_t_GSZ,sig_GSZ_el_max);
    PrepareModelHistograms(str_models[8],n_GSZ,abs_t_GSZ,sig_GSZ_el_min);

    // load the histograms
    for(Int_t i = 0; i < 9; i++)
    {
        h_models[i] = LoadHistograms("continuous/" + str_models[i]);
        FindAverageAbsTPerBin(h_models[i],str_models[i]);
        h_models_binned[i] = LoadHistograms("binned/" + str_models[i]);
    } 

    // save the histograms with continuous model predictions
    TList *l_continuous = new TList();
    for(Int_t i = 0; i < 9; i++) l_continuous->Add(h_models[i]);
    path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/continuous/all.root";
    TFile *f_continuous = new TFile(path.Data(),"RECREATE");
    l_continuous->Write("HistList", TObject::kSingleKey);
    f_continuous->ls();
    f_continuous->Close();
    
    // save the histograms with binned model predictions
    TList *l_binned = new TList();
    for(Int_t i = 0; i < 9; i++) l_binned->Add(h_models_binned[i]);
    path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/all.root";
    TFile *f_binned = new TFile(path.Data(),"RECREATE");
    l_binned->Write("HistList", TObject::kSingleKey);
    f_binned->ls();
    f_binned->Close();

    return;
}

void PrepareModelHistograms(TString str_name, Int_t n_data, Double_t *abs_t_val, Double_t *sig_val)
{
    TCanvas *c_tot = new TCanvas("c_tot","c_tot",900,600);
    // integrate the model over the whole range: 0.04 to 1.00 and save the value
    // prepare histogram with the model curve and store it in a root file
    Double_t int_tot = GraphIntegral_Calculate(c_tot,"continuous/" + str_name,n_data,abs_t_val,sig_val,0.04,1.00,kTRUE); // in mili barns
    PrintIntegralValue(int_tot,str_name);
    // integrate the model in the pT bins
    // prepare the histogram with the pT binning from the analysis
    int_tot = 0.;
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/" + str_name + "/");
    TH1D *h_binned = new TH1D("hist", "hist", nPtBins, tBoundaries);
    TCanvas *c[5] = { NULL };
    Double_t integrals[5] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++)
    {
        c[iBin] = new TCanvas(Form("uptobin%i",iBin+1),Form("uptobin%i",iBin+1),900,600);
        Double_t t_low = tBoundaries[iBin];
        Double_t t_upp = tBoundaries[iBin+1];
        integrals[iBin] = GraphIntegral_Calculate(c[iBin],str_name,n_data,abs_t_val,sig_val,0.04,t_upp); // in mili barns
        path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/" + str_name + Form("/uptobin%i.pdf",iBin+1); 
        c[iBin]->Print(path.Data());
        integrals[iBin] -= GraphIntegral_Calculate(c[iBin],str_name,n_data,abs_t_val,sig_val,0.04,t_low); // in mili barns
        delete c[iBin];
        h_binned->SetBinContent(iBin+1,integrals[iBin]);
        int_tot += integrals[iBin];
    }
    Printf("%s: sum of integrals over bins: %.2f", str_name.Data(), int_tot * 1e3);
    h_binned->Scale(1.,"width"); // normalize by the bin widths
    TCanvas *c_hist = new TCanvas("c_hist","c_hist",900,600);
    c_hist->SetLogy();
    h_binned->Draw("HIST");
    path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/histogram_plots/" + str_name + ".pdf";
    c_hist->Print(path.Data());
    // save the binned histogram to a root file
    TList *l = new TList();
    l->Add(h_binned);
    path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/binned/" + str_name + ".root";
    TFile *file = new TFile(path.Data(),"RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    file->ls();
    file->Close();

    delete c_tot;
    delete c_hist;
    return;
}

void FindAverageAbsTPerBin(TH1D *h, TString str)
{
    // print original histogram: edges and bin contents
    PrintHistogramContent(h,str);
    vector<Double_t> edges;
    vector<Double_t> sigma;
    Int_t iPtBinCurrent = 1; // index of a bin inside of which we currently are
    for(Int_t i = 1; i <= h->GetNbinsX(); i++)
    {
        Double_t tLowCurrent = h->GetBinLowEdge(i);
        if(tLowCurrent > tBoundaries[iPtBinCurrent]){
            edges.push_back(tBoundaries[iPtBinCurrent]);
            sigma.push_back(h->GetBinContent(i-1));
            iPtBinCurrent++;
        }
        edges.push_back(tLowCurrent);
        sigma.push_back(h->GetBinContent(i));
    }
    edges.push_back(h->GetBinLowEdge(h->GetNbinsX() + 1));
    Double_t *edges_arr = &edges[0];
    TH1D *h_new = new TH1D("h_new","h_new",sigma.size(),edges_arr);
    for(Int_t i = 1; i <= h_new->GetNbinsX(); i++) h_new->SetBinContent(i,sigma[i-1]);
    PrintHistogramContent(h_new,str + "_new");
    // find the average |t| per bin according to the models
    iPtBinCurrent = 1;
    Double_t t_max_bin = tBoundaries[iPtBinCurrent];
    Double_t bin_int = 0;
    Double_t bin_weighted_int = 0;
    vector<Double_t> bins_int;
    vector<Double_t> bins_weighted_int;
    vector<Double_t> bins_avgT;
    for(Int_t i = 1; i <= h_new->GetNbinsX(); i++) 
    {
        if(h_new->GetBinCenter(i) > t_max_bin)
        {
            iPtBinCurrent++;
            t_max_bin = tBoundaries[iPtBinCurrent];
            bins_int.push_back(bin_int);
            bins_weighted_int.push_back(bin_weighted_int);
            bins_avgT.push_back(bin_weighted_int / bin_int);
            bin_int = 0;
            bin_weighted_int = 0;
        }
        bin_int += h_new->GetBinContent(i) * (h_new->GetBinLowEdge(i+1) - h_new->GetBinLowEdge(i));
        bin_weighted_int += h_new->GetBinCenter(i) * h_new->GetBinContent(i) * (h_new->GetBinLowEdge(i+1) - h_new->GetBinLowEdge(i));
    }
    // save the values from the last bin
    bins_int.push_back(bin_int);
    bins_weighted_int.push_back(bin_weighted_int);
    bins_avgT.push_back(bin_weighted_int / bin_int);
    // print the results
    Double_t sum_int(0.), sum_weighted_int(0.);
    TString path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/avg_t/info_" + str + ".txt";
    ofstream outfile(path.Data());
    outfile << "[in micro barns]\n";
    outfile << "bin\twInt\tint\tavgT\n";
    outfile << std::fixed << std::setprecision(3);
    for(Int_t i = 1; i <= nPtBins; i++)
    {
        sum_int += bins_int[i-1] * 1e3; // in micro barns
        sum_weighted_int += bins_weighted_int[i-1] * 1e3; // in micro barns
        outfile << i << "\t" << bins_weighted_int[i-1] * 1e3 << "\t" << bins_int[i-1] * 1e3 << "\t" << bins_avgT[i-1] << "\n"; 
    }
    outfile << "sum:\t" << sum_weighted_int << "\t" << sum_int << "\n";
    outfile.close();
    // print the values of avg |t| only (to be read by PhotoCrossSec_PlotWithRatios.C ...)
    bins_avgT.push_back(1.);
    bins_avgT.push_back(1.);
    bins_avgT.push_back(1.);
    bins_avgT.push_back(1.);
    bins_avgT.push_back(1.);
    path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/avg_t/" + str + ".txt";
    outfile.open(path.Data());
    outfile << std::fixed << std::setprecision(3);
    for(Int_t i = 1; i <= nPtBins; i++)
    {
        outfile << bins_avgT[i-1] << "\n"; 
    }
    outfile.close();

    //delete h_new;
    return;
}

void PrintIntegralValue(Double_t integral, TString str) // print the integral value in micro barns
{
    TString path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/integrals/" + str + ".txt";
    ofstream outfile(path.Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << integral * 1e3; 
    outfile.close();
    return;
}

void PrintHistogramContent(TH1D *h, TString str)
{
    TString path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/print_histograms/" + str + ".txt";
    ofstream outfile(path.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("This histogram has %i bins.\n", h->GetNbinsX());
    outfile << "bin\tt_low\tt_upp\tsigma\n";
    for(Int_t i = 1; i <= h->GetNbinsX(); i++)
    {
        outfile << i << "\t" << h->GetBinLowEdge(i) << "\t" << h->GetBinLowEdge(i+1) << "\t" << h->GetBinContent(i) << "\n";
    }
    outfile.close();
    return;
}

TH1D* LoadHistograms(TString str)
{
    TString path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/" + str + ".root";
    TFile *f = TFile::Open(path.Data(),"read");
    if(f) Printf("Input file %s loaded.", f->GetName()); 
    TList *l = (TList*) f->Get("HistList");
    if(l) Printf("List %s loaded.", l->GetName()); 
    TH1D *h = (TH1D*)l->FindObject("hist");
    h->SetNameTitle(str,str);
    if(h) Printf("Histogram %s loaded.", h->GetName());
    Printf("Integral value: %.2f microb", h->Integral("width") * 1e3);
    return h;
}