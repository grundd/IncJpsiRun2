// CrossSec_PrepareHistosAndGraphs.C
// David Grund, Sep 4, 2022

#include "CrossSec_Utilities.h"

TH1D *h_models_binned[7] = { NULL };
TGraph *gr_models_binned[7] = { NULL };
TGraph* gr_GSZ_err_binned[2] = { NULL };

void CrossSec_PrepareHistosAndGraphs(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    InitObjects();

    // load graphs
    LoadGraphs_data();
    LoadGraphs_SL();
    LoadGraphs_CCK();
    LoadGraphs_MS();
    LoadGraphs_GSZ();

    // integrate the data
    Printf("Data integral is: %.3f", IntegrateData()*1e3);

    // create histograms from the graphs and save them to a file
    TList *lh = new TList(); // list of histograms
    TList *lg = new TList(); // list of graphs
    lg->Add(gr_data_uncr);
    lg->Add(gr_data_corr);
    for(Int_t i = 0; i < 7; i++)
    {
        CreateHistogramFromGraph(i);
        lh->Add(h_models[i]);
        lg->Add(gr_models[i]);
    }
    lg->Add(gr_GSZ_err[0]);
    lg->Add(gr_GSZ_err[1]);

    // find the average value of |t| in each pT bin
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/");
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/AverageT/");

    for(Int_t i = 0; i < 7; i++)
    {
        // prepare the output file for the average values of |t|
        ofstream os;
        os.open("Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/AverageT/" + str_models[i] + ".txt");
        os << std::fixed << std::setprecision(4);
        // prepare binned graph and histogram
        h_models_binned[i] = new TH1D("hBinned_" + str_models[i], "hBinned_" + str_models[i], nPtBins, tBoundaries); 
        gr_models_binned[i] = new TGraph(nPtBins);
        gr_models_binned[i]->SetName("grBinned_" + str_models[i]);
        // go over pT bins
        for(Int_t iBin = 0; iBin < nPtBins; iBin++)
        {
            Double_t integral(0.), avgt(0.);
            IntegrateModel(i,tBoundaries[iBin],tBoundaries[iBin+1],integral,avgt);
            integral = integral / (tBoundaries[iBin+1] - tBoundaries[iBin]); // normalize by the bin-widths
            h_models_binned[i]->SetBinContent(iBin+1,integral);
            gr_models_binned[i]->SetPoint(iBin,avgt,integral);
            // print the results
            os << tBoundaries[iBin] << "\t" << tBoundaries[iBin+1] << "\t" << avgt << "\n";
        }
        os.close();
        lh->Add(h_models_binned[i]);
        lg->Add(gr_models_binned[i]);
    } 
    // prepare binned histograms with the errors of GSZ models
    gr_GSZ_err_binned[0] = new TGraph(2*2*nPtBins);
    gr_GSZ_err_binned[0]->SetName("grBinned_err_GSZ-el+diss");
    gr_GSZ_err_binned[1] = new TGraph(2*2*nPtBins);
    gr_GSZ_err_binned[1]->SetName("grBinned_err_GSZ-el");
    // fill graphs showing the error bands
    for (Int_t i = 0; i < nPtBins; i++)
    {
        gr_GSZ_err_binned[0]->SetPoint(2*i,   tBoundaries[i],   gr_models_binned[5]->GetPointY(i) * GSZ_err_scale_upp[0]);
        gr_GSZ_err_binned[0]->SetPoint(2*i+1, tBoundaries[i+1], gr_models_binned[5]->GetPointY(i) * GSZ_err_scale_upp[0]);
        gr_GSZ_err_binned[0]->SetPoint(2*nPtBins+2*i,   tBoundaries[nPtBins-i],   gr_models_binned[5]->GetPointY(nPtBins-i-1) * GSZ_err_scale_low[0]);
        gr_GSZ_err_binned[0]->SetPoint(2*nPtBins+2*i+1, tBoundaries[nPtBins-i-1], gr_models_binned[5]->GetPointY(nPtBins-i-1) * GSZ_err_scale_low[0]);
        gr_GSZ_err_binned[1]->SetPoint(2*i,   tBoundaries[i],   gr_models_binned[6]->GetPointY(i) * GSZ_err_scale_upp[1]);
        gr_GSZ_err_binned[1]->SetPoint(2*i+1, tBoundaries[i+1], gr_models_binned[6]->GetPointY(i) * GSZ_err_scale_upp[1]);
        gr_GSZ_err_binned[1]->SetPoint(2*nPtBins+2*i,   tBoundaries[nPtBins-i],   gr_models_binned[6]->GetPointY(nPtBins-i-1) * GSZ_err_scale_low[1]);
        gr_GSZ_err_binned[1]->SetPoint(2*nPtBins+2*i+1, tBoundaries[nPtBins-i-1], gr_models_binned[6]->GetPointY(nPtBins-i-1) * GSZ_err_scale_low[1]);
    }
    lg->Add(gr_GSZ_err_binned[0]);
    lg->Add(gr_GSZ_err_binned[1]);

    TFile *f = new TFile("Results/" + str_subfolder + "CrossSec/PrepareHistosAndGraphs/histograms_and_graphs.root","RECREATE");    
    lh->Write("histograms", TObject::kSingleKey);
    lg->Write("graphs", TObject::kSingleKey);
    f->ls();
    f->Close();

    return;
}