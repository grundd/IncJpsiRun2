// PhotoCrossSec_ExpFits.C
// David Grund, Aug 27, 2022

// root headers
#include "TSystem.h"
// my headers
#include "PhotoCrossSec_Utilities.h"

TH1D *h_models[9] = { NULL };
Double_t avg_t_models[9][5] = { 0 };
TGraph *gr_models[9] = { NULL };
TGraphErrors *gr_data = NULL;
Double_t par_pure_val[2];
Double_t par_pure_err[2];
Double_t par_quad_val[3];
Double_t par_quad_err[3];
TString binning;

void DoExpFit(Int_t iModel);
// iModel == 0 ... 8 => models, == 9 => data

void PhotoCrossSec_ExpFits(Int_t iAnalysis, Int_t iBinning)
{
    InitAnalysis(iAnalysis);

    if(iBinning == 0)      binning = "continuous/";
    else if(iBinning == 1) binning = "binned/";
    else return;
    gSystem->Exec("mkdir -p Results/" + str_subfolder + "PhotoCrossSec/ExpFits/" + binning);

    // load data
    ReadInput_data();
    Double_t sig_err_tot[5] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        // from microbarns to milibarns
        sig_val[i] = sig_val[i] / 1e3;
        sig_err_stat[i] = sig_err_stat[i] / 1e3;
        sig_err_syst[i] = sig_err_syst[i] / 1e3;
        sig_err_tot[i] = TMath::Sqrt(TMath::Power(sig_err_stat[i],2) + TMath::Power(sig_err_syst[i],2));
        Printf("Cross section in bin %i:\n", i+1);
        Printf("%.5f pm %.5f(stat.) pm %.5f(syst.).", sig_val[i], sig_err_stat[i], sig_err_syst[i]);
        Printf("%.5f pm %.5f(tot.).", sig_val[i], sig_err_tot[i]);
    }
    gr_data = new TGraphErrors(nPtBins,abs_t_val,sig_val,NULL,sig_err_tot);

    // load histograms
    TString path = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/" + binning + "all.root";
    TFile *f = TFile::Open(path.Data(),"read");
    if(f) Printf("Input file %s loaded.", f->GetName()); 
    TList *l = (TList*) f->Get("HistList");
    if(l) Printf("List %s loaded.", l->GetName()); 
    for(Int_t i = 0; i < 9; i++) h_models[i] = (TH1D*)l->FindObject(binning + str_models[i]);
    // load average values of |t|
    ifstream ifs;
    for(Int_t i = 0; i < 9; i++)
    {
        TString str = "Results/" + str_subfolder + "PhotoCrossSec/PrepareHistograms/avg_t/" + str_models[i] + ".txt";
        ifs.open(str.Data());
        if(!ifs.fail())
        {
            Printf("Loaded avg |t| in pT bins for %s:", str_models[i].Data());
            for(Int_t iBin = 0; iBin < nPtBins; iBin++)
            {  
                ifs >> avg_t_models[i][iBin];
                Printf("bin %i: %.3f", iBin+1, avg_t_models[i][iBin]);
            } 
            ifs.close();
        }
        else 
        {
            Printf("File %s not found. Terminating...", path.Data());
            return;
        }
    }    
    // create graphs
    for(Int_t i = 0; i < 9; i++)
    {
        if(iBinning == 0) gr_models[i] = new TGraph(h_models[i]);
        else if(iBinning == 1)
        {
            gr_models[i] = new TGraph(nPtBins); // a new TGraph with nPtBins points
            for(Int_t iBin = 0; iBin < nPtBins; iBin++)
            {
                gr_models[i]->SetPoint(iBin,avg_t_models[i][iBin],h_models[i]->GetBinContent(iBin+1));
            }
        }
        else return;
        
    } 

    // do the fits and print the results
    ofstream outfile;
    outfile.open("Results/" + str_subfolder + "PhotoCrossSec/ExpFits/" + binning + "#parameters.txt");
    outfile << std::fixed << std::setprecision(2);
    outfile << "\t\tpure exponential \t\texp with a quadratic term\n";
    outfile << "model \t\tA \terr \tb \terr \tA \terr \tb \terr \tc \terr \n";
    DoExpFit(9);
    outfile << "data\t\t" << par_pure_val[0] * 1e3 << "\t" << par_pure_err[0] * 1e3 << "\t"
                        << par_pure_val[1] << "\t" << par_pure_err[1] << "\t"
                        << par_quad_val[0] * 1e3 << "\t" << par_quad_err[0] * 1e3 << "\t"
                        << par_quad_val[1] << "\t" << par_quad_err[1] << "\t"
                        << par_quad_val[2] << "\t" << par_quad_err[2] << "\n";
    // do exp fits
    for(Int_t i = 0; i < 9; i++){
        DoExpFit(i);
        outfile << Form("%s\t", str_models[i].Data());
        if(i == 3 || i == 4) outfile << "\t";
        outfile << par_pure_val[0] * 1e3 << "\t" << par_pure_err[0] * 1e3 << "\t"
                << par_pure_val[1] << "\t" << par_pure_err[1] << "\t"
                << par_quad_val[0] * 1e3 << "\t" << par_quad_err[0] * 1e3 << "\t"
                << par_quad_val[1] << "\t" << par_quad_err[1] << "\t"
                << par_quad_val[2] << "\t" << par_quad_err[2] << "\n";
    } 
    outfile.close();

    return;
}

void DoExpFit(Int_t iModel)
{
    TGraph *gr = NULL;
    if(iModel == 9) gr = gr_data;
    else            gr = gr_models[iModel];

    TF1 *f_exp_pure = new TF1("f_exp_pure", "[0] * exp(-[1] * x)");
    TF1 *f_exp_quad = new TF1("f_exp_quad", "[0] * exp(-[1] * x + [2] * x * x)");
    // set function properties
    f_exp_pure->SetLineStyle(1);
    f_exp_pure->SetLineColor(kRed);
    f_exp_pure->SetLineWidth(2);
    f_exp_quad->SetLineStyle(1);
    f_exp_quad->SetLineColor(215);
    f_exp_quad->SetLineWidth(2);
    gr->Fit(f_exp_pure);
    gr->Fit(f_exp_quad);
    for(Int_t ipar = 0; ipar < 3; ipar++){
        if(ipar != 2){
            par_pure_val[ipar] = f_exp_pure->GetParameter(ipar);
            par_pure_err[ipar] = f_exp_pure->GetParError(ipar);
        }
        par_quad_val[ipar] = f_exp_quad->GetParameter(ipar);
        par_quad_err[ipar] = f_exp_quad->GetParError(ipar);
    }
    // set data properties 
    gStyle->SetEndErrorSize(4);         
    gr->SetMarkerStyle(kFullCross);
    gr->SetMarkerSize(0.7);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kBlack);
    // canvas
    TCanvas *cFit = new TCanvas("cFit","cFit",900,600);
    // margins
    cFit->SetTopMargin(0.03);
    cFit->SetBottomMargin(0.14);
    cFit->SetRightMargin(0.03);
    cFit->SetLeftMargin(0.11);
    cFit->SetLogy();
    // plot the graph
    TH1 *h = (TH1*) gr->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
    // vertical axis
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.15);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetDecimals(1);
    //h->SetMinimum(0.00002);
    //h->SetMaximum(0.04);
    // horizontal axis
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.045);
    //h->GetXaxis()->SetRangeUser(0.04,1.00);
    // draw graph and curves
    gr_models[iModel]->Draw("AP SAME");
    f_exp_pure->Draw("AL SAME");
    f_exp_quad->Draw("AL SAME");

    TString name;
    if(iModel == 9) name = "data";
    else            name = str_models[iModel];

    TLegend *l = new TLegend(0.12,0.18,0.40,0.50);
    l->AddEntry((TObject*)0,name,"");
    l->AddEntry(f_exp_quad,"A exp(-bt + ct^{2})","L");
    l->AddEntry((TObject*)0,Form("A = (%.2f #pm %.2f) #mub #it{c}^{2} GeV^{-2}", par_quad_val[0] * 1e3, par_quad_err[0] * 1e3),"");
    l->AddEntry((TObject*)0,Form("b = (%.2f #pm %.2f) #it{c}^{2} GeV^{-2}", par_quad_val[1], par_quad_err[1]),"");
    l->AddEntry((TObject*)0,Form("c = (%.2f #pm %.2f) #it{c}^{4} GeV^{-4}", par_quad_val[2], par_quad_err[2]),"");
    l->AddEntry((TObject*)0,Form("integral 0.04 < |t| < 1.00 GeV^{2} #it{c}^{-2}: %.2f #mub", f_exp_quad->Integral(0.04,1.00) * 1e3),"");
    l->AddEntry(f_exp_pure,"A exp(-bt)","L");
    l->AddEntry((TObject*)0,Form("A = (%.2f #pm %.2f) #mub #it{c}^{2} GeV^{-2}", par_pure_val[0] * 1e3, par_pure_err[0] * 1e3),"");
    l->AddEntry((TObject*)0,Form("b = (%.2f #pm %.2f) #it{c}^{2} GeV^{-2}", par_pure_val[1], par_pure_err[1]),"");
    l->AddEntry((TObject*)0,Form("integral 0.04 < |t| < 1.00 GeV^{2} #it{c}^{-2}: %.2f #mub", f_exp_pure->Integral(0.04,1.00) * 1e3),"");
    l->SetTextSize(0.032);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    if(iModel == 9) cFit->Print("Results/" + str_subfolder + "PhotoCrossSec/ExpFits/" + name + ".pdf");
    else cFit->Print("Results/" + str_subfolder + "PhotoCrossSec/ExpFits/" + binning + name + ".pdf");

    return;
}