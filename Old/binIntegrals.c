//-------------------------------------------------------
// macro to test other normalisations of models to data
// for the cross section of incoherent j/psi prod. at y=0
//-------------------------------------------------------

//-------------------------------------------------------
// headers
//-------------------------------------------------------
// c++ headers
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


// root headers
#include <TGraph.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>

//-------------------------------------------------------
// info about the models
// reserve space for data to be read
//-------------------------------------------------------

// HS predictions
const Int_t nData_HS = 76;
Double_t abs_t_HS[nData_HS];
Double_t sig_HS_coh_n[nData_HS];
Double_t sig_HS_coh_n_err[nData_HS];
Double_t sig_HS_inc_n[nData_HS];
Double_t sig_HS_inc_n_err[nData_HS];
Double_t sig_HS_coh_hs[nData_HS];
Double_t sig_HS_coh_hs_err[nData_HS];
Double_t sig_HS_inc_hs[nData_HS];
Double_t sig_HS_inc_hs_err[nData_HS];

// Guzey predictions
const Int_t nData_GZ = 100;
Double_t abs_t_GZ[nData_GZ];
Double_t sig_GZ_el_min[nData_GZ];
Double_t sig_GZ_el_max[nData_GZ];
Double_t sig_GZ_diss_min[nData_GZ];
Double_t sig_GZ_diss_max[nData_GZ];
Double_t sig_GZ_tot_min[nData_GZ];
Double_t sig_GZ_tot_max[nData_GZ];
Double_t sig_GZ_tot_avg[nData_GZ];
Double_t sig_GZ_el_avg[nData_GZ];

// Heikki predictions
const Int_t nData_HM = 183;
Double_t abs_t_HM[nData_HM];
Double_t sig_HM_fluct[nData_HM];
Double_t sig_HM_noflu[nData_HM];

// STARlight predictions
const Int_t nData_SL = 125;
Double_t abs_t_SL[nData_SL];
Double_t sig_SL[nData_SL];

// bins in t
const int nBins = 5;
double min_t[nBins] = {0.0400,0.0801,0.1521,0.2581,0.4775};
double max_t[nBins] = {0.0801,0.1521,0.2581,0.4775,1.0};

//-------------------------------------------------------
// functions to read in the models
//-------------------------------------------------------

//-------------------------------------------------------
void ReadInput_HSModel()
{
    // read the input file for hot-spot model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/HSModel/data-dtdy-y_0.6-Run1.txt");
    for(Int_t i = 0; i < nData_HS; i++){
        Double_t x;
        ifs >> x;
        Double_t tmp;
        ifs >> tmp;
        ifs >> tmp;
        ifs >> abs_t_HS[i];
        ifs >> sig_HS_coh_n[i];
        ifs >> sig_HS_coh_n_err[i];        
        ifs >> sig_HS_inc_n[i];
        ifs >> sig_HS_inc_n_err[i];        
        ifs >> sig_HS_coh_hs[i];
        ifs >> sig_HS_coh_hs_err[i];      
        ifs >> sig_HS_inc_hs[i];
        ifs >> sig_HS_inc_hs_err[i];
        //std::cout << i << " " << abs_t[i] << " " <<sig_HS_inc_n[i]<< " " <<sig_HS_inc_hs[i] << endl;
    }
    ifs.close();
    Printf("Predictions of the HS model loaded.");

    return;
}

//-------------------------------------------------------
void ReadInput_Guzey()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Guzey/incoh_tdep_nuc_run2.dat");
    for(Int_t i = 0; i < nData_GZ; i++){
        ifs >> abs_t_GZ[i];
        ifs >> sig_GZ_el_min[i];
        ifs >> sig_GZ_el_max[i];
        ifs >> sig_GZ_diss_min[i];
        ifs >> sig_GZ_diss_max[i];
        ifs >> sig_GZ_tot_min[i];
        ifs >> sig_GZ_tot_max[i];
        //std::cout << i << " " << abs_t_GZ[i] << " " << sig_GZ_diss_min[i]<< " " << sig_GZ_diss_max[i] << endl;
	sig_GZ_tot_avg[i] = 0.5*(sig_GZ_tot_min[i]+sig_GZ_tot_max[i])*1e-6;
	sig_GZ_el_avg[i] = 0.5*(sig_GZ_el_min[i]+sig_GZ_el_max[i])*1e-6;	
    }
    ifs.close();
    Printf("Predictions of Guzey's model loaded.");

    return;
}

//-------------------------------------------------------
void ReadInput_Heikki()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_fluct");
    for(Int_t i = 0; i < nData_HM; i++){
        ifs >> abs_t_HM[i];
        ifs >> sig_HM_fluct[i];
        //std::cout << i << " " << abs_t_HM[i] << " " << sig_HM_fluct[i] << endl;
    }
    ifs.close();
    ifs.open("Trees/PhotoCrossSec/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_nofluct");
    for(Int_t i = 0; i < nData_HM; i++){
        ifs >> abs_t_HM[i];
        ifs >> sig_HM_noflu[i];
        //std::cout << i << " " << abs_t_HM[i] << " " << sig_HM_noflu[i] << endl;
    }
    ifs.close();
    Printf("Predictions of Heikki's model loaded.");

    return;
}

//-------------------------------------------------------
void ReadInput_STARlight()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("Trees/PhotoCrossSec/STARlight/IncJ_tDep_0.00-2.50.txt");
    for(Int_t i = 0; i < nData_SL; i++){
        ifs >> abs_t_SL[i];
        ifs >> sig_SL[i];
        //std::cout << i << " " << abs_t_SL[i] << " " << sig_SL[i] << endl;
    }
    ifs.close();
    Printf("STARlight predictions loaded.");

    return;
}

//-------------------------------------------------------
// functions to integrate the models
//-------------------------------------------------------

void doOneBin(double tmin, double tmax, int nData, double *abs_t, double *sig, TGraph *gr)
{

  // some variables to store the info
  double tot = 0; // total cross section in bin
  double avgt = 0; // average |t| in bin
  int iFirst = -1; // index to first data point inside bin
  int iLast = -1; // index to last data point inside bin  

  // get the average cross section and average t
  cout << " +++++++++++++++++++++++++++++++++++ " << endl;
  for(int i=0;i<nData;i++) {
    // select range
    if (abs_t[i] < tmin) continue;
    if (abs_t[i] > tmin && iFirst == -1) iFirst = i;
    if (abs_t[i] < tmax) iLast = i;
    if (abs_t[i] > tmax) continue;
    // compute width of the bin (half distances to previous/next bin)
    double Dt1 = 0.5*(abs_t[i+1]-abs_t[i]);
    double Dt2 = 0.5*(abs_t[i]-abs_t[i-1]);
    //cout << "i-1:" << abs_t[i-1] << " i:" << abs_t[i] << " i+1:" << abs_t[i+1] << endl;
    double Dt = Dt1+Dt2;
    cout << "integrating from " << abs_t[i]-Dt2 << " to " << abs_t[i]+Dt1 << endl;
    // update total cross section
    tot += (Dt*sig[i]);
    // update average <t>
    avgt += (Dt*sig[i]*abs_t[i]);
    //    cout << abs_t[i] << " " << sig[i] << " " << Dt << endl;
  }
  // compute borders of integral
  double Dtmin = 0.5*(abs_t[iFirst]-abs_t[iFirst-1]);
  double tminDt = abs_t[iFirst]-Dtmin;
  double Dtmax = 0.5*(abs_t[iLast+1]-abs_t[iLast]);
  double tmaxDt = abs_t[iLast]+Dtmax;
  // correct for borther effects at tmin
  double dt1 = 0.5*(abs_t[iFirst+1]-abs_t[iFirst-1]);
  double dt2 = dt1+(tminDt-tmin);
  double sig1 = dt1*sig[iFirst];
  double sig2 = dt2*sig[iFirst];
  tot += (sig2-sig1);
  avgt += (abs_t[iFirst]*(sig2-sig1));
  // correct for borther effects at tmax
  double dt3 = 0.5*(abs_t[iLast+1]-abs_t[iLast-1]);
  double dt4 = dt3-(tmaxDt-tmax);
  double sig3 = dt3*sig[iLast];
  double sig4 = dt4*sig[iLast];
  tot += (sig4-sig3);
  avgt += (abs_t[iLast]*(sig4-sig3));
  // print borders
  cout << " Limits: (" << tmin << "," << tmax << ") ... first and last = ("<< tminDt <<","<<tmaxDt<<")" << endl;
  cout << " Integral:" << (tot * 1e3) << " mub,  <t>" << (avgt/tot) << endl;
  // store point
  gr->SetPoint(gr->GetN(),(avgt/tot),(tot/(tmax-tmin)));
}

//-------------------------------------------------------
void doOneModel(int nData, double *abs_t, double *sig, TGraph *gr)
{
  // compute cross section in each bin and put it in the graph

  // loop over bins
  for(int iBin=0;iBin<nBins;iBin++) {
    doOneBin(min_t[iBin],max_t[iBin],nData,abs_t,sig,gr);
  }

}

/*
void doOneBin(double tmin, double tmax, int nData, double *abs_t, double *sig, TGraph *gr)
{
  // get the average cross section and average t
  double n = 0;
  double tot = 0;
  double avgt = 0;
  for(int i=0;i<nData;i++) {
    if (abs_t[i] < tmin) continue;
    if (abs_t[i] > tmax) continue;
    tot += sig[i];
    avgt += (sig[i]*abs_t[i]);
    n++;
  }
  // cout << " cs = " << (tot/n) << " <t> = " << (avgt/tot) << endl;
  gr->SetPoint(gr->GetN(),(avgt/tot),(tot/n));
}
*/

//-------------------------------------------------------
// K-S test auxiliary functions
//-------------------------------------------------------

void getCum(int n, TGraph *g, double *dt, double *cum)
{
  double tot = 0;
  for(int i=0;i<n;i++) {
    cum[i]=dt[i]*(g->GetPointY(i));
    tot+=cum[i];
  }
  for(int i=0;i<n;i++) {
    cum[i]/=tot;
    cout << cum[i] << " " << tot << endl;
  }
}

//-------------------------------------------------------
double doKS(int n, double *d, double *t)
{
  double max = 0;
  for(int i=0;i<n;i++) {
    double m = abs(d[i]-t[i]);
    if (m>max) max = m;
  }
  return max;
}

//-------------------------------------------------------
// main entry point
//-------------------------------------------------------

void binIntegrals()
{

  // read models
  ReadInput_STARlight();
  ReadInput_HSModel();
  ReadInput_Heikki();
  ReadInput_Guzey();


  // define graphs for models
  TGraph *gr_inc_hs = new TGraph();
  TGraph *gr_inc_n = new TGraph();
  TGraph *gr_hm_nof = new TGraph();
  TGraph *gr_hm_f = new TGraph();
  TGraph *gr_gs_el = new TGraph();
  TGraph *gr_gs_tot = new TGraph();
  TGraph *gr_sl = new TGraph();

  cout << " HS " << endl;
  doOneBin(0.04, 1.0, nData_HS,abs_t_HS,sig_HS_inc_hs,gr_inc_hs);
  /*
  cout << " n " << endl;
  doOneBin(0.04, 1.0, nData_HS,abs_t_HS,sig_HS_inc_n,gr_inc_n);	   
  cout << " GZ_el " << endl;
  doOneBin(0.04, 1.0, nData_GZ,abs_t_GZ,sig_GZ_el_avg,gr_gs_el);
  cout << " GZ_tot " << endl;
  doOneBin(0.04, 1.0, nData_GZ,abs_t_GZ,sig_GZ_tot_avg,gr_gs_tot);	   
  cout << " HM_no " << endl;
  doOneBin(0.04, 1.0, nData_HM,abs_t_HM,sig_HM_noflu,gr_hm_nof);	   
  cout << " HM_f " << endl;
  doOneBin(0.04, 1.0, nData_HM,abs_t_HM,sig_HM_fluct,gr_hm_f);
  cout << " SL " << endl;
  doOneBin(0.04, 1.0, nData_SL,abs_t_SL,sig_SL,gr_sl);
  */

  return;
  
  doOneModel(nData_HS,abs_t_HS,sig_HS_inc_hs,gr_inc_hs);
  doOneModel(nData_HS,abs_t_HS,sig_HS_inc_n,gr_inc_n);
  doOneModel(nData_GZ,abs_t_GZ,sig_GZ_el_avg,gr_gs_el);
  doOneModel(nData_GZ,abs_t_GZ,sig_GZ_tot_avg,gr_gs_tot);
  doOneModel(nData_HM,abs_t_HM,sig_HM_noflu,gr_hm_nof);
  doOneModel(nData_HM,abs_t_HM,sig_HM_fluct,gr_hm_f);
  doOneModel(nData_SL,abs_t_SL,sig_SL,gr_sl);
  
  // set cosmetics
  gr_inc_hs->SetMarkerStyle(20);gr_inc_hs->SetMarkerColor(kRed);gr_inc_hs->SetMarkerSize(1.25);
  gr_inc_n->SetMarkerStyle(23);gr_inc_n->SetMarkerColor(kRed);gr_inc_n->SetMarkerSize(1.25);
  gr_gs_tot->SetMarkerStyle(20);gr_gs_tot->SetMarkerColor(kGreen);gr_gs_tot->SetMarkerSize(1.25);
  gr_gs_el->SetMarkerStyle(23);gr_gs_el->SetMarkerColor(kGreen);  gr_gs_el->SetMarkerSize(1.25);
  gr_hm_f->SetMarkerStyle(20);gr_hm_f->SetMarkerColor(kCyan);gr_hm_f->SetMarkerSize(1.25);
  gr_hm_nof->SetMarkerStyle(23);gr_hm_nof->SetMarkerColor(kCyan);gr_hm_nof->SetMarkerSize(1.25);
  gr_sl->SetMarkerStyle(20);gr_sl->SetMarkerColor(kBlue);gr_sl->SetMarkerSize(1.25);

  ///////////////////////////////
  // now data
  ///////////////////////////////
  const int nData = 5;
  double incCS[nData][4] = {
    {21.8,2.4,2.1},
    {19.1,1.9,1.5},
    {13.1,1.6,1.0},
    {8.1,1.1,0.5},
    {4.4,0.6,0.3},
  };
  double inct[nData] = {0.06,0.114,0.201,0.351,0.652};
  double Deltat[nData] = {0.04,0.072,0.106,0.219,0.523};
  TGraphErrors *incGE = new  TGraphErrors();
  for(int i=0;i<nData;i++) {
    double err = std::sqrt(incCS[i][1]*incCS[i][1]+incCS[i][2]*incCS[i][2]);
    incGE->SetPoint(i,inct[i],incCS[i][0]/1000.);
    incGE->SetPointError(i,0,err/1000);    
  }

  ///////////////////////////////
  // Kolmogorov–Smirnov test
  ///////////////////////////////

  // get cumulants
  double cumData[nData];
  getCum(nData,incGE,Deltat,cumData);
  double cumHS[nData];
  getCum(nData,gr_inc_hs,Deltat,cumHS);
  double cumN[nData];
  getCum(nData,gr_inc_n,Deltat,cumN);
  double cumHMF[nData];
  getCum(nData,gr_hm_f,Deltat,cumHMF);
  double cumHMNF[nData];
  getCum(nData,gr_hm_nof,Deltat,cumHMNF);
  double cumEL[nData];
  getCum(nData,gr_gs_el,Deltat,cumEL);
  double cumELD[nData];
  getCum(nData,gr_gs_tot,Deltat,cumELD);
  double cumSL[nData];
  getCum(nData,gr_sl,Deltat,cumSL);

  // perform test
  double ksHS = doKS(nData,cumData,cumHS);
  double ksN = doKS(nData,cumData,cumN);
  double ksHMF = doKS(nData,cumData,cumHMF);
  double ksHMNF = doKS(nData,cumData,cumHMNF);
  double ksEL = doKS(nData,cumData,cumEL);
  double ksELD = doKS(nData,cumData,cumELD);  
  double ksSL = doKS(nData,cumData,cumSL);

  // print results
  cout << endl << " ===> Results of KS test: " << endl;
  cout << "     HMF  = " << ksHMF << endl;  
  cout << "     EL+D = " << ksELD << endl;    
  cout << "     HS   = " << ksHS << endl;
  cout << "     SL   = " << ksSL << endl;    
  cout << "     HMNF = " << ksHMNF << endl;  
  cout << "     EL   = " << ksEL << endl;
  cout << "     N    = " << ksN << endl;

  return;
  
  ///////////////////////////////
  // normalise to first bin
  ///////////////////////////////
  double sc_hs = incGE->GetPointY(0)/gr_inc_hs->GetPointY(0);
  double sc_n = incGE->GetPointY(0)/gr_inc_n->GetPointY(0);
  double sc_hmf = incGE->GetPointY(0)/gr_hm_f->GetPointY(0);
  double sc_hmnf = incGE->GetPointY(0)/gr_hm_nof->GetPointY(0);
  double sc_gs_el = incGE->GetPointY(0)/gr_gs_el->GetPointY(0);
  double sc_gs_tot = incGE->GetPointY(0)/gr_gs_tot->GetPointY(0);
  double sc_sl = incGE->GetPointY(0)/gr_sl->GetPointY(0);
  cout << incGE->GetPointY(0) << " " << gr_sl->GetPointY(0) << " " << sc_sl << endl;
  cout << incGE->GetPointY(0) << " gr_gs_tot " << gr_gs_tot->GetPointY(0) << " " << sc_gs_tot << endl;
  cout << incGE->GetPointY(0) << " gr_hm_f " << gr_hm_f->GetPointY(0) << " " << sc_hmf << endl;
  cout << incGE->GetPointY(0) << " gr_inc_hs " << gr_inc_hs->GetPointY(0) << " " << sc_hs  << endl;  

  
   for(int i=0;i<5;i++) {
     gr_inc_hs->SetPoint(i,gr_inc_hs->GetPointX(i),sc_hs*gr_inc_hs->GetPointY(i));
     gr_inc_n->SetPoint(i,gr_inc_n->GetPointX(i),sc_n*gr_inc_n->GetPointY(i));
     gr_hm_f->SetPoint(i,gr_hm_f->GetPointX(i),sc_hmf*gr_hm_f->GetPointY(i));
     gr_hm_nof->SetPoint(i,gr_hm_nof->GetPointX(i),sc_hmnf*gr_hm_nof->GetPointY(i));
     gr_gs_el->SetPoint(i,gr_gs_el->GetPointX(i),sc_gs_el*gr_gs_el->GetPointY(i));
     gr_gs_tot->SetPoint(i,gr_gs_tot->GetPointX(i),sc_gs_tot*gr_gs_tot->GetPointY(i));
     gr_sl->SetPoint(i,gr_sl->GetPointX(i),sc_sl*gr_sl->GetPointY(i));
  } 
  
  // plot data
  ////////////////
  TCanvas *c = new TCanvas("c","c",800,600);
  TH1F* f = gPad->DrawFrame(0.02,0.0005,1.02,.10);
  f->SetTitle("Cross section dependence on |#it{t}|;|#it{t}| (GeV^{2} #it{c}^{-2});d#sigma_{#gammaPb}/d|#it{t}| (mb #it{c}^{2} GeV^{-2})");
  f->Draw("AXIS");
  gPad->SetLogy();

  // add models
  gr_inc_hs->Draw("p,same");
  gr_inc_n->Draw("p,same");
  gr_gs_tot->Draw("p,same");
  gr_gs_el->Draw("p,same");
  gr_hm_f->Draw("p,same");
  gr_hm_nof->Draw("p,same");
  gr_sl->Draw("p,same");

  // add data
  incGE->SetMarkerStyle(20);
  incGE->SetMarkerSize(0.8);
  incGE->SetMarkerColor(kBlack);
  incGE->Draw("p,same");


  // add legend
  TLegend *leg1 = new TLegend(0.65,0.55,0.85,0.85);
  leg1->SetTextSize(0.025);
  leg1->AddEntry(gr_sl,Form("STARlight * %4.2f",sc_sl), "p");
  leg1->AddEntry(gr_hm_f,Form("MS: flu. * %4.2f",sc_hmf), "p");
  leg1->AddEntry(gr_hm_nof,Form("MS: no-flu. * %4.2f",sc_hmnf), "p");
  leg1->AddEntry(gr_gs_tot,Form("GSZ:el+diss * %4.2f",sc_gs_tot), "p");
  leg1->AddEntry(gr_gs_el,Form("GSZ:el * %4.2f",sc_gs_el), "p");
  leg1->AddEntry(gr_inc_n,Form("CCK: n * %4.2f",sc_n), "p");
  leg1->AddEntry(gr_inc_hs,Form("CCK: hs * %4.2f",sc_hs), "p");      
  leg1->Draw();

    TLatex* latex = new TLatex();
    latex->SetTextSize(0.035);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.93,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
  
}

