// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include <chrono> // sleep_for, sleep_until
#include <thread> // nanoseconds, system_clock, seconds
// https://stackoverflow.com/questions/158585/how-do-you-add-a-timed-delay-to-a-c-program 
#include <vector>
// root headers
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
// roofit headers
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooRealSumPdf.h"
// ???
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooTFnBinding.h"

using namespace RooFit;