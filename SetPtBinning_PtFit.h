// SetPtBinning_PtFit.h
// David Grund, Mar 20, 2022

#include <fstream>
#include <sstream> 
#include "RooBinning.h"

vector<Double_t> ptBoundaries_PtFit_vec;
vector<Double_t> tBoundaries_PtFit_vec;
Double_t fPtCutLow_PtFit = 0.0;
Double_t fPtCutUpp_PtFit = 2.0;
RooBinning fPtBins_PtFit(fPtCutLow_PtFit, fPtCutUpp_PtFit);
RooBinning fTBins_PtFit(fPtCutLow_PtFit*fPtCutLow_PtFit, fPtCutUpp_PtFit*fPtCutUpp_PtFit);

void SetPtBinning_PtFit()
{
    ifstream ifs;
    ifs.open("Results/" + str_subfolder + "PtFit_SubtractBkg/bins_defined.txt");
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            //Printf("Reading bin boundary no. %i: %s", i, str.data());
            istringstream istr(str);
            Double_t binBoundary;
            istr >> binBoundary;
            ptBoundaries_PtFit_vec.push_back(binBoundary);
            tBoundaries_PtFit_vec.push_back(binBoundary*binBoundary);
            i++;   
        }
    }  
    ifs.close();
    // set the number of bins
    nPtBins_PtFit = ptBoundaries_PtFit_vec.size() - 1;
    // set the values of boundaries
    ptBoundaries_PtFit = &ptBoundaries_PtFit_vec[0];
    tBoundaries_PtFit  = &tBoundaries_PtFit_vec[0];
    // set RooBinning
    for(Int_t i = 0; i < nPtBins_PtFit - 1; i++)
    {
        fPtBins_PtFit.addBoundary(ptBoundaries_PtFit[i + 1]);
        fTBins_PtFit.addBoundary(tBoundaries_PtFit[i + 1]);
    }

    // print the loaded values
    Printf("*** Bin boundaries loaded: ***");
    for(Int_t i = 0; i < nPtBins_PtFit; i++){
        Printf("Bin %i: pT_low = %.3f, pT_upp = %.3f",i+1, ptBoundaries_PtFit[i], ptBoundaries_PtFit[i+1]);
    }
    Printf("*** Bin boundaries loaded: ***");
    for(Int_t i = 0; i < nPtBins_PtFit; i++){
        Printf("Bin %i: t_low = %.4f, t_upp = %.4f",i+1, tBoundaries_PtFit[i], tBoundaries_PtFit[i+1]);
    }

    return;
}