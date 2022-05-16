// SetPtBinning_PtFit.h
// David Grund, Mar 20, 2022

#include <fstream>
#include <sstream> 
#include "RooBinning.h"

vector<Double_t> BinsBoundaries_PtFit;
Double_t fPtCutLow_PtFit = 0.0;
Double_t fPtCutUpp_PtFit = 2.0;
RooBinning fPtBins_PtFit(fPtCutLow_PtFit, fPtCutUpp_PtFit);

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
            BinsBoundaries_PtFit.push_back(binBoundary);
            i++;   
        }
    }  
    ifs.close();
    // set the number of bins
    nPtBins_PtFit = BinsBoundaries_PtFit.size() - 1;
    // set the values of boundaries
    ptBoundaries_PtFit = &BinsBoundaries_PtFit[0];
    // set RooBinning
    for(Int_t i = 0; i < nPtBins_PtFit - 1; i++) fPtBins_PtFit.addBoundary(ptBoundaries_PtFit[i + 1]);

    // print the loaded values
    Printf("*** Bin boundaries loaded: ***");
    for(Int_t i = 0; i < nPtBins_PtFit; i++){
        Printf("Bin %i: pT_low = %.3f, pT_upp = %.3f",i+1, ptBoundaries_PtFit[i], ptBoundaries_PtFit[i+1]);
    }

    return;
}