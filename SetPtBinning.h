// SetPtBinning.h
// David Grund, Feb 27, 2022

Double_t *ptBoundaries = NULL;
Double_t PtBins_4bins[5] = { 0 };
Double_t PtBins_5bins[6] = { 0 };

void SetPtBinning_main(){

    Double_t YieldPerBin_val[5] = { 0 };
    Double_t YieldPerBin_err[5] = { 0 };    

    if(nPtBins == 4) ptBoundaries = PtBins_4bins;
    if(nPtBins == 5) ptBoundaries = PtBins_5bins;

    ifstream ifs;
    ifs.open(Form("Results/" + str_subfolder + "BinsThroughMassFit/%ibins_defined.txt", nPtBins));
    for(Int_t i = 0; i < nPtBins; i++){
        // load the boundaries and yields with errors
        Double_t bin_boundary_low;
        Double_t bin_yield_val;
        Double_t bin_yield_err;
        ifs >> bin_boundary_low >> bin_yield_val >> bin_yield_err;
        ptBoundaries[i] = bin_boundary_low;
        YieldPerBin_val[i] = bin_yield_val;
        YieldPerBin_err[i] = bin_yield_err;
    }
    // load the last boundary
    ifs >> ptBoundaries[nPtBins];
    // close the file
    ifs.close();
    // print the loaded values
    Printf("*** Bin boundaries loaded: ***");
    for(Int_t i = 0; i < nPtBins; i++){
        Printf("Bin %i: pT_low = %.3f, pT_upp = %.3f, yield_val = %.3f, yield_err = %.3f",
            i+1, ptBoundaries[i], ptBoundaries[i+1], YieldPerBin_val[i], YieldPerBin_err[i]);
    }

    return;
}