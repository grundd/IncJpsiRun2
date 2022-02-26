// SetPtBinning.h
// David Grund, Feb 26, 2022

Double_t *ptBoundaries = NULL;
Double_t pt4bins_p1[5] = {0.200, 0.280, 0.375, 0.569, 1.000};
Double_t pt5bins_p1[6] = {0.200, 0.264, 0.335, 0.445, 0.658, 1.000};
Double_t pt4bins_p3[5] = {0.200, 0.283, 0.390, 0.572, 1.000};
Double_t pt5bins_p3[6] = {0.200, 0.265, 0.336, 0.453, 0.659, 1.000};

void SetPtBinning(Bool_t isPass3){

    if(!pass3){
        // PtBinning "Method 3": BinsThroughMassFit.c
        if(nPtBins == 4) ptBoundaries = pt4bins_p1;
        if(nPtBins == 5) ptBoundaries = pt5bins_p1;
    } else {
        if(nPtBins == 4) ptBoundaries = pt4bins_p3;
        if(nPtBins == 5) ptBoundaries = pt5bins_p3;            
    }
}