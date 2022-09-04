// CrossSec_main.C
// David Grund, Sep 3, 2022

#include "CrossSec_Utilities.h"

void CrossSec_main(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);
    InitObjects();

    // load graphs
    LoadGraphs_data();
    LoadGraphs_SL();
    LoadGraphs_CCK();
    LoadGraphs_MS();
    LoadGraphs_GSZ();

    Printf("Data integral is: %.3f", IntegrateData()*1e3);

    // create histograms from graphs
    for(Int_t i = 0; i < 7; i++)
    {
        CreateHistogramFromGraph(i);
        Double_t integral(0.), avgt(0.);
        IntegrateModel(i,tBoundaries[4],tBoundaries[5],integral,avgt);
    } 

    return;
}