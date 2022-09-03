// CrossSec_main.C
// David Grund, Sep 3, 2022

#include "CrossSec_Utilities.h"

void CrossSec_main(Int_t iAnalysis)
{
    InitAnalysis(iAnalysis);

    InitObjects();

    // load graphs
    LoadGraphs_SL();
    LoadGraphs_CCK();
    LoadGraphs_MS();
    LoadGraphs_GSZ();

    // create histograms from graphs
    for(Int_t i = 0; i < 7; i++)
    {
        CreateHistogramFromGraph(i);
        IntegrateModel(i,0.04,1.00);
    } 

    return;
}