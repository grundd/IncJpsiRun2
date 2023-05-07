# shell script must be first allowed: chmod +x RunAnalysis.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./RunAnalysis.sh

# define the type of the analysis (see AnalysisConfig.h)
declare -i iAnalysis=3
# define if compile each macro
declare -i compile=0
# define which macros to run
declare -a arr=("0" "1" "2" "3" "4" "5" "6y" "7" "8y" "9y" "10")
#declare -a arr=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#declare -a arr=("0" "1y" "2y" "3y" "4y" "5y" "6y" "7y" "8y" "9y" "10y")

# 1) count events (data & MC) and do run list check
if [ "${arr[1]}" = "1y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q CountEvents.C\($iAnalysis\)
        root -q CountEvents_MC.C\($iAnalysis\)
        root -q RunListCheck.C\($iAnalysis\)
    else 
        root -q CountEvents.C+\($iAnalysis\)
        root -q CountEvents_MC.C+\($iAnalysis\)
        root -q RunListCheck.C+\($iAnalysis\)
    fi
fi

# 2) Integrated luminosity:
#    - get trigger counters for both periods
#    - calculate the lumi
if [ "${arr[2]}" = "2y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q GetTriggerCounters.C\($iAnalysis\)
        root -q IntegratedLuminosity.C\($iAnalysis\)
    else 
        root -q GetTriggerCounters.C+\($iAnalysis\)
        root -q IntegratedLuminosity.C+\($iAnalysis\)
    fi
fi

# 3) Invariant mass fits of coh, inc, all and allbins
#    - MC
#    - data
if [ "${arr[3]}" = "3y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q InvMassFit_MC.C\($iAnalysis,0\)
        root -q InvMassFit.C\($iAnalysis,0\)
    else 
        root -q InvMassFit_MC.C+\($iAnalysis,0\)
        root -q InvMassFit.C+\($iAnalysis,0\)
    fi
fi

# 4) Set pT binning via the invariant mass fitting
if [ "${arr[4]}" = "4y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q BinsThroughMassFit.C\($iAnalysis\)
    else root -q BinsThroughMassFit.C+\($iAnalysis\)
    fi
fi

# 5) Invariant mass fits in pT bins:
#    - MC 
#    - data
if [ "${arr[5]}" = "5y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q InvMassFit_MC.C\($iAnalysis,1\)
        root -q InvMassFit.C\($iAnalysis,1\)
    else 
        root -q InvMassFit_MC.C+\($iAnalysis,1\)
        root -q InvMassFit.C+\($iAnalysis,1\)
    fi
fi

# 6) AxE in pT bins and veto efficiency:
if [ "${arr[6]}" = "6y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q AxE_Dissociative.cxx\($iAnalysis\)
        root -q AxE_PtBins.cxx\($iAnalysis\)
        root -q AxE_PtDep.cxx\($iAnalysis\)
        root -q VetoEfficiency.C\($iAnalysis\)
    else 
        root -q AxE_PtBins.cxx+\($iAnalysis\)
        root -q AxE_Dissociative.cxx+\($iAnalysis\)
        root -q AxE_PtDep.cxx+\($iAnalysis\)
        root -q VetoEfficiency.C+\($iAnalysis\)
    fi
fi

# 7) Fits of the transverse momentum distribution:
#    - create pT bins for the fit of the pT distribution and subtract background in these bins
#    - prepare MC templates (PDFs) to be used in pT fits
#    - calculate normalizations of feed-down curves in pT fits
#    - run all pT fits
#    - find the optimal RA (for which there is a minimum in chi2 of the pT fit)
if [ "${arr[7]}" = "7y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q PtFit_SubtractBkg.C\($iAnalysis\)
        root -q PtFit_PrepareMCTemplates.C\($iAnalysis\)
        root -q PtFit_FeedDownNormalization.C\($iAnalysis\)
        root -q PtFit_NoBkg.C\($iAnalysis\)
        root -q STARlight_OptimalRA.C\($iAnalysis\)
    else 
        root -q PtFit_SubtractBkg.C+\($iAnalysis\)
        root -q PtFit_PrepareMCTemplates.C+\($iAnalysis\)
        root -q PtFit_FeedDownNormalization.C+\($iAnalysis\)
        root -q PtFit_NoBkg.C+\($iAnalysis\)
        root -q STARlight_OptimalRA.C+\($iAnalysis\)
    fi
fi

# 8) Systematic uncertainties:
#    - signal extraction
#    - fD and fC (via changes in R)
#    - Z_vertex selection
if [ "${arr[8]}" = "8y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q InvMassFit_SystUncertainties.C\($iAnalysis\)
        root -q VertexZ_SystUncertainties.C\($iAnalysis\)
        root -q PtFit_SystUncertainties.C\($iAnalysis\)
    else 
        root -q InvMassFit_SystUncertainties.C+\($iAnalysis\)
        root -q VertexZ_SystUncertainties.C+\($iAnalysis\)
        root -q PtFit_SystUncertainties.C+\($iAnalysis\)
    fi
fi

# 9) Photonuclear cross section:
#    - calculate the average |t| per bin based on STARlight predictions
#    - calculate the photonuclear cross section
#    - prepare the histograms and graphs for the data and all the models
#    - plot the photonuclear cross section: measurement and my plot
#    - plot the photonuclear cross section: paper plot with ratios
#    - plot the fiducial cross section (integrated over |t| within 0.04 and 1.00)
#    - do exponential fits of the data and all the models
if [ "${arr[9]}" = "9y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q STARlight_tVsPt.C\($iAnalysis\)
        root -q CrossSec_Calculate.C\($iAnalysis\)
        root -q CrossSec_PrepareHistosAndGraphs.C\($iAnalysis\)
        root -q CrossSec_Plot.C\($iAnalysis\)
        root -q CrossSec_PlotWithRatios.C\($iAnalysis\)
        root -q CrossSec_Fiducial.C\($iAnalysis\)
        root -q CrossSec_ExpFits.cxx\($iAnalysis\)
    else 
        root -q STARlight_tVsPt.C+\($iAnalysis\)
        root -q CrossSec_Calculate.C+\($iAnalysis\)
        root -q CrossSec_PrepareHistosAndGraphs.C+\($iAnalysis\)
        root -q CrossSec_Plot.C+\($iAnalysis\)
        root -q CrossSec_PlotWithRatios.C+\($iAnalysis\)
        root -q CrossSec_Fiducial.C+\($iAnalysis\)
        root -q CrossSec_ExpFits.cxx+\($iAnalysis\)
    fi
fi

# 10) Extra macros:
#     - calculate pT resolution in bins as FWHM values
#     - migration of events between pT bins
#     - make PID plots (electron-muon separation)
#     - make plots to answers comments from Raphaelle
#     - detailed pT dependence of AxE
if [ "${arr[10]}" = "10y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then 
        root -q ResolutionPt.C\($iAnalysis\)
        root -q MigrationPtRecGen.C\($iAnalysis\)
        root -q ElectronsMuonsPID.C\($iAnalysis\)
        root -q RaphaelleComments.C\($iAnalysis\)
    else 
        root -q ResolutionPt.C+\($iAnalysis\)
        root -q MigrationPtRecGen.C+\($iAnalysis\)
        root -q ElectronsMuonsPID.C+\($iAnalysis\)
        root -q RaphaelleComments.C+\($iAnalysis\)
    fi
fi