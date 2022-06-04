# shell script must be first allowed: chmod +x RunAnalysis.sh
#!/bin/bash
# to run it do (inside Ali shell):
# ./RunAnalysis.sh

# define the type of the analysis (see AnalysisConfig.h)
declare -i iAnalysis=14
# define if compile each macro
declare -i compile=0
# define which macros to run
#declare -a arr=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25")
declare -a arr=("0y" "1y" "2y" "3y" "4y" "5y" "6y" "7y" "8y" "9y" "10y" "11" "12y" "13y" "14y" "15y" "16y" "17y" "18y" "19y" "20y" "21y" "22y" "23y" "24y" "25y")

# 0) Count events (data)
if [ "${arr[0]}" = "0y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q CountEvents.C\($iAnalysis\)
    else root -q CountEvents.C+\($iAnalysis\)
    fi
fi

# 1) Count events (MC)
if [ "${arr[1]}" = "1y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q CountEvents_MC.C\($iAnalysis\)
    else root -q CountEvents_MC.C+\($iAnalysis\)
    fi
fi

# 2) Run list check
if [ "${arr[2]}" = "2y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q RunListCheck.C\($iAnalysis\)
    else root -q RunListCheck.C+\($iAnalysis\)
    fi
fi

# 3) Get trigger counters for both periods
if [ "${arr[3]}" = "3y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q GetTriggerCounters.C\($iAnalysis\)
    else root -q GetTriggerCounters.C+\($iAnalysis\)
    fi
fi

# 4) Calculate the integrated luminosity
if [ "${arr[4]}" = "4y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q IntegratedLuminosity.C\($iAnalysis\)
    else root -q IntegratedLuminosity.C+\($iAnalysis\)
    fi
fi

# 5) MC invariant mass fits of coh, inc, all and allbins
if [ "${arr[5]}" = "5y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q InvMassFit_MC.C\($iAnalysis,0\)
    else root -q InvMassFit_MC.C+\($iAnalysis,0\)
    fi
fi

# 6) Invariant mass fits of coh, inc, all and allbins
if [ "${arr[6]}" = "6y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q InvMassFit.C\($iAnalysis,0\)
    else root -q InvMassFit.C+\($iAnalysis,0\)
    fi
fi

# 7) Set pT binning via the invariant mass fitting
if [ "${arr[7]}" = "7y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q BinsThroughMassFit.C\($iAnalysis\)
    else root -q BinsThroughMassFit.C+\($iAnalysis\)
    fi
fi

# 8) MC invariant mass fits in pT bins
if [ "${arr[8]}" = "8y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q InvMassFit_MC.C\($iAnalysis,1\)
    else root -q InvMassFit_MC.C+\($iAnalysis,1\)
    fi
fi

# 9) Invariant mass fits in pT bins
if [ "${arr[9]}" = "9y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q InvMassFit.C\($iAnalysis,1\)
    else root -q InvMassFit.C+\($iAnalysis,1\)
    fi
fi

# 10) AxE in pT bins
if [ "${arr[10]}" = "10y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q AxE_PtBins.C\($iAnalysis\)
    else root -q AxE_PtBins.C+\($iAnalysis\)
    fi
fi

# 11) AxE: dependence on pT
if [ "${arr[11]}" = "11y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q AxE_PtDep.C\($iAnalysis\)
    else root -q AxE_PtDep.C+\($iAnalysis\)
    fi
fi

# 12) Create pT bins for the fit of the pT distribution and subtract background in these bins
if [ "${arr[12]}" = "12y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PtFit_SubtractBkg.C\($iAnalysis\)
    else root -q PtFit_SubtractBkg.C+\($iAnalysis\)
    fi
fi

# 13) Prepare MC templates (PDFs) to be used in pT fits
if [ "${arr[13]}" = "13y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PtFit_PrepareMCTemplates.C\($iAnalysis\)
    else root -q PtFit_PrepareMCTemplates.C+\($iAnalysis\)
    fi
fi

# 14) Various pT fits with background subtracted
if [ "${arr[14]}" = "14y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PtFit_NoBkg.C\($iAnalysis\)
    else root -q PtFit_NoBkg.C+\($iAnalysis\)
    fi
fi

# 15) Find the optimal RA (for which there is a minimum in chi2 of the pT fit)
if [ "${arr[15]}" = "15y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q STARlight_OptimalRA.C\($iAnalysis\)
    else root -q STARlight_OptimalRA.C+\($iAnalysis\)
    fi
fi

# 16) Systematic uncertainties corresponding to signal extraction (invariant mass fitting)
if [ "${arr[16]}" = "16y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q InvMassFit_SystUncertainties.C\($iAnalysis\)
    else root -q InvMassFit_SystUncertainties.C+\($iAnalysis\)
    fi
fi

# 17) Calculate the average |t| per bin based on STARlight predictions
if [ "${arr[17]}" = "17y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q STARlight_tVsPt.C\($iAnalysis\)
    else root -q STARlight_tVsPt.C+\($iAnalysis\)
    fi
fi

# 18) Calculate the photonuclear cross section
if [ "${arr[18]}" = "18y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PhotoCrossSec_Calculate.C\($iAnalysis\)
    else root -q PhotoCrossSec_Calculate.C+\($iAnalysis\)
    fi
fi

# 19) Plot the photonuclear cross section
if [ "${arr[19]}" = "19y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PhotoCrossSec_Plot.C\($iAnalysis\)
    else root -q PhotoCrossSec_Plot.C+\($iAnalysis\)
    fi
fi

# 20) Plot the photonuclear cross section with ratios
if [ "${arr[20]}" = "20y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PhotoCrossSec_PlotWithRatios.C\($iAnalysis\)
    else root -q PhotoCrossSec_PlotWithRatios.C+\($iAnalysis\)
    fi
fi

# 21) Plot the total cross section (integrated over |t|)
if [ "${arr[21]}" = "21y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q PhotoCrossSec_Total.C\($iAnalysis\)
    else root -q PhotoCrossSec_Total.C+\($iAnalysis\)
    fi
fi

# 22) Migration of events between pT bins
if [ "${arr[22]}" = "22y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q MigrationPtReecGen.C\($iAnalysis\)
    else root -q MigrationPtReecGen.C+\($iAnalysis\)
    fi
fi

# 23) Calculate pT resolution as a FWHM value
if [ "${arr[23]}" = "23y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q ResolutionPt.C\($iAnalysis\)
    else root -q ResolutionPt.C+\($iAnalysis\)
    fi
fi

# 24) Make PID plots (muon-electron separation)
if [ "${arr[24]}" = "24y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q ElectronsMuonsPID.C\($iAnalysis\)
    else root -q ElectronsMuonsPID.C+\($iAnalysis\)
    fi
fi

# 25) Answers to comments from Raphaelle
if [ "${arr[25]}" = "25y" ] 
then
    if [[ "$compile" -eq 0 ]]
    then root -q RaphaelleComments.C\($iAnalysis\)
    else root -q RaphaelleComments.C+\($iAnalysis\)
    fi
fi