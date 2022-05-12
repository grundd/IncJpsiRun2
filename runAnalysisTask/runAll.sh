# shell script must be first allowed: chmod +x runAll.sh
#!/bin/bash

# run over LHC data, kCohJpsiToEl and kCohJpsiToMu
for (( i = 0; i <= 2; i++ ))
do
    root -q runAnalysis.C\($i\)
done

# run kCohPsi2sToElPi and kCohPsi2sToMuPi (charged and neutral)
for (( i = 3; i <= 4; i++ ))
do
    root -q runAnalysis.C\($i,0\)
    root -q runAnalysis.C\($i,1\)
done

# run over kIncohJpsiToEl and kIncohJpsiToMu
for (( i = 5; i <= 6; i++ ))
do
    root -q runAnalysis.C\($i\)
done

# run kCohPsi2sToElPi and kIncohPsi2sToMuPi (charged and neutral)
for (( i = 7; i <= 8; i++ ))
do
    root -q runAnalysis.C\($i,0\)
    root -q runAnalysis.C\($i,1\)
done

# run over kTwoGammaToElMedium and kTwoGammaToMuMedium
for (( i = 9; i <= 10; i++ ))
do
    root -q runAnalysis.C\($i\)
done
