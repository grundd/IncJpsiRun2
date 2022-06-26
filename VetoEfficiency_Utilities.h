// VetoEfficiency_Utilities.h
// David Grund, June 26, 2022

// cpp headers
#include <stdio.h> // printf
#include <iostream> // cout, cin
// root headers
#include "TString.h"

// neutron bins:
const Int_t nBinsN = 5; 
Double_t fNumberOfN[nBinsN+1] = {0.0, 1.5, 5.5, 10.5, 20.5, 50.5};
TString  sNumberOfN[nBinsN] = {"0-1", "2-5", "6-10", "11-20", "21-50"};
Double_t fVetoIneff_A[nBinsN+1] = {0.0, 0.085, 0.157, 0.265, 0.413, 0.579};
Double_t fVetoIneff_C[nBinsN+1] = {0.0, 0.146, 0.333, 0.425, 0.542, 0.844};
Double_t fVetoEff_A[nBinsN+1] = { 0 };
Double_t fVetoEff_C[nBinsN+1] = { 0 };

class NeutronMatrix
{
    // first index (rows) = neutron bin in A
    // second index (columns) = neutron bin in C
    // first bin = no neutrons, then five neutron bins
    public:
        NeutronMatrix();
        ~NeutronMatrix();
        void     AddEvent(Int_t iBinA, Int_t iBinC) {fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] + 1;}
        Double_t GetBinContent(Int_t iBinA, Int_t iBinC) {return fEv_neutron_bins[iBinA][iBinC];}
        Double_t CountEvents_tot();
        Double_t CountEvents_0n0n() {return fEv_neutron_bins[0][0];}
        Double_t CountEvents_Xn0n();
        Double_t CountEvents_0nXn();
        Double_t CountEvents_XnXn();
        void     Multiply(Double_t x);
        void     SubtractMatrix(NeutronMatrix *nm);
        void     ApplyEfficiencies();
        void     PrintToConsole();
        void     PrintToFile(TString name, Int_t precision = 0);
        void     LoadFromFile(TString name);
    private:
        Double_t fEv_neutron_bins[nBinsN+1][nBinsN+1];
};

NeutronMatrix::NeutronMatrix()
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) fEv_neutron_bins[iBinA][iBinC] = 0;
    }
    Printf("Neutron matrix created.");
}

NeutronMatrix::~NeutronMatrix(){}

Double_t NeutronMatrix::CountEvents_tot()
{
    Double_t sum = 0;
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) sum += fEv_neutron_bins[iBinA][iBinC];
    }
    return sum;
}

Double_t NeutronMatrix::CountEvents_Xn0n()
{
    Double_t sum = 0;
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++) sum += fEv_neutron_bins[iBinA][0];
    return sum;
}

Double_t NeutronMatrix::CountEvents_0nXn()
{
    Double_t sum = 0;
    for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) sum += fEv_neutron_bins[0][iBinC];
    return sum;
}

Double_t NeutronMatrix::CountEvents_XnXn()
{
    Double_t sum = 0;
    for(Int_t iBinA = 1; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 1; iBinC < nBinsN+1; iBinC++) sum += fEv_neutron_bins[iBinA][iBinC];
    }
    return sum;
}

void NeutronMatrix::Multiply(Double_t x)
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] * x;
    }
    return;
}

void NeutronMatrix::SubtractMatrix(NeutronMatrix *nm)
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) 
            fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] - nm->GetBinContent(iBinA,iBinC);
    }
    return;
}

void NeutronMatrix::ApplyEfficiencies()
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) 
            fEv_neutron_bins[iBinA][iBinC] = fEv_neutron_bins[iBinA][iBinC] / fVetoEff_A[iBinA] / fVetoEff_C[iBinC];
    }
    return;
}

void NeutronMatrix::PrintToConsole()
{
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) cout << fEv_neutron_bins[iBinA][iBinC] << "\t";
        cout << endl;
    }
    return;
}

void NeutronMatrix::PrintToFile(TString name, Int_t precision = 0)
{
    ofstream outfile;
    outfile.open(name.Data());
    outfile << std::fixed << std::setprecision(precision);
    for(Int_t iBinA = 0; iBinA < nBinsN+1; iBinA++){
        for(Int_t iBinC = 0; iBinC < nBinsN+1; iBinC++) outfile << fEv_neutron_bins[iBinA][iBinC] << "\t";
        outfile << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", name.Data());
    return;
}

void NeutronMatrix::LoadFromFile(TString name)
{
    ifstream ifs;
    ifs.open(name.Data());
    if(!ifs.fail())
    {
        // Read data from the file
        for(Int_t iRow = 0; iRow < nBinsN+1; iRow++){
            for(Int_t iCol = 0; iCol < nBinsN+1; iCol++){
                ifs >> fEv_neutron_bins[iRow][iCol];
            }
        }
        Printf("Neutron matrix loaded from %s.", name.Data());
    } else {
        Printf("Cannot open the file %s.", name.Data());
    }
    ifs.close();
    return;
}