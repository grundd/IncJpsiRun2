// VetoEff.C
// David Grund, June 26, 2022

#include "VetoEfficiency_Utilities.h"

void VetoEff()
{
    NeutronMatrix *nm = new NeutronMatrix();
    nm->LoadFromFile("Results/5bins_pass3/VetoEfficiency/mass_3.00to3.20/nEv_PtAll.txt");
    nm->PrintToConsole();

    return;
}