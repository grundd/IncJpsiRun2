// _STARlight_Utilities.h
// David Grund, May 24, 2022

// cpp headers
#include <iostream>
#include <fstream> 
#include <sstream> 
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
// needed by STARlight macros
#include "TLorentzVector.h"
#include "TClonesArray.h"

Double_t fPtGm, fPtVM, fPtPm;
TLorentzVector *parent;
TClonesArray *daughters;

//###############################################################################
// To connect branch addresses of STARlight and GammaVMPom trees:

void ConnectTreeVariables_tPtGammaVMPom(TTree *tPtGammaVMPom)
{
    tPtGammaVMPom->SetBranchAddress("fPtGm", &fPtGm);
    tPtGammaVMPom->SetBranchAddress("fPtVM", &fPtVM);
    tPtGammaVMPom->SetBranchAddress("fPtPm", &fPtPm);

    Printf("Variables from %s connected.", tPtGammaVMPom->GetName());

    return;
}

void ConnectTreeVariables_tSL(TTree *tSL)
{
    tSL->SetBranchAddress("parent", &parent);
    tSL->SetBranchAddress("daughters", &daughters);

    Printf("Variables from %s connected.", tSL->GetName());

    return;
}

//###############################################################################
// To create tree from the file PtGammaVMPom.txt

void PrepareTreesPtGammaVMPom(Int_t nGenEv, TString folder_in, TString folder_out)
{
	TString name_out = folder_out + "tree_tPtGammaVMPom.root";
    TFile *f_out = TFile::Open(name_out.Data(),"read");
    if(f_out){
        Printf("Tree %s already created.", name_out.Data());
        return;

    } else {  

		Printf("Tree %s will be created.", name_out.Data());

		// create the output file and tree
		f_out = new TFile(name_out.Data(), "RECREATE");
		if(!f_out){
			Printf("Could not create output file %s. Terminating...", name_out.Data());
			return;
		}

		TTree *tPtGammaVMPom = new TTree("tPtGammaVMPom", "tPtGammaVMPom");
		tPtGammaVMPom->Branch("fPtGm", &fPtGm, "fPtGm/D");
		tPtGammaVMPom->Branch("fPtVM", &fPtVM, "fPtVM/D");
		tPtGammaVMPom->Branch("fPtPm", &fPtPm, "fPtPm/D");

		Int_t nEntriesAnalysed = 0;
		Int_t nEntriesProgress = (Double_t)nGenEv / 20.;
		Int_t nPercent = 0;

		ifstream ifs;
		ifs.open((folder_in + "PtGammaVMPom.txt").Data());
		if(!ifs.fail()){
			for(Int_t i = 0; i < nGenEv; i++){
				ifs >> fPtGm >> fPtVM >> fPtPm;
				tPtGammaVMPom->Fill();

				if((i+1) % nEntriesProgress == 0){
				nPercent += 5;
				nEntriesAnalysed += nEntriesProgress;
				Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
				}
			}
			ifs.close();
		} else {
			Printf("File %s missing. Terminating.", (folder_in + "PtGammaVMPom.txt").Data());
			return;
		}

		tPtGammaVMPom->Write("",TObject::kWriteDelete);
		if(f_out){
			f_out->Close();
			delete f_out;
		}

		Printf("*****");
		Printf("Done.");
		Printf("*****");
		Printf("\n\n");

		return;
	}
}

//###############################################################################
// From STARlight:

double IDtoMass(int particleCode)
{
    double mass;
    if (particleCode == 2 || particleCode==3) {mass = 0.00051099907;} // electron
    else if (particleCode == 5 || particleCode==6) {mass = 0.105658389;} // muon
    else if (particleCode == 8 || particleCode==9)  {mass = 0.13956995;} // charged pion
    else if (particleCode == 7) {mass = 0.1345766;} // neutral pion
    else if (particleCode == 11|| particleCode==12) {mass = 0.493677;} // charged kaon
    else if (particleCode == 10 || particleCode == 16)  {mass = 0.497614;} // neutral kaon
    else if (particleCode == 14)	{mass = 0.93827231;} // proton
    else {
        std::cout << "unknown daughter particle (ID = " << particleCode << "), please modify code to accomodate" << std::endl;
        mass = -1.0;
        //exit(0); 
    } 

    return mass;
}

void ConvertStarlightAsciiToTree(Int_t nGenEv, TString folder_in, TString folder_out)
{
	TString name_out = folder_out + "tree_STARlight.root";
    TFile *f_out = TFile::Open(name_out.Data(),"read");
    if(f_out){
        Printf("STARlight tree %s already created.", name_out.Data());
        return;

    } else {   

        Printf("STARlight tree %s will be created.", name_out.Data());

		// create the output file and tree
		f_out = new TFile(name_out.Data(), "RECREATE");
		if(!f_out){
			Printf("Could not create output file %s. Terminating...", name_out.Data());
			return;
		}

		TTree*          outTree           = new TTree("starlightTree", "starlightTree");
		TLorentzVector* parentParticle    = new TLorentzVector();
		TClonesArray*   daughterParticles = new TClonesArray("TLorentzVector");
		outTree->Branch("parent",    "TLorentzVector", &parentParticle,    32000, -1);
		outTree->Branch("daughters", "TClonesArray",   &daughterParticles, 32000, -1);

		Int_t nEntriesAnalysed = 0;
		Int_t nEntriesProgress = (Double_t)nGenEv / 100.;
		Int_t nPercent = 0;
		Int_t i = 0;

		ifstream ifs;
		ifs.open((folder_in + "slight.out").Data());
		unsigned int countLines = 0;
		while (ifs.good()) {
			string       line;
			stringstream lineStream;
			
			// read EVENT
			string label;
			int    eventNmb, nmbTracks;
			// no more lines => end
			if (!getline(ifs, line))
				break;
			++countLines;
			lineStream.str(line);
			lineStream >> label >> eventNmb >> nmbTracks;
			//Printf("%s", line.data()); // DGRUND
			//cout << countLines << "\t" << eventNmb << "\t" << nmbTracks << endl; // DGRUND
			if (!(label == "EVENT:"))
				continue;
			
			// read VERTEX
			// no more lines => end
			if (!getline(ifs, line))
				break;
			++countLines;
			lineStream.str(line);
			lineStream >> label;
			//Printf("%s", line.data()); // DGRUND
			assert(label == "VERTEX:");
				
			*parentParticle = TLorentzVector(0, 0, 0, 0);
			for (int i = 0; i < nmbTracks; ++i) {
				// read tracks
				int    particleCode;
				double momentum[3];
				// no more lines => end
				if (!getline(ifs, line))
					break;
				++countLines;
				lineStream.str(line);
				lineStream >> label >> particleCode >> momentum[0] >> momentum[1] >> momentum[2];
				//Printf("%s", line.data()); // DGRUND
				assert(label == "TRACK:");
				Double_t daughterMass = IDtoMass(particleCode);
				if (daughterMass < 0) {break;}
				const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
									+ momentum[2] * momentum[2] + daughterMass * daughterMass);
				new ( (*daughterParticles)[i] ) TLorentzVector(momentum[0], momentum[1], momentum[2], E);
				*parentParticle += *(static_cast<TLorentzVector*>(daughterParticles->At(i)));
			}
			daughterParticles->Compress();
			outTree->Fill();

			if((i+1) % nEntriesProgress == 0){
			nPercent += 1;
			nEntriesAnalysed += nEntriesProgress;
			Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
			}
			i++;
		}

		outTree->Write("",TObject::kWriteDelete);
		if(f_out){
			f_out->Close();
			delete f_out;
		}

		Printf("*****");
		Printf("Done.");
		Printf("*****");
		Printf("\n\n");

		return;
	} 
}
//###############################################################################