//09/10/2017
//Program inspired by gs_project1.cc and main92.cc
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-saves some features of the final particles of each event in n-tuples;
//		-restricts decays for eta: eta->mu anti-mu
//		-tries to reconstruct invariant mass of eta from decay products by pairing up the decay particles

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

double invmass(Vec4 p1, Vec4 p2){

	double inv_mass=(p1 + p2).mCalc();
	return inv_mass;

}

int main() {

	HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	HepMC::IO_GenEvent ascii_io("gs_project2.dat", std::ios::out);

	// Generator. Process selection. LHC initialization. Histogram.
    	Pythia pythia;
    	pythia.readString("Random:setSeed = on");
    	pythia.readString("Random:seed = 0");
    	//pythia.readString("HardQCD:all = on");
    	pythia.readString("SoftQCD:all = on");
	pythia.readString("PhaseSpace:pTHatMin = 20.");

	//Set omega to decay to 2 electrons and gamma
	//pythia.readString("221:onMode = 1 1 11 22 11 -11");

	//Set omega to decay to 2 electrons
	pythia.readString("223:onMode = 0");
	pythia.readString("223:onIfMatch = 11 -11");
	//pythia.readString("333:onMode = 0");
	//pythia.readString("333:onIfMatch = 11 -11");

	//Initilise for p (2212) and  p (2212) collisins at 13 TeV
    	//initialization of LHC environment
    	pythia.readString("Beams:idA = 2212");
    	pythia.readString("Beams:idB = 2212");

    	//Start at 13 TeV
    	double myEcm =13000.;
    	pythia.settings.parm("Beams:eCM", myEcm);

	pythia.init();

	// Set up the ROOT TFile and TTree.
	TFile *file = TFile::Open("gs_project2.root","recreate");

	//Create event
	Event *event = &pythia.event;

	//Define a histograms into which we accumulate the invariant mass distribution for electron combinations
 
    	TH1F *omega_invmass = new TH1F("omega_invmass","Reconstructed omega invariant mass from electron pairs", 100, 0.0, 5);
    	omega_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *electron_number_event = new TH1F("electron_number_event","Electron-antielectron number per event", 10000, 0.0, 10000.0);
    	electron_number_event -> GetXaxis()-> SetTitle("event index");
	electron_number_event -> GetYaxis()-> SetTitle("number of #e^{#pm}");

	//Create TTree
	TTree *T = new TTree("T","ev1 Tree");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, mother2_var, motherid1_var, motherid2_var;

	TBranch *index = T->Branch("index", &index_var);
	TBranch *id = T->Branch("id", &id_var);
	TBranch *energy = T->Branch("energy", &energy_var);
	TBranch *mass = T->Branch("mass", &mass_var);
	//TBranch *p = T->Branch("p", &p_var);
	TBranch *px = T->Branch("px", &px_var);
	TBranch *py = T->Branch("py", &py_var);
	TBranch *pz = T->Branch("pz", &pz_var);
	TBranch *mother1 = T->Branch("mother1", &mother1_var);
	TBranch *mother2 = T->Branch("mother2", &mother2_var);
	TBranch *motherid1 = T->Branch("motherid1", &motherid1_var);
	TBranch *motherid2 = T->Branch("motherid2", &motherid2_var);

	//How many events shall we generate?
    	Long64_t nEvents=100000;

	double inv_mass;
	const int nPrint=5;

	Long64_t total_final_state = 0;
	Long64_t total_electron_antielectron = 0;

	bool antielectron;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia.next()) continue;

		//Print out entire event contents and decays for first five events.
		if (iEvent < nPrint) {pythia.info.list(); pythia.event.list();}

		Long64_t final_state = 0;
		Long64_t electron_number = 0;
		Long64_t antielectron_number = 0;

		antielectron = false; //sets anti-electrons status to false, as in they haven't been saved in the tree

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia.event.size(); ++i) {	

			//but only consider those that are "final state", i.e. ignore
			//intermediate ones that have decayed.
			if (pythia.event[i].isFinal()){

				++final_state;

				//Looking for electron pairs: first find a electron
				if (pythia.event[i].id()==11){

					++electron_number;

					//adding electron to the tree
					index_var = iEvent;
					id_var = pythia.event[i].id();
					energy_var = pythia.event[i].e();
					mass_var = pythia.event[i].m();
					px_var = pythia.event[i].px();
					py_var = pythia.event[i].py();
					pz_var = pythia.event[i].pz();
					mother1_var = pythia.event[i].mother1();
					mother2_var = pythia.event[i].mother2();
					motherid1_var = pythia.event[pythia.event[i].mother1()].id();
					motherid2_var = pythia.event[pythia.event[i].mother2()].id();

					T->Fill();
					//std::cout<<"iEvent = "<<iEvent<<"    "<<" id = "<<id_var<<"    "<< " energy = "<<energy_var<<std::endl;

				   	for (Long64_t j = 0; j < pythia.event.size(); ++j){

						//If particle is anti-electron
						if (pythia.event[j].id()==-11){

							//check if the anti-electron has not been filled in the tree (prevent multiple saving of the same anti-electron)
							if(antielectron == false){

								++antielectron_number;

								//adding anti-electron to the tree
								index_var = iEvent;
								id_var = pythia.event[j].id();
								energy_var = pythia.event[j].e();
								mass_var = pythia.event[j].m();
								px_var = pythia.event[j].px();
								py_var = pythia.event[j].py();
								pz_var = pythia.event[j].pz();
								mother1_var = pythia.event[j].mother1();
								mother2_var = pythia.event[j].mother2();
								motherid1_var = pythia.event[pythia.event[j].mother1()].id();
								motherid2_var = pythia.event[pythia.event[j].mother2()].id();

								T->Fill();
								//std::cout<<"i = "<<i<<"    "<<"j = "<<j<<"    "<<"iEvent = "<<iEvent<<"    "<<"id = "<<id_var<<"    "<< "energy = "<<energy_var<<std::endl;
								
							}

							//Reconstructing invariant mass of all electron-anti-electron combinations
							inv_mass = invmass(pythia.event[i].p(), pythia.event[j].p());

							//Plot histogram of all electron-anti-electron combinations
							omega_invmass->Fill(inv_mass);	

						}

					}
					antielectron = true;//all the anti-electrons have been saved in the event

				}

			} 

		}

		//std::cout<<"final state = "<<final_state<<setw(5)<<" electron number = "<<electron_number<<setw(5)<<" anti electron = "<<antielectron_number<<setw(5)<<" difference = "<<final_state - electron_number - antielectron_number<<std::endl;

		total_final_state = total_final_state + final_state;
		total_electron_antielectron = total_electron_antielectron + electron_number + antielectron_number;
		electron_number_event->Fill(total_electron_antielectron);
	
	}// End event loop.

	//std::cout<<"total final state = "<<total_final_state<<setw(5)<<" electron+antielectron number = "<<total_electron_antielectron<<std::endl;

	// Statistics on event generation.
	//pythia.stat();

	//  Print and Write tree.
	T->Print();
	T->Write();

	file->Write();
	delete file;

	// Done.
	return 0;
}
