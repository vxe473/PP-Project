//12/10/2017
//Program reading in the tree generated in gs_project2.cc and performing analysis on the reconstructed invariant mass
//credits to Phil for the reading the tree from file part
//credits to Nigel for analysis procedure suggestions
//#include "treeio.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "cmath"
#include <iostream>
#include <vector>

struct vect{
	
	std::vector<double> index, id, energy, mass, px, py, pz, mother1, mother2, motherid1, motherid2;
	
};

double invmass(vect *v, Long64_t i, Long64_t j, Long64_t k){

	//std::cout<<"i = "<<i<<"    "<<"j = "<<j<<"    "<<"k = "<<k<<"    "<< "energy i = "<<v->energy[i]<<"    "<<"energy j = "<<v->energy[j]<<"energy k = "<<v->energy[k]<<std::endl;

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j]+v->energy[k],2)-pow(v->px[i]+v->px[j]+v->px[k],2)-pow(v->py[i]+v->py[j]+v->py[k],2)-pow(v->pz[i]+v->pz[j]+v->pz[k],2));

	//std::cout<<"invariant mass = "<<setprecision(15)<<inv_mass<<std::endl;

	return inv_mass;

}

double invmass2(vect *v, Long64_t i, Long64_t j){

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j],2)-pow(v->px[i]+v->px[j],2)-pow(v->py[i]+v->py[j],2)-pow(v->pz[i]+v->pz[j],2));

	return inv_mass;

}

double p(vect *v, Long64_t i){

	double p1=sqrt(pow(v->px[i],2)+pow(v->py[i],2)+pow(v->pz[i],2));

	return p1;

}

double pt(vect *v, Long64_t i){

	double pt1=sqrt(pow(v->px[i],2)+pow(v->py[i],2));

	return pt1;

}

void analyze_event(vect *v, TH1F *omega_invmass, TH1F *electronpi_invmass_mother, TH1F *omega_invmass_mother, TH1F *electron_invmass_mother, TH1F *electron_pt_mother, TH1F *electron_ptcut_mother, TH1F *pi_pt_mother, TH1F *electron_number_event, TH1F *pi_number_event, TH1F *total_number_event){

	double inv_mass;

	Long64_t nentries = v->index.size();

//std::cout<<"event size = "<<nentries<<std::endl;

	Long64_t electron_number = 0;
	Long64_t antielectron_number = 0;
	Long64_t pi_number = 0;

	bool antielectron = false; 
	bool pi = false; 

	for(Long64_t i=0;i<nentries;i++){

		//if entry is electron
		if(v->id[i] == 11){

			++electron_number;

			//Loop through the event and find an anti-electron
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-electron
				if(v->id[j] == -11){

					if(antielectron == false){
						++antielectron_number;
					}

					for(Long64_t k=0;k<nentries;k++){

						//if entry is pi-0
						if(v->id[k] == 111){

							if(pi == false){
								++pi_number;
							}

							//Reconstructing invariant mass of all electron-antielectron-pi-0 combinations
							inv_mass = invmass(v, i, j, k);

							//Plot histogram of all electron-antielectron-pi-0 combinations
							omega_invmass->Fill(inv_mass);

							if((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k])){

								//If the electron-antielectron-pi-0 pair come from the same mother particle plot them in a separate histogram

								if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j])){

									//If the electron-antielectron pair come from the same mother particle plot them in a separate histogram
									electronpi_invmass_mother->Fill(inv_mass);

									if(v->motherid1[i] == 223){

										omega_invmass_mother->Fill(inv_mass);

										double electron_invmass = invmass2(v, i, j);
										electron_invmass_mother->Fill(electron_invmass);

										double pt_mother1 = pt(v, i);
										electron_pt_mother->Fill(pt_mother1);

										double pt_mother2 = pt(v, j);
										electron_pt_mother->Fill(pt_mother2);

										double p_mother1 = p(v, i);
										//electron_p_mother->Fill(p_mother1);

										if((pt_mother1>0.5) && (p_mother1>10.0))
										{
										electron_ptcut_mother->Fill(pt_mother1);
										}

										double p_mother2 = p(v, j);
										//electron_p_mother->Fill(p_mother2);

										if((pt_mother2>0.5) && (p_mother2>10.0))
										{
										electron_ptcut_mother->Fill(pt_mother2);
										}

										double pt_mother3 = pt(v, k);
										pi_pt_mother->Fill(pt_mother3);

									}

								}

							}

						}
					}
					pi = true;

				}

			}
			antielectron = true;

		}

	}

	//if(int(v->index[0]) % 1000 == 0){std::cout<<"electron number = "<<electron_number<<setw(5)<<"anti electron = "<<antielectron_number<<setw(5)<<"pi-0 = "<<pi_number<<std::endl;}

	if (v->index.size() > 0) 
	{
		double total_electron_antielectron = electron_number+antielectron_number;
		double total_number = electron_number+antielectron_number+pi_number;
		electron_number_event->SetBinContent((v->index[1])+1,total_electron_antielectron);
		pi_number_event->SetBinContent((v->index[1])+1,pi_number);
		total_number_event->SetBinContent((v->index[1])+1,total_number);

	}
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

}

int gs_project5() {

        gROOT->ProcessLine(".x lhcbStyle.C");

	// Open the input TFile 
	TFile  *file_in = new TFile("gs_project3.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project5.root", "recreate");

	//Get tree from file
	TTree *T = (TTree*)file_in->Get("T");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, mother2_var, motherid1_var, motherid2_var;

	//Get branches of the tree
	T->SetBranchAddress("index",&index_var);
   	T->SetBranchAddress("id",&id_var);
	T->SetBranchAddress("energy",&energy_var);
   	T->SetBranchAddress("mass",&mass_var);
	T->SetBranchAddress("px",&px_var);
   	T->SetBranchAddress("py",&py_var);
	T->SetBranchAddress("pz",&pz_var);
	T->SetBranchAddress("mother1",&mother1_var);
	T->SetBranchAddress("mother2",&mother2_var);
	T->SetBranchAddress("motherid1",&motherid1_var);
	T->SetBranchAddress("motherid2",&motherid2_var);

	TH1F *omega_invmass = new TH1F("omega_invmass","Reconstructed omega invariant mass from electron pairs and pi", 400, 0.0, 40.0);
    	omega_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *electronpi_invmass_mother = new TH1F("electronpi_invmass_mother","Reconstructed electron pair + #pi invariant mass with same mother", 100, 0.0, 1.0);
    	electronpi_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *omega_invmass_mother = new TH1F("omega_invmass_mother","Reconstructed invariant mass from electron pairs + #pi with same mother (#omega)", 100, 0.0, 1.0);
    	omega_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *electron_invmass_mother = new TH1F("electron_invmass_mother","Reconstructed electron pair invariant mass with same mother(#omega)", 2000, 0.0, 1.0);
    	electron_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");;
    	electron_invmass_mother -> GetYaxis()-> SetTitle("Number of Entries");
    	electron_invmass_mother -> SetLineColor(2);

	TH1F *electron_pt_mother = new TH1F("electron_pt_mother","#e^{-} and #e^{+} p_{t} distribution with same mother (#omega)", 100, 0.0, 1.5);
    	electron_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	//electron_pt_mother -> GetYaxis()-> SetLog();

	TH1F *electron_ptcut_mother = new TH1F("electron_ptcut_mother","#e^{-} and #e^{+} p_{t} distribution after cuts (with same mother (#omega))", 150, 0.0, 1.5);
    	electron_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	
	TH1F *pi_pt_mother = new TH1F("pi_pt_mother","#pi p_{t} distribution as #omega decay product", 150, 0.0, 1.5);
    	electron_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	//electron_pt_mother -> GetYaxis()-> SetLog();

	TH1F *electron_number_event = new TH1F("electron_number_event","Electron-antielectron number per event", 100000, 0.0, 100000.0);
    	electron_number_event -> GetXaxis()-> SetTitle("event index");
	electron_number_event -> GetYaxis()-> SetTitle("number of #e^{#pm}");

	TH1F *pi_number_event = new TH1F("pi_number_event","#pi number per event", 100000, 0.0, 100000.0);
    	pi_number_event -> GetXaxis()-> SetTitle("event index");
	pi_number_event -> GetYaxis()-> SetTitle("number of #pi");

	TH1F *total_number_event = new TH1F("total_number_event","Total number of particles per event", 100000, 0.0, 100000.0);
    	total_number_event -> GetXaxis()-> SetTitle("event index");
	total_number_event -> GetYaxis()-> SetTitle("number of particles");

	vect* v= new vect;//declares a vector in which we are going to store all the particles from the same event
	double prev_index = 0.0;//saving the event index of the previous entry of the tree

	//Loop through the entries of the tree
	Long64_t nentries = T->GetEntries();
std::cout<<"Total number of entries: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T->GetEntry(i);

		//checks if the entry is the very last one in the tree
		if(i+1==nentries){

			//add last particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			analyze_event(v, omega_invmass, electronpi_invmass_mother, omega_invmass_mother, electron_invmass_mother, electron_pt_mother, electron_ptcut_mother, pi_pt_mother, electron_number_event, pi_number_event, total_number_event);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}

		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){

			analyze_event(v, omega_invmass, electronpi_invmass_mother, omega_invmass_mother, electron_invmass_mother, electron_pt_mother, electron_ptcut_mother, pi_pt_mother, electron_number_event, pi_number_event, total_number_event);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			prev_index = index_var;//setting the index for current event

		}
		else{

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}
	      	
   	}
    
	TCanvas *c1=new TCanvas("c1","",600,600);

	omega_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_omegainvmass.pdf","pdf");

	electron_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_electroninvmass_mother.pdf","pdf");

	omega_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_omegainvmass_mother.pdf","pdf");

	electron_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_electron_pt_mother.pdf","pdf");

	electron_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_electron_ptcut_mother.pdf","pdf");

	pi_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_pi_pt_mother.pdf","pdf");

	electron_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_electron_number_event.pdf","pdf");

	pi_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_pi_number_event.pdf","pdf");

	total_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_total_number_event.pdf","pdf");

	file_out->Write();

	delete file_in;
	delete file_out;

	// Done.
	return 0;
}


