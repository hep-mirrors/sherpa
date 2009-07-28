/*
This is an example setup for an analysis of multiplicities and inclusive
quantities in hadron decays.
It has to include Main_FullDecay.C, and a FullDecay has to be
started with the following specialities:

  - ./FullDecay DECAYER=<kfcode> EVENTS=100000 ANALYSIS=1
*/

#include "Main_FullDecay.C"
#include "ATOOLS/Org/Shell_Tools.H"

#ifdef USING__ROOT
static TFile* rootfile;
static map<Flavour,TH1D*> multiplicities;
static map<Flavour,TH1D*> energies;
static map<Flavour,TH1D*> directenergies;
typedef map<Flavour,TH1D*>::iterator FTMapIt;
static map<Flavour,int> int_multiplicities;
typedef map<Flavour,int>::iterator FIMapIt;
#endif

void InitialiseAnalysis()
{
#ifdef USING__ROOT
  std::string adir = "";
  ATOOLS::MakeDir(adir,493);
  rootfile = new TFile(string(adir+"/Multiplicities_"+
                              mother_flav.ShellName()+".root").c_str(), "RECREATE");
  list<Flavour> flavours;
  flavours.push_back(Flavour(kf_pi_plus));
  flavours.push_back(Flavour(kf_pi_plus).Bar());
  flavours.push_back(Flavour(kf_K_plus));
  flavours.push_back(Flavour(kf_K_plus).Bar());
  flavours.push_back(Flavour(kf_K_L));
  flavours.push_back(Flavour(kf_e));
  flavours.push_back(Flavour(kf_e).Bar());
  flavours.push_back(Flavour(kf_nue));
  flavours.push_back(Flavour(kf_nue).Bar());
  flavours.push_back(Flavour(kf_numu));
  flavours.push_back(Flavour(kf_numu).Bar());
  flavours.push_back(Flavour(kf_mu));
  flavours.push_back(Flavour(kf_mu).Bar());
  flavours.push_back(Flavour(kf_tau));
  flavours.push_back(Flavour(kf_tau).Bar());
  flavours.push_back(Flavour(kf_nutau));
  flavours.push_back(Flavour(kf_nutau).Bar());
  flavours.push_back(Flavour(kf_n));
  flavours.push_back(Flavour(kf_n).Bar());
  flavours.push_back(Flavour(kf_p_plus));
  flavours.push_back(Flavour(kf_p_plus).Bar());
  flavours.push_back(Flavour(kf_photon));

  for(list<Flavour>::iterator it=flavours.begin(); it!=flavours.end(); it++) {
    multiplicities[*it] =    makeTH1D(it->ShellName()+"_multi", "",
                                      9, -0.5, 8.5,
                                      it->RootName()+" multiplicity","Events");
    energies[*it]       =    makeTH1D(it->ShellName()+"_energy", "",
                                      25, 0.0, mother_flav.HadMass()/2.0,
                                      "E_{"+it->RootName()+"}",
                                      "#frac{1}{#Gamma} #frac{d#Gamma}{dE}");
    directenergies[*it] =    makeTH1D(it->ShellName()+"_directenergy", "",
                                      25, 0.0, mother_flav.HadMass()/2.0,
                                      "E_{"+it->RootName()+"}",
                                      "#frac{1}{#Gamma} #frac{d#Gamma}{dE}");
    int_multiplicities[*it] = 0;
  }
#endif
}


void AnalyseEvent(Blob_List* blobs)
{
#ifdef USING__ROOT
  int outgoing=1;
  Particle_List outparts = blobs->ExtractParticles(part_status::active, outgoing);
  for(Particle_List::iterator part=outparts.begin(); part!=outparts.end(); part++) {
    if(multiplicities.find((*part)->Flav())!=multiplicities.end()) {
      int_multiplicities[(*part)->Flav()]++;
      energies[(*part)->Flav()]->Fill((*part)->Momentum()[0]);
      if((*part)->ProductionBlob()->InParticle(0)->Flav()==mother_flav) {
        directenergies[(*part)->Flav()]->Fill((*part)->Momentum()[0]);
      }
    }
    else {
      PRINT_INFO("stable but unanalysed outparticle encountered: "
                 <<**part);
    }
  }
  for(FIMapIt it=int_multiplicities.begin(); it!=int_multiplicities.end(); it++) {
    Flavour flav = it->first;
    multiplicities[flav]->Fill(it->second);
    it->second = 0;
  }
#endif
}


void FinishAnalysis()
{
#ifdef USING__ROOT
  for(FTMapIt it=multiplicities.begin(); it!=multiplicities.end(); it++) {
    it->second->Write();
  }
  for(FTMapIt it=energies.begin(); it!=energies.end(); it++) {
    it->second->Write();
  }
  for(FTMapIt it=directenergies.begin(); it!=directenergies.end(); it++) {
    it->second->Write();
  }
#endif
}
