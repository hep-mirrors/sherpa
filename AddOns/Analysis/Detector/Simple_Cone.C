#include "AddOns/Analysis/Detector/Simple_Cone.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include <algorithm>

#ifdef USING__ROOT
#include "ATOOLS/Math/Scaling.H"
#include "TH1D.h"
#include "TH2D.h"
#endif 

using namespace ANALYSIS;
using namespace ATOOLS;


Simple_Cone::Simple_Cone(const double Etcut,const double Etmin, double sep) : 
  m_Etcut(Etcut), m_Etmin(Etmin), m_dR(sep), m_dR2(sep*sep),
  m_flav(Flavour(kf_jet))
{ 
#ifdef USING__ROOT
  std::string name("JES");
  (*MYROOT::myroot)(new TH2D(name.c_str(),name.c_str(),45,20.,200.,50,0.2,2.0),name);
#endif 
}

Simple_Cone::~Simple_Cone() {}

bool Simple_Cone::ConstructJets(ObjectList * jets) {
  CalcJets(jets);
  if (jets->size()>0) return true;
  return false;
}


void Simple_Cone::CalcJets(ObjectList * jets)
{
  bool newjet(false);
  std::list<Cell *> * ecells(p_ECal->GetHitCells()), * hcells(p_HCal->GetHitCells());
  if ((!ecells || ecells->size()==0) && (!hcells || hcells->size()==0)) return;
  
  std::list<Track *> mctracks;
  p_MC->GetTracks(mctracks);

  std::vector<Cell *> cone;

  std::list<Cell *>::iterator cit;
  std::list<Track *>::iterator trit;
  
  double maxet(m_Etmin),eta(0.),phi(0.);
  
  Cell * seed(NULL);
  std::set<Cell *> badseeds;
  
  Reconstructed_Object * jet(NULL);
  msg_Debugging()<<METHOD<<" :"<<std::endl;
  do {
    newjet = false;
    maxet  = m_Etmin;
    if (hcells && hcells->size()>0) {
      for (cit=hcells->begin();cit!=hcells->end();cit++) {
	(*cit)->Centroid(eta,phi);
	if (!(*cit)->Used() && (*cit)->EPerp()>maxet &&
	    badseeds.find((*cit))==badseeds.end()) {
	  seed  = (*cit);
	  maxet = (*cit)->EPerp();
	}
      }
    }
    if (ecells&&ecells->size()>0) {
      for (cit=ecells->begin();cit!=ecells->end();cit++) {
	(*cit)->Centroid(eta,phi);
	if (!(*cit)->Used() && (*cit)->EPerp()>maxet &&
	    badseeds.find((*cit))==badseeds.end()) {
	  seed  = (*cit);
	  maxet = (*cit)->EPerp();
	}
      }
    }

    if (maxet>m_Etmin) {
      msg_Debugging()<<"   Found potential cone-seed ("<<maxet<<"): "
	       <<seed->ParticleEntries()->begin()->first->Flav()<<" with "
	       <<seed->ParticleEntries()->begin()->first->Momentum()<<" --> ";

      newjet = true;
      cone.clear();
      seed->Centroid(eta,phi);
      badseeds.insert(seed);
      cone.push_back(seed);

      msg_Debugging()<<" adding in cells:"<<std::endl;
      for (cit=hcells->begin();cit!=hcells->end();cit++) {
	if ((*cit)!=seed && !(*cit)->Used() &&  
	    (*cit)->R2(eta,phi)<m_dR2) {
	  msg_Debugging()<<"      add hcal :"
		   <<(*cit)->ParticleEntries()->begin()->first->Flav()<<" with "
		   <<(*cit)->ParticleEntries()->begin()->first->Momentum()<<std::endl;
	  cone.push_back((*cit));
	  (*cit)->SetUsed(true);  
	}
      }
      for (cit=ecells->begin();cit!=ecells->end();cit++) {
	if ((*cit)!=seed && !(*cit)->Used() && 
	    (*cit)->R2(eta,phi)<m_dR2) {
	  msg_Debugging()<<"      add ecal :"
		   <<(*cit)->ParticleEntries()->begin()->first->Flav()<<" with "
		   <<(*cit)->ParticleEntries()->begin()->first->Momentum()<<std::endl;
	  cone.push_back((*cit));
	  (*cit)->SetUsed(true);	
	}
      }
      jet = new Reconstructed_Object(m_flav,cone);
      for (trit=mctracks.begin();trit!=mctracks.end();trit++) {
	if (!(*trit)->used && (sqr(eta-(*trit)->eta)+sqr(phi-(*trit)->phi))<m_dR2) {
	  jet->AddTrack((*trit));
	  (*trit)->used = true;
	}
      }
      jet->SetIncludeTracks(true);
      jet->Update();
      
      msg_Debugging()<<"   Cone momentum = "<<jet->Mom()<<" --> ";
      if (jet->Mom().EPerp()<m_Etcut) {
	msg_Debugging()<<"too little energy in cone ("<<jet->Mom().EPerp()<<") delete it."<<std::endl;
	jet->SetUsed(false);
	delete jet;
      }
      else {
	jets->push_back(jet);
	msg_Debugging()<<"enough energy in cone ("<<jet->Mom().EPerp()<<") take it."<<std::endl;
#ifdef USING__ROOT
	std::string name("JES");
	((TH2D*)(*MYROOT::myroot)[name])->Fill(jet->Mom().EPerp(),
					       jet->Mom().EPerp()/jet->TrueMom().EPerp(),1.);
#endif 
      }
    }
  } while (newjet);
}


void Simple_Cone::FillShape(const int jetno,ATOOLS::Histogram * histo,
			    const double weight,const double ncount)
{
}
