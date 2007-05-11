#include "Simple_Cone.H"
#include "Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;


Simple_Cone::Simple_Cone(const double Etcut,const double Etmin, double sep) : 
  m_Etcut(Etcut), m_Etmin(Etmin), m_dR(sep), m_dR2(sep*sep),
  m_flav(Flavour(kf::jet))
{ }

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
  std::list<Cell *>::iterator cit;
  
  double maxet,eta,phi;
  
  Cell * seed(NULL);
  std::set<Cell *> badseeds;
  
  Reconstructed_Object * jet(NULL);
  do {
    newjet = false;
    maxet  = m_Etmin;
    double eta,phi;
    if (hcells && hcells->size()>0) {
      for (cit=hcells->begin();cit!=hcells->end();cit++) {
	(*cit)->Centroid(eta,phi);
	//std::cout<<"Test hcal-cell:"<<(*cit)<<", E = "<<(*cit)->TotalDeposit()
	//	 <<" for "<<(*cit)->ParticleEntries()->begin()->first->Flav()
	//	 <<" at ("<<eta<<","<<phi<<")."<<std::endl;
	if (!(*cit)->Used() && (*cit)->EPerp()>maxet &&
	    badseeds.find((*cit))==badseeds.end()) {
	  //std::cout<<"       ....... take it."<<std::endl;
	  seed  = (*cit);
	  maxet = (*cit)->EPerp();
	}
      }
    }
    if (ecells&&ecells->size()>0) {
      for (cit=ecells->begin();cit!=ecells->end();cit++) {
	(*cit)->Centroid(eta,phi);
	//std::cout<<"Test ecal-cell: E = "<<(*cit)->TotalDeposit()
	//	 <<" for "<<(*cit)->ParticleEntries()->begin()->first->Flav()
	//	 <<" at ("<<eta<<","<<phi<<")."<<std::endl;
	if (!(*cit)->Used() && (*cit)->EPerp()>maxet &&
	    badseeds.find((*cit))==badseeds.end()) {
	  seed  = (*cit);
	  maxet = (*cit)->EPerp();
	}
      }
    }
    if (maxet>m_Etmin) {
      newjet = true;
      seed->Centroid(eta,phi);
      badseeds.insert(seed);
      //std::cout<<METHOD<<" try "<<maxet<<" at "<<eta<<", "<<phi<<" for "<<seed<<std::endl;
      std::vector<Cell *> * p_cone(new std::vector<Cell *>);
      p_cone->push_back(seed);
      for (cit=hcells->begin();cit!=hcells->end();cit++) {
	if ((*cit)!=seed && !(*cit)->Used() &&  
	    (*cit)->R2(eta,phi)<m_dR2) {
	  p_cone->push_back((*cit));
	  (*cit)->SetUsed(true);  
	}
      }
      for (cit=ecells->begin();cit!=ecells->end();cit++) {
	if ((*cit)!=seed && !(*cit)->Used() && 
	    (*cit)->R2(eta,phi)<m_dR2) {
	  p_cone->push_back((*cit));
	  (*cit)->SetUsed(true);	
	}
      }
      jet = new Reconstructed_Object(m_flav,p_cone);
      std::cout<<METHOD<<" : New jet : "<<jet<<" with "<<p_cone->size()<<" cells."<<std::endl;
      if (jet->Mom().EPerp()<m_Etcut) {
	std::cout<<"        ... delete it."<<std::endl;
	for (std::vector<Cell *>::iterator cit=p_cone->begin();
	     cit!=p_cone->end(); cit++) (*cit)->SetUsed(false);
	delete jet;
      }
      else {
	jets->push_back(jet);
      }
    }
  } while (newjet);
}


void Simple_Cone::FillShape(const int jetno,ATOOLS::Histogram * histo,
			    const double weight,const int ncount)
{
}
