#include "AddOns/Analysis/Triggers/MySISCone.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Analysis/Triggers/momentum.h"
#include "AddOns/Analysis/Triggers/siscone.h"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;
using namespace siscone;


SISCone::SISCone(ATOOLS::Particle_Qualifier_Base * const qualifier,double f) : 
  Jet_Algorithm_Base(qualifier), m_f(f)
{
  p_ConeCluster = new Csiscone();
}

SISCone::~SISCone() 
{
  delete p_ConeCluster;
}

Flavour SISCone::GetBFlavour(const Vec4D & mom) {
  Particle_List particles;
  for (Particle_List::const_iterator it=p_orig->begin(); it!=p_orig->end();++it) {
    if ((*p_qualifier)(*it)) {
      if ((*it)->Momentum()==mom) {
	if (((*it)->Flav()).Kfcode()==kf_b || ((*it)->Flav()).Kfcode()==kf_bjet) 
	  return Flavour(kf_bjet);
	return Flavour(kf_jet);
      }
      particles.push_back(*it);
    }
  }
  if (particles.size()>31) {
    msg_Error()<<METHOD<<": too many particles to determine jet flavour"<<endl;
    return Flavour(kf_jet);
  }
  size_t nb=particles.size();
  size_t cnt=1<<nb;
  for (size_t i=3;i<cnt;i++) {
    Vec4D sum(0.,0.,0.,0.);
    int bf(0);
    for (size_t j=0;j<nb;j++) if (i&(1<<j)) {
      sum+=particles[j]->Momentum();
      if (m_bflag==0) bf+=((particles[j]->Flav()).Kfcode()==kf_b || 
			   (particles[j]->Flav()).Kfcode()==kf_bjet);
      else if (m_bflag==-1) bf+=(1-2*particles[j]->Flav().IsAnti())*
	((particles[j]->Flav()).Kfcode()==kf_b || (particles[j]->Flav()).Kfcode()==kf_bjet);      
    }
    if (sum==mom) {
      if (bf!=0) return Flavour(kf_bjet);
      return Flavour(kf_jet);
    }
  }
  msg_Error()<<METHOD<<": could not determine jet flavour"<<endl;
  msg_Error()<<" mom="<<mom<<endl;
  for (size_t j=0;j<nb;j++) msg_Error()<<j<<" "<<particles[j]->Momentum()<<" "<<particles[j]->Flav()<<endl;
  return Flavour(kf_jet);  
}

void SISCone::AddToKtlist(double kt2) {
  if (p_kts) {
    p_kts->push_back(kt2);
  }
}

void SISCone::AddToJetlist(const Vec4D & mom) {
  if (p_jets) {
    Flavour jf=(m_bflag==1)?Flavour(kf_jet):GetBFlavour(mom);
    p_jets->push_back(new Particle(p_jets->size(),jf,mom));
  }
}

bool SISCone::ConstructJets(const Particle_List *pl, Particle_List * jets,
				     std::vector<double> * kt2, double R)
{
  // assume empty containers
  p_orig = pl;
  p_jets = jets;
  p_kts  = kt2;

  // create vector list from particle list
  vector<Cmomentum> particles;
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    if ((*p_qualifier)(*it)) {
      particles.push_back(Cmomentum((*it)->Momentum()[1],(*it)->Momentum()[2],(*it)->Momentum()[3],
				    (*it)->Momentum()[0]));
    }
  }

  // cluster
  p_ConeCluster->compute_jets(particles, R, m_f);
  for (vector<Cjet>::iterator it_j= p_ConeCluster->jets.begin();
       it_j!= p_ConeCluster->jets.end();it_j++) {
    AddToJetlist(Vec4D(it_j->v.E,it_j->v.px,it_j->v.py,it_j->v.pz));
    AddToKtlist(it_j->v.perp());
  }

  p_orig = 0;
  p_jets=0;
  p_kts =0;

  return true;
}


