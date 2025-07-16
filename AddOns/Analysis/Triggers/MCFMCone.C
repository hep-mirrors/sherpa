#include "AddOns/Analysis/Triggers/MCFMCone.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Triggers/ConeMCFM.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;


MCFMCone::MCFMCone(ATOOLS::Particle_Qualifier_Base * const qualifier,double rsep) : 
  Jet_Algorithm_Base(qualifier), m_rsep(rsep)
{
  p_ConeCluster = 0;
}

MCFMCone::~MCFMCone() 
{
  if (p_ConeCluster) delete p_ConeCluster;
}


void MCFMCone::AddToKtlist(double kt2) {
  if (p_kts) {
    p_kts->push_back(kt2);
  }
}

void MCFMCone::AddToJetlist(const Vec4D & mom) {
  if (p_jets) {
    p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_jet),mom));
  }
}

bool MCFMCone::ConstructJets(const Particle_List *pl, Particle_List * jets,
				     std::vector<double> * kt2, double R)
{
  if (!p_ConeCluster) p_ConeCluster=new ConeMCFM(m_rsep,R);

  // assume empty containers
  p_jets = jets;
  p_kts  = kt2;

  // create vector list from particle list
  vector<Vec4D> particles;
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    if ((*p_qualifier)(*it)) {
      particles.push_back((*it)->Momentum());
    }
  }

  // cluster
  p_ConeCluster->ConstructJets(particles);

  for (vector<Vec4D>::iterator it_j= p_ConeCluster->m_pjets.begin();
       it_j!= p_ConeCluster->m_pjets.end();it_j++) {
    AddToJetlist((*it_j));
    AddToKtlist((*it_j).PPerp());
  }

  p_jets=0;
  p_kts =0;

  return true;
}


