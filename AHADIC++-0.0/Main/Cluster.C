#include "Cluster.H"
#include "Message.H"
#include "Poincare.H"

using namespace AHADIC;
using namespace ATOOLS;

namespace AHADIC {
  std::ostream & operator<<(std::ostream & s, const Part_List & pl) {
    s<<"Particle List with "<<pl.size()<<" elements"<<std::endl;
    for (Part_Const_Iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
      s<<**pit<<std::endl;
    }
    return s;
  };
}

Cluster::Cluster() :
  p_fpair(NULL), m_momentum(Vec4D(0.,0.,0.,0.)), m_hasboost(false),
  p_left(NULL), p_right(NULL)
{
  p_parts[0] = p_parts[1] = NULL;
}

Cluster::Cluster(ATOOLS::Particle * trip, ATOOLS::Particle * antitrip) :
  p_fpair(new FlavourPair), m_hasboost(false),
  p_left(NULL), p_right(NULL)
{
  Flavour qflav   = trip->Flav();
  Flavour qbflav  = antitrip->Flav();
  if (((qflav.IsQuark() && !qflav.IsAnti()) ||
       (qflav.IsDiQuark() && qflav.IsAnti())) &&
      ((qbflav.IsQuark() && qbflav.IsAnti()) ||
       (qbflav.IsDiQuark() && !qbflav.IsAnti()))) {
    p_parts[0]      = trip;
    p_parts[1]      = antitrip;
    p_fpair->first  = qflav;
    p_fpair->second = qbflav;
  }
  else if (((qflav.IsQuark() && qflav.IsAnti()) ||
	    (qflav.IsDiQuark() && !qflav.IsAnti())) &&
	   ((qbflav.IsQuark() && !qbflav.IsAnti()) ||
	    (qbflav.IsDiQuark() && qbflav.IsAnti()))) {
    p_parts[0]      = antitrip;
    p_parts[1]      = trip;
    p_fpair->first  = qbflav;
    p_fpair->second = qflav;
  }
  else {
    msg.Error()<<"Error in Cluster::Cluster("<<trip<<","<<antitrip<<") :"<<std::endl
	       <<"   Cannot handle this colour structure."<<std::endl
	       <<"   Abort the run."<<std::endl;
    abort();
  }
  if (qflav.IsQuark() && qbflav.IsQuark())     m_type = 11;
  if (qflav.IsQuark() && qbflav.IsDiQuark())   m_type = 12;
  if (qflav.IsDiQuark() && qbflav.IsQuark())   m_type = 21;
  if (qflav.IsDiQuark() && qbflav.IsDiQuark()) m_type = 22;

  m_momentum = trip->Momentum()+antitrip->Momentum();
  //std::cout<<"New cluster "<<m_momentum<<std::endl;
}

Cluster::~Cluster() 
{
  if (p_left)  { delete p_left;  p_left=NULL;  }
  if (p_right) { delete p_right; p_right=NULL; }
  if (p_fpair) { delete p_fpair; p_fpair=NULL; }
  if (p_parts) {
    for (int i=0;i<2;i++) {
      if (p_parts[i] && p_parts[i]->ProductionBlob()==NULL) delete p_parts[i];
    }
  }
}


void Cluster::RescaleMomentum(ATOOLS::Vec4D newmom)
{
  Poincare rest(m_momentum);
  Poincare back(newmom);
  Vec4D help1 = p_parts[0]->Momentum();
  rest.Boost(help1);
  back.BoostBack(help1);
  p_parts[0]->SetMomentum(help1);
  Vec4D help2 = p_parts[1]->Momentum();
  rest.Boost(help2);
  back.BoostBack(help2);
  p_parts[1]->SetMomentum(help2);
  Vec4D help = m_momentum;
  m_momentum = newmom;


  Vec4D testmom = m_momentum-Momentum(0)-Momentum(1);
  if (dabs(testmom.Abs2())>1.e-6 || testmom[0]>1.e-6) {
    std::cout<<"Maybe error in RescaleMomentum: "<<testmom<<std::endl<<(*this)
	     <<m_momentum<<" "<<Momentum(0)<<" "<<Momentum(1)<<" :: "
	     <<m_momentum.Abs2()<<" <---> "<<newmom.Abs2()<<std::endl
	     <<help<<"  "<<help1<<" "<<help2<<std::endl;
  }
}

ATOOLS::Particle * Cluster::GetParticle(const int i) {
  if (i>-1&&i<2) return p_parts[i];
  ATOOLS::msg.Error()<<"Error in Cluster::GetParticle("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return NULL."<<std::endl;
  return NULL;
}

ATOOLS::Vec4D Cluster::Momentum(const int i) { 
  if (i==-1) return m_momentum;        
  if (i<2)   return p_parts[i]->Momentum();
  ATOOLS::msg.Error()<<"Error in Cluster::Momentum("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster momentum."<<std::endl;
  return m_momentum;
}

double Cluster::Mass(const int i) { 
  if (i==-1) return sqrt(m_momentum.Abs2());        
  if (i<2)   return p_parts[i]->Flav().Mass();
  ATOOLS::msg.Error()<<"Error in Cluster::Mass("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster mass."<<std::endl;
  return sqrt(m_momentum.Abs2());
}


void Cluster::BoostInCMS() {
  if (m_hasboost) return;
  m_boost    = ATOOLS::Poincare(m_momentum);
  m_boost.Boost(m_momentum);
  ATOOLS::Vec4D mom;
  for (int i=0;i<2;i++) {
    mom = p_parts[i]->Momentum();
    m_boost.Boost(mom);
    p_parts[i]->SetMomentum(mom);
  }
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);
  m_hasboost = true;
}

void Cluster::BoostBack() {
  if (!m_hasboost) return;
  m_boost.BoostBack(m_momentum);
  ATOOLS::Vec4D mom;
  for (int i=0;i<2;i++) {
    mom = p_parts[i]->Momentum();
    m_boost.BoostBack(mom);
    p_parts[i]->SetMomentum(mom);
  }
  if (p_left!=NULL)  p_left->BoostBack(m_boost);
  if (p_right!=NULL) p_right->BoostBack(m_boost);
  m_hasboost = false;
}

void Cluster::Boost(Poincare & boost) {
  boost.Boost(m_momentum);
  ATOOLS::Vec4D mom;
  for (int i=0;i<2;i++) {
    mom = p_parts[i]->Momentum();
    boost.Boost(mom);
    p_parts[i]->SetMomentum(mom);
  }
  if (p_left!=NULL)  p_left->Boost(boost);
  if (p_right!=NULL) p_right->Boost(boost);
}

void Cluster::BoostBack(Poincare & boost) {
  boost.BoostBack(m_momentum);
  ATOOLS::Vec4D mom;
  for (int i=0;i<2;i++) {
    mom = p_parts[i]->Momentum();
    boost.BoostBack(mom);
    p_parts[i]->SetMomentum(mom);
  }
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);
}

std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"-------------------------------------------------------------"<<std::endl
     <<"Cluster ("<<cluster.p_fpair->first<<" "<<cluster.p_fpair->second
     <<", "<<cluster.m_momentum<<","<<cluster.m_momentum.Abs2()<<" ): "<<std::endl
     <<(*cluster.p_parts[0])<<std::endl<<(*cluster.p_parts[1])<<std::endl;
  return str;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Cluster_List & cl) {
  s<<"Cluster List with "<<cl.size()<<" elements"<<std::endl;
  for (Cluster_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

