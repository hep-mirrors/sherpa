#include "Cluster.H"
#include "Hadronisation_Parameters.H"
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
  m_type(ctp::no),
  m_momentum(Vec4D(0.,0.,0.,0.)),
  m_hasboost(false), p_left(NULL), p_right(NULL), p_prev(NULL)
{
  m_momenta[0]  = m_momenta[1]  = Vec4D(0.,0.,0.,0.);  
  m_flavours[0] = m_flavours[1] = Flavour(kf::none);
}

Cluster::Cluster(const ATOOLS::Flavour & flav1,const ATOOLS::Vec4D & mom1,
		 const ATOOLS::Flavour & flav2,const ATOOLS::Vec4D & mom2) :
  m_hasboost(false), p_left(NULL), p_right(NULL), p_prev(NULL)
{
  m_momentum = mom1+mom2;
  if (flav1.IsQuark()) {
    if (!flav1.IsAnti()) {
      if (flav2.IsQuark() && flav2.IsAnti()) {
	m_flavours[0] = flav1; m_flavours[1] = flav2; 
	m_momenta[0]  = mom1;  m_momenta[1]  = mom2; 
	m_type = ctp::qq; return;
      }
      if (flav2.IsDiQuark() && !flav2.IsAnti()) {
	m_flavours[0] = flav1; m_flavours[1] = flav2; 
	m_momenta[0]  = mom1;  m_momenta[1]  = mom2; 
	m_type = ctp::qd; return;
      }
    }
    else {
      if (flav2.IsQuark() && !flav2.IsAnti()) {
	m_flavours[0] = flav2; m_flavours[1] = flav1; 
	m_momenta[0]  = mom2;  m_momenta[1]  = mom1; 
	m_type = ctp::qq; return;
      }
      if (flav2.IsDiQuark() && flav2.IsAnti()) {
	m_flavours[0] = flav2; m_flavours[1] = flav1; 
	m_momenta[0]  = mom2;  m_momenta[1]  = mom1; 
	m_type = ctp::dq; return;
      }
    }
  }
  else if (flav1.IsDiQuark()) {
    if (flav1.IsAnti()) {
      if (flav2.IsQuark() && flav2.IsAnti()) {
	m_flavours[0] = flav1; m_flavours[1] = flav2; 
	m_momenta[0]  = mom1;  m_momenta[1]  = mom2; 
	m_type = ctp::dq; return;
      }
      if (flav2.IsDiQuark() && !flav2.IsAnti()) {
	m_flavours[0] = flav1; m_flavours[1] = flav2; 
	m_momenta[0]  = mom1;  m_momenta[1]  = mom2; 
	m_type = ctp::dd; return;
      }
    }
    else {
      if (flav2.IsQuark() && !flav2.IsAnti()) {
	m_flavours[0] = flav2; m_flavours[1] = flav1; 
	m_momenta[0]  = mom2;  m_momenta[1]  = mom1; 
	m_type = ctp::qd; return;
      }
      if (flav2.IsDiQuark() && flav2.IsAnti()) {
	m_flavours[0] = flav2; m_flavours[1] = flav1; 
	m_momenta[0]  = mom2;  m_momenta[1]  = mom1; 
	m_type = ctp::dd; return;
      }
    }
  }

  msg.Error()<<"Error in Cluster::Cluster("<<flav1<<","<<flav2<<") :"<<std::endl
	     <<"   Cannot handle this colour structure."<<std::endl
	     <<"   Abort the run."<<std::endl;
  abort();
}

Cluster::~Cluster() 
{
  if (p_left)  { delete p_left;  p_left=NULL;  }
  if (p_right) { delete p_right; p_right=NULL; }
}

void Cluster::Update()
{
  m_momentum = m_momenta[0]+m_momenta[1];
  if (m_flavours[0].IsQuark()) {
    if (!m_flavours[0].IsAnti()) {
      if (m_flavours[1].IsQuark() && m_flavours[1].IsAnti())    { m_type = ctp::qq; return; }
      if (m_flavours[1].IsDiQuark() && !m_flavours[1].IsAnti()) { m_type = ctp::qd; return; }
    }
    else {
      if (m_flavours[1].IsQuark() && !m_flavours[1].IsAnti())   { m_type = ctp::qq; return; }
      if (m_flavours[1].IsDiQuark() && m_flavours[1].IsAnti())  { m_type = ctp::dq; return; }
    }
  }
  else if (m_flavours[0].IsDiQuark()) {
    if (m_flavours[0].IsAnti()) {
      if (m_flavours[1].IsQuark() && m_flavours[1].IsAnti())    { m_type = ctp::dq; return; }
      if (m_flavours[1].IsDiQuark() && !m_flavours[1].IsAnti()) { m_type = ctp::dd; return; }
    }
    else {
      if (m_flavours[1].IsQuark() && !m_flavours[1].IsAnti())   { m_type = ctp::qd; return; }
      if (m_flavours[1].IsDiQuark() && m_flavours[1].IsAnti())  { m_type = ctp::dd; return; }
    }  
  }
  msg.Error()<<"ERROR in Cluster::Update():"<<std::endl
	     <<"   Funny flavour/colour structure : "<<m_flavours[0]<<", "<<m_flavours[1]<<std::endl
	     <<"   Will abort."<<std::endl;
  abort();
}

void Cluster::RescaleMomentum(ATOOLS::Vec4D newmom)
{
  Poincare rest(m_momentum);
  Poincare back(newmom);

  rest.Boost(m_momenta[0]);
  back.BoostBack(m_momenta[0]);
  rest.Boost(m_momenta[1]);
  back.BoostBack(m_momenta[1]);

  Vec4D help = m_momentum;
  m_momentum = newmom;

  Vec4D testmom = m_momentum-m_momenta[0]-m_momenta[1];
  if (dabs(testmom.Abs2())>1.e-6 || testmom[0]>1.e-6) {
    std::cout<<"Maybe error in RescaleMomentum("<<help<<") : "<<testmom<<std::endl<<(*this)
	     <<m_momentum<<" "<<Momentum(1)<<" "<<Momentum(2)<<" :: "
	     <<help.Abs2()<<" <---> "<<newmom.Abs2()<<std::endl;
  }
}


ATOOLS::Vec4D Cluster::Momentum(const int i) const { 
  if (i==0) return m_momentum;        
  if (i<3)  return m_momenta[i-1];
  ATOOLS::msg.Error()<<"Error in Cluster::Momentum("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster momentum."<<std::endl;
  abort();
  return m_momentum;
}

void Cluster::SetMomentum(const int i,const ATOOLS::Vec4D & mom) { 
  if (i==0) { m_momentum     = mom; return; }       
  if (i<3)  { m_momenta[i-1] = mom; return; }       
  ATOOLS::msg.Error()<<"Error in Cluster::SetMomentum("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster momentum."<<std::endl;
  abort();
}

double Cluster::Mass(const int i) const { 
  if (i==0) return sqrt(m_momentum.Abs2());        
  Flavour flav = m_flavours[i-1];
  if (i<3)  return  hadpars.GetConstituents()->Mass(flav);
  ATOOLS::msg.Error()<<"Error in Cluster::Mass("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster mass."<<std::endl;
  abort();
  return sqrt(m_momentum.Abs2());
}

void Cluster::SetFlav(const int i,const ATOOLS::Flavour & flav) { 
  if (i==1) { m_flavours[0] = flav; return; } 
  if (i==2) { m_flavours[1] = flav; return; } 
  ATOOLS::msg.Error()<<"Error in Cluster::SetFlav("<<i<<")"<<std::endl
		     <<"   Out of bounds. Do nothing and hope for the best."<<std::endl;  
  abort();
}

Flavour Cluster::GetFlav(const int i) const {
  if (i==1) { return m_flavours[0]; }
  if (i==2) { return m_flavours[1]; }
  ATOOLS::msg.Error()<<"Error in Cluster::GetFlav("<<i<<")"<<std::endl
		     <<"   Out of bounds. Do nothing and hope for the best."<<std::endl;  
  abort();
}

void Cluster::BoostInCMS() {
  if (m_hasboost) return;
  m_boost = ATOOLS::Poincare(m_momentum);
  m_boost.Boost(m_momentum);
  for (int i=0;i<2;i++) m_boost.Boost(m_momenta[i]);
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);
  m_hasboost = true;
}

void Cluster::BoostBack() {
  if (!m_hasboost) return;
  m_boost.BoostBack(m_momentum);
  for (int i=0;i<2;i++) m_boost.BoostBack(m_momenta[i]);
  if (p_left!=NULL)  p_left->BoostBack(m_boost);
  if (p_right!=NULL) p_right->BoostBack(m_boost);
  m_hasboost = false;
}

void Cluster::Boost(Poincare & boost) {
  boost.Boost(m_momentum);
  for (int i=0;i<2;i++) boost.Boost(m_momenta[i]);
  if (p_left!=NULL)  p_left->Boost(boost);
  if (p_right!=NULL) p_right->Boost(boost);
}

void Cluster::BoostBack(Poincare & boost) {
  boost.BoostBack(m_momentum);
  for (int i=0;i<2;i++) boost.BoostBack(m_momenta[i]);
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);
}

void Cluster::BoostBack(ATOOLS::Vec4D & mom) {
  if (!m_hasboost) return;
  m_boost.BoostBack(mom);
}

std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"-------------------------------------------------------------"<<std::endl
     <<"Cluster ("<<cluster.m_flavours[0]<<" "<<cluster.m_flavours[1]
     <<", "<<cluster.m_momentum<<","<<cluster.m_momentum.Abs2()<<" ): "<<std::endl
     <<cluster.m_momenta[0]<<" "<<cluster.m_momenta[1]<<std::endl;
  if (cluster.p_left)  str<<"  ---- Left  ----> "<<std::endl<<(*cluster.p_left);
  if (cluster.p_right) str<<"  ---- Right ----> "<<std::endl<<(*cluster.p_right);
  return str;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Cluster_List & cl) {
  s<<"Cluster List with "<<cl.size()<<" elements"<<std::endl;
  for (Cluster_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

