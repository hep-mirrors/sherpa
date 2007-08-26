#include "Cluster.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"
#include "Poincare.H"

using namespace AHADIC;
using namespace ATOOLS;

namespace AHADIC {
  long int Cluster::s_cluster_count=0;
  long int Cluster::s_cluster_number=0;

  long int control::s_AHAparticles=0;
  long int control::s_AHAprotoparticles=0;
  long int control::s_AHAblobs=0;  
}


Cluster::Cluster() :
  m_type(ctp::no),m_leads(ltp::no),
  m_momentum(Vec4D(0.,0.,0.,0.)),
  m_hasboost(false), m_hasrotate(false), 
  p_left(NULL), p_right(NULL), p_prev(NULL), p_self(NULL),
  m_number(s_cluster_number++)
{
  //std::cout<<METHOD<<" for "<<m_number<<"."<<std::endl;
  m_momenta[0]  = m_momenta[1]  = Vec4D(0.,0.,0.,0.);  
  m_flavours[0] = m_flavours[1] = Flavour(kf::none);

  s_cluster_count++;
}

Cluster::Cluster(const ATOOLS::Flavour & flav1,const ATOOLS::Vec4D mom1,
		 const ATOOLS::Flavour & flav2,const ATOOLS::Vec4D mom2) :
  m_hasboost(false), m_hasrotate(false), 
  m_leads(ltp::no), p_left(NULL), p_right(NULL), p_prev(NULL), p_self(NULL),
  m_number(s_cluster_number++)
{
  m_momentum = mom1+mom2;

  //std::cout<<METHOD<<" for "<<m_number<<" & "<<m_momentum
  //	   <<"  ("<<sqrt(m_momentum.Abs2())<<")"<<std::endl;
  s_cluster_count++;
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

  msg_Error()<<"Error in Cluster::Cluster("<<flav1<<","<<flav2<<") :"<<std::endl
	     <<"   Cannot handle this colour structure."<<std::endl
	     <<"   Abort the run."<<std::endl;
  abort();
}

Cluster::~Cluster() 
{
  // std::cout<<METHOD<<" delete : "<<m_number<<" "<<Mass()<<std::endl;
  s_cluster_count--;
}

Blob * Cluster::CHHDecayBlob()
{
  Blob * blob = new Blob();
  control::s_AHAblobs++;
  blob->SetType(btp::Cluster_Decay);
  blob->SetTypeSpec("AHADIC-1.0");
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetId();
  blob->AddToInParticles(p_self);
  p_self->SetStatus(part_status::decayed);
  p_self->ProductionBlob()->UnsetStatus(blob_status::needs_hadrondecays);

  if (p_left && p_left->p_self) {
    blob->AddToOutParticles(p_left->p_self);
    if (p_left->p_self->Flav()!=Flavour(kf::cluster)) {
      delete p_left; p_left=NULL;
    }
  }
  if (p_right && p_right->p_self) {
    blob->AddToOutParticles(p_right->p_self);
    if (p_right->p_self->Flav()!=Flavour(kf::cluster)) {
      delete p_right; p_right=NULL;
    }
  }

  return blob;
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
  msg_Error()<<"ERROR in Cluster::Update():"<<std::endl
	     <<"   Funny flavour/colour structure : "<<m_flavours[0]<<", "<<m_flavours[1]<<std::endl
	     <<"   Will abort."<<std::endl;
  abort();
}

void Cluster::RescaleMomentum(ATOOLS::Vec4D newmom)
{
  Poincare rest(m_momentum);
  Poincare back(newmom);

  Vec4D save[3];
  save[0] = m_momentum;
  save[1] = m_momenta[0];
  save[2] = m_momenta[1];

  rest.Boost(m_momenta[0]);
  back.BoostBack(m_momenta[0]);
  rest.Boost(m_momenta[1]);
  back.BoostBack(m_momenta[1]);
  m_momentum = newmom;

  Vec4D testmom = m_momentum-m_momenta[0]-m_momenta[1];
  if (dabs(testmom.Abs2())>1.e-6 || testmom[0]>1.e-6) {
    msg_Error()<<"Maybe error in RescaleMomentum("<<save[0]<<" -> "<<m_momentum<<")"<<std::endl
	       <<" Was : "<<save[0]<<" ("<<save[0].Abs2()<<")"<<" -> "
	       <<save[1]<<" ("<<save[1].Abs2()<<") + "
	       <<save[2]<<" ("<<save[2].Abs2()<<")"<<std::endl
	       <<" Is  : "<<Momentum(0)<<" ("<<Momentum(0).Abs2()<<") -> "
	       <<Momentum(1)<<" ("<<Momentum(1).Abs2()<<")"
	       <<Momentum(2)<<" ("<<Momentum(2).Abs2()<<")"<<std::endl
	       <<"   diff: "<<testmom;
    rest.Boost(m_momentum); back.BoostBack(m_momentum);
    msg_Error()<<" from "<<newmom<<" = "<<m_momentum<<"."<<std::endl;
  }
}


ATOOLS::Vec4D Cluster::Momentum(const int i) const { 
  if (i==0) return m_momentum;        
  if (i<3)  return m_momenta[i-1];
  msg_Error()<<"Error in Cluster::Momentum("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster momentum."<<std::endl;
  abort();
  return m_momentum;
}

void Cluster::SetMomentum(const int i,const ATOOLS::Vec4D mom) { 
  if (i==0) { m_momentum     = mom; return; }       
  if (i<3)  { m_momenta[i-1] = mom; return; }       
  msg_Error()<<"Error in Cluster::SetMomentum("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster momentum."<<std::endl;
  abort();
}

double Cluster::Mass(const int i) const { 
  if (i==0) return sqrt(m_momentum.Abs2());        
  Flavour flav = m_flavours[i-1];
  if (i<3)  return  hadpars.GetConstituents()->Mass(flav);
  msg_Error()<<"Error in Cluster::Mass("<<i<<")"<<std::endl
		     <<"   Out of bounds. Return cluster mass."<<std::endl;
  abort();
  return sqrt(m_momentum.Abs2());
}

void Cluster::SetFlav(const int i,const ATOOLS::Flavour & flav) { 
  if (i==1) { m_flavours[0] = flav; return; } 
  if (i==2) { m_flavours[1] = flav; return; } 
  msg_Error()<<"Error in Cluster::SetFlav("<<i<<")"<<std::endl
		     <<"   Out of bounds. Do nothing and hope for the best."<<std::endl;  
  abort();
}

Flavour Cluster::GetFlav(const int i) const {
  if (i==1) { return m_flavours[0]; }
  if (i==2) { return m_flavours[1]; }
  msg_Error()<<"Error in Cluster::GetFlav("<<i<<")"<<std::endl
		     <<"   Out of bounds. Do nothing and hope for the best."<<std::endl;  
  abort();
}


void Cluster::BoostInCMSAndRotateOnZ() {
  BoostInCMS();

  m_rotate = ATOOLS::Poincare(m_momenta[0],ATOOLS::Vec4D(1.,ATOOLS::Vec3D::ZVEC));
  ATOOLS::Vec4D copy0(m_momenta[0]), copy1(m_momenta[1]);
  m_rotate.Rotate(copy0);
  m_rotate.Rotate(copy1);
  if (copy0[3]<copy1[3]) {
    m_rotate = ATOOLS::Poincare(m_momenta[0],ATOOLS::Vec4D(1.,(-1.)*ATOOLS::Vec3D::ZVEC));
  }
  m_hasrotate = true;
  Rotate(m_rotate);
}

void Cluster::RotateAndBoostBack() {
  if (!m_hasboost || !m_hasrotate) return;

  RotateBack(m_rotate);
  m_hasrotate = false;
  BoostBack();
  m_hasboost = false;
}

void Cluster::BoostInCMS() {
  if (m_hasboost || m_hasrotate) return;
  m_boost = ATOOLS::Poincare(m_momentum);
  m_boost.Boost(m_momentum);
  for (int i=0;i<2;i++) m_boost.Boost(m_momenta[i]);
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);

  if (p_self) p_self->SetMomentum(m_momentum);
  m_hasboost = true;
}

void Cluster::BoostBack() {
  if (!m_hasboost) return;
  m_boost.BoostBack(m_momentum);
  for (int i=0;i<2;i++) m_boost.BoostBack(m_momenta[i]);
  if (p_left!=NULL)  p_left->BoostBack(m_boost);
  if (p_right!=NULL) p_right->BoostBack(m_boost);

  if (p_self) p_self->SetMomentum(m_momentum);
  m_hasboost = false;
}

void Cluster::Boost(Poincare & boost) {
  boost.Boost(m_momentum);
  for (int i=0;i<2;i++) boost.Boost(m_momenta[i]);
  if (p_left!=NULL)  p_left->Boost(boost);
  if (p_right!=NULL) p_right->Boost(boost);

  if (p_self) p_self->SetMomentum(m_momentum);
}

void Cluster::BoostBack(Poincare & boost) {
  boost.BoostBack(m_momentum);
  for (int i=0;i<2;i++) boost.BoostBack(m_momenta[i]);
  if (p_left!=NULL)  p_left->BoostBack(boost);
  if (p_right!=NULL) p_right->BoostBack(boost);

  if (p_self) p_self->SetMomentum(m_momentum);
}

void Cluster::Rotate(Poincare & rotate) {
  rotate.Rotate(m_momentum);
  for (int i=0;i<2;i++) rotate.Rotate(m_momenta[i]);
  if (p_left!=NULL)  p_left->Rotate(rotate);
  if (p_right!=NULL) p_right->Rotate(rotate);

  if (p_self) p_self->SetMomentum(m_momentum);
}

void Cluster::RotateBack(Poincare & rotate) {
  rotate.RotateBack(m_momentum);
  for (int i=0;i<2;i++) rotate.RotateBack(m_momenta[i]);
  if (p_left!=NULL)  p_left->RotateBack(rotate);
  if (p_right!=NULL) p_right->RotateBack(rotate);

  if (p_self) p_self->SetMomentum(m_momentum);
}

void Cluster::BoostBack(ATOOLS::Vec4D & mom) {
  if (!m_hasboost) return;
  m_boost.BoostBack(mom);
}

void Cluster::RotateAndBoostBack(ATOOLS::Vec4D & mom) {
  if (!m_hasboost || !m_hasrotate) return;
  m_rotate.RotateBack(mom);
  m_boost.BoostBack(mom);
}

std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"-------------------------------------------------------------"<<std::endl
     <<"Cluster ["<<cluster.m_number<<"] ("<<cluster.m_flavours[0]<<" "<<cluster.m_flavours[1]
     <<", "<<cluster.m_momentum<<","<<sqrt(cluster.m_momentum.Abs2())<<" ): "<<std::endl
     <<cluster.m_momenta[0]<<" ("<<sqrt(cluster.m_momenta[0].Abs2())<<"), "
     <<cluster.m_momenta[1]<<" ("<<sqrt(cluster.m_momenta[1].Abs2())<<")";
  if (cluster.m_leads==ltp::leadingtrip) str<<" L_trip";
  if (cluster.m_leads==ltp::leadinganti) str<<" L_anti";
  if (cluster.m_leads==ltp::leadingboth) str<<" L_both";
  str<<std::endl;
  if (cluster.p_left)  str<<"  ---- Left  ----> "<<std::endl<<(*cluster.p_left);
  if (cluster.p_right) str<<"  ---- Right ----> "<<std::endl<<(*cluster.p_right);
  return str;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Cluster_List & cl) {
  s<<"Cluster List with "<<cl.size()<<" elements:"<<std::endl;
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Proto_Particle_List & pl) {
  s<<"Proto_Particle_List with "<<pl.size()<<" elements:"<<std::endl;
  for (PPL_Const_Iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
    s<<"   "<<(*pit)->m_info<<" : "<<(*pit)->m_flav<<" "<<(*pit)->m_mom<<std::endl;
  }
  return s;
}
