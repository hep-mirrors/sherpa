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


Proto_Particle::Proto_Particle(const Proto_Particle & pp) :
  m_flav(pp.m_flav), m_mom(pp.m_mom), m_info(pp.m_info), m_mass(pp.m_mass)
{ 
  control::s_AHAprotoparticles++; 
}

Proto_Particle::Proto_Particle(ATOOLS::Flavour flav,ATOOLS::Vec4D mom,char info) :
  m_flav(flav), m_mom(mom), m_info(info), 
  m_mass(hadpars.GetConstituents()->Mass(flav)) 
{ 
  control::s_AHAprotoparticles++; 
}


Proto_Particle::~Proto_Particle()
{ 
  control::s_AHAprotoparticles--; 
}


std::ostream & AHADIC::operator<<(std::ostream & s, const Proto_Particle& proto) {
  s<<"   "<<proto.m_info<<" : "<<proto.m_flav<<" "<<proto.m_mom
   <<" "<<sqrt(ATOOLS::Max(0.,proto.m_mom.Abs2()))<<std::endl;
  return s;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Proto_Particle_List & pl) {
  s<<"Proto_Particle_List with "<<pl.size()<<" elements:"<<std::endl;
  for (PPL_Const_Iterator pit=pl.begin(); pit!=pl.end(); ++pit) s<<(**pit)<<std::endl;
  return s;
}

Cluster::Cluster(Vec4D mom,Flavour flav,bool active) :
  m_active(active), p_trip(NULL), p_anti(NULL), 
  m_momentum(mom), m_flav(flav),
  m_hasboost(false), m_hasrotate(false), 
  p_left(NULL), p_right(NULL), p_prev(NULL),
  m_number(s_cluster_number++)
{
  s_cluster_count++;
}

Cluster::Cluster(Proto_Particle * trip,Proto_Particle * anti) :
  m_active(true), p_trip(trip), p_anti(anti), 
  m_momentum(p_trip->m_mom+p_anti->m_mom), m_flav(Flavour(kf_cluster)),
  m_hasboost(false), m_hasrotate(false), 
  p_left(NULL), p_right(NULL), p_prev(NULL),
  m_number(s_cluster_number++)
{
  ////PRINT_VAR(m_momentum);
  s_cluster_count++;
  if (((p_trip->m_flav.IsQuark() && !p_trip->m_flav.IsAnti()) || 
       (p_trip->m_flav.IsDiQuark() && p_trip->m_flav.IsAnti())) &&
      ((p_anti->m_flav.IsQuark() && p_anti->m_flav.IsAnti()) || 
       (p_anti->m_flav.IsDiQuark() && !p_anti->m_flav.IsAnti()))) return;

  msg_Error()<<"Error in Cluster::Cluster("<<p_trip->m_flav<<","<<p_anti->m_flav<<") :"<<std::endl
	     <<"   Cannot handle this colour structure, will abort the run."<<std::endl
	     <<"   Please contact the Sherpa group for further assistance."<<std::endl;
  abort();
}

Cluster::~Cluster() 
{
  if (p_trip) { 
#ifdef memchecker
    std::cout<<"### Delete Proto_Particle ("<<p_trip
    	     <<"/"<<p_trip->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
    delete p_trip; p_trip=NULL; 
  }
  if (p_anti) { 
#ifdef memchecker
    std::cout<<"### Delete Proto_Particle ("<<p_anti
    	     <<"/"<<p_anti->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
    delete p_anti; p_anti=NULL; 
  }
  s_cluster_count--;
}

void Cluster::Update()
{
  m_momentum = p_trip->m_mom + p_anti->m_mom;
  //PRINT_VAR(m_momentum);
  if (!p_trip && !p_anti) return;
  if (((p_trip->m_flav.IsQuark() && !p_trip->m_flav.IsAnti()) || 
       (p_trip->m_flav.IsDiQuark() && p_trip->m_flav.IsAnti())) &&
      ((p_anti->m_flav.IsQuark() && p_anti->m_flav.IsAnti()) || 
       (p_anti->m_flav.IsDiQuark() && !p_anti->m_flav.IsAnti()))) return;

  msg_Error()<<"Error in Cluster::Cluster("<<p_trip->m_flav<<","<<p_anti->m_flav<<") :"<<std::endl
	     <<"   Cannot handle this colour structure, will abort the run."<<std::endl
	     <<"   Please contact the Sherpa group for further assistance."<<std::endl;
  exit(0);
}

Particle * Cluster::GetSelf() const { 
  Particle * part(new Particle(-1,m_flav,m_momentum));
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(m_flav.PSMass());
  control::s_AHAparticles++;
  return part;
}

Blob * Cluster::ConstructDecayBlob()
{
  Blob * blob = new Blob();
  control::s_AHAblobs++;
  blob->SetType(btp::Cluster_Decay);
  blob->SetTypeSpec("AHADIC-1.0");
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetId();
  Particle * part(GetSelf());
  blob->AddToInParticles(part);
  part->SetStatus(part_status::decayed);
  part->ProductionBlob()->UnsetStatus(blob_status::needs_hadrondecays);

  if (p_left) {
    part = p_left->GetSelf();
    blob->AddToOutParticles(part);
    if (part->Flav()!=Flavour(kf_cluster)) p_left->m_active=false;
  }
  if (p_right) {
    part = p_right->GetSelf();
    blob->AddToOutParticles(part);
    if (part->Flav()!=Flavour(kf_cluster)) p_right->m_active=false;
  }

  return blob;
}

void Cluster::RescaleMomentum(ATOOLS::Vec4D newmom)
{
  Poincare rest(m_momentum);
  Poincare back(newmom);

  Vec4D save[3+int(p_left!=NULL)+int(p_right!=NULL)];
  save[0] = m_momentum;
  if (p_trip)  save[1] = p_trip->m_mom;
  if (p_anti)  save[2] = p_anti->m_mom;
  if (p_left)  save[3] = p_trip->m_mom;
  if (p_right) save[4] = p_anti->m_mom;

  if (p_trip)  rest.Boost(p_trip->m_mom);
  if (p_trip)  back.BoostBack(p_trip->m_mom);
  if (p_anti)  rest.Boost(p_anti->m_mom);
  if (p_anti)  back.BoostBack(p_anti->m_mom);
  if (p_left)  p_left->Boost(rest);
  if (p_left)  p_left->BoostBack(back);
  if (p_right)  p_right->Boost(rest);
  if (p_right)  p_right->BoostBack(back);
  m_momentum = newmom;
  //PRINT_VAR(m_momentum);
  //PRINT_VAR(newmom);
  //PRINT_VAR(p_trip->m_mom);
  //PRINT_VAR(p_anti->m_mom);

  Vec4D testmom = m_momentum-p_trip->m_mom-p_anti->m_mom;
  if (dabs(testmom.Abs2())>1.e-6 || testmom[0]>1.e-6) {
    msg_Error()<<"Maybe error in RescaleMomentum("<<save[0]<<" -> "<<m_momentum<<")"<<std::endl
	       <<" Was : "<<save[0]<<" ("<<save[0].Abs2()<<")"<<" -> "
	       <<save[1]<<" ("<<save[1].Abs2()<<") + "
	       <<save[2]<<" ("<<save[2].Abs2()<<")"<<std::endl
	       <<" Is  : "<<m_momentum<<" ("<<m_momentum.Abs2()<<") -> ";
    if (p_trip) msg_Error()<<p_trip->m_mom<<" ("<<p_trip->m_mom.Abs2()<<")";
    else msg_Error()<<"No triplet: "<<p_trip<<" ";
    if (p_anti) msg_Error()<<p_anti->m_mom<<" ("<<p_anti->m_mom.Abs2()<<")"<<std::endl;
    else msg_Error()<<"No antitriplet: "<<p_anti<<" ";
    msg_Error()<<"   diff: "<<testmom;
    rest.Boost(m_momentum); back.BoostBack(m_momentum);
    //PRINT_VAR(m_momentum);
    msg_Error()<<" from "<<newmom<<" = "<<m_momentum<<"."<<std::endl;
  }
  if (p_left!=NULL) {
    msg_Error()<<"Maybe error in RescaleMomentum("<<save[0]<<" -> "<<m_momentum<<")"<<std::endl
	       <<"   How about the left/right offsprings?"<<std::endl;
  }
}

void Cluster::BoostInCMSAndRotateOnZ() {
  //PRINT_INFO(*this);
  if (!p_trip) return;
  BoostInCMS();
  //PRINT_INFO(*this);

  m_rotate = ATOOLS::Poincare(p_trip->m_mom,ATOOLS::Vec4D(1.,ATOOLS::Vec3D::ZVEC));
  ATOOLS::Vec4D copy0(p_trip->m_mom), copy1(p_anti->m_mom);
  m_rotate.Rotate(copy0);
  m_rotate.Rotate(copy1);
  //PRINT_INFO(*this);
  if (copy0[3]<copy1[3]) {
    m_rotate = ATOOLS::Poincare(p_trip->m_mom,ATOOLS::Vec4D(1.,(-1.)*ATOOLS::Vec3D::ZVEC));
  }
  m_hasrotate = true;
  Rotate(m_rotate);
  //PRINT_INFO(*this);
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
  //PRINT_INFO(this<<" "<<m_momentum);
  m_boost = ATOOLS::Poincare(m_momentum);
  m_boost.Boost(m_momentum);
  //PRINT_INFO(this<<" "<<m_momentum);
  if (p_trip) m_boost.Boost(p_trip->m_mom);
  if (p_anti) m_boost.Boost(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);

  m_hasboost = true;
}

void Cluster::BoostBack() {
  if (!m_hasboost) return;
  m_boost.BoostBack(m_momentum);
  if (p_trip) m_boost.BoostBack(p_trip->m_mom);
  if (p_anti) m_boost.BoostBack(p_anti->m_mom);
  if (p_left!=NULL)  p_left->BoostBack(m_boost);
  if (p_right!=NULL) p_right->BoostBack(m_boost);

  m_hasboost = false;
}

void Cluster::Boost(Poincare & boost) {
  boost.Boost(m_momentum);
  if (p_trip) boost.Boost(p_trip->m_mom);
  if (p_anti) boost.Boost(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Boost(boost);
  if (p_right!=NULL) p_right->Boost(boost);
}

void Cluster::BoostBack(Poincare & boost) {
  boost.BoostBack(m_momentum);
  if (p_trip) boost.BoostBack(p_trip->m_mom);
  if (p_anti) boost.BoostBack(p_anti->m_mom);
  if (p_left!=NULL)  p_left->BoostBack(boost);
  if (p_right!=NULL) p_right->BoostBack(boost);
}

void Cluster::Rotate(Poincare & rotate) {
  rotate.Rotate(m_momentum);
  if (p_trip) rotate.Rotate(p_trip->m_mom);
  if (p_anti) rotate.Rotate(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Rotate(rotate);
  if (p_right!=NULL) p_right->Rotate(rotate);
}

void Cluster::RotateBack(Poincare & rotate) {
  rotate.RotateBack(m_momentum);
  if (p_trip) rotate.RotateBack(p_trip->m_mom);
  if (p_anti) rotate.RotateBack(p_anti->m_mom);
  if (p_left!=NULL)  p_left->RotateBack(rotate);
  if (p_right!=NULL) p_right->RotateBack(rotate);
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

void Cluster::Print() {
  msg_Out()<<"   Cluster ["<<m_active<<", "<<m_number<<", "<<size()<<"] "
	   <<"("<<m_flav<<" with "<<m_momentum<<") ---> ";
  if (p_left)  msg_Out()<<" (l) "<<p_left->m_number<<" ";
  if (p_right) msg_Out()<<" (r) "<<p_right->m_number;
  msg_Out()<<std::endl;
  if (p_left) p_left->Print(); 
  if (p_right) p_right->Print(); 
}

void Cluster::Delete() {
  if (p_left)  p_left->Delete();
  if (p_right) p_right->Delete();

  delete this;
}

std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"-------------------------------------------------------------"<<std::endl
     <<"Cluster ["<<cluster.m_flav<<", "<<cluster.m_number<<", "<<cluster.size()<<"] "
     <<"("<<cluster.m_momentum<<","<<sqrt(cluster.m_momentum.Abs2())<<" ): "<<std::endl;
  if (cluster.p_trip) str<<"  "<<(*cluster.p_trip);
  if (cluster.p_anti) str<<"  "<<(*cluster.p_anti);
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

