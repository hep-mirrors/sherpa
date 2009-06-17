#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include <algorithm>

using namespace AHADIC;
using namespace ATOOLS;

namespace AHADIC {
  long int Cluster::s_cluster_count=0;
  long int Cluster::s_cluster_number=0;

  long int control::s_AHAparticles=0;
  long int control::s_AHAprotoparticles=0;
  long int control::s_AHAblobs=0;

  std::list<Proto_Particle *>      Proto_Particle::s_actives;
  std::list<Proto_Particle_List *> Proto_Particle_List::s_actives;
  std::list<ListOfPPLs *>          ListOfPPLs::s_actives;
  std::list<Cluster *>             Cluster::s_actives;
  std::list<Cluster_List *>        Cluster_List::s_actives;
}

Proto_Particle::Proto_Particle(const Proto_Particle & pp) :
  m_flav(pp.m_flav), m_mom(pp.m_mom), m_info(pp.m_info), 
  m_mass(pp.m_mass), m_kt2max(pp.m_kt2max),
  p_partner(pp.p_partner)
{ 
  control::s_AHAprotoparticles++; 
  s_actives.push_back(this);
}

Proto_Particle::Proto_Particle(ATOOLS::Flavour flav,ATOOLS::Vec4D mom,char info) :
  m_flav(flav), m_mom(mom), m_info(info), 
  m_mass(hadpars.GetConstituents()->Mass(flav)), m_kt2max(0.), 
  p_partner(NULL)
{ 
  control::s_AHAprotoparticles++; 
  s_actives.push_back(this);
}


Proto_Particle::~Proto_Particle()
{ 
#ifdef memchecker
  std::cout<<"### delete Proto_Particle: ("<<m_flav<<"/"<<this<<")."<<std::endl;
#endif
  control::s_AHAprotoparticles--; 
  s_actives.remove(this);
}

bool Proto_Particle::CheckConsistency(std::ostream & s,std::string method) {
  if (dabs(m_mass-hadpars.GetConstituents()->Mass(m_flav))>1.e-6 ||
      dabs(m_mass-sqrt(m_mom.Abs2()))>1.e-6 ||
      dabs(sqrt(m_mom.Abs2())-hadpars.GetConstituents()->Mass(m_flav))>1.e-6) {
    s<<"Error in "<<METHOD<<" called by "<<method<<":"<<std::endl
	 <<"   Masses and momenta not consistent for "<<m_flav<<"("<<m_mass<<"),"
     <<" sqrt(mom^2) = "<<sqrt(m_mom.Abs2())
     <<" & constituent mass = "<<hadpars.GetConstituents()->Mass(m_flav)<<"."<<std::endl;
    return false;
  }
  return true;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Proto_Particle& proto) {
  s<<"   "<<proto.m_info<<" : "<<proto.m_flav<<" "<<proto.m_mom
   <<" "<<sqrt(ATOOLS::Max(0.,proto.m_mom.Abs2()))
   <<", kt2max = "<<sqrt(ATOOLS::Max(0.,proto.m_kt2max))<<std::endl;
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
  p_left(NULL), p_right(NULL), p_prev(NULL), p_nbtrip(NULL), p_nbanti(NULL),
  m_number(++s_cluster_number)
{
  s_cluster_count++;
  s_actives.push_back(this);
}

Cluster::Cluster(Proto_Particle * trip,Proto_Particle * anti) :
  m_active(true), p_trip(trip), p_anti(anti), 
  m_momentum(p_trip->m_mom+p_anti->m_mom), m_flav(Flavour(kf_cluster)),
  m_hasboost(false), m_hasrotate(false), 
  p_left(NULL), p_right(NULL), p_prev(NULL), p_nbtrip(NULL), p_nbanti(NULL),
  m_number(++s_cluster_number)
{
  ////PRINT_VAR(m_momentum);
  s_cluster_count++;
  s_actives.push_back(this);
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
#ifdef memchecker
  std::cout<<"*** delete Cluster: "<<m_number<<" with "
	   <<" pps: ("<<p_trip<<"/"<<p_anti<<")."<<std::endl;
#endif
  s_cluster_count--;
  s_actives.remove(this);
}

void Cluster::Update()
{
  m_momentum = p_trip->m_mom + p_anti->m_mom;
  if (p_trip==NULL && p_anti==NULL) return;
  if (((p_trip->m_flav.IsQuark() && !p_trip->m_flav.IsAnti()) || 
       (p_trip->m_flav.IsDiQuark() && p_trip->m_flav.IsAnti())) &&
      ((p_anti->m_flav.IsQuark() && p_anti->m_flav.IsAnti()) || 
       (p_anti->m_flav.IsDiQuark() && !p_anti->m_flav.IsAnti()))) return;

  msg_Error()<<"Error in Cluster::Cluster("<<p_trip->m_flav<<","<<p_anti->m_flav<<") :"<<std::endl
	     <<"   Cannot handle this colour structure, will abort the run."<<std::endl
	     <<"   Please contact the Sherpa group for further assistance."<<std::endl;
  exit(0);
}

bool Cluster::CheckConsistency(std::ostream & s,std::string method) {
  bool passed(dabs(Mass2()-Momentum().Abs2())<1.e-8);
  if (p_trip)  passed = passed && p_trip->CheckConsistency(s,method);
  if (p_anti)  passed = passed && p_anti->CheckConsistency(s,method);
  if (!passed) {
    s<<"Error in "<<METHOD<<" called by "<<method<<":"<<std::endl
     <<"   Masses and momenta not consistent for cluster "<<m_number<<": "
     <<Mass2()<<" vs. "<<Momentum()<<" ("<<Momentum().Abs2()<<")"<<std::endl;
  }
  if (p_left)  passed = passed && p_left->CheckConsistency(s,method);
  if (p_right) passed = passed && p_right->CheckConsistency(s,method);
  if (p_left && p_right) {
    Vec4D check(Momentum()-p_left->Momentum()-p_right->Momentum());
    if (!IsZero(check.Abs2()) || !IsZero(check[0]/1.e6)) {
      s<<"Error in "<<METHOD<<" called by "<<method<<":"<<std::endl
       <<"   Four-momentum not conserved: "<<check<<" ("<<check.Abs2()<<") for "<<std::endl
       <<"   "<<Momentum()<<"  ---> "<<std::endl
       <<"   "<<p_left->Momentum()<<" + "<<p_right->Momentum()<<"."<<std::endl;
    }
  }
  return passed;   
}

Particle * Cluster::GetSelf() const { 
  Particle * part(new Particle(-1,size()==1?m_decayproducts[0]:m_flav,m_momentum));
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(m_flav.HadMass());
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

  if (p_left!=NULL) {
    part = p_left->GetSelf();
    blob->AddToOutParticles(part);
    if (part->Flav()!=Flavour(kf_cluster)) p_left->m_active=false;
  }
  if (p_right!=NULL) {
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

  Vec4D_Vector save(3+int(p_left!=NULL)+int(p_right!=NULL));
  save[0] = m_momentum;
  if (p_trip!=NULL)  save[1] = p_trip->m_mom;
  if (p_anti!=NULL)  save[2] = p_anti->m_mom;
  if (p_left!=NULL)  save[3] = p_trip->m_mom;
  if (p_right!=NULL) save[4] = p_anti->m_mom;

  if (p_trip!=NULL)  rest.Boost(p_trip->m_mom);
  if (p_trip!=NULL)  back.BoostBack(p_trip->m_mom);
  if (p_anti!=NULL)  rest.Boost(p_anti->m_mom);
  if (p_anti!=NULL)  back.BoostBack(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Boost(rest);
  if (p_left!=NULL)  p_left->BoostBack(back);
  if (p_right!=NULL)  p_right->Boost(rest);
  if (p_right!=NULL)  p_right->BoostBack(back);
  m_momentum = newmom;
  //PRINT_VAR(m_momentum);
  //PRINT_VAR(newmom);
  //PRINT_VAR(p_trip->m_mom);
  //PRINT_VAR(p_anti->m_mom);

  Vec4D testmom = m_momentum-p_trip->m_mom-p_anti->m_mom;
  if (dabs(testmom.Abs2()/save[0][0])>1.e-6 || testmom[0]/save[0][0]>1.e-6) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   From "<<save[0]<<"("<<sqrt(Max(0.,save[0].Abs2()))<<") to "<<m_momentum<<" with "<<std::endl
	       <<"   "<<save[1]<<" ("<<save[1].Abs2()<<") + "
	       <<save[2]<<" ("<<save[2].Abs2()<<")"<<std::endl;
    if (p_trip!=NULL) msg_Error()<<p_trip->m_mom<<" ("<<p_trip->m_mom.Abs2()<<")";
           else msg_Error()<<"No triplet: "<<p_trip<<" ";
    if (p_anti!=NULL) msg_Error()<<p_anti->m_mom<<" ("<<p_anti->m_mom.Abs2()<<")"<<std::endl;
           else msg_Error()<<"No antitriplet: "<<p_anti<<" ";
    msg_Error()<<"   diff: "<<testmom;
    rest.Boost(m_momentum); back.BoostBack(m_momentum);
    //PRINT_VAR(m_momentum);
    msg_Error()<<" from "<<newmom<<" = "<<m_momentum<<"."<<std::endl;
    exit(1);
  }
  if (p_left!=NULL) {
    msg_Error()<<"Maybe error in RescaleMomentum("<<save[0]<<" -> "<<m_momentum<<")"<<std::endl
	       <<"   How about the left/right offsprings?"<<std::endl;
  }
}

void Cluster::BoostInCMSAndRotateOnZ() {
  //PRINT_INFO(*this);
  if (p_trip==NULL) return;
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
  if (p_trip!=NULL) m_boost.Boost(p_trip->m_mom);
  if (p_anti!=NULL) m_boost.Boost(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Boost(m_boost);
  if (p_right!=NULL) p_right->Boost(m_boost);

  m_hasboost = true;
}

void Cluster::BoostBack() {
  if (!m_hasboost) return;
  m_boost.BoostBack(m_momentum);
  if (p_trip!=NULL) m_boost.BoostBack(p_trip->m_mom);
  if (p_anti!=NULL) m_boost.BoostBack(p_anti->m_mom);
  if (p_left!=NULL)  p_left->BoostBack(m_boost);
  if (p_right!=NULL) p_right->BoostBack(m_boost);

  m_hasboost = false;
}

void Cluster::Boost(Poincare & boost) {
  boost.Boost(m_momentum);
  if (p_trip!=NULL) boost.Boost(p_trip->m_mom);
  if (p_anti!=NULL) boost.Boost(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Boost(boost);
  if (p_right!=NULL) p_right->Boost(boost);
}

void Cluster::BoostBack(Poincare & boost) {
  boost.BoostBack(m_momentum);
  if (p_trip!=NULL) boost.BoostBack(p_trip->m_mom);
  if (p_anti!=NULL) boost.BoostBack(p_anti->m_mom);
  if (p_left!=NULL)  p_left->BoostBack(boost);
  if (p_right!=NULL) p_right->BoostBack(boost);
}

void Cluster::Rotate(Poincare & rotate) {
  rotate.Rotate(m_momentum);
  if (p_trip!=NULL) rotate.Rotate(p_trip->m_mom);
  if (p_anti!=NULL) rotate.Rotate(p_anti->m_mom);
  if (p_left!=NULL)  p_left->Rotate(rotate);
  if (p_right!=NULL) p_right->Rotate(rotate);
}

void Cluster::RotateBack(Poincare & rotate) {
  rotate.RotateBack(m_momentum);
  if (p_trip!=NULL) rotate.RotateBack(p_trip->m_mom);
  if (p_anti!=NULL) rotate.RotateBack(p_anti->m_mom);
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
  msg_Out()<<"   Cluster [active = "<<m_active<<", number = "<<m_number<<", size = "<<size()<<"], "
	   <<"constituents = "<<p_trip->m_flav<<" & "<<p_anti->m_flav<<std::endl
	   <<"      flavour = "<<m_flav<<" with "<<m_momentum<<"), "<<m_momentum.Abs2()<<" ---> ";
  if (m_decayproducts.size()>0) {
    for (size_t i=0;i<m_decayproducts.size();i++) 
      msg_Out()<<m_decayproducts[i]<<" ";
    msg_Out()<<"."<<std::endl;
    return;
  }
  if (p_left!=NULL)  msg_Out()<<" (l) "<<p_left->m_number<<" ";
  if (p_right!=NULL) msg_Out()<<" (r) "<<p_right->m_number;
  msg_Out()<<std::endl;
  if (p_left!=NULL) p_left->Print(); 
  if (p_right!=NULL) p_right->Print(); 
}

void Cluster::Delete() {
  if (p_left!=NULL)  p_left->Delete();
  if (p_right!=NULL) p_right->Delete();

}


std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"-------------------------------------------------------------"<<std::endl
     <<"Cluster ["<<cluster.m_flav<<", "<<cluster.m_number<<", "<<cluster.size()<<"] "
     <<"("<<cluster.m_momentum<<","<<sqrt(cluster.m_momentum.Abs2())<<" ) ";
  if (cluster.p_nbtrip!=NULL) str<<" [> "<<cluster.p_nbtrip->Number()<<"] ";
  if (cluster.p_nbanti!=NULL) str<<" [< "<<cluster.p_nbanti->Number()<<"] ";
  str<<":"<<std::endl;
  if (cluster.p_trip!=NULL) str<<"  "<<(*cluster.p_trip);
  if (cluster.p_anti!=NULL) str<<"  "<<(*cluster.p_anti);
  if (cluster.p_left!=NULL)  str<<"  ---- Left  ----> "<<std::endl<<(*cluster.p_left);
  if (cluster.p_right!=NULL) str<<"  ---- Right ----> "<<std::endl<<(*cluster.p_right);
  return str;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Cluster_List & cl) {
  s<<"Cluster List with "<<cl.size()<<" elements:"<<std::endl;
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

