#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Blob.H"
#include <algorithm>

using namespace AHADIC;
using namespace ATOOLS;

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
  if (p_trip && p_anti &&
      ((p_trip->m_flav.IsQuark() && !p_trip->m_flav.IsAnti()) || 
       (p_trip->m_flav.IsDiQuark() && p_trip->m_flav.IsAnti())) &&
      ((p_anti->m_flav.IsQuark() && p_anti->m_flav.IsAnti()) || 
       (p_anti->m_flav.IsDiQuark() && !p_anti->m_flav.IsAnti()))) return;

  msg_Error()<<"Error in Cluster::Cluster"
	     <<"("<<p_trip->m_flav<<","<<p_anti->m_flav<<") :\n"
	     <<"   Cannot handle this colour structure, will ignore it."
	     <<std::endl;
}

Cluster::~Cluster() 
{
#ifdef memchecker
  std::cout<<"*** delete Cluster: "<<m_number<<" with "
	   <<" pps: ("<<p_trip<<"/"<<p_anti<<").\n";
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

  msg_Error()
    <<"Error in Cluster::Cluster("<<p_trip->m_flav<<","<<p_anti->m_flav<<"):\n"
    <<"   Cannot handle this colour structure, will abort the run.\n"
    <<"   Please contact the Sherpa group for further assistance.";
  exit(0);
}

bool Cluster::CheckConsistency(std::ostream & s,std::string method) {
  bool passed(dabs(Mass2()-Momentum().Abs2())<1.e-8);
  if (p_trip)  passed = passed && p_trip->CheckConsistency(s,method);
  if (p_anti)  passed = passed && p_anti->CheckConsistency(s,method);
  if (!passed) {
    s<<"Error in "<<METHOD<<" called by "<<method<<":\n"
     <<"   Masses and momenta not consistent for cluster "<<m_number<<": "
     <<Mass2()<<" vs. "<<Momentum()<<" ("<<Momentum().Abs2()<<")\n";
  }
  if (!m_clusters.empty()) {
    Vec4D check(Momentum());
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) {
      passed = passed && (*cit)->CheckConsistency(s,method);
      check -= (*cit)->Momentum();
    }
    if (!IsZero(check.Abs2()) || !IsZero(check[0]/1.e6)) {
      s<<"Error in "<<METHOD<<" called by "<<method<<":\n"
       <<"   Four-momentum not conserved: "<<check<<" ("<<check.Abs2()<<") "
       <<"for "<<Momentum()<<"  ---> \n"
       <<"   "<<p_left->Momentum()<<" + "<<p_right->Momentum()<<".\n";
    }
  }
  return passed;   
}

Particle * Cluster::GetSelf() const { 
  Particle * part(new Particle(-1,size()==1?
			       m_decayproducts[0]:m_flav,m_momentum));
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

void Cluster::RescaleMomentum(Vec4D newmom)
{
  Poincare rest(m_momentum);
  Poincare back(newmom);

  Vec4D_Vector save(3+m_clusters.size());
  save[0] = m_momentum;
  if (p_trip!=NULL)   save[1] = p_trip->m_mom;
  if (p_anti!=NULL)   save[2] = p_anti->m_mom;
  if (p_trip!=NULL)   rest.Boost(p_trip->m_mom);
  if (p_trip!=NULL)   back.BoostBack(p_trip->m_mom);
  if (p_anti!=NULL)   rest.Boost(p_anti->m_mom);
  if (p_anti!=NULL)   back.BoostBack(p_anti->m_mom);
  size_t i(3);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) {
      save[i++] = (*cit)->Momentum();
      (*cit)->Boost(rest);
      (*cit)->BoostBack(back);
    }
  }
  m_momentum = newmom;
  Vec4D testmom(m_momentum-p_trip->m_mom-p_anti->m_mom);
  if (dabs(testmom.Abs2()/save[0][0])>1.e-6 || testmom[0]/save[0][0]>1.e-6) {
    msg_Debugging()<<"Error in "<<METHOD<<":\n"
	       <<"   From "<<save[0]<<" ("
		   <<sqrt(Max(0.,save[0].Abs2()))<<") to "
	       <<m_momentum<<" with \n"
	       <<"   "<<save[1]<<" ("<<save[1].Abs2()<<") + "
	       <<save[2]<<" ("<<save[2].Abs2()<<")\n";
    if (p_trip!=NULL) {
      msg_Debugging()<<"  Trip: "<<p_trip->m_mom
		     <<" ("<<p_trip->m_mom.Abs2()<<")";
    }
    else { 
      msg_Debugging()<<"No triplet: "<<p_trip<<" ";
    }
    if (p_anti!=NULL) {
      msg_Debugging()<<" Anti: "<<p_anti->m_mom
		     <<" ("<<p_anti->m_mom.Abs2()<<")\n";
    }
    else { 
      msg_Debugging()<<"No antitriplet: "<<p_anti<<" ";
    }
    msg_Debugging()<<"   diff: "<<testmom;
    rest.Boost(m_momentum); 
    back.BoostBack(m_momentum);
    DEBUG_VAR(m_momentum);
    msg_Debugging()<<" from "<<newmom<<" --> "<<m_momentum<<".\n";
  }
  if (p_left!=NULL) {
    msg_Error()<<"Maybe error in RescaleMomentum("
	       <<save[0]<<" -> "<<m_momentum<<")\n"
	       <<"   How about the left/right offsprings?\n";
  }
}

void Cluster::BoostInCMSAndRotateOnZ() {
  if (p_trip==NULL) return;
  BoostInCMS();

  m_rotate = Poincare(p_trip->m_mom,
			      Vec4D(1.,Vec3D::ZVEC));
  Vec4D copy0(p_trip->m_mom), copy1(p_anti->m_mom);
  m_rotate.Rotate(copy0);
  m_rotate.Rotate(copy1);
  if (copy0[3]<copy1[3]) {
    m_rotate = Poincare(p_trip->m_mom,Vec4D(1.,(-1.)*Vec3D::ZVEC));
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
  m_boost = Poincare(m_momentum);
  m_boost.Boost(m_momentum);
  if (p_trip!=NULL) m_boost.Boost(p_trip->m_mom);
  if (p_anti!=NULL) m_boost.Boost(p_anti->m_mom);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      (*cit)->Boost(m_boost);
  }
  m_hasboost = true;
}

void Cluster::BoostBack() {
  if (!m_hasboost) return;
  m_boost.BoostBack(m_momentum);
  if (p_trip!=NULL) m_boost.BoostBack(p_trip->m_mom);
  if (p_anti!=NULL) m_boost.BoostBack(p_anti->m_mom);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      (*cit)->BoostBack(m_boost);
  }
  m_hasboost = false;
}

void Cluster::Boost(Poincare & boost) {
  boost.Boost(m_momentum);
  if (p_trip!=NULL) boost.Boost(p_trip->m_mom);
  if (p_anti!=NULL) boost.Boost(p_anti->m_mom);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      (*cit)->Boost(boost);
  }
}

void Cluster::BoostBack(Poincare & boost) {
  boost.BoostBack(m_momentum);
  if (p_trip!=NULL) boost.BoostBack(p_trip->m_mom);
  if (p_anti!=NULL) boost.BoostBack(p_anti->m_mom);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      (*cit)->BoostBack(boost);
  }
}

void Cluster::Rotate(Poincare & rotate) {
  rotate.Rotate(m_momentum);
  if (p_trip!=NULL) rotate.Rotate(p_trip->m_mom);
  if (p_anti!=NULL) rotate.Rotate(p_anti->m_mom);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      (*cit)->Rotate(rotate);
  }
}

void Cluster::RotateBack(Poincare & rotate) {
  rotate.RotateBack(m_momentum);
  if (p_trip!=NULL) rotate.RotateBack(p_trip->m_mom);
  if (p_anti!=NULL) rotate.RotateBack(p_anti->m_mom);
  if (!m_clusters.empty()) {
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      (*cit)->RotateBack(rotate);
  }
}

void Cluster::BoostBack(Vec4D & mom) {
  if (!m_hasboost) return;
  m_boost.BoostBack(mom);
}

void Cluster::RotateAndBoostBack(Vec4D & mom) {
  if (!m_hasboost || !m_hasrotate) return;
  m_rotate.RotateBack(mom);
  m_boost.BoostBack(mom);
}

void Cluster::Print() {
  msg_Out()<<"   Cluster [active = "<<m_active<<", number = "
	   <<m_number<<", size = "<<size()<<"], "
	   <<"constituents = "<<p_trip->m_flav<<" & "
	   <<p_anti->m_flav<<std::endl
	   <<"      flavour = "<<m_flav<<" with "
	   <<m_momentum<<"), "<<m_momentum.Abs2()<<" ---> ";
  if (m_decayproducts.size()>0) {
    for (size_t i=0;i<m_decayproducts.size();i++) 
      msg_Out()<<m_decayproducts[i]<<" ";
    msg_Out()<<".\n";
    return;
  }
  if (!m_clusters.empty()) {
    msg_Out()<<" ("<<m_clusters.size()<<"): { ";
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      msg_Out()<<(*cit)->m_number<<" ";
    msg_Out()<<"}\n";
    for (Cluster_Iterator cit=m_clusters.begin();cit!=m_clusters.end();cit++) 
      msg_Out()<<(**cit)<<"\n";
  }
}

void Cluster::Delete() {
  while (!m_clusters.empty()) {
    delete (m_clusters.front());
    m_clusters.pop_front();
  }
}

bool Cluster::EnsureMomentum() {
  Vec4D check(Momentum());
  for (Cluster_Iterator cit(m_clusters.begin());cit!=m_clusters.end();cit++)
    check -= (*cit)->Momentum();
  if (dabs(check.Abs2())>1.e-10 || dabs(check[0])>1.e-10 || 
      dabs(check.PSpat())>1.e-10) {
    // msg_Out()<<"*** "<<METHOD<<" yields "<<check<<" ("<<check.Abs2()<<") "
    // 	     <<"for cluster "<<Number()<<"\n"<<(*this)
    // 	     <<"********************************************************\n";
    return false;
  }
  return true;
}


std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"-------------------------------------------------------------\n"
     <<"Cluster ["<<cluster.m_flav<<", "<<cluster.m_number<<", "
     <<cluster.size()<<"] "
     <<"("<<cluster.m_momentum<<", "
     <<"mass = "<<sqrt(cluster.m_momentum.Abs2())<<", "
     <<"y = "<<cluster.m_momentum.Y()<<") ";
  if (cluster.p_nbtrip!=NULL) str<<" [> "<<cluster.p_nbtrip->Number()<<"] ";
  if (cluster.p_nbanti!=NULL) str<<" [< "<<cluster.p_nbanti->Number()<<"] ";
  str<<":\n";
  if (cluster.p_trip!=NULL) str<<"  "<<(*cluster.p_trip);
  if (cluster.p_anti!=NULL) str<<"  "<<(*cluster.p_anti);
  if (!cluster.m_clusters.empty()) {
    msg_Out()<<" ("<<cluster.m_clusters.size()<<"): { ";
    for (Cluster_Const_Iterator cit=cluster.m_clusters.begin();
	 cit!=cluster.m_clusters.end();cit++) 
      msg_Out()<<(*cit)->m_number<<" ";
    msg_Out()<<"}\n";
  }
  else msg_Out()<<"\n";
  return str;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Cluster_List & cl) {
  Vec4D totmom(0.,0.,0.,0.);
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) 
    totmom += (*cit)->Momentum();
  s<<"Cluster List with "<<cl.size()<<" elements, mom = "<<totmom<<":\n";
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

