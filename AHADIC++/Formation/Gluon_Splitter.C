#include "AHADIC++/Formation/Gluon_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

bool Gluon_Splitter::MakeLongitudinalMomenta() {
  //msg_Out()<<"   "<<METHOD<<" called.\n";
  CalculateLimits();
  m_z1 = m_zselector(m_z1min,m_z1max,p_part1->Flavour());
  m_z2 = m_zselector(m_z2min,m_z2max,p_part2->Flavour());
  double R2 = 0.;
  if (CalculateXY()) {
    R2 = (1.-m_y)*(1.-m_z2)*m_Q2-m_kt2;  
    //msg_Out()<<"   "<<METHOD<<": "<<sqrt(R2)<<" > "<<m_minQ_1<<".\n";
    //if (R2<0.) {
    //  msg_Out()<<"   negative cluster mass: z1 = "<<m_z1<<", z2 = "<<m_z2
    //	       <<" and x = "<<m_x<<" from z1 in ["<<m_z1min<<", "<<m_z1max<<"]"
    //	       <<".\n";
    //}
    //else if (R2>sqr(m_minQ_1)) {
    //  msg_Out()<<"*** works:  z1 = "<<m_z1<<", z2 = "<<m_z2<<".\n";
    //}
  }
  return (R2>sqr(m_minQ_1));
}

void Gluon_Splitter::CalculateLimits() {
  m_z1min   = m_m12/m_Q2;
  m_z1max   = 1.;
  double r2 = (m_popped_mass2+m_kt2)/m_Q2;
  double R2 = (m_m12+m_kt2)/m_Q2;
  double ce = (1. + r2 - R2)/2.;
  double la = Lambda(1,r2,R2);
  m_z2min   = Max(ce - la,r2);
  m_z2max   = ce + la;
  if (m_z2max<m_z2min) abort();
}

bool Gluon_Splitter::CalculateXY() {
  m_x = m_m12/(m_z1*m_Q2);
  m_y = (m_popped_mass2+m_kt2)/(m_z2*m_Q2);
  return (m_x<=1. && m_y<=1.);
}

bool Gluon_Splitter::CheckKinematics() {
  //msg_Out()<<METHOD<<"("<<p_part3<<").\n";
  if ((1.-m_x-m_z2)*(1.-m_y-m_z1)*m_Q2<m_kt2) return false;
  if (p_part3==0) return true;
  Vec4D newmom2 = m_E*(m_y*s_AxisP+m_z2*s_AxisM)-m_ktvec;
  m_rotat.RotateBack(newmom2);
  m_boost.BoostBack(newmom2);
  //msg_Out()<<"   "<<METHOD<<"("<<m_newflav2<<" "<<p_part3->Flavour()<<"):\n"
  //	   <<"   "<<newmom2<<" + "<<p_part3->Momentum()
  //	   <<"   -> "<<(newmom2+p_part3->Momentum()).Abs2()
  //	   <<" vs. "<<sqr(m_minQ_2)<<".\n";
  return ((newmom2+p_part3->Momentum()).Abs2()>sqr(m_minQ_2));
}

void Gluon_Splitter::FillParticlesInLists() {
  Cluster * cluster = MakeCluster();
  if (p_softclusters->Treat(cluster)) delete cluster;
  else p_cluster_list->push_back(cluster);
  UpdateSpectator();
}

void Gluon_Splitter::UpdateSpectator() {
  // Replace splitted gluon with (anti-)(di-)quark and correct momentum
  Vec4D newmom2 = m_E*(m_y*s_AxisP+m_z2*s_AxisM)-m_ktvec;
  m_rotat.RotateBack(newmom2);
  m_boost.BoostBack(newmom2);
  p_part2->SetFlavour(m_newflav2);
  p_part2->SetMomentum(newmom2);
}

Cluster * Gluon_Splitter::MakeCluster() {
  // If kinematically allowed, i.e. if the emerging cluster is heavy enough,
  // we calculate the split of cluster momentum into momenta for its
  // constituents; otherwise we just use some ad-hoc kinematics with one of
  // the light cluster constituents being off-shell.

  // This is the overall cluster momentum - we do not need it - and its
  // individual components, i.e. the momenta of the Proto_Particles
  // it is made of --- transverse momentum goes to new particle
  //msg_Out()<<METHOD<<".\n";
  Vec4D newmom11 = m_E*(         m_z1*s_AxisP+          m_x*s_AxisM);
  Vec4D newmom12 = m_E*((1.-m_y-m_z1)*s_AxisP+(1.-m_x-m_z2)*s_AxisM)+m_ktvec;
  // back into lab system
  m_rotat.RotateBack(newmom11);
  m_boost.BoostBack(newmom11);
  m_rotat.RotateBack(newmom12);
  m_boost.BoostBack(newmom12);
  if (newmom11.Abs2()-m_m12 > 1.e-8)          abort();
  // Update momentum of original (anti-)(di-) quark after gluon splitting
  //if (p_part1->Flavour()==Flavour(kf_b) ||
  //    p_part1->Flavour()==Flavour(kf_b).Bar())
  //  msg_Out()<<METHOD<<": "<<p_part1->Momentum()<<" --> "<<newmom11<<".\n";
  p_part1->SetMomentum(newmom11);
  // Make new particle
  Proto_Particle * newp12 = new Proto_Particle(m_newflav1,newmom12,'l');
  // Take care of sequence in cluster = triplet + anti-triplet
  //msg_Out()<<METHOD<<" constructs new cluster with mass = "
  //	   <<sqrt((newmom11+newmom12).Abs2())<<" for "
  //	   <<"["<<newp12->Flavour()<<" "<<p_part1->Flavour()<<"].\n";
  return (m_barrd?new Cluster(newp12,p_part1):new Cluster(p_part1,newp12));
}
