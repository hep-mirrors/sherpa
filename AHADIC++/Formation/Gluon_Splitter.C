#include "AHADIC++/Formation/Gluon_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;

bool Gluon_Splitter::MakeLongitudinalMomenta() {
  if (4.*(m_popped_mass2+m_kt2)/sqr(m_Q-m_mass1)>1.) {
    msg_Error()<<METHOD<<" throws error: impossible kt^2 = "<<m_kt2<<"\n"
	       <<"   for Q = "<<m_Q<<" and M = "<<m_mass1<<" "
	       <<"with m = "<<m_popped_mass<<"\n";
    return false;
  }
  CalculateLimits();
  do {
    m_z = m_zselector(m_zmin,m_zmax);
  } while (!CalculateXY());
  return true;
}

void Gluon_Splitter::CalculateLimits() {
  double arg = sqrt(1.-4.*(m_popped_mass2+m_kt2)/sqr(m_Q-m_mass1));
  m_zmin = (1.-arg)/2.;
  m_zmax = (1.+arg)/2.;
}

bool Gluon_Splitter::CalculateXY() {
  double R2  = (m_popped_mass2+m_kt2)/(m_z*(1.-m_z));
  double arg = sqr(m_Q2+R2-m_m12)-4.*m_Q2*R2;
  if (arg<0.) return false;
  m_y = (m_Q2-m_m12+R2-sqrt(arg))/(2.*m_Q2);
  m_x = m_m12/((1.-m_y)*m_Q2);
  //msg_Out()<<METHOD<<":\n"
  //	   <<"*** from {x, y, z, kt = "<<m_x<<", "<<m_y<<", "<<m_z<<", "
  //	   <<m_kt<<"}\n"
  //	   <<"*** m2 = "<<m_popped_mass2<<", kt2 = "<<m_kt2<<", "
  //	   <<"M2 = "<<m_m12<<"}.\n";
  return (m_x<=1. && m_x>=0. && m_y<=1. && m_y>=0.);
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

bool Gluon_Splitter::FillParticlesInLists() {
  Cluster * cluster = MakeCluster();
  switch (p_softclusters->Treat(cluster)) {
  case 1:
    delete cluster;
    break;
  case -1:
    return false;
  default:
    p_cluster_list->push_back(cluster);
    break;
  }
  UpdateSpectator();
  return true;
}

void Gluon_Splitter::UpdateSpectator() {
  // Replace splitted gluon with (anti-)(di-)quark and correct momentum
  Vec4D newmom2 = m_E*(m_y*(1.-m_z)*s_AxisP + m_z*(1.-m_x)*s_AxisM)-m_ktvec;
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
  Vec4D newmom11 = m_E*((1.-m_y)*s_AxisP +               m_x*s_AxisM);
  Vec4D newmom12 = m_E*( m_y*m_z*s_AxisP + (1.-m_z)*(1.-m_x)*s_AxisM)+m_ktvec;
  Vec4D newmom2 = m_E*(m_y*(1.-m_z)*s_AxisP + m_z*(1.-m_x)*s_AxisM)-m_ktvec;
  // back into lab system
  m_rotat.RotateBack(newmom11);
  m_boost.BoostBack(newmom11);
  m_rotat.RotateBack(newmom12);
  m_boost.BoostBack(newmom12);
  if (dabs(newmom11.Abs2()-m_m12)          > 1.e-8 ||
      dabs(newmom12.Abs2()-m_popped_mass2) > 1.e-8) {
    msg_Error()<<"Error in "<<METHOD<<": masses not respected.\n"
	       <<newmom11<<" -> "<<sqrt(dabs(newmom11.Abs2()))
	       <<" vs. "<<m_mass1<<"\n"
	       <<newmom12<<" -> "<<sqrt(dabs(newmom12.Abs2()))
	       <<" vs. "<<m_popped_mass<<" from "<<m_newflav1<<"\n"
      	       <<newmom2<<" -> "<<sqrt(dabs(newmom2.Abs2()))
	       <<" vs. "<<m_popped_mass<<" from "<<m_newflav1<<"\n"
	       <<"*** from {x, y, z, kt} = "
	       <<"{"<<m_x<<", "<<m_y<<", "<<m_z<<", "<<m_kt<<"}, "
	       <<" Q = "<<m_Q<<"}.\n"
	       <<"*** mom = "<<p_part1->Momentum()<<" and "
	       <<p_part2->Momentum()<<".\n";
    exit(1);
  }
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
