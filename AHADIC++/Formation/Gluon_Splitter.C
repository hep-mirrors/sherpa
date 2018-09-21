#include "AHADIC++/Formation/Gluon_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;

Gluon_Splitter::~Gluon_Splitter() {
  msg_Out()<<METHOD<<" with "<<m_kin_fails<<" kinematic fails.\n";
}


void Gluon_Splitter::Init(const bool & isgluon) {
  Splitter_Base::Init(true);
  m_alpha = hadpars->Get("alphaG");
  m_beta  = hadpars->Get("betaG");
  m_gamma = hadpars->Get("gammaG");
}
  
bool Gluon_Splitter::MakeLongitudinalMomenta() {
  m_arg = (sqr(m_Q2-m_minQ_12-m_popped_mass2)-
	   4.*(m_Q2*m_kt2 + m_minQ_12*m_popped_mass2));
  if (m_arg<0.) return false;
  CalculateLimits();
  do { m_z2 = m_zselector(m_z2min,m_z2max); } while (!CalculateXY());
  return true;
}

void Gluon_Splitter::CalculateLimits() {
  double mean1 = (m_Q2+m_minQ_12-m_popped_mass2)/(2.*m_Q2);
  double delta = sqrt(m_arg)/(2.*m_Q2);
  m_z1min = Max(0.0,mean1-delta);
  m_z1max = Min(1.0,mean1+delta);
  double mean2 = (m_Q2-m_minQ_12+m_popped_mass2)/(2.*m_Q2);
  m_z2min = Max(0.0,mean2-delta/2.);
  m_z2max = Min(1.0,mean2+delta);
}

bool Gluon_Splitter::CalculateXY() {
  m_z1 = 1.-(m_popped_mass2+m_kt2)/(m_z2*m_Q2);
  double M2  = m_z1*(1.-m_z2)*m_Q2;
  double arg = sqr(M2-m_kt2-m_m12)-4*m_m12*m_kt2;
  if (arg<0.) return false;
  m_x = ((M2-m_kt2+m_m12)+sqrt(arg))/(2.*M2);
  m_y = m_kt2>1.e-6?m_kt2/M2/(1.-m_x):1.;
  return (!(m_x>1.) && !(m_x<0.) && !(m_y>1.) && !(m_y<0.));
}

double Gluon_Splitter::
WeightFunction(const double & z,const double & zmin,const double & zmax) {
  double norm = 1.;
  if (m_alpha<=0. && m_beta<=0.)
    norm *= Max(pow(zmin,m_alpha), pow(1.-zmax,m_beta));
  else {
    if (m_alpha<=0.) norm *= pow(zmin,m_alpha);
    if (m_beta<=0.)  norm *= pow(1.-zmax,m_beta);
  }
  return (pow(z,m_alpha)+pow(1.-z,m_beta))/norm;
}

bool Gluon_Splitter::CheckKinematics() {
  // check if:
  // 1. new cluster mass larger than minimal mass
  // 2. spectator still on its mass shell
  // 3. new cluster particle with mass 0
  // 4. new particle after gluon splitting on its mass-shell.
  double M2 = m_z1*(1.-m_z2)*m_Q2;
  if (M2-m_kt2-m_minQ_12 < 1.e-6*m_Q2 ||
      dabs(m_x*(1.-m_y)*M2-m_m12) > 1.e-6*m_Q2 ||
      dabs((1.-m_x)*m_y*M2-m_kt2) > 1.e-6*m_Q2 ||
      dabs((1.-m_z1)*m_z2*m_Q2-m_kt2-m_popped_mass2) > 1.e-6*m_Q2) {
    msg_Tracking()<<"Error in "<<METHOD<<": failed to reconstruct masses.\n"
		  <<"   cluster mass:"<<(m_z1*(1.-m_z2)*m_Q2-m_kt2)<<" > "
		  <<m_minQ_12<<",\n"
		  <<"   spectator mass:"<<(m_x*(1.-m_y)*m_z1*(1.-m_z2)*m_Q2)
		  <<" vs. "<<m_m12<<" ("<<p_part2->Flavour()<<"),\n"
		  <<"   new in-quark:"<<((1.-m_x)*m_y*m_z1*(1.-m_z2)*m_Q2-m_kt2)
		  <<" should be 0 for ("<<m_newflav1<<")\n"
		  <<"   new out-quark:"<<((1.-m_z1)*m_z2*m_Q2-m_kt2)<<" vs. "
		  <<m_popped_mass2<<".\n";
    m_kin_fails++;
    return false;
  }
  if (p_part3==0) return true;
  Vec4D newmom2 = m_E*((1.-m_z1)*s_AxisP+m_z2*s_AxisM)-m_ktvec;
  m_rotat.RotateBack(newmom2);
  m_boost.BoostBack(newmom2);
  return ((newmom2+p_part3->Momentum()).Abs2()>sqr(m_minQ_2));
}

bool Gluon_Splitter::FillParticlesInLists() {
  Cluster * cluster = MakeCluster();
  switch (p_softclusters->Treat(cluster)) {
  case 1:
    delete cluster;
    break;
  case -1:
    delete cluster;
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
  Vec4D newmom2 = m_E*((1.-m_z1)*s_AxisP+m_z2*s_AxisM)-m_ktvec;
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
  Vec4D newmom11 = (m_E*(m_z1*     m_x*s_AxisP + (1.-m_z2)*(1.-m_y)*s_AxisM));
  Vec4D newmom12 = (m_E*(m_z1*(1.-m_x)*s_AxisP + (1.-m_z2)*     m_y*s_AxisM) +
		    m_ktvec);
  // back into lab system
  m_rotat.RotateBack(newmom11);
  m_boost.BoostBack(newmom11);
  m_rotat.RotateBack(newmom12);
  m_boost.BoostBack(newmom12);
  if (dabs(newmom11.Abs2()-m_m12) > 1.e-3*m_Q2 ||
      dabs(newmom12.Abs2())       > 1.e-3*m_Q2) {
    Vec4D newmom2  = m_E*((1.-m_z1)*s_AxisP+m_z2*s_AxisM) - m_ktvec;
    m_rotat.RotateBack(newmom2);
    m_boost.BoostBack(newmom2);
    msg_Error()<<"Error in "<<METHOD<<": masses not respected.\n"
	       <<newmom11<<" -> "<<sqrt(dabs(newmom11.Abs2()))
	       <<" vs. "<<m_mass1<<"\n"
	       <<newmom12<<" -> "<<sqrt(dabs(newmom12.Abs2()))
	       <<" vs. "<<m_popped_mass<<" from "<<m_newflav1<<"\n"
      	       <<newmom2<<" -> "<<sqrt(dabs(newmom2.Abs2()))
	       <<" vs. "<<m_popped_mass<<" from "<<m_newflav1<<"\n"
	       <<"*** from {x, y, z1, z2, kt} = "
	       <<"{"<<m_x<<", "<<m_y<<", "<<m_z1<<", "<<m_z2<<", "<<m_kt<<"}, "
	       <<" Q = "<<m_Q<<", M = "<<sqrt(m_Q2*m_z1*(1.-m_z2)-m_kt2)<<", "
	       <<"ktvec = "<<m_ktvec<<"("<<m_ktvec.Abs()<<").\n"
	       <<"*** mom = "
	       <<p_part1->Momentum()<<"("<<p_part1->Flavour()<<") and "
	       <<p_part2->Momentum()<<"("<<p_part2->Flavour()<<").\n";
    m_kin_fails++;
    return NULL;
  }
  // Update momentum of original (anti-)(di-) quark after gluon splitting
  p_part1->SetMomentum(newmom11);
  // Make new particle
  Proto_Particle * newp12 = new Proto_Particle(m_newflav1,newmom12,'l');
  // Take care of sequence in cluster = triplet + anti-triplet
  Cluster * cluster(m_barrd?new Cluster(newp12,p_part1):new Cluster(p_part1,newp12));
  return cluster;
}
