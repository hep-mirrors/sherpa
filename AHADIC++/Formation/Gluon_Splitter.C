#include "AHADIC++/Formation/Gluon_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;

Gluon_Splitter::~Gluon_Splitter() {
  msg_Debugging()<<METHOD<<" with "<<m_kin_fails<<" kinematic fails.\n";
}


void Gluon_Splitter::Init(const bool & isgluon) {
  Splitter_Base::Init(true);
  // Gluon Decay Form: 1 = default
  // 0: z ~ z^alpha * (1-z)^alpha
  // 1: z ~ z^alpha + (1-z)^alpha
  m_mode  = hadpars->Switch("GluonDecayForm");
  m_alpha = hadpars->Get("alphaG");
  m_analyse = true;
  if (m_analyse) {
    m_histograms[std::string("Yasym_frag_2")] = new Histogram(0,0.,8.,32);
  }
}
  
bool Gluon_Splitter::MakeLongitudinalMomenta() {
  m_arg = (sqr(m_Q2-m_minQ2[0]-m_popped_mass2)-
	   4.*(m_Q2*m_kt2 + m_minQ2[0]*m_popped_mass2));
  if (m_arg<0.) return false;
  CalculateLimits();
  do { m_z[1] = m_zselector(m_zmin[1],m_zmax[1]); } while (!CalculateXY());
  return true;
}

void Gluon_Splitter::CalculateLimits() {
  double mean1 = (m_Q2+m_minQ2[0]-m_popped_mass2)/(2.*m_Q2);
  double delta = sqrt(m_arg)/(2.*m_Q2);
  m_zmin[0] = Max(0.0,mean1-delta);
  m_zmax[0] = Min(1.0,mean1+delta);
  double mean2 = (m_Q2-m_minQ2[0]+m_popped_mass2)/(2.*m_Q2);
  m_zmin[1] = Max(0.0,mean2-delta/2.);
  m_zmax[1] = Min(1.0,mean2+delta);
}

bool Gluon_Splitter::CalculateXY() {
  m_z[0] = 1.-(m_popped_mass2+m_kt2)/(m_z[1]*m_Q2);
  double M2 = m_z[0]*(1.-m_z[1])*m_Q2;
  //This is a new addition w.r.t. original master
  double R2 = M2 - m_kt2;
  if (R2 < m_mdec[0]) {
    M2     = m_mdec2[0];
    m_z[0] = M2/((1.-m_z[1])*m_Q2);
  }
  if ((M2/m_m2[0] > 1e6 && M2/m_kt2 > 1e6) || m_kt2<1.e-12) {
    // Use Taylor expansion to first order in m12/M2 and kt2/M2 to avoid
    // numerical instability for x, y -> 1.0
    m_x = 1.0 - m_kt2/M2;
    m_y = 1.0;
  } else {
    double arg = sqr(M2-m_kt2-m_m2[0])-4*m_m2[0]*m_kt2;
    if (arg<0.) return false;
    m_x = ((M2-m_kt2+m_m2[0])+sqrt(arg))/(2.*M2);
    m_y = m_kt2/M2/(1.-m_x);
  }
  return (!(m_x>1.) && !(m_x<0.) && !(m_y>1.) && !(m_y<0.));
}

double Gluon_Splitter::
WeightFunction(const double & z,const double & zmin,const double & zmax,
	       const unsigned int & cnt) {
  double norm = 1.;
  switch (m_mode) {
  case 1:
    norm = pow(0.5,2*m_alpha);
    return pow(z*(1.-z),m_alpha)/norm;
  case 0:
  default:
    break;
  }
  if (m_alpha<=0.) norm = pow(zmin,m_alpha) + pow(1.-zmax,m_alpha);
  return (pow(z,m_alpha)+pow(1.-z,m_alpha))/norm;
}


bool Gluon_Splitter::CheckKinematics() {
  // check if:
  // 1. new cluster mass larger than minimal mass
  // 2. spectator still on its mass shell
  // 3. new cluster particle with mass 0
  // 4. new particle after gluon splitting on its mass-shell.
  double M2 = m_z[0]*(1.-m_z[1])*m_Q2;
  if (M2-m_kt2-m_minQ2[0] < 1.e-6*m_Q2 ||
      dabs(m_x*(1.-m_y)*M2-m_m2[0]) > 1.e-6*m_Q2 ||
      dabs((1.-m_x)*m_y*M2-m_kt2) > 1.e-6*m_Q2 ||
      dabs((1.-m_z[0])*m_z[1]*m_Q2-m_kt2-m_popped_mass2) > 1.e-6*m_Q2) {
    msg_Tracking()<<"Error in "<<METHOD<<": failed to reconstruct masses.\n"
		  <<"   cluster mass:"<<(m_z[0]*(1.-m_z[1])*m_Q2-m_kt2)<<" > "
		  <<m_minQ2[0]<<",\n"
		  <<"   spectator mass:"<<(m_x*(1.-m_y)*m_z[0]*(1.-m_z[1])*m_Q2)
		  <<" vs. "<<m_m2[0]<<" ("<<p_part[1]->Flavour()<<"),\n"
		  <<"   new in-quark:"<<((1.-m_x)*m_y*m_z[0]*(1.-m_z[1])*m_Q2-m_kt2)
		  <<" should be 0 for ("<<m_newflav[0]<<")\n"
		  <<"   new out-quark:"<<((1.-m_z[0])*m_z[1]*m_Q2-m_kt2)<<" vs. "
		  <<m_popped_mass2<<".\n";
    m_kin_fails++;
    return false;
  }
  if (p_part[2]==0) return true;
  Vec4D newmom2 = m_E*((1.-m_z[0])*s_AxisP+m_z[1]*s_AxisM)-m_ktvec;
  m_rotat.RotateBack(newmom2);
  m_boost.BoostBack(newmom2);
  return ((newmom2+p_part[2]->Momentum()).Abs2()>sqr(m_minQ[1]));
}

bool Gluon_Splitter::FillParticlesInLists() {
  Cluster * cluster = MakeCluster();
  if (cluster==NULL) return false;
  Vec4D  mom = cluster->Momentum();
  Flavour fl = Flavour(kf_none);
  if (p_softclusters->PromptTransit(cluster,fl)) {
    ReplaceClusterWithHadron(fl,mom);
    delete cluster;
  }
  else {
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
  }
  UpdateSpectator(mom);
  return true;
}

void Gluon_Splitter::ReplaceClusterWithHadron(const Flavour & fl,Vec4D & mom) {
  double M2 = m_Q2, mt12 = sqr(fl.Mass())+m_kt2, mt22 = m_m2[1]+m_kt2; 
  double alpha1 = ((M2+mt12-mt22)+sqrt(sqr(M2+mt12-mt22)-4.*M2*mt12))/(2.*M2);
  double beta1  = mt12/(M2*alpha1);
  mom = m_E*(alpha1*s_AxisP + beta1*s_AxisM)+m_ktvec;
  m_rotat.RotateBack(mom);
  m_boost.BoostBack(mom);
  p_softclusters->GetHadrons()->push_back(new Proto_Particle(fl,mom,false));
}


void Gluon_Splitter::UpdateSpectator(const Vec4D & clumom) {
  // Replace splitted gluon with (anti-)(di-)quark and correct momentum
  p_part[1]->SetFlavour(m_newflav[1]);
  p_part[1]->SetMomentum(m_Qvec-clumom);
  p_part[1]->SetKT2_Max(m_kt2);
}

Cluster * Gluon_Splitter::MakeCluster() {
  // If kinematically allowed, i.e. if the emerging cluster is heavy enough,
  // we calculate the split of cluster momentum into momenta for its
  // constituents; otherwise we just use some ad-hoc kinematics with one of
  // the light cluster constituents being off-shell.

  // This is the overall cluster momentum - we do not need it - and its
  // individual components, i.e. the momenta of the Proto_Particles
  // it is made of --- transverse momentum goes to new particle
  Vec4D newmom11 = (m_E*(m_z[0]*     m_x*s_AxisP + (1.-m_z[1])*(1.-m_y)*s_AxisM));
  Vec4D newmom12 = (m_E*(m_z[0]*(1.-m_x)*s_AxisP + (1.-m_z[1])*     m_y*s_AxisM) +
		    m_ktvec);
  // back into lab system
  m_rotat.RotateBack(newmom11);
  m_boost.BoostBack(newmom11);
  m_rotat.RotateBack(newmom12);
  m_boost.BoostBack(newmom12);
  if (!CheckConstituentKinematics(newmom11,newmom12)) {
    m_kin_fails++;
    return NULL;
  }
  // Update momentum of original (anti-)(di-) quark after gluon splitting
  p_part[0]->SetMomentum(newmom11);
  // Make new particle
  Proto_Particle * newp12 = new Proto_Particle(m_newflav[0],newmom12,false,
					       p_part[0]->IsBeam() || p_part[1]->IsBeam());
  newp12->SetKT2_Max(m_kt2);
  // Take care of sequence in cluster = triplet + anti-triplet
  Cluster * cluster(m_barrd?
		    new Cluster(newp12,p_part[0]):
		    new Cluster(p_part[0],newp12));
  // this is for a simple analysis only
  if (m_analyse) {
    m_lastmass = sqrt(dabs(cluster->Momentum().Abs2()));
    m_lastB    = (newp12->Flavour()==Flavour(kf_b) ||
		  newp12->Flavour()==Flavour(kf_b).Bar() ||
		  p_part[0]->Flavour()==Flavour(kf_b) ||
		  p_part[0]->Flavour()==Flavour(kf_b).Bar());
    m_lastC    = (!m_lastB &&
		  (newp12->Flavour()==Flavour(kf_b) ||
		   newp12->Flavour()==Flavour(kf_b).Bar() ||
		   p_part[0]->Flavour()==Flavour(kf_b) ||
		   p_part[0]->Flavour()==Flavour(kf_b).Bar()));
    double y = cluster->Momentum().Y();
    m_histograms[std::string("Yasym_frag_2")]->Insert(dabs(y),(y>0.?1.:-1.));
  }
  return cluster;
}

bool Gluon_Splitter::CheckConstituentKinematics(const ATOOLS::Vec4D & newmom11,
						const ATOOLS::Vec4D & newmom12) {
  if (dabs(newmom11.Abs2()-m_m2[0]) < 1.e-3*m_Q2 &&
      dabs(newmom12.Abs2())         < 1.e-3*m_Q2) return true;
  Vec4D newmom2  = m_E*((1.-m_z[0])*s_AxisP+m_z[1]*s_AxisM) - m_ktvec;
  m_rotat.RotateBack(newmom2);
  m_boost.BoostBack(newmom2);
  msg_Error()<<"Error in "<<METHOD<<": masses not respected.\n"
	     <<newmom11<<" -> "<<sqrt(dabs(newmom11.Abs2()))
	     <<" vs. "<<m_mass[0]<<"\n"
	     <<newmom12<<" -> "<<sqrt(dabs(newmom12.Abs2()))
	     <<" vs. "<<m_popped_mass<<" from "<<m_newflav[0]<<"\n"
	     <<newmom2<<" -> "<<sqrt(dabs(newmom2.Abs2()))
	     <<" vs. "<<m_popped_mass<<" from "<<m_newflav[0]<<"\n"
	     <<"*** from {x, y, z1, z2, kt} = "
	     <<"{"<<m_x<<", "<<m_y<<", "<<m_z[0]<<", "<<m_z[1]<<", "<<m_kt<<"}, "
	     <<" Q = "<<m_Q<<", M = "<<sqrt(m_Q2*m_z[0]*(1.-m_z[1])-m_kt2)<<", "
	     <<"ktvec = "<<m_ktvec<<"("<<m_ktvec.Abs()<<").\n"
	     <<"*** mom = "
	     <<p_part[0]->Momentum()<<"("<<p_part[0]->Flavour()<<") and "
	     <<p_part[1]->Momentum()<<"("<<p_part[1]->Flavour()<<").\n";
  return false;
}
