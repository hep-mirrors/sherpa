#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;

Cluster_Splitter::Cluster_Splitter(std::list<Cluster *> * cluster_list,
				   Soft_Cluster_Handler * softclusters) :
  Splitter_Base(cluster_list,softclusters) {}  

void Cluster_Splitter::Init() {
  Splitter_Base::Init();
  m_alpha[0] = hadpars->Get("alphaL");
  m_beta[0]  = hadpars->Get("betaL");
  m_gamma[0] = hadpars->Get("gammaL");
  m_alpha[1] = hadpars->Get("alphaH");
  m_beta[1]  = hadpars->Get("betaH");
  m_gamma[1] = hadpars->Get("gammaH");
  m_alpha[2] = hadpars->Get("alphaD");
  m_beta[2]  = hadpars->Get("betaD");
  m_gamma[2] = hadpars->Get("gammaD");
}

bool Cluster_Splitter::MakeLongitudinalMomenta() {
  CalculateLimits();
  double ran1, ran2, dM12, dM22, arg1, arg2;
  double mt2 = m_popped_mass2+m_kt2; 
  double Mup2 = sqr(sqrt(m_Q2)-m_mass1-m_mass2);
  size_t trials(0);
  do {
    dM12 = DeltaM2(m_Q2,mt2,m_flavs1.first);
    dM22 = DeltaM2(m_Q2,mt2,m_flavs2.second);
    m_R12 = m_m12min+dM12-m_kt2;
    m_R21 = m_m22min+dM22-m_kt2;
    arg1  = sqr(m_Q2-m_R21+m_R12)-4.*m_Q2*(m_R12+m_kt2);
    arg2  = sqr(m_Q2-m_R12+m_R21)-4.*m_Q2*(m_R21+m_kt2);
    if (arg1<0. || arg2<0.) {
      m_z1 = m_z2 = -1.;
      continue;
    }
    m_z1  = ((m_Q2-m_R21+m_R12+sqrt(arg1))/(2.*m_Q2));
    m_z2  = ((m_Q2-m_R12+m_R21+sqrt(arg2))/(2.*m_Q2));
    //msg_Out()<<"\n\n\n"
    //         <<METHOD<<"("<<m_flavs1.first<<" "<<m_flavs2.second<<",
    //	     <<"kt = "<<m_kt<<", mass = "<<m_popped_mass<<"): "
    // 	     <<"M1 = "<<sqrt(dM12)<<", M2 = "<<sqrt(dM12)<<".\n"
    //	     <<"   --> "<<m_z1<<" and "<<m_z2
    //	     <<" from minimals = "<<m_m12min<<" "<<m_m22min
    //	     <<" --> "<<m_R12<<" and "<<m_R21<<".\n";
  } while (trials++<1000 && (m_z1>1. || m_z2>1. || m_z1<0. || m_z2<0.));
  if (trials>=999) {
    msg_Out()<<METHOD<<" failed for cluster with mass = "
	     <<sqrt(m_Q2)<<", popped = "<<m_newflav2
	     <<" and kt = "<<sqrt(m_kt2)<<":\n"
	     <<(*p_part1)<<(*p_part2)<<"\n";
    return false;
  }
  return CheckIfAllowed();
}

double Cluster_Splitter::DeltaM2(const double & MC2,const double & mt2,
				 const Flavour & flav) {
  m_flcnt = 0;
  if (flav==Flavour(kf_b) || flav==Flavour(kf_b).Bar() ||
      flav==Flavour(kf_c) || flav==Flavour(kf_c).Bar()) {
    m_flcnt = 1;
  }
  else if (flav.IsDiQuark()) {
    m_flcnt = 2;
  }
  double m2, M2 = (m_flcnt==1?sqr(flav.HadMass()):MC2);
  double max = (pow(m_alpha[m_flcnt]*m_gamma[m_flcnt],-m_alpha[m_flcnt]) *
		exp(-1./m_alpha[m_flcnt]));
  do {
    m2 = ran->Get()*MC2;
  } while (pow(m2/M2,m_alpha[m_flcnt])*
	   exp(-m_gamma[m_flcnt]*m2/M2)<ran->Get()*max);
  return m2;
}

bool Cluster_Splitter::CheckIfAllowed() {
  //msg_Out()<<METHOD<<"(z1 = "<<m_z1<<", z2 = "<<m_z2<<"):\n"
  //	   <<"   ("<<p_part1->Flavour()<<"+"<<m_newflav1<<"), "
  //	   <<"m = "<<sqrt(m_R12)<<" > "<<m_minQ_1<<"\n"
  //	   <<"   ("<<m_newflav2<<"+"<<p_part2->Flavour()<<"), "
  //	   <<"m = "<<sqrt(m_R21)<<" > "<<m_minQ_2<<"\n";
  return (m_R12>sqr(m_minQ_1) && m_R21>sqr(m_minQ_2));
}

void Cluster_Splitter::FillParticlesInLists() {
  for (size_t i=0;i<2;i++) {
    Cluster * cluster = MakeCluster(i);
    if (p_softclusters->Treat(cluster)) {
      //msg_Out()<<"   ---> decays to hadrons, erase it.\n";
      delete cluster;
    }
    else {
      //msg_Out()<<"   ---> keep it.\n";
      p_cluster_list->push_back(cluster);
    }
  }
}

Cluster * Cluster_Splitter::MakeCluster(size_t i) {
  double alpha = (i==0? m_z1  : 1.-m_z1 );
  double beta  = (i==0? m_z2  : 1.-m_z2 );
  double mass2 = (i==0? m_m12 : m_m22);
  double sign  = (i==0?    1. : -1.);
  double R2    = (i==0? m_R12 : m_R21);
  double R02   = mass2+(m_popped_mass2+m_kt2);
  double ab    = 4.*mass2*(m_popped_mass2+m_kt2);
  double x = 1;
  if (sqr(R2-R02)>ab) {
    double centre = (R2+mass2-(m_popped_mass2+m_kt2))/(2.*R2);
    double lambda = Lambda(R2,mass2,m_popped_mass2+m_kt2);
    x = centre+lambda;
  }
  double y      = mass2/(x*R2);
  // This is the overall cluster momentum - we do not need it - and its
  // individual components, i.e. the momenta of the Proto_Particles
  // it is made of.
  Vec4D newmom11 = (m_E*(     x*alpha*s_AxisP+     y*(1.-beta)*s_AxisM));
  Vec4D newmom12 = (m_E*((1.-x)*alpha*s_AxisP+(1.-y)*(1.-beta)*s_AxisM) +
		    sign * m_ktvec);
  Vec4D clumom = m_E*(alpha*s_AxisP + (1.-beta)*s_AxisM) + sign * m_ktvec;
  
  // back into lab system
  m_rotat.RotateBack(newmom11);
  m_boost.BoostBack(newmom11);
  m_rotat.RotateBack(newmom12);
  m_boost.BoostBack(newmom12);
  (i==0?p_part1:p_part2)->SetMomentum(newmom11);
  Proto_Particle * newp =
    new Proto_Particle((i==0?m_newflav1:m_newflav2),newmom12,'l');

  //msg_Out()<<METHOD<<"["<<p_part1->Flavour()<<", "<<p_part2->Flavour()<<"]"
  //	   <<" makes cluster [i = "<<i<<": "
  //	   <<(i==0?p_part1->Flavour():m_newflav2)<<", "
  //	   <<(i==0?m_newflav1:p_part2->Flavour())<<"]\n";

  return (i==0?new Cluster(p_part1,newp):new Cluster(newp,p_part2));
}

void Cluster_Splitter::CalculateLimits() {
  m_z1min = m_z2min = 0.;
  m_z1max = m_z2max = 1.;
  // sum of constituent masses squared + transverse momentum squared 
  m_m12min = sqr(m_minQ_1); 
  m_m22min = sqr(m_minQ_2);
}

