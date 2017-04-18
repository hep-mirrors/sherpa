#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

bool Cluster_Splitter::MakeLongitudinalMomenta() {
  CalculateLimits();
  m_z1  = m_zselector(m_z1min,m_z1max,p_part1->Flavour());
  m_z2  = m_zselector(m_z2min,m_z2max,p_part2->Flavour());
  m_R12 = m_z1*(1.-m_z2)*m_Q2-m_kt2;
  m_R21 = (1.-m_z1)*m_z2*m_Q2-m_kt2;
  return CheckIfAllowed();
}

bool Cluster_Splitter::CheckIfAllowed() {
  //msg_Out()<<METHOD<<":\n"
  //	   <<"   ("<<p_part1->Flavour()<<"+"<<m_newflav1<<"), "
  //	   <<"m = "<<sqrt(m_R12)<<" > "<<m_minQ_1<<"\n"
  //	   <<"   ("<<m_newflav2<<"+"<<p_part2->Flavour()<<"), "
  //	   <<"m = "<<sqrt(m_R21)<<" > "<<m_minQ_2<<"\n";
  return (m_R12>sqr(m_minQ_1) && m_R21>sqr(m_minQ_2));
}

void Cluster_Splitter::FillParticlesInLists() {
  for (size_t i=0;i<2;i++) {
    Cluster * cluster = MakeCluster(i);
    //msg_Out()<<METHOD<<" for "<<(*cluster);
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
  //Vec4D newmom1= m_E* (       alpha*s_AxisP+       (1.-beta)*s_AxisM) +
  //               sign * m_ktvec;
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
}

