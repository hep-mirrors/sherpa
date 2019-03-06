#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Splitter::Cluster_Splitter(list<Cluster *> * cluster_list,
				   Soft_Cluster_Handler * softclusters) :
  Splitter_Base(cluster_list,softclusters),
  m_mode(1),
  m_output(false) 
{
  m_analyse = false;
  if (m_analyse) {
    m_histograms[string("kt")]    = new Histogram(0,0.,5.,100);
    m_histograms[string("z1")]    = new Histogram(0,0.,1.,100);
    m_histograms[string("z2")]    = new Histogram(0,0.,1.,100);
    m_histograms[string("mass")]  = new Histogram(0,0.,100.,200);
    m_histograms[string("Rmass")] = new Histogram(0,0.,2.,100);
  }
}

void Cluster_Splitter::Init(const bool & isgluon) {
  Splitter_Base::Init(false);
  m_alpha[0] = hadpars->Get("alphaL");
  m_beta[0]  = hadpars->Get("betaL");
  m_gamma[0] = hadpars->Get("gammaL");
  m_alpha[1] = hadpars->Get("alphaH");
  m_beta[1]  = hadpars->Get("betaH");
  m_gamma[1] = hadpars->Get("gammaH");
  m_alpha[2] = hadpars->Get("alphaD");
  m_beta[2]  = hadpars->Get("betaD");
  m_gamma[2] = hadpars->Get("gammaD");
  m_alpha[3] = hadpars->Get("alphaB");
  m_beta[3]  = hadpars->Get("betaB");
  m_gamma[3] = hadpars->Get("gammaB");
  m_kt02     = sqr(hadpars->Get("kt_o"));
}

bool Cluster_Splitter::MakeLongitudinalMomenta() {
  CalculateLimits();
  if (m_mode==0) {
    FixCoefficients(p_part1->Flavour(),p_part2->Flavour());
    m_z1  = m_zselector(m_z1min,m_z1max);
    FixCoefficients(p_part2->Flavour(),p_part1->Flavour());
    m_z2  = m_zselector(m_z2min,m_z2max);
    m_R12 = m_z1*(1.-m_z2)*m_Q2-m_kt2;
    m_R21 = (1.-m_z1)*m_z2*m_Q2-m_kt2;
  }
  else if (m_mode==1) {
    do {
      FixCoefficients(p_part1->Flavour(),p_part2->Flavour());
      m_R12   = m_m12min + DeltaM2(0); 
      FixCoefficients(p_part2->Flavour(),p_part1->Flavour());
      m_R21   = m_m22min + DeltaM2(1);
    } while (m_R12+m_R21>m_Q2);
    double e12  = (m_R12+m_kt2)/m_Q2, e21 = (m_R21+m_kt2)/m_Q2;
    double disc = sqr(1-e12-e21)-4.*e12*e21;
    if (disc<0.) return false;
    disc = sqrt(disc);
    m_z1 = (1.+e12-e21+disc)/2.;
    m_z2 = (1.-e12+e21+disc)/2.;
  }
  return CheckIfAllowed();
}

void Cluster_Splitter::
FixCoefficients(const Flavour & flav1,const Flavour & flav2) {
  m_flcnt   = 0;
  m_masses2 = sqr(Max(1.,(p_constituents->Mass(flav1) +
			  p_constituents->Mass(flav2))));
  if (flav1==Flavour(kf_b) || flav1==Flavour(kf_b).Bar() ||
      flav1==Flavour(kf_c) || flav1==Flavour(kf_c).Bar()) 
    m_flcnt = 1;
  else if (p_part1->IsLeading() ||
	   ((m_mode==0) && p_part2->IsLeading())) {
    m_flcnt = 1;
    m_masses2 *= 4.;
  }
  else if (flav1.IsDiQuark())
    m_flcnt = 2;
  else if (p_part1->IsBeam() || p_part2->IsBeam())
    m_flcnt = 3;
  m_a = m_alpha[m_flcnt];
  m_b = m_beta[m_flcnt];
  m_c = m_gamma[m_flcnt];
}

void Cluster_Splitter::CalculateLimits() {
  m_m12min =
    sqr(Min(p_softclusters->MinSingleMass(m_flavs1.first,m_newflav1),
	    p_softclusters->MinDoubleMass(m_flavs1.first,m_newflav1)));
  m_m22min =
    sqr(Min(p_softclusters->MinSingleMass(m_newflav2,m_flavs2.second),
	    p_softclusters->MinDoubleMass(m_newflav2,m_flavs2.second)));
  double lambda  = sqrt(sqr(m_Q2-m_m12min-m_m22min)-
			4.*(m_m12min+m_kt2)*(m_m22min+m_kt2));
  double centre1 = m_Q2-m_m22min+m_m12min;
  double centre2 = m_Q2-m_m12min+m_m22min;
  m_z1min   = (centre1-lambda)/(2.*m_Q2);
  m_z1max   = (centre1+lambda)/(2.*m_Q2);
  m_z2min   = (centre2-lambda)/(2.*m_Q2);
  m_z2max   = (centre2+lambda)/(2.*m_Q2);
  m_mean12  = (m_minQ_1+m_maxQ_1)/2.;
  m_mean21  = (m_minQ_2+m_maxQ_2)/2.;
  m_sigma12 = sqr(Max(1.,m_maxQ_1-m_minQ_1));
  m_sigma21 = sqr(Max(1.,m_maxQ_2-m_minQ_2));
}

double Cluster_Splitter::DeltaM2(const size_t & cl) {
  double deltaMmax = Min(sqrt(m_Q2-m_m12min-m_m22min),10.+(cl==0?m_maxQ_1:m_maxQ_2));
  double sigmaM    = m_c * (cl==0?m_sigma12:m_sigma21);
  double wtmax     = pow(deltaMmax/sigmaM,m_a);
  double deltaM, weight;
  do {
    deltaM = ran->Get()*deltaMmax;
    weight = pow(deltaM/sigmaM,m_a)*exp(-sqr(deltaM)/sigmaM);
  } while (weight<ran->Get()*wtmax);
  return sqr(deltaM);
}

bool Cluster_Splitter::CheckIfAllowed() {
  bool allowed = (m_R12>m_minQ_12 && m_R21>m_minQ_22);
  //if (p_part1->IsLeading() && m_R12>sqr(m_maxQ_1))
  //  allowed = allowed &&
  //    (exp(-sqr(sqrt(m_R12)-m_maxQ_1)/(2.*m_sigma12))>ran->Get());
  //if (p_part2->IsLeading() && m_R21>sqr(m_maxQ_2))
  //  allowed = allowed &&
  //    (exp(-sqr(sqrt(m_R21)-m_maxQ_2)/(2.*m_sigma21))>ran->Get()); 
  return allowed;
}

double Cluster_Splitter::
WeightFunction(const double & z,const double & zmin,const double & zmax) {
  double arg  = dabs(m_c)>1.e-2? m_c*(m_kt2+m_masses2)/m_kt02 : 0. ;
  double norm = exp(-arg/zmax); //*(1.-zmax)); 
  if (m_a<=0. && m_b<=0.)
    norm *= Max(pow(zmin,m_a), pow(1.-zmax,m_b));
  else {
    if (m_a<=0.) norm *= pow(zmin,m_a);
    if (m_b<=0.) norm *= pow(1.-zmax,m_b);
  }
  double wt = pow(z,m_a) * pow(1.-z,m_b) * exp(-arg/z); //*(1.-z));
  if (wt>norm) {
    msg_Error()<<"Error in "<<METHOD<<": wt(z) = "<<wt<<"("<<z<<") "
	       <<"for wtmax = "<<norm<<" "
	       <<"[a, b, c = "<<m_a<<", "<<m_b<<", "<<m_c<<"].\n";
    exit(1);
  }
  return wt / norm;
}

bool Cluster_Splitter::FillParticlesInLists() {
  for (size_t i=0;i<2;i++) p_cluster_list->push_back(MakeCluster(i));
  return true;
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
    x = (i==0)? centre+lambda : centre-lambda;
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
    new Proto_Particle((i==0?m_newflav1:m_newflav2),newmom12,
		       false,p_part1->IsBeam()||p_part2->IsBeam());
  Cluster * cluster(i==0?new Cluster(p_part1,newp):new Cluster(newp,p_part2));
  /*if (m_output) {
  double mass = sqrt((newmom11+newmom12).Abs2());
  if (mass<m_popped_mass+(i==0?m_mass1:m_mass2)) {
    msg_Out()<<"Small masses(i="<<i<<") --> [";
    if (i==0) msg_Out()<<p_part1->Flavour()<<"/"<<m_newflav1;
    if (i==1) msg_Out()<<m_newflav2<<"/"<<p_part2->Flavour();
    msg_Out()<<"], m = "<<mass<<" vs. ";
    if (i==0) msg_Out()<<m_mass1<<"/"<<m_popped_mass;
    if (i==1) msg_Out()<<m_popped_mass<<"/"<<m_mass2;
    msg_Out()<<"\n"<<"    "<<(*cluster);
  }
  }*/
  if (m_analyse && m_Q>91.) {
    if (i==1) {
      m_histograms[string("kt")]->Insert(sqrt(m_kt2));
      m_histograms[string("z1")]->Insert(m_z1);
      m_histograms[string("z2")]->Insert(m_z2);
    }
    m_histograms[string("mass")]->Insert(sqrt(R2));
    m_histograms[string("Rmass")]->Insert(2.*sqrt(R2/m_Q2));
  }
  return cluster;
}
