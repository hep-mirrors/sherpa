#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Splitter::Cluster_Splitter(list<Cluster *> * cluster_list,
				   Soft_Cluster_Handler * softclusters) :
  Splitter_Base(cluster_list,softclusters),
  m_output(false) 
{
}

void Cluster_Splitter::Init(const bool & isgluon) {
  Splitter_Base::Init(false);
  m_defmode  = hadpars->Switch("ClusterSplittingForm");
  m_beammode = hadpars->Switch("RemnantSplittingForm");
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
  m_kt02     = sqr(hadpars->Get("kT_0"));
  m_analyse  = false; //hadpars->Switch("Analysis");
  if (m_analyse) {
    m_histograms[string("kt")]      = new Histogram(0,0.,5.,100);
    m_histograms[string("z1")]      = new Histogram(0,0.,1.,100);
    m_histograms[string("z2")]      = new Histogram(0,0.,1.,100);
    m_histograms[string("mass")]    = new Histogram(0,0.,100.,200);
    m_histograms[string("Rmass")]   = new Histogram(0,0.,2.,100);
    m_histograms[string("kt_0")]    = new Histogram(0,0.,5.,100);
    m_histograms[string("z1_0")]    = new Histogram(0,0.,1.,100);
    m_histograms[string("z2_0")]    = new Histogram(0,0.,1.,100);
    m_histograms[string("mass_0")]  = new Histogram(0,0.,100.,200);
    m_histograms[string("Rmass_0")] = new Histogram(0,0.,2.,100);
  }
}

bool Cluster_Splitter::MakeLongitudinalMomenta() {
  CalculateLimits();
  FixCoefficients();
  switch (m_mode) {
  case 3:
    return MakeLongitudinalMomentaZ();
  case 2:
    return MakeLongitudinalMomentaZSimple();
  case 1:
    return MakeLongitudinalMomentaMassSimple();
  case 0:
  default:
    return MakeLongitudinalMomentaMass();
  }
  return false;
}

void Cluster_Splitter::FixCoefficients() {
  // this is where the magic happens.  
  m_mode = m_defmode;
  double sum_mass = 0, massfac;
  double threshold = p_softclusters->DecayThreshold(p_part[0]->Flavour(),p_part[1]->Flavour()); 
  for (size_t i=0;i<2;i++) {
    Proto_Particle * part = p_part[i];
    Flavour flav = part->Flavour();
    massfac      = 1.;
    size_t flcnt = 0;
    if (p_part[i]->IsLeading() ||
	(m_mode==0 && p_part[1-i]->IsLeading())) {
      flcnt   = 1;
      massfac = 2.;
    }
    else if (flav.IsDiQuark())
      flcnt = 2;
    if (part->IsBeam()) {
      flcnt  = 3;
      m_mode = m_beammode;
    }
    m_a[i] = m_alpha[flcnt]; // * m_Q/threshold
    m_b[i] = m_beta[flcnt]  * threshold/m_Q;
    m_c[i] = m_gamma[flcnt];
    sum_mass += massfac * p_constituents->Mass(flav);
  }
  m_masses = Max(1.,sum_mass);
}

void Cluster_Splitter::CalculateLimits() {
  // Masses from Splitter_Base:
  // - constitutents:
  //   m_mass[0,1] and m_m2[0,1] = sqr(m_mass[0,1]), m_popped_mass, m_popped_mass2
  //   m_msum[0,1] = mass+mass for the pairs, m_msum2[0,1] = sqr(m_msum[0,1])
  // - hadrons:
  //   m_minQ[0,1] is lightest single or double transition (double for di-di pairs)
  //   m_mdec is lightest decay transition
  for (size_t i=0;i<2;i++) m_m2min[i] = Min(m_minQ2[i],m_mdec2[i]);
  double lambda  = sqrt(sqr(m_Q2-m_m2min[0]-m_m2min[1])-
			4.*(m_m2min[0]+m_kt2)*(m_m2min[1]+m_kt2));
  for (size_t i=0;i<2;i++) {
    double centre = m_Q2-m_m2min[1-i]+m_m2min[i];
    m_zmin[i]  = (centre-lambda)/(2.*m_Q2);
    m_zmax[i]  = (centre+lambda)/(2.*m_Q2);
    m_mean[i]  = sqrt(m_kt02);
    m_sigma[i] = sqrt(m_kt02);
  }
}

bool Cluster_Splitter::MakeLongitudinalMomentaZ() {
  size_t maxcounts=1000;
  while ((maxcounts--)>0) {
    if (MakeLongitudinalMomentaZSimple()) {
      double weight=1.;
      for (size_t i=0;i<2;i++) {
	if (m_gamma[i]>1.e-4) {
	  double DeltaM2 = m_R2[i]-m_minQ2[i];
	  weight *= DeltaM2>0.?exp(-m_gamma[i]*DeltaM2/m_sigma[i]):0.;
	}
      }
      if (weight>=ran->Get()) return true;
    }
  }
  return false;
}

bool Cluster_Splitter::MakeLongitudinalMomentaZSimple() {
  bool mustrecalc = false;
  for (size_t i=0;i<2;i++) m_z[i]  = m_zselector(m_zmin[i],m_zmax[i],i);
  for (size_t i=0;i<2;i++) {
    m_R2[i] = m_z[i]*(1.-m_z[1-i])*m_Q2-m_kt2;
    //This is a difference w.r.t. master
    if (m_R2[i]<m_mdec2[i]+m_kt2) {
      m_R2[i] = m_mdec2[i]+m_kt2;  //(ran->Get()>0.5?m_mdec2[i]:m_m2min[i])+m_kt2;
      mustrecalc = true;
    }
  }
  // another check: bool allowed = (m_R12>m_minQ_12 && m_R21>m_minQ_22);
  bool ok = (m_R2[0]>m_mdec2[0]+m_kt2) && (m_R2[1]>m_mdec2[1]+m_kt2);
  return (ok && (mustrecalc?RecalculateZs():true));
}

double Cluster_Splitter::
WeightFunction(const double & z,const double & zmin,const double & zmax,
	       const unsigned int & cnt) {
  // identical, just have to check the m_a, m_b, m_c
  double norm = 1., arg;
  double value = 1.;
  if (m_a[cnt]>=0.) norm *= pow(zmax,m_a[cnt]);
               else norm *= pow(zmin,m_a[cnt]);
  if (m_b[cnt]>=0.) norm *= pow(1.-zmin,m_b[cnt]);
               else norm *= pow(1.-zmax,m_b[cnt]);
  double wt = pow(z,m_a[cnt]) * pow(1.-z,m_b[cnt]);

  value = wt/norm;
  
  if (m_mode==2) {
    arg   = dabs(m_c[cnt])>1.e-2 ? m_c[cnt]*(m_kt2+m_masses*m_masses)/m_kt02 : 0.;
    value *= exp(-arg*((zmax-z)/(z*zmax)));
    norm *= exp(-arg/zmax);
    wt   *= exp(-arg/z);
  }

  if (wt>norm) THROW(fatal_error,"wt is larger than assumed maxwt - this should never happen.");
  return value;
}

bool Cluster_Splitter::RecalculateZs() {
  double e12  = (m_R2[0]+m_kt2)/m_Q2, e21 = (m_R2[1]+m_kt2)/m_Q2;
  double disc = sqr(1-e12-e21)-4.*e12*e21;
  if (disc<0.) return false;
  disc = sqrt(disc);
  m_z[0] = (1.+e12-e21+disc)/2.;
  m_z[1] = (1.-e12+e21+disc)/2.;
  return true;
}

bool Cluster_Splitter::MakeLongitudinalMomentaMassSimple() {
  bool success;
  long int trials = 1000;
  do {
    for (size_t i=0;i<2;i++) {
      m_R2[i] = sqr(m_minQ[i] + DeltaM(i));
      if (m_R2[i]<=m_mdec2[i]+m_kt2) {
	m_R2[i] = m_minQ2[i]+m_kt2; //Min(m_minQ2[i],m_mdec2[i])+m_kt2;
      }
    }
    success = m_R2[0]+m_R2[1]<m_Q2 && RecalculateZs();
  } while ((trials--)>0 && !success);
  return trials>0;
}

bool Cluster_Splitter::MakeLongitudinalMomentaMass() {
  size_t maxcounts=1000;
  while ((maxcounts--)>0) {
    if (MakeLongitudinalMomentaMassSimple()) {
      double weight=1.;
      for (size_t i=0;i<2;i++) {
	if (m_alpha[i]>1.e-4) weight *= pow(m_z[i],m_alpha[i]);
	if (m_beta[i]>1.e-4)  weight *= pow(1.-m_z[i],m_beta[i]);
      }
      if (weight>=ran->Get()) return true;
    }
  }
  return false;
}

double Cluster_Splitter::DeltaM(const size_t & cl) {
  double deltaM, deltaMmax = m_Q-sqrt(m_m2min[0])-sqrt(m_m2min[1]);
  double mean =  m_mean[cl], sigma = 1./(m_c[cl] * sqrt(m_kt02));
  double arg  =  1.-exp(-sigma * deltaMmax);
  size_t trials = 1000;
  do {
    // Weibull distribution
    //deltaM = sqrt(offset+pow(-log(ran->Get()),1./m_a[cl])*lambda);
    // Normal distribution
    //deltaM = mean + sigma * ran->GetGaussian();
    // Log-Normal distribution
    //deltaM = exp(log(mean)+log(sigma)*ran->GetGaussian());
    // simple exponential
    deltaM = -1./sigma*log(1.-ran->Get()*arg);
  } while ((deltaM>deltaMmax) && (trials--)>1000);
  return trials>0?deltaM:0.;
}


bool Cluster_Splitter::FillParticlesInLists() {
  size_t shuffle = MakeAndCheckClusters();
  if (shuffle) MakeNewMomenta(shuffle);
  for (size_t i=0;i<2;i++) {
    if (shuffle&(i+1)) FillHadronAndDeleteCluster(i);
    else if (shuffle)  UpdateAndFillCluster(i);
    else p_cluster_list->push_back(p_out[i]);
  }
  /*
  if (shuffle>0) 
    msg_Out()<<METHOD<<" shuffled momenta:\n"
	     <<m_cms<<" -> "<<(m_newmom[0]+m_newmom[1])<<"\n = "<<m_newmom[0]<<" + "<<m_newmom[1]<<"\n";
  else {
    msg_Out()<<METHOD<<" didn't shuffle momenta:\n"
	     <<m_cms<<"("<<sqrt(m_cms.Abs2())<<") -> \n"<<(*p_out[0])<<(*p_out[1]);
    double mass1_0 = sqrt(((*p_out[0])[0]->Momentum()+(*p_out[0])[1]->Momentum()).Abs2());
    double mass2_0 = sqrt(((*p_out[1])[0]->Momentum()+(*p_out[1])[1]->Momentum()).Abs2());
    double mass1_1 = sqrt(((*p_out[0])[0]->Momentum()+(*p_out[1])[0]->Momentum()).Abs2());
    double mass2_1 = sqrt(((*p_out[0])[1]->Momentum()+(*p_out[1])[1]->Momentum()).Abs2());
    msg_Out()<<"--> mass shuffle "<<mass1_0<<" + "<<mass2_0<<" --> "
	     <<mass1_1<<" + "<<mass2_1<<"\n\n";
  }
  */
  return true;
}

size_t Cluster_Splitter::MakeAndCheckClusters() {
  size_t  shuffle = 0;
  for (size_t i=0;i<2;i++) {
    p_out[i]     = MakeCluster(i);
    m_cms       += m_mom[i] = p_out[i]->Momentum();
    m_mass2[i]   = m_mom[i].Abs2();
    if (p_softclusters->PromptTransit(p_out[i],m_fl[i])) shuffle += (i+1);
    else m_fl[i] = Flavour(kf_none);
  }
  return shuffle;
}

void Cluster_Splitter::MakeNewMomenta(size_t shuffle) {    
  double mt2[2], alpha[2], beta[2];
  for (size_t i=0;i<2;i++) {
    mt2[i]    = (shuffle&(i+1) ? sqr(m_fl[i].Mass()) : m_mass2[i] ) + m_kt2;
  }
  alpha[0]    = ((m_Q2+mt2[0]-mt2[1])+sqrt(sqr(m_Q2+mt2[0]-mt2[1])-4.*m_Q2*mt2[0]))/(2.*m_Q2); 
  beta[0]     = mt2[0]/(m_Q2*alpha[0]);
  alpha[1]    = 1.-alpha[0];
  beta[1]     = 1.-beta[0];
  m_newmom[0] = m_E*(alpha[0]*s_AxisP + beta[0]*s_AxisM)+m_ktvec;
  m_newmom[1] = Vec4D(m_Q,0.,0.,0.)-m_newmom[0];
}

void Cluster_Splitter::FillHadronAndDeleteCluster(size_t i) {
  delete p_out[i];
  m_rotat.RotateBack(m_newmom[i]);
  m_boost.BoostBack(m_newmom[i]);
  p_softclusters->GetHadrons()->push_back(new Proto_Particle(m_fl[i],m_newmom[i],false));
}

void Cluster_Splitter::UpdateAndFillCluster(size_t i) {
  Poincare BoostIn(m_mom[i]);
  Poincare BoostOut(m_newmom[i]);
  //Vec4D check(0.,0.,0.,0.);
  for (size_t j=0;j<2;j++) {
    Vec4D partmom = (*p_out[i])[j]->Momentum(); 
    BoostIn.Boost(partmom);
    BoostOut.BoostBack(partmom);
    m_rotat.RotateBack(partmom);
    m_boost.BoostBack(partmom);
    //check += partmom;
    (*p_out[i])[j]->SetMomentum(partmom);
  }
  m_rotat.RotateBack(m_newmom[i]);
  m_boost.BoostBack(m_newmom[i]);
  p_out[i]->SetMomentum(m_newmom[i]);
  p_cluster_list->push_back(p_out[i]);
}

Cluster * Cluster_Splitter::MakeCluster(size_t i) {
  double lca   = (i==0? m_z[0]  : 1.-m_z[0] );
  double lcb   = (i==0? m_z[1]  : 1.-m_z[1] );
  double sign  = (i==0?    1. : -1.);
  double R02   = m_m2[i]+(m_popped_mass2+m_kt2);
  double ab    = 4.*m_m2[i]*(m_popped_mass2+m_kt2);
  double x = 1.;
  if (sqr(m_R2[i]-R02)>ab) {
    double centre = (m_R2[i]+m_m2[i]-(m_popped_mass2+m_kt2))/(2.*m_R2[i]);
    double lambda = Lambda(m_R2[i],m_m2[i],m_popped_mass2+m_kt2);
    x = (i==0)? centre+lambda : centre-lambda;
  }
  double y      = m_m2[i]/(x*m_R2[i]);
  // This is the overall cluster momentum - we do not need it - and its
  // individual components, i.e. the momenta of the Proto_Particles
  // it is made of.
  Vec4D newmom11 = (m_E*(     x*lca*s_AxisP+     y*(1.-lcb)*s_AxisM));
  Vec4D newmom12 = (m_E*((1.-x)*lca*s_AxisP+(1.-y)*(1.-lcb)*s_AxisM) +
		    sign * m_ktvec);
  Vec4D clumom = m_E*(lca*s_AxisP + (1.-lcb)*s_AxisM) + sign * m_ktvec;

  // back into lab system
  m_rotat.RotateBack(newmom11);
  m_boost.BoostBack(newmom11);
  m_rotat.RotateBack(newmom12);
  m_boost.BoostBack(newmom12);
  p_part[i]->SetMomentum(newmom11);

  Proto_Particle * newp =
    new Proto_Particle(m_newflav[i],newmom12,false,
		       p_part[0]->IsBeam()||p_part[1]->IsBeam());
  newp->SetKT2_Max(m_kt2);
  Cluster * cluster;
  if (i==0) cluster = new Cluster(p_part[0],newp);
  if (i==1) cluster = new Cluster(newp,p_part[1]);
  newp->SetGeneration(p_part[i]->Generation()+1);
  p_part[i]->SetGeneration(p_part[i]->Generation()+1);
  if (m_analyse) {
    if (m_Q>91.) {
      if (i==1) {
	m_histograms[string("kt_0")]->Insert(sqrt(m_kt2));
	m_histograms[string("z1_0")]->Insert(m_z[0]);
	m_histograms[string("z2_0")]->Insert(m_z[1]);
      }
      m_histograms[string("mass_0")]->Insert(sqrt(m_R2[i]));
      m_histograms[string("Rmass_0")]->Insert(2.*sqrt(m_R2[i]/m_Q2));
    }
    else {
      if (i==1) {
	m_histograms[string("kt")]->Insert(sqrt(m_kt2));
	m_histograms[string("z1")]->Insert(m_z[0]);
	m_histograms[string("z2")]->Insert(m_z[1]);
      }
      m_histograms[string("mass")]->Insert(sqrt(m_R2[i]));
      m_histograms[string("Rmass")]->Insert(2.*sqrt(m_R2[i]/m_Q2));
    }
  }
  return cluster;
}

