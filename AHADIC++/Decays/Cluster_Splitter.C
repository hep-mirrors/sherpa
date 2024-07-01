#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


// mode 0: old mode
// mode 1: new mode, equivalent functionality to mode 0
// mode 2: new mode, not using z boundaries computed before (most likely broken)
#define AHADIC_CLUSTER_SPLITTER_MODE 1

// mode 0: usual fragmentation function
// mode 1: simplified fragmentation function, that can be integrated
// mode 2: ...
#define AHADIC_FRAGMENTATION_FUNCTION 0

Cluster_Splitter::Cluster_Splitter(list<Cluster *> * cluster_list,
				   Soft_Cluster_Handler * softclusters,
				   Flavour_Selector     * flavourselector,
				   KT_Selector          * ktselector) :
  Splitter_Base(cluster_list,softclusters,flavourselector,ktselector),
  m_output(false)
{
}

void Cluster_Splitter::Init() {
  Splitter_Base::Init();
  m_defmode  = hadpars->Switch("ClusterSplittingForm");
  m_beammode = hadpars->Switch("RemnantSplittingForm");

  m_alpha[0] = hadpars->GetVec("alphaL");
  m_beta[0]  = hadpars->GetVec("betaL");
  m_gamma[0] = hadpars->GetVec("gammaL");

  m_alpha[1] = hadpars->GetVec("alphaH");
  m_beta[1]  = hadpars->GetVec("betaH");
  m_gamma[1] = hadpars->GetVec("gammaH");

  m_alpha[2] = hadpars->GetVec("alphaD");
  m_beta[2]  = hadpars->GetVec("betaD");
  m_gamma[2] = hadpars->GetVec("gammaD");

  m_alpha[3] = hadpars->GetVec("alphaB");
  m_beta[3]  = hadpars->GetVec("betaB");
  m_gamma[3] = hadpars->GetVec("gammaB");

  const std::vector<double> _kt0s = hadpars->GetVec("kT_0");
  for (auto _kt0 : _kt0s)
    m_kt02.push_back(sqr(_kt0));

  // for the reweighting we need to find the min/max values of each
  // of the parameters
  for(int i{0}; i<4; ++i) {
    m_alpha_max[i] = *std::max_element(m_alpha[i].begin(), m_alpha[i].end());
    m_alpha_min[i] = *std::min_element(m_alpha[i].begin(), m_alpha[i].end());

    m_beta_max[i] = *std::max_element(m_beta[i].begin(), m_beta[i].end());
    m_beta_min[i] = *std::min_element(m_beta[i].begin(), m_beta[i].end());

    m_gamma_max[i] = *std::max_element(m_gamma[i].begin(), m_gamma[i].end());
    m_gamma_min[i] = *std::min_element(m_gamma[i].begin(), m_gamma[i].end());
  }


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
    m_type[i] = flcnt;
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
  for (size_t i=0;i<2;i++)
    m_m2min[i] = Min(m_minQ2[i],m_mdec2[i]);
  const double lambda = sqrt(sqr(m_Q2-m_m2min[0]-m_m2min[1])-
		       4.*(m_m2min[0]+m_kt2)*(m_m2min[1]+m_kt2));
  for (size_t i=0;i<2;i++) {
    const double centre = m_Q2-m_m2min[1-i]+m_m2min[i];
    m_zmin[i] = (centre-lambda)/(2.*m_Q2);
    m_zmax[i] = (centre+lambda)/(2.*m_Q2);
    m_mean[i]  = sqrt(m_kt02[0]);
    m_sigma[i] = sqrt(m_kt02[0]);
  }
}

bool Cluster_Splitter::MakeLongitudinalMomentaZ() {
  msg_Out() << "Got to a non-ported place" << std::endl;
  msg_Out() << "bool Cluster_Splitter::MakeLongitudinalMomentaZ()\n";
  size_t maxcounts=1000;
  while ((maxcounts--)>0) {
    if (MakeLongitudinalMomentaZSimple()) {
      double weight=1.;
      for (size_t i=0;i<2;i++) {
	if (m_gamma[i][0]>1.e-4) {
	  double DeltaM2 = m_R2[i]-m_minQ2[i];
	  weight *= DeltaM2>0.?exp(-m_gamma[i][0]*DeltaM2/m_sigma[i]):0.;
	}
      }
      if (weight>=ran->Get()) return true;
    }
  }
  return false;
}

bool Cluster_Splitter::MakeLongitudinalMomentaZSimple() {
  // todo: remove old cluster mode, and m_R2
  bool mustrecalc = false;

#if AHADIC_CLUSTER_SPLITTER_MODE == 0
  for (size_t i=0;i<2;i++) m_z[i]  = select_z(m_zmin[i],m_zmax[i],i);
  for (size_t i=0;i<2;i++) {
    m_R2[i] = m_z[i]*(1.-m_z[1-i])*m_Q2-m_kt2;
    if (m_R2[i]<m_mdec2[i]+m_kt2) {
      m_R2[i] = m_mdec2[i]+m_kt2;
      mustrecalc = true;
    }
  }
  bool ok = (m_R2[0]>m_mdec2[0]+m_kt2) && (m_R2[1]>m_mdec2[1]+m_kt2);
  return (ok && (mustrecalc?RecalculateZs():true));
#endif

  // order : lead > beam > rest
  //     -> 1 > 3 > rest
  //     -> stored in m_a[i] = flcnt;
  const int p0 = m_type[0];
  const int p1 = m_type[1];
  int i1{0}, i2{1};

  if(false) {
    if(p0 == p1) {
      if(ran->Get() < 0.5) {
	i1 = 0;
	i2 = 1;
      } else {
	i1 = 1;
	i2 = 0;
      }
    } else if (p0 == 1) {
      i1 = 0;
      i2 = 1;
    } else if (p1 == 1) {
      i1 = 1;
      i2 = 0;
    } else if (p0 == 3) {
      i1 = 0;
      i2 = 1;
    } else if (p1 == 3) {
      i1 = 1;
      i2 = 0;
    } else {
      if(ran->Get() < 0.5) {
	i1 = 0;
	i2 = 1;
      } else {
	i1 = 1;
	i2 = 0;
      }
    }
  }

  const double a0 = (m_mdec2[0]+2*m_kt2) / m_Q2;
  const double a1 = (m_mdec2[1]+2*m_kt2) / m_Q2;

  double p,q;
  if (i1 == 0) {
    p = (a1-a0-1);
    q = a0;
  } else if (i1 == 1) {
    p = (a0-a1-1);
    q = a1;
  }
  double _sqrt {p*p/4 - q};
  if(_sqrt < 0)
    return false;

  double lower = -p/2 - sqrt(p*p/4 - q);
  double upper = -p/2 + sqrt(p*p/4 - q);
  if(lower > upper)
    std::cout << "Needs fixing" << std::endl;

#if AHADIC_CLUSTER_SPLITTER_MODE == 1
  m_z[i1]  = select_z(std::max(m_zmin[i1],lower),
			 std::min(m_zmax[i1],upper),
			 i1);
  std::cout << "DEBUG: ARGS: "
	    << 0 << " "
	    << m_type[i1] << " "
	    << m_z[i1] << " "
	    << m_nsplit << " "
	    << m_zmin[i1] << " "
	    << m_zmax[i1] << " "
	    << lower << " " << upper << " "
	    << std::endl;
#endif
#if AHADIC_CLUSTER_SPLITTER_MODE == 2
  m_z[i1]  = select_z(0.,1.,i1);
#endif

  if(i1 == 0) {
    lower = a1/(1-m_z[i1]);
    upper = 1-a0/m_z[i1];
  } else {
    lower = a0/(1-m_z[i1]);
    upper = 1-a1/m_z[i1];
  }

#if AHADIC_CLUSTER_SPLITTER_MODE == 1
  m_z[i2]  = select_z(std::max(m_zmin[i2],lower),
			std::min(m_zmax[i2],upper),
			i2);
  std::cout << "DEBUG: ARGS: "
	    << 1 << " "
	    << m_type[i2] << " "
    	    << m_z[i2] << " "
    	    << m_nsplit << " "
	    << m_zmin[i2] << " " << m_zmax[i2] << " "
	    << lower << " " << upper << " "
	    << std::endl;
#endif
#if AHADIC_CLUSTER_SPLITTER_MODE == 2
  m_z[i2]  = select_z(0.,1.,i2);
#endif

  return true;
}

bool Cluster_Splitter::CheckKinematics() {
  for (size_t i=0;i<2;i++) {
    if(m_z[i] < m_zmin[i] || m_zmax[i] < m_z[i])
      return false;
    m_R2[i] = m_z[i]*(1.-m_z[1-i])*m_Q2-m_kt2;
    if (m_R2[i]<m_mdec2[i]+m_kt2)
      return false;
  }
  return true;
}

double Cluster_Splitter::FragmentationFunctionProb(double z, double zmin, double zmax,
						   double gamma, double kt02) {
  if(zmax - zmin < 0.05) {
    // std::cout << "Rejected: " << zmin << " " << zmax << " " << zmax - zmin
    // 	      << " " << m_accepted << " " << m_rejected
    // 	      << std::endl;
    m_rejected++;
    return 1;
  } else {
    m_accepted++;
  }
  // We just need to reweight the Fragmentation-Function witht the corresponding
  // integral
  auto f = [](double _z, double arg, double zmax) -> double {
    return (1-_z) * exp(-arg/_z);
  };
  auto F = [](double _z, double arg) -> double {
    // 1/2 (-e^(-g/x) x (-2 - g + x) + g (2 + g) Ei(-g/x))
    return 0.5*(-exp(-arg/_z)*_z*(-2-arg+_z) + arg*(2+arg)*std::expint(-arg/_z));
  };
  double arg    = gamma*(m_kt2+m_masses*m_masses)/kt02;

  //double arg    = gamma*(m_kt2)/kt02*10;
  //double arg    = gamma;
  const auto integral = (F(zmax,arg) - F(zmin,arg));
  return f(z,arg,zmax) / integral;
}

double Cluster_Splitter::FragmentationFunction(double z, double zmin, double zmax,
					       int cnt, int i_var) {
  const auto type {m_type[cnt]};
  const double alpha = m_alpha[type][i_var];
  const double beta  = m_beta [type][i_var];
  const double gamma = m_gamma[type][i_var];
  const double kt02  = m_kt02[i_var];

#if AHADIC_FRAGMENTATION_FUNCTION == 0
  double arg, norm {1.}, value {1.};
  // if (alpha>=0.)
  //   norm *= pow(zmax,m_alpha_max[type]);
  // else
  //   norm *= pow(zmin,m_alpha_max[type]);
  if (alpha>=0.)
    norm *= pow(zmax,alpha);
  else
    norm *= pow(zmin,alpha);

  // if (beta>=0.)
  //   norm *= pow(1.-zmin,m_beta_max[type]);
  // else
  //   norm *= pow(1.-zmax,m_beta_max[type]);
  if (beta>=0.)
    norm *= pow(1.-zmin,beta);
  else
    norm *= pow(1.-zmax,beta);

  double wt = pow(z,alpha) * pow(1.-z,beta);
  value = wt/norm;

  if (m_mode==2) {
    const auto gamma_min {gamma};
    arg    = dabs(gamma)>5.e-3 ? (m_kt2+m_masses*m_masses)/kt02 : 0.;
    value *= exp(-arg*((zmax*gamma-z*gamma_min)/(z*zmax)));

    // irrelevant
    norm  *= exp(-arg*gamma_min/zmax);
    wt    *= exp(-arg*gamma_min/z);
  }

  if (wt>norm) {
    std::cout << "z = " << z << std::endl;
    std::cout << "Something is wrong in the Cluster Fragmentation function"
	      << std::endl;
    exit(1);
  }
  //value /= 10.;
  return value;
# endif

#if AHADIC_FRAGMENTATION_FUNCTION == 1
  auto f = [](double _z, double arg, double zmax) -> double {
    return (1-_z) * exp(-arg*((zmax-_z)/(_z*zmax)));
  };

  double arg    = gamma*(m_kt2+m_masses*m_masses)/kt02;
  //double arg    = gamma*(m_kt2)/kt02*10;
  //double arg    = gamma;
  // compute max
  double norm {1-zmin};
  return f(z, arg,zmax) / norm;
#endif
}

double Cluster_Splitter::
WeightFunction(const double & z,const double & zmin,const double & zmax,
	       const unsigned int & cnt) {
  const auto value = FragmentationFunction(z,zmin,zmax,cnt,0);
  return value;
}

void Cluster_Splitter::z_rejected(const double wgt, const double & z,
				  const double & zmin,const double & zmax,
				  const unsigned int & cnt) {
#if AHADIC_FRAGMENTATION_FUNCTION == 1
  return;
#endif
  const auto type = m_type[cnt];
  for (int i{0}; i<m_alpha[0].size(); i++) {
    const auto a    = m_alpha[type][i];
    const auto b    = m_beta [type][i];
    const auto c    = m_gamma[type][i];
    const auto kt   = m_kt02[i];
    const auto wgt_new = FragmentationFunction(z,zmin,zmax,cnt,i);
    tmp_variation_weights[i] *= (1.-wgt_new) / (1.-wgt);
  }
}

void Cluster_Splitter::z_accepted(const double wgt, const double & z,
				  const double & zmin,const double & zmax,
				  const unsigned int & cnt) {
  const auto type = m_type[cnt];
#if AHADIC_FRAGMENTATION_FUNCTION == 1
  const double wgt_old = FragmentationFunctionProb(z,zmin,zmax,m_gamma[type][0],m_kt02[0]);
#else
  const double wgt_old = wgt;
#endif
  for (int i{0}; i<m_alpha[0].size(); i++) {
    const auto a    = m_alpha[type][i];
    const auto b    = m_beta [type][i];
    const auto c    = m_gamma[type][i];
    const auto kt   = m_kt02[i];
#if AHADIC_FRAGMENTATION_FUNCTION == 1
    const auto wgt_new = FragmentationFunctionProb(z,zmin,zmax,c,kt);
#else
    const auto wgt_new = FragmentationFunction(z,zmin,zmax,cnt,i);
#endif
    // TODO: find out why this becomes non from time to time
    const auto frac = wgt_new / wgt_old;
    std::cout << "DEBUG: WSLS: " << type << " "
	      << m_nsplit << " "
	      << i << " "
	      << wgt_new << " "
	      << wgt_old << " "
	      << frac << std::endl;

    if(!std::isnan(frac))
      tmp_variation_weights[i] *= frac;
  }
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
  msg_Out() << "Got to a non-ported place" << std::endl;
  msg_Out() << "Cluster_Splitter::MakeLongitudinalMomentaMassSimple()" << std::endl;
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
  msg_Out() << "Got to a non-ported place" << std::endl;
  msg_Out() << "bool Cluster_Splitter::MakeLongitudinalMomentaMass()" << std::endl;
  size_t maxcounts=1000;
  while ((maxcounts--)>0) {
    if (MakeLongitudinalMomentaMassSimple()) {
      double weight=1.;
      for (size_t i=0;i<2;i++) {
	if (m_alpha[i][0]>1.e-4) weight *= pow(m_z[i],m_alpha[i][0]);
	if (m_beta[i][0]>1.e-4)  weight *= pow(1.-m_z[i],m_beta[i][0]);
      }
      if (weight>=ran->Get()) return true;
    }
  }
  return false;
}

double Cluster_Splitter::DeltaM(const size_t & cl) {
  msg_Out() << "Got to a non-ported place" << std::endl;
  msg_Out() << "Cluster_Splitter::DeltaM\n";
  double deltaM, deltaMmax = m_Q-sqrt(m_m2min[0])-sqrt(m_m2min[1]);
  double mean =  m_mean[cl], sigma = 1./(m_type[cl] * sqrt(m_kt02[0]));
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
  cluster->m_nsplit = m_nsplit + 1;
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

