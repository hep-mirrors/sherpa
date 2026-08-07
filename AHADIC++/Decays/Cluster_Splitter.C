#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <limits>
#include <string>

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


// mode 0: old mode
// mode 1: new mode, equivalent functionality to mode 0
// mode 2: new mode, not using z boundaries computed before (most likely broken)
#define AHADIC_CLUSTER_SPLITTER_MODE 1

Cluster_Splitter::Cluster_Splitter(list<Cluster *> * cluster_list,
				   Soft_Cluster_Handler * softclusters,
				   Flavour_Selector     * flavourselector,
				   KT_Selector          * ktselector,
				   Hadronisation_Reweighting   * reweighting) :
  Splitter_Base(cluster_list,softclusters,flavourselector,ktselector,
		reweighting)
{
}

void Cluster_Splitter::Init() {
  Splitter_Base::Init();
  m_defmode  = hadpars->Switch("ClusterSplittingForm");
  m_beammode = hadpars->Switch("RemnantSplittingForm");
  if (p_reweighting->Active() && (m_defmode != 2 || m_beammode != 2)) {
    THROW(fatal_error, std::string("Reweighting of AHADIC only ported for cluster splitting mode 2.\n")
                      + "Found CLUSTER_SPLITTING_MODE = " + std::to_string(m_defmode)
                      + ", REMNANT_CLUSTER_MODE = " + std::to_string(m_beammode) + ".\n"
                      + "Please adjust your settings.");
  }
  m_n_cluster_variations = p_reweighting->NumberOfClusterVariations();

  m_alpha[0] = p_reweighting->GetVariationVector("alphaL");
  m_beta[0]  = p_reweighting->GetVariationVector("betaL");
  m_gamma[0] = p_reweighting->GetVariationVector("gammaL");

  m_alpha[1] = p_reweighting->GetVariationVector("alphaH");
  m_beta[1]  = p_reweighting->GetVariationVector("betaH");
  m_gamma[1] = p_reweighting->GetVariationVector("gammaH");

  m_alpha[2] = p_reweighting->GetVariationVector("alphaD");
  m_beta[2]  = p_reweighting->GetVariationVector("betaD");
  m_gamma[2] = p_reweighting->GetVariationVector("gammaD");

  m_alpha[3] = p_reweighting->GetVariationVector("alphaB");
  m_beta[3]  = p_reweighting->GetVariationVector("betaB");
  m_gamma[3] = p_reweighting->GetVariationVector("gammaB");

  const std::vector<double> _kt0s = p_reweighting->GetVariationVector("kT_0");
  m_kt02.clear();
  m_kt02.reserve(_kt0s.size());
  for (auto _kt0 : _kt0s)
    m_kt02.push_back(sqr(_kt0));

  m_cvals.resize(m_n_cluster_variations);
  m_logprobs.resize(m_n_cluster_variations);
  m_probs.resize(m_n_cluster_variations);

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
  if (!CalculateLimits()) return false;
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

bool Cluster_Splitter::CalculateLimits() {
  // Masses from Splitter_Base:
  // - constitutents:
  //   m_mass[0,1] and m_m2[0,1] = sqr(m_mass[0,1]), m_popped_mass, m_popped_mass2
  //   m_msum[0,1] = mass+mass for the pairs, m_msum2[0,1] = sqr(m_msum[0,1])
  // - hadrons:
  //   m_minQ[0,1] is lightest single or double transition (double for di-di pairs)
  //   m_mdec is lightest decay transition
  for (size_t i=0;i<2;i++)
    m_m2min[i] = Min(m_minQ2[i],m_mdec2[i]);
  const double arg = sqr(m_Q2-m_m2min[0]-m_m2min[1])-
                     4.*(m_m2min[0]+m_kt2)*(m_m2min[1]+m_kt2);
  if (arg<0.) return false;
  const double lambda = sqrt(arg);
  for (size_t i=0;i<2;i++) {
    const double centre = m_Q2-m_m2min[1-i]+m_m2min[i];
    m_zmin[i] = (centre-lambda)/(2.*m_Q2);
    m_zmax[i] = (centre+lambda)/(2.*m_Q2);
    m_mean[i]  = sqrt(m_kt02[0]);
    m_sigma[i] = sqrt(m_kt02[0]);
  }
  return true;
}

bool Cluster_Splitter::MakeLongitudinalMomentaZ() {
  msg_Error() << "Got to a non-ported place: "
              << "bool Cluster_Splitter::MakeLongitudinalMomentaZ()\n";
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

  for (size_t i=0;i<2;i++) {
    m_z[i] = SelectZ(m_zmin[i],m_zmax[i],i);
    if (m_z[i] < 0.) return false;
  }
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
    msg_Error() << "Inconsistent z bounds: lower > upper in MakeLongitudinalMomentaZSimple\n";
#if AHADIC_CLUSTER_SPLITTER_MODE == 1
  m_z[i1] = SelectZ(std::max(m_zmin[i1],lower), std::min(m_zmax[i1],upper), i1);
#endif
#if AHADIC_CLUSTER_SPLITTER_MODE == 2
  m_z[i1] = SelectZ(0.,1.,i1);
#endif
  if (m_z[i1] < 0.) return false;

  if(i1 == 0) {
    lower = a1/(1-m_z[i1]);
    upper = 1-a0/m_z[i1];
  } else {
    lower = a0/(1-m_z[i1]);
    upper = 1-a1/m_z[i1];
  }

#if AHADIC_CLUSTER_SPLITTER_MODE == 1
  m_z[i2] = SelectZ(std::max(m_zmin[i2],lower), std::min(m_zmax[i2],upper), i2);
#endif
#if AHADIC_CLUSTER_SPLITTER_MODE == 2
  m_z[i2] = SelectZ(0.,1.,i2);
#endif
  if (m_z[i2] < 0.) return false;

  p_reweighting->RecordClusterZ(m_nsplit, m_type[0], m_z[0],
                                m_type[1], m_z[1], m_Q); // OUTPUT
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

double Cluster_Splitter::FragExponent(const double gamma, const double kt02,
				      const double scale) const {
  return Frag_Norm::Exponent(gamma,kt02,scale);
}

bool Cluster_Splitter::FillLogDensities(const double z, const double zmin,
					const double zmax,
					const unsigned int cnt) {
  // The density of the z accepted by SelectZ() is f(z)/int f, irrespective of
  // the flat proposal and of the rejected trials, see ZAccepted.  Fill the log
  // of it for every variation.
  if (!(zmin>0.) || !(zmax<1.) || !(z>zmin) || !(z<zmax)) return false;
  const int t = m_type[cnt];
  const double scale = m_kt2 + m_masses*m_masses;
  // The exponents of all variations are needed before the nodes can be placed:
  // the rule has to resolve exp(-c/z) for the most strongly peaked of them
  // while staying common to all of them, see Frag_Norm.
  double cmin = std::numeric_limits<double>::max(), cmax = 0.;
  size_t ipeak = 0;
  for (size_t ivar=0; ivar<m_n_cluster_variations; ++ivar) {
    m_cvals[ivar] = FragExponent(m_gamma[t][ivar],m_kt02[ivar],scale);
    const double ac = dabs(m_cvals[ivar]);
    if (ac<cmin) cmin = ac;
    if (ac>cmax) { cmax = ac; ipeak = ivar; }
  }
  m_fragnorm.SetRange(zmin,zmax,m_alpha[t][ipeak],m_beta[t][ipeak],cmin,cmax);
  const double logz = std::log(z), log1mz = std::log1p(-z), invz = 1./z;
  for (size_t ivar=0; ivar<m_n_cluster_variations; ++ivar) {
    const double alpha = m_alpha[t][ivar], beta = m_beta[t][ivar];
    const double c     = m_cvals[ivar];
    m_logprobs[ivar] = (alpha*logz + beta*log1mz - c*invz)
                     - m_fragnorm(alpha,beta,c);
  }
  return true;
}

void Cluster_Splitter::FillProbs(const double wgt, const double z,
				 const double zmin, const double zmax,
				 const unsigned int cnt) {
  m_probs[0] = wgt;
  for (size_t ivar=1; ivar<m_n_cluster_variations; ++ivar)
    m_probs[ivar] = FragmentationFunction(z,zmin,zmax,cnt,ivar);
}

double Cluster_Splitter::FragmentationFunction(double z, double zmin, double zmax,
					       int cnt, int i_var) {
  const auto type {m_type[cnt]};
  return FragmentationFunction(z, zmin, zmax,
			       m_alpha[type][i_var], m_beta[type][i_var],
			       m_gamma[type][i_var], m_kt02[i_var]);
}

double Cluster_Splitter::FragmentationFunction(double z, double zmin, double zmax,
					       double alpha, double beta,
					       double gamma, double kt02) {
  // This is the sampling path: it defines the accept/reject decisions of
  // SelectZ and hence the random number sequence.  Do not rewrite the
  // expressions below - not even into the mathematically equivalent
  // exp(alpha*log(z)+...) of LogNorm() - or the nominal event sample changes.
  // The one exception is the guarded branch at the end, which is only taken
  // where the expressions below have no meaning at all; see there.
  if (m_mode == 2) {
    const double c = dabs(gamma) > 5.e-3
      ? gamma * (m_kt2 + m_masses*m_masses) / kt02
      : 0.;
    // The maximum of the fragmentation function on [zmin,zmax] is attained at
    // one of the interval ends or at an interior stationary point of
    // alpha ln z + beta ln(1-z) - c/z, i.e. a root of
    // (alpha+beta) z^2 - (alpha-c) z - c = 0. Collect those candidates once so
    // that both evaluations below use exactly the same set.
    double probe[4] = { zmin, zmax, 0., 0. };
    size_t nprobe = 2;
    const double A = alpha + beta;
    if (std::abs(A) > 1e-10) {
      const double disc = sqr(alpha - c) + 4. * A * c;
      if (disc >= 0.) {
        for (const double sign : {1., -1.}) {
          const double z_crit = ((alpha - c) + sign * sqrt(disc)) / (2. * A);
          if (z_crit > zmin && z_crit < zmax)
            probe[nprobe++] = z_crit;
        }
      }
    } else if (std::abs(alpha - c) > 1e-10) {
      const double z_crit = c / (c - alpha);
      if (z_crit > zmin && z_crit < zmax)
        probe[nprobe++] = z_crit;
    }
    auto g = [&](double _z) {
      return pow(_z, alpha) * pow(1.-_z, beta) * exp(-c / _z);
    };
    double norm = std::max(g(probe[0]), g(probe[1]));
    for (size_t i=2; i<nprobe; ++i) norm = std::max(norm, g(probe[i]));
    // Take the direct route only while norm is a normal, finite, positive
    // number, i.e. while g(z)/norm carries full precision. Everything else -
    // zero, subnormal, infinite, NaN - goes through the log-space branch.
    if (norm >= std::numeric_limits<double>::min() && std::isfinite(norm))
      return std::min(1.0, g(z) / norm);

    ///////////////////////////////////////////////////////////////////////////
    // g is built as a product, so once c/z exceeds about 745 its exp(-c/z)
    // factor underflows to zero at every probe point and norm becomes zero (or
    // subnormal, which is just as useless). g(z)/norm is then 0/0 or inf, and
    // because std::min(1.0,NaN) and std::min(1.0,inf) both return 1.0 every
    // trial would be accepted: SelectZ would hand back a z drawn uniformly on
    // [zmin,zmax] while the reweighting keeps using the true, sharply peaked
    // density. That mismatch is not a small effect - a single such draw
    // produced an event weight of 2e15 - and it also means the nominal z was
    // not distributed according to the fragmentation function.
    //
    // The same degeneracy exists in the opposite direction: for gamma < 0 the
    // factor becomes exp(+|c|/z) and overflows, giving norm = inf and again
    // inf/inf = NaN. Both are covered by the guard above.
    //
    // Repeat the identical comparison with the difference taken in log space,
    // where nothing can underflow. This branch is reached only when the
    // expression above is meaningless, so the accept/reject decisions, and
    // with them the random number sequence, are untouched wherever it was well
    // defined. Measured incidence: 1 in 6.7e6 z draws.
    ///////////////////////////////////////////////////////////////////////////
    auto lg = [&](double _z) {
      return alpha*std::log(_z) + beta*std::log1p(-_z) - c/_z;
    };
    double lognorm = std::max(lg(probe[0]), lg(probe[1]));
    for (size_t i=2; i<nprobe; ++i) lognorm = std::max(lognorm, lg(probe[i]));
    return std::min(1.0, std::exp(lg(z) - lognorm));
  }

  // Note: the branch below has no exp(-c/z) factor, so f can only underflow
  // for z below ~1e-123, which the kinematics cannot reach. It therefore needs
  // no equivalent guard.

  // f(z) = z^alpha * (1-z)^beta
  // interior mode from d/dz[log f] = alpha/z - beta/(1-z) = 0 => z* = alpha/(alpha+beta)
  auto f = [&](double _z) {
    return pow(_z, alpha) * pow(1.-_z, beta);
  };
  double norm = std::max(f(zmin), f(zmax));
  if (std::abs(alpha + beta) > 1e-10) {
    const double z_mode = alpha / (alpha + beta);
    if (z_mode > zmin && z_mode < zmax)
      norm = std::max(norm, f(z_mode));
  }
  return std::min(1.0, f(z) / norm);
}

double Cluster_Splitter::
WeightFunction(const double & z,const double & zmin,const double & zmax,
	       const unsigned int & cnt) {
  return FragmentationFunction(z, zmin, zmax, cnt, 0);
}

void Cluster_Splitter::ZAccepted(const double wgt, const double & z,
				 const double & zmin,const double & zmax,
				 const unsigned int & cnt) {
  // The value accepted by an accept/reject loop is distributed according to
  // f(z)/int f, independent of the proposal and of the number of rejected
  // trials.  Averaging the trial-by-trial weight of the whole accept/reject
  // chain over the rejections gives exactly the ratio of these normalised
  // densities, so using it directly is unbiased and has a strictly smaller
  // variance.  This mirrors what KT_Selector and Gluon_Splitter do.
  if(!p_reweighting->DoClusterSplittingReweighting(m_nsplit)) return;
  if (p_reweighting->ARFragReweighting()) {
    FillProbs(wgt,z,zmin,zmax,cnt);
    p_reweighting->ClusterSplittingReweightingAR(true,m_probs);
    return;
  }
  if(!FillLogDensities(z,zmin,zmax,cnt)) return;
  p_reweighting->ClusterSplittingReweighting(m_logprobs);
}

void Cluster_Splitter::ZRejected(const double wgt, const double & z,
				 const double & zmin,const double & zmax,
				 const unsigned int & cnt) {
  if(!p_reweighting->ARFragReweighting()) return;
  if(!p_reweighting->DoClusterSplittingReweighting(m_nsplit)) return;
  FillProbs(wgt,z,zmin,zmax,cnt);
  p_reweighting->ClusterSplittingReweightingAR(false,m_probs);
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
  msg_Error() << "Got to a non-ported place: "
              << "Cluster_Splitter::MakeLongitudinalMomentaMassSimple()\n";
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
  msg_Error() << "Got to a non-ported place: "
              << "bool Cluster_Splitter::MakeLongitudinalMomentaMass()\n";
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
  msg_Error() << "Got to a non-ported place: Cluster_Splitter::DeltaM\n";
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
  for (size_t j=0;j<2;j++) {
    Vec4D partmom = (*p_out[i])[j]->Momentum();
    BoostIn.Boost(partmom);
    BoostOut.BoostBack(partmom);
    m_rotat.RotateBack(partmom);
    m_boost.BoostBack(partmom);
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

