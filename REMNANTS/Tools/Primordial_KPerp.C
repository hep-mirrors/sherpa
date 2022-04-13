#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

Primordial_KPerp::Primordial_KPerp() :
  m_defform("gauss_limited"),
  m_defmean(0.0), m_defsigma(1.5), m_refE(7000.), m_scaleexpo(0.08),
  m_defQ2(0.77), m_defktmax(3.0), m_defeta(5.),
  m_analysis(false)
{ }
    
Primordial_KPerp::~Primordial_KPerp() {
  if (m_analysis) FinishAnalysis();
}

void Primordial_KPerp::Initialize() {
  // Setting the mean and width of the Gaussian - we could discuss
  // more functional forms of the distribution here.  Default for leptons
  // is zero transverse momentum, only hadrons have a kT distribution of
  // remnants - again, this may have to change or become more
  // sophisticated for photons with hadronic structure.
  // Note: To play it safe, we have set the maximal intrinsic kperp that we
  // use for the hard-wired limitation of Gauss and Dipole factor to at
  // least 1 GeV for each (hadronic) beam.
  auto s = Settings::GetMainSettings()["INTRINSIC_KPERP"];
  auto forms   = s["FORM"]      .SetDefault(m_defform)  .GetTwoVector<string>();
  auto means   = s["MEAN"]      .SetDefault(m_defmean)  .GetTwoVector<double>();
  auto sigmas  = s["SIGMA"]     .SetDefault(m_defsigma) .GetTwoVector<double>();
  auto Q2s     = s["Q2"]        .SetDefault(m_defQ2)    .GetTwoVector<double>();
  auto ktmaxs  = s["MAX"]       .SetDefault(m_defktmax) .GetTwoVector<double>();
  auto refEs   = s["REFE"]      .SetDefault(m_refE)     .GetTwoVector<double>();
  auto expos   = s["SCALE_EXPO"].SetDefault(m_scaleexpo).GetTwoVector<double>();
  auto ktexpos = s["CUT_EXPO"]  .SetDefault(m_defeta)   .GetTwoVector<double>();
  for (size_t beam=0;beam<2;beam++) {
    const auto escale = pow(rpa->gen.Ecms()/refEs[beam], expos[beam]);
    m_form[beam]      = SelectForm(forms[beam]);
    m_mean[beam]      = means[beam];
    m_sigma[beam]     = sigmas[beam] * escale;
    m_Q2[beam]        = Q2s[beam] * escale;
    m_ktmax[beam]     = Max(1.0, ktmaxs[beam] * escale);
    m_eta[beam]       = ktexpos[beam];
  }
  if (m_analysis) InitAnalysis();
}

prim_kperp_form::code Primordial_KPerp::SelectForm(const std::string & form) {
  prim_kperp_form::code pkf = prim_kperp_form::undefined;
  if (form=="None" || form=="none")pkf = prim_kperp_form::none;
  else if (form=="gauss")          pkf = prim_kperp_form::gauss;
  else if (form=="gauss_limited")  pkf = prim_kperp_form::gauss_limited;
  else if (form=="dipole")         pkf = prim_kperp_form::dipole;
  else if (form=="dipole_limited") pkf = prim_kperp_form::dipole_limited;
  else THROW(not_implemented,"Intrinsic KPerp model not implemented.");
  return pkf;
}

bool Primordial_KPerp::
CreateBreakupKinematics(const size_t & beam,ParticleMomMap * ktmap,const double & scale) {
  m_beam       = beam;
  p_ktmap      = ktmap;
  Vec4D kt_tot = Vec4D(0.,0.,0.,0.);
  double E_tot = 0.;
  // harvesting particles from the beam blob of beam "beam" and
  for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
       pmmit!=p_ktmap->end();pmmit++) {
    if (pmmit->first->Momentum()[0]<0.) {
      msg_Error()<<"Error in "<<METHOD<<" found particle with negative energy:\n"
		 <<(*pmmit->first)<<"\n"
		 <<"Will exit the run, please notify the authors.\n";
      exit(1);
    }
    kt_tot += pmmit->second =
      scale * KT(Min(m_ktmax[m_beam],pmmit->first->Momentum()[0]));
    E_tot  += pmmit->first->Momentum()[0];
    if (m_analysis) m_histos[string("KT_before")]->Insert(sqrt(dabs(pmmit->second.Abs2())));
  }
  // making sure the transverse momenta add up to 0.
  BalanceKT(kt_tot,E_tot);
  return true;
}

void Primordial_KPerp::BalanceKT(const Vec4D & kt_tot,const double & E_tot) {
  // Taking the net kT in a single blob/beam break-up, kt_tot, and
  // distributing it in proportion to the particle energies.
  for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
       pmmit!=p_ktmap->end();pmmit++) {
    pmmit->second = pmmit->second - pmmit->first->Momentum()[0]/E_tot * kt_tot;
    if (m_analysis) m_histos[string("KT_after")]->Insert(sqrt(dabs(pmmit->second.Abs2())));
  }
}

Vec4D Primordial_KPerp::KT(const double & ktmax) {
  double kt=0.;
  do {
    switch (m_form[m_beam]) {
    case prim_kperp_form::none:
      kt = 0.;
      break;
    case prim_kperp_form::gauss:
      kt = KT_Gauss(ktmax);
      break;
    case prim_kperp_form::gauss_limited:
      kt = KT_Gauss_Limited(ktmax);
      break;
    case prim_kperp_form::dipole:
      kt = KT_Dipole(ktmax);
      break;
    case prim_kperp_form::dipole_limited:
      kt = KT_Dipole_Limited(ktmax);
      break;
    default:
      THROW(fatal_error,"Unknown KPerp form.");
    }
  } while (kt<0. || kt>ktmax);
  // Add angle and construct the vector
  if (kt==0.) return Vec4D(0.,0.,0.,0.);
  const auto phi = 2*M_PI*ran->Get();
  return kt * Vec4D(0., cos(phi), sin(phi), 0.);
}

double Primordial_KPerp::KT_Gauss(const double & ktmax) const {
  double kt(0.);
  if (ktmax<1.e-3) return kt; // save to return no kt for small ktmax
  // too small ktmax can lead to an infinite loop due to the low prob. of generating small values
  if ((ktmax<m_mean[m_beam] - 2.*m_sigma[m_beam])) return ktmax*ran->Get();
  // if ktmax is too small wrt. to sigma, the kt range is too narrow
  if (ktmax<0.1 * m_sigma[m_beam]) return ktmax*ran->Get();
  do { kt = abs(m_mean[m_beam]+Sign(0.5-ran->Get())*m_sigma[m_beam]*sqrt(-log(std::max(1.e-5,ran->Get()))));}
  while (kt>ktmax);
  return kt;
}

double Primordial_KPerp::KT_Gauss_Limited(const double & ktmax) const {
  // Generate normalised Gaussian random numbers
  // with an additional polynomial limitation
  double kt;
  do {
    kt = KT_Gauss(ktmax);
  } while (LimitedWeight(kt)<ran->Get());
  return kt;
}

double Primordial_KPerp::KT_Dipole(const double & ktmax) const {
  // Dipole form factor
  double kt;
  do {
    kt = ran->Get()*ktmax;
  } while (DipoleWeight(kt)<ran->Get());
  return kt;
}

double Primordial_KPerp::KT_Dipole_Limited(const double & ktmax) const {
  // Dipole form factor, with an additional polynomial limitation
  double kt;
  do {
    kt = ran->Get()*ktmax;
  } while (DipoleWeight(kt)*LimitedWeight(kt)<ran->Get());
  return kt;
}

double Primordial_KPerp::DipoleWeight(const double & kt) const {
  return 1./(1.+sqr(kt)/m_Q2[m_beam]);
}

double Primordial_KPerp::LimitedWeight(const double & kt) const {
  if (kt > m_ktmax[m_beam]) return 0.0;
  return 1.0 - pow(kt/m_ktmax[m_beam], m_eta[m_beam]);
}

void Primordial_KPerp::InitAnalysis() {
  m_histos[string("KT_before")]  = new Histogram(0,0,20,200);
  m_histos[string("KT_after")]   = new Histogram(0,0,20,200);
}

void Primordial_KPerp::FinishAnalysis() {
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator 
	 hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = string("Remnant_Analysis/")+hit->first+string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}

