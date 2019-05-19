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
    const auto e = pow(rpa->gen.Ecms()/refEs[beam], expos[beam]);
    m_form[beam]     = SelectForm(forms[beam]);
    m_mean[beam]     = means[beam];
    m_sigma[beam]    = sigmas[beam] * e;
    m_Q2[beam]       = Q2s[beam] * e;
    m_ktmax[beam]
      = Max(1.0, ktmaxs[beam] * pow((rpa->gen.Ecms()/refEs[beam]), expos[beam]));
    m_eta[beam]      = ktexpos[beam];
  }
  if (m_analysis) InitAnalysis();
}

prim_kperp_form::code Primordial_KPerp::SelectForm(const std::string & form) {
  prim_kperp_form::code pkf = prim_kperp_form::undefined;
  if      (form=="gauss")          pkf = prim_kperp_form::gauss;
  else if (form=="gauss_limited")  pkf = prim_kperp_form::gauss_limited;
  else if (form=="dipole")         pkf = prim_kperp_form::dipole;
  else if (form=="dipole_limited") pkf = prim_kperp_form::dipole_limited;
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
      msg_Out()<<(*pmmit->first)<<"\n";
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
      msg_Error()<<METHOD<<": kperp form undefined.  Exit the run.\n";
      exit(1);
    }
  } while (kt<0. || kt>ktmax);
  // Add angles and construct the vector
  double cosphi = cos(2.*M_PI*ran->Get()), sinphi = sqrt(1.-cosphi*cosphi);
  return kt * Vec4D(0.,cosphi,sinphi,0.);
}

double Primordial_KPerp::KT_Gauss(const double & ktmax) const {
//  double kt=m_sigma[m_beam]*sqrt(log(1-ran->Get()*(1-exp(-pow(m_ktmax[m_beam]/m_sigma[m_beam],2)))));
  double ran1=std::max(0.00001,ran->Get());
  double kt=m_sigma[m_beam]*sqrt(-log(ran1));
  return kt;
}

double Primordial_KPerp::KT_Gauss_Limited(const double & ktmax) const {
  // Generate normalised Gaussian random numbers
  // with an additional polynomial limitation
  double ran1, ran2, R, kt;
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
  if (kt < m_ktmax[m_beam])
    return 1.0 - pow(kt/m_ktmax[m_beam], m_eta[m_beam]);
  else
    return 0.0;
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

