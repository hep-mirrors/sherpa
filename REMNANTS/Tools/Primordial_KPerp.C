#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

Primordial_KPerp::Primordial_KPerp() : m_analysis(false) { }

Primordial_KPerp::~Primordial_KPerp() {
  if (m_analysis) FinishAnalysis();
}

void Primordial_KPerp::Initialize(Remnant_Handler * rhandler) {
  msg_Out()<<"Initializing primordial transverse momentum ...\n";
  m_on = Settings::GetMainSettings()["INTRINSIC_KPERP"].Get<bool>();
  if (!m_on) {
    m_form = {primkT_form::none, primkT_form::none};
    return;
  }
  for (size_t beam=0;beam<2;beam++) {
    m_beamflav          = rhandler->GetRemnant(beam)->Flav();
    msg_Debugging()<<METHOD<<"("<<beam<<"): flav = "<<m_beamflav<<", "<<m_beamflav.Kfcode()<<"\n";
    m_form[beam]        = rempars->KT_Form(m_beamflav);
    m_recoil[beam]      = rempars->KT_Recoil(m_beamflav);
    if (m_form[beam]==primkT_form::none) continue;
    double escale       = pow(rpa->gen.Ecms()/rempars->Get(m_beamflav,"REFERENCE_ENERGY"),
			   rempars->Get(m_beamflav,"ENERGY_SCALING_EXPO"));
    if (m_form[beam]==primkT_form::gauss || m_form[beam]==primkT_form::gauss_limited) {
      m_SIMean[beam]    = rempars->Get(m_beamflav,"SHOWER_INITIATOR_MEAN") * escale;
      m_SISigma[beam]   = rempars->Get(m_beamflav,"SHOWER_INITIATOR_SIGMA") * escale;
      m_SpecMean[beam]  = rempars->Get(m_beamflav,"BEAM_SPECTATOR_MEAN");
      m_SpecSigma[beam] = rempars->Get(m_beamflav,"BEAM_SPECTATOR_SIGMA");
    } else if (m_form[beam]==primkT_form::dipole || m_form[beam]==primkT_form::dipole_limited) {
      m_SIQ2[beam]      = rempars->Get(m_beamflav,"SHOWER_INITIATOR_Q2") * escale;
      m_SpecQ2[beam]    = rempars->Get(m_beamflav,"BEAM_SPECTATOR_Q2");
    }
    if (m_form[beam]==primkT_form::gauss_limited || m_form[beam]==primkT_form::dipole_limited) {
      m_SIKtmax[beam]   = Max(0.0, rempars->Get(m_beamflav,"SHOWER_INITIATOR_KTMAX") * escale);
      m_SIEta[beam]     = rempars->Get(m_beamflav,"SHOWER_INITIATOR_KTEXPO");
      m_SpecKtmax[beam] = Max(0.0, rempars->Get(m_beamflav,"BEAM_SPECTATOR_KTMAX"));
      m_SpecEta[beam]   = rempars->Get(m_beamflav,"BEAM_SPECTATOR_KTEXPO");
    } else {
	m_SIKtmax[beam]    = m_SpecKtmax[beam] = 1000.;
	m_SIEta[beam]      = m_SpecEta[beam]   = 0.;
    }
  }
  if (m_analysis) InitAnalysis();
}

bool Primordial_KPerp::
CreateBreakupKinematics(const size_t & beam,ParticleMomMap * ktmap,const double & scale) {
  if (!m_on) return true;
  m_beam  = beam;
  p_ktmap = ktmap;
  Vec4D  kt_Show = Vec4D(0.,0.,0.,0.), kt_Spec = Vec4D(0.,0.,0.,0.);
  double E_Show  = 0., E_Spec = 0.;
  // harvesting particles from the beam blob of beam "beam" and
  for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
       pmmit!=p_ktmap->end();pmmit++) {
    if (pmmit->first->Momentum()[0]<0.)
      THROW(fatal_error,"particle with negative energy");
    pmmit->second = scale * KT(pmmit->first);
    if (pmmit->first->Info()=='I') {
      kt_Show += pmmit->second;
      E_Show  += pmmit->first->Momentum()[0];
    }
    else {
      kt_Spec += pmmit->second;
      E_Spec  += pmmit->first->Momentum()[0];
    }
  }
  if (m_analysis) {
    for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
	 pmmit!=p_ktmap->end();pmmit++)
      m_histos[std::string("KT_before")]->Insert(sqrt(dabs(pmmit->second.Abs2())));
  }
  // making sure the transverse momenta add up to 0.
  BalanceKT(kt_Show,E_Show,kt_Spec,E_Spec);
  return true;
}

void Primordial_KPerp::BalanceKT(const Vec4D & kt_Show,const double & E_Show,
				 const Vec4D & kt_Spec,const double & E_Spec) {
  // Taking the net kT in a single blob/beam break-up, kt_tot, and
  // distributing it in proportion to the particle energies.
  if (m_recoil[m_beam]==primkT_recoil::democratic) {
    Vec4D kt_tot = kt_Show + kt_Spec;
    double E_tot = E_Show + E_Spec;
    for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
	 pmmit!=p_ktmap->end();pmmit++) {
      pmmit->second = pmmit->second - pmmit->first->Momentum()[0]/E_tot * kt_tot;
    }
  }
  else if (m_recoil[m_beam]==primkT_recoil::beam_vs_shower) {
    for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
	 pmmit!=p_ktmap->end();pmmit++) {
      if (pmmit->first->Info()=='I') {
	pmmit->second = pmmit->second - pmmit->first->Momentum()[0]/E_Show * kt_Spec;
      }
      else {
	pmmit->second = pmmit->second - pmmit->first->Momentum()[0]/E_Spec * kt_Show;
      }
    }
  }
  if (m_analysis) {
    for (ParticleMomMap::iterator pmmit=p_ktmap->begin();
	 pmmit!=p_ktmap->end();pmmit++) {
      m_histos[std::string("KT_after")]->Insert(sqrt(dabs(pmmit->second.Abs2())));
    }
  }
}

Vec4D Primordial_KPerp::KT(const Particle * part,const double & ktext) {
  if (!m_on || m_form[m_beam]==primkT_form::none) return Vec4D(0.,0.,0.,0.);
  if (part->Info()=='I') {
    m_mean  = m_SIMean[m_beam];  m_sigma = m_SISigma[m_beam]; m_Q2 = m_SIQ2[m_beam];
    m_ktmax = m_SIKtmax[m_beam]; m_eta   = m_SIEta[m_beam];
  }
  else {
    m_mean  = m_SpecMean[m_beam];  m_sigma = m_SpecSigma[m_beam]; m_Q2 = m_SpecQ2[m_beam];
    m_ktmax = m_SpecKtmax[m_beam]; m_eta   = m_SpecEta[m_beam];
  }
  if (m_ktmax<=0.) return Vec4D(0.,0.,0.,0.);
  double ktmax = Min(m_ktmax,part->Momentum()[0]), kt = 0.;
  if (ktext>0.) ktmax = Min(ktmax,ktext);
  do {
    switch (m_form[m_beam]) {
    case primkT_form::none:           kt = 0.;                       break;
    case primkT_form::gauss:          kt = KT_Gauss(ktmax);          break;
    case primkT_form::gauss_limited:  kt = KT_Gauss_Limited(ktmax);  break;
    case primkT_form::dipole:         kt = KT_Dipole(ktmax);         break;
    case primkT_form::dipole_limited: kt = KT_Dipole_Limited(ktmax); break;
    default:
      THROW(fatal_error,"Unknown KPerp form.");
    }
  } while (kt<0. || kt>ktmax);
  // Add angle and construct the vector
  if (kt==0.) return Vec4D(0.,0.,0.,0.);
  const double phi = 2.*M_PI*ran->Get();
  return kt * Vec4D(0., cos(phi), sin(phi), 0.);
}

double Primordial_KPerp::KT_Gauss(const double & ktmax) const {
  double kt(0.);
  if (ktmax<1.e-3) return kt; // save to return no kt for small ktmax
  // too small ktmax can lead to an infinite loop due to the low prob. of generating small values
  // if ktmax is too small wrt. to sigma, the kt range is too narrow
  if (ktmax<(m_mean-3.*m_sigma) || ktmax<0.1 * m_sigma) return ktmax*ran->Get();
  do {
    kt = abs(m_mean + Sign(0.5-ran->Get())*m_sigma*sqrt(-log(std::max(1.e-5,ran->Get()))));
  } while (kt>ktmax);
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
  return 1./(1.+sqr(kt)/m_Q2);
}

double Primordial_KPerp::LimitedWeight(const double & kt) const {
  if (kt > m_ktmax) return 0.0;
  return 1.0 - pow(kt/m_ktmax, m_eta);
}

void Primordial_KPerp::InitAnalysis() {
  m_histos[std::string("KT_before")]  = new Histogram(0,0,20,200);
  m_histos[std::string("KT_after")]   = new Histogram(0,0,20,200);
}

void Primordial_KPerp::FinishAnalysis() {
  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator
	 hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = std::string("Remnant_Analysis/")+hit->first+std::string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}

