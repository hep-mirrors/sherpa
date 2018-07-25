#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

Primordial_KPerp::Primordial_KPerp() :
  m_analysis(false)
{ }
    
Primordial_KPerp::~Primordial_KPerp() {
  if (m_analysis) FinishAnalysis();
}

void Primordial_KPerp::Initialize() {
  for (size_t beam=0;beam<2;beam++) {
    string beamtag   = to_string(beam);
    m_form[beam]     = repars->GetKTForm(beam);
    m_mean[beam]     = (*repars)("kt_mean_"+beamtag);
    m_sigma[beam]    = (*repars)("kt_sigma_"+beamtag);
    m_Q2[beam]       = (*repars)("Q2_"+beamtag);
    m_ktmax[beam]    = (*repars)("kt_max_"+beamtag);
    m_eta[beam]      = (*repars)("eta_"+beamtag);
  }
  if (m_analysis) InitAnalysis();
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
    kt_tot += pmmit->second = scale * KT(pmmit->first->Momentum());
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

Vec4D Primordial_KPerp::KT(const Vec4D & mom) {
  double kt=0.,ktmax=mom[0];
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
  double cosphi = 2.*ran->Get()-1, sinphi = sqrt(1.-cosphi*cosphi);
  return kt * Vec4D(0.,cosphi,sinphi,0.);
}

double Primordial_KPerp::KT_Gauss(const double & ktmax) const {
  // Generate normalised Gaussian random numbers according to the
  // Marsaglia method
  double ran1, ran2, R;
  do {
    ran1 = 2.*ran->Get()-1.;
    ran2 = 2.*ran->Get()-1.;
    R    = ran1*ran1+ran2*ran2;
  } while (R>1. || R==0.);
  R  = sqrt(-2.*log(R)/R);
  // shift the Gaussian random numbers
  return m_mean[m_beam] + R * m_sigma[m_beam];
}

double Primordial_KPerp::KT_Gauss_Limited(const double & ktmax) const {
  // Generate normalised Gaussian random numbers according to the
  // Marsaglia method, with an additional polynomial limitation
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
  return Max(0.,pow(m_ktmax[m_beam],m_eta[m_beam])-pow(kt,m_eta[m_beam]));
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

