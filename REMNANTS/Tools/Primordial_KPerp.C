#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

Primordial_KPerp::Primordial_KPerp() :
  m_defform("gauss_limited"),
  m_defmean(0.0), m_defsigma(1.0), m_refE(7000.), m_scaleexpo(0.55),
  m_defQ2(0.77), m_defktmax(5.), m_defeta(5.),
  m_analysis(false)
{ }
    
Primordial_KPerp::~Primordial_KPerp() {
  if (m_analysis) FinishAnalysis();
}

void Primordial_KPerp::Initialize(const string path,const string file) {
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(path);
  dataread.SetInputFile(file);
  // Setting the mean and width of the Gaussian - we could discuss
  // more functional forms of the distribution here.  Default for leptons
  // is zero transverse momentum, only hadrons have a kT distribution of
  // remnants - again, this may have to change or become more
  // sophisticated for photons with hadronic structure.
  for (size_t beam=0;beam<2;beam++) {
    string formtag   = "K_PERP_FORM_"+to_string(beam);
    string meantag   = "K_PERP_MEAN_"+to_string(beam);
    string sigmatag  = "K_PERP_SIGMA_"+to_string(beam);
    string Q2tag     = "K_PERP_Q2_"+to_string(beam);
    string refEtag   = "K_PERP_REFE_"+to_string(beam);
    string expotag   = "K_PERP_SCALE_EXPO_"+to_string(beam);
    string ktmaxtag  = "K_PERP_MAX_"+to_string(beam);
    string ktexpotag = "K_PERP_CUT_EXPO_"+to_string(beam);
    double refE      = dataread.GetValue<double>(refEtag,m_refE);
    double expo      = dataread.GetValue<double>(expotag,m_scaleexpo);
    double mean      = dataread.GetValue<double>(meantag,m_defmean);
    double sigma     = dataread.GetValue<double>(sigmatag,m_defsigma);
    double Q2        = dataread.GetValue<double>(Q2tag,m_defQ2);
    m_form[beam]     = SelectForm(dataread.GetValue<string>(formtag,m_defform));
    m_mean[beam]     = mean;
    m_sigma[beam]    = sigma * pow((rpa->gen.Ecms()/refE),expo);
    m_Q2[beam]       = Q2 * pow((rpa->gen.Ecms()/refE),expo);
    m_ktmax[beam]    = dataread.GetValue<double>(ktmaxtag,m_defktmax);
    m_eta[beam]      = dataread.GetValue<double>(ktexpotag,m_defeta);
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

