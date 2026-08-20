#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/FormFactors/Propagator.H"
#include "METOOLS/FormFactors/Line_Shapes.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/SM/FormFactorParameters.H"   // new parameter system

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include <set>
#include <sstream>
#include <iomanip>

namespace METOOLS {
  enum FF_0_PP_mode {
    pipi_plus       = 1,
    KK_plus         = 2,
    Kpi_plus        = 11,
    pipi_0          = 101,
    KK_0            = 102,
    pipipi_0        = 201,
    V_pi            = 301,
    unknown         = 999
  };

  inline std::ostream & operator<<(std::ostream & str,const FF_0_PP_mode & m) {
    switch (m) {
    case pipi_plus: return str<<"pipi_plus";
    case KK_plus  : return str<<"KK_plus";
    case Kpi_plus : return str<<"Kpi_plus";
    case pipi_0   : return str<<"pipi_0 (gamma* -> pi+ pi-)";
    case KK_0     : return str<<"KK_0 (gamma* -> K+ K-)";
    case pipipi_0 : return str<<"pipipi_0 (gamma* -> pi+ pi- pi0)";
    case V_pi     : return str<<"V_pi (anomalous gamma rho pi)";
    case unknown  : return str<<"unknown";
    }
    return str<<"invalid("<<int(m)<<")";
  }

  inline bool IsChargeFormFactor(const FF_0_PP_mode & m) {
    return (m==pipi_plus || m==KK_plus || m==Kpi_plus ||
	    m==pipi_0    || m==KK_0);
  }

  class FFVMD: public Form_Factor {
  private:
    METOOLS::Summed_Propagator m_props;
    FF_0_PP_mode m_mode;
    int          m_pos;
    void FixMode(const Vertex_Key &key);
    void Construct();
    void ConstructPionFormFactor();
    void ConstructKaonFormFactor();
    void ConstructThreePionFormFactor();
    void ConstructVectorPionFormFactor();
  public:
    FFVMD(const Vertex_Key &key);
    Complex FF();
  };

  // Transition form factor (gamma* gamma* -> P)
  class FFTFF: public Form_Factor {
  private:
    METOOLS::Summed_Propagator m_prop;
    void Construct(const kf_code & ps);
  public:
    FFTFF(const Vertex_Key &key);
    Complex FF();
  };
} // end namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

// ----------------------------------------------------------------------
// FFVMD implementation
// ----------------------------------------------------------------------

FFVMD::FFVMD(const Vertex_Key &key):
  Form_Factor("VMD",key),
  m_props(), m_mode(FF_0_PP_mode::unknown),
  m_pos(-1) {
  set<kf_code> kfs;
  msg_Debugging()<<METHOD<<"("<<key.m_j.size()<<"): "<<key.m_j[0]->Flav()<<", "
                 <<key.m_j[1]->Flav()<<"\n";
  msg_Debugging()<<key.p_mv->in[0]<<", "<<key.p_mv->in[1]<<", "
                 <<key.p_mv->in[2]<<"\n";
  for (size_t i(0); i<key.m_j.size(); ++i) {
    if (!key.m_j[i]->Flav().IsHadron()) m_pos=i;
    else kfs.insert(key.m_j[i]->Flav().Kfcode());
  }
  FixMode(key);
  static set<int> reported;
  const bool report(reported.insert(int(m_mode)).second);
  ostringstream tag;
  tag<<"VMD form factor, mode "<<m_mode;
  if (report) FFTableHead(tag.str());   // from FormFactorParameters
  Construct();
  const Complex f0(m_props(0.));
  if (report) FFTableFoot("F(0)", f0);
  CheckFFNormalisation(tag.str(), f0, IsChargeFormFactor(m_mode));
  msg_Debugging()<<METHOD<<":\n";
}

void FFVMD::FixMode(const Vertex_Key &key) {
  multiset<kf_code> kfs;
  for (size_t i(0); i<key.p_mv->in.size(); ++i)
    if (key.p_mv->in[i].IsHadron()) kfs.insert(key.p_mv->in[i].Kfcode());

  if (kfs == multiset<kf_code>{kf_pi, kf_pi_plus, kf_pi_plus}) {
    m_mode = FF_0_PP_mode::pipipi_0;
    return;
  }
  if (kfs == multiset<kf_code>{kf_rho_770, kf_pi} ||
      kfs == multiset<kf_code>{kf_rho_770_plus, kf_pi_plus}) {
    m_mode = FF_0_PP_mode::V_pi;
    return;
  }
  if (kfs == multiset<kf_code>{kf_pi_plus, kf_pi_plus})
    m_mode = FF_0_PP_mode::pipi_0;
  else if (kfs == multiset<kf_code>{kf_K_plus, kf_K_plus})
    m_mode = FF_0_PP_mode::KK_0;
}

void FFVMD::Construct() {
  switch (int(m_mode)) {
  case int(FF_0_PP_mode::pipi_0):     ConstructPionFormFactor(); break;
  case int(FF_0_PP_mode::pipipi_0):   ConstructThreePionFormFactor(); break;
  case int(FF_0_PP_mode::V_pi):       ConstructVectorPionFormFactor(); break;
  case int(FF_0_PP_mode::KK_0):       ConstructKaonFormFactor(); break;
  case int(FF_0_PP_mode::unknown):
  default:
    msg_Out()<<METHOD<<" yields no form factor.\n";
    break;
  }
}

void FFVMD::ConstructPionFormFactor() {
  msg_Debugging()<<METHOD<<": Gounaris-Sakurai parameters\n";
  // Read all pion parameters from settings
  PionFF::ReadAll();

  // Build the omega/phi mixture (fixed widths)
  METOOLS::Summed_Propagator * rhofac = new METOOLS::Summed_Propagator(false);
  rhofac->Add(new METOOLS::Unity(), Complex(1.,0.));
  rhofac->Add(new METOOLS::WeightedBreitWigner(
                LineShapes->Get(Flavour(kf_omega_782)),
                resonance_type::fixed,
                PionFF::Omega782().Mass(),
                PionFF::Omega782().Width()),
              PionFF::Omega782().Weight());
  rhofac->Add(new METOOLS::WeightedBreitWigner(
                NULL,
                resonance_type::fixed,
                PionFF::Phi1020().Mass(),
                PionFF::Phi1020().Width()),
              PionFF::Phi1020().Weight());

  // rho(770) with omega/phi interference
  METOOLS::Multiplied_Propagator * rho = new METOOLS::Multiplied_Propagator();
  rho->Add(new METOOLS::GounarisSakurai(
             LineShapes->Get(Flavour(kf_rho_770)),
             resonance_type::GS,
             PionFF::Rho().Mass(),
             PionFF::Rho().Width()));
  rho->Add(rhofac);
  m_props.Add(rho, PionFF::Rho().Weight());

  // rho(1450), rho(1700), rho(2150) – GS
  m_props.Add(new METOOLS::GounarisSakurai(
                LineShapes->Get(Flavour(kf_rho_1450)),
                resonance_type::GS,
                PionFF::Rho1450().Mass(),
                PionFF::Rho1450().Width()),
              PionFF::Rho1450().Weight());
  m_props.Add(new METOOLS::GounarisSakurai(
                LineShapes->Get(Flavour(kf_rho_1700)),
                resonance_type::GS,
                PionFF::Rho1700().Mass(),
                PionFF::Rho1700().Width()),
              PionFF::Rho1700().Weight());
  m_props.Add(new METOOLS::GounarisSakurai(
                NULL,
                resonance_type::GS,
                PionFF::Rho2150().Mass(),
                PionFF::Rho2150().Width()),
              PionFF::Rho2150().Weight());
}

void FFVMD::ConstructKaonFormFactor() {
  msg_Debugging()<<METHOD<<": Gounaris-Sakurai parameters\n";
  KaonFF::ReadAll();   // read all kaon parameters

  // Three independent SU(3) channels – summed with correct weights
  METOOLS::Summed_Propagator * kaon = new METOOLS::Summed_Propagator(false);

  // rho family (isovector, weight 1/2)
  kaon->Add(new METOOLS::GounarisSakurai(
              LineShapes->Get(Flavour(kf_rho_770)),
              resonance_type::GS,
              KaonFF::Rho().Mass(),
              KaonFF::Rho().Width()),
            0.5 * KaonFF::Rho().Weight());
  kaon->Add(new METOOLS::GounarisSakurai(
              LineShapes->Get(Flavour(kf_rho_1450)),
              resonance_type::GS,
              KaonFF::Rho1450().Mass(),
              KaonFF::Rho1450().Width()),
            0.5 * KaonFF::Rho1450().Weight());
  kaon->Add(new METOOLS::GounarisSakurai(
              LineShapes->Get(Flavour(kf_rho_1700)),
              resonance_type::GS,
              KaonFF::Rho1700().Mass(),
              KaonFF::Rho1700().Width()),
            0.5 * KaonFF::Rho1700().Weight());
  // rho(2150) commented out in original – keep as is

  // omega family (isoscalar, weight 1/6) – fixed width
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,
              KaonFF::Omega().Mass(),
              KaonFF::Omega().Width()),
            KaonFF::Omega().Weight() / 6.);
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,
              KaonFF::Omega1420().Mass(),
              KaonFF::Omega1420().Width()),
            KaonFF::Omega1420().Weight() / 6.);
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,
              KaonFF::Omega1650().Mass(),
              KaonFF::Omega1650().Width()),
            KaonFF::Omega1650().Weight() / 6.);

  // phi family (isoscalar, weight 1/3) – fixed width
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,
              KaonFF::Phi().Mass(),
              KaonFF::Phi().Width()),
            KaonFF::Phi().Weight() / 3.);
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,
              KaonFF::Phi1680().Mass(),
              KaonFF::Phi1680().Width()),
            KaonFF::Phi1680().Weight() / 3.);

  m_props.Add(kaon, Complex(1.,0.));
}

void FFVMD::ConstructThreePionFormFactor() {
  msg_Debugging()<<METHOD<<": Gounaris-Sakurai parameters\n";
  PionFF::ReadAll();   // we reuse the pion rho parameters

  // rho family for the two-pion sub-channel
  m_props.Add(new METOOLS::GounarisSakurai(
                LineShapes->Get(Flavour(kf_rho_770)),
                resonance_type::GS,
                PionFF::Rho().Mass(),
                PionFF::Rho().Width()),
              PionFF::Rho().Weight());
  m_props.Add(new METOOLS::GounarisSakurai(
                LineShapes->Get(Flavour(kf_rho_1450)),
                resonance_type::GS,
                PionFF::Rho1450().Mass(),
                PionFF::Rho1450().Width()),
              PionFF::Rho1450().Weight());
  m_props.Add(new METOOLS::GounarisSakurai(
                LineShapes->Get(Flavour(kf_rho_1700)),
                resonance_type::GS,
                PionFF::Rho1700().Mass(),
                PionFF::Rho1700().Width()),
              PionFF::Rho1700().Weight());
}

void FFVMD::ConstructVectorPionFormFactor() {
  msg_Debugging()<<METHOD<<": three-pion isoscalar parameters\n";
  ThreePionFF::ReadAll();   // omega(782), phi(1020), J/psi

  // gamma -> omega/phi/Jpsi -> rho pi – fixed width
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,
              ThreePionFF::Omega782().Mass(),
              ThreePionFF::Omega782().Width()),
            ThreePionFF::Omega782().Weight());
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,
              ThreePionFF::Phi1020().Mass(),
              ThreePionFF::Phi1020().Width()),
            ThreePionFF::Phi1020().Weight());
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,
              ThreePionFF::Jpsi().Mass(),
              ThreePionFF::Jpsi().Width()),
            ThreePionFF::Jpsi().Weight());
}

Complex FFVMD::FF() {
  double Q2(0.);
  bool found(false);
  for (size_t i=0; i<p_v->J().size(); ++i) {
    if (p_v->J(i)->Flav().IsHadron()) continue;
    const double q2 = p_v->J(i)->P().Abs2();
    if (!found || std::abs(q2)>std::abs(Q2)) { Q2=q2; found=true; }
  }
  if (!p_v->JC()->Flav().IsHadron()) {
    const double q2 = p_v->JC()->P().Abs2();
    if (!found || std::abs(q2)>std::abs(Q2)) { Q2=q2; found=true; }
  }
  if (!found) Q2 = p_v->JC()->P().Abs2();
  return m_props(Q2);
}

// ----------------------------------------------------------------------
// FFTFF implementation (transition form factor)
// ----------------------------------------------------------------------

FFTFF::FFTFF(const Vertex_Key &key):
  Form_Factor("TFF",key), m_prop() {
  kf_code ps = kf_none;
  for (size_t i=0; i<key.p_mv->in.size(); ++i)
    if (key.p_mv->in[i].IsHadron()) ps = key.p_mv->in[i].Kfcode();

  static set<kf_code> reported;
  const bool report = reported.insert(ps).second;
  ostringstream title;
  title << "Transition form factor, " << Flavour(ps);
  if (report) FFTableHead(title.str());

  Construct(ps);
  const Complex f0 = m_prop(0.) * m_prop(0.);
  if (report) FFTableFoot("F(0,0)", f0);
  CheckFFNormalisation("FFTFF(kf = " + std::to_string(ps) + ")", f0);
}

void FFTFF::Construct(const kf_code & ps) {
  msg_Debugging()<<METHOD<<": VMD poles for kf = "<<ps<<"\n";
  // Read all two‑photon parameters (they are all in the same sector)
  TwoPhotonFF::ReadAll();

  // Select the appropriate parameter set for this pseudoscalar
  ResonanceParameters *rho(nullptr), *omega(nullptr), *phi(nullptr);
  switch (ps) {
  case kf_pi:
    rho   = &TwoPhotonFF::Pi0Rho();
    omega = &TwoPhotonFF::Pi0Omega();
    phi   = &TwoPhotonFF::Pi0Phi();
    break;
  case kf_eta:
    rho   = &TwoPhotonFF::EtaRho();
    omega = &TwoPhotonFF::EtaOmega();
    phi   = &TwoPhotonFF::EtaPhi();
    break;
  case kf_eta_prime_958:
    rho   = &TwoPhotonFF::EtapRho();
    omega = &TwoPhotonFF::EtapOmega();
    phi   = &TwoPhotonFF::EtapPhi();
    break;
  default:
    msg_Out()<<METHOD<<": no transition form factor for kf = "<<ps
             <<", falling back to F = 1.\n";
    return;
  }

  // All fixed‑width propagators – normalisation enforced by Summed_Propagator
  m_prop.Add(new METOOLS::FixedBreitWigner(NULL, rho->Mass(), rho->Width()),
             rho->Weight());
  m_prop.Add(new METOOLS::FixedBreitWigner(NULL, omega->Mass(), omega->Width()),
             omega->Weight());
  m_prop.Add(new METOOLS::FixedBreitWigner(NULL, phi->Mass(), phi->Width()),
             phi->Weight());
}

Complex FFTFF::FF() {
  double q2[2];
  size_t n=0;
  for (size_t i=0; i<p_v->J().size() && n<2; ++i) {
    Current * j = p_v->J(i);
    if (!j->Flav().IsHadron()) q2[n++] = j->P().Abs2();
  }
  if (n<2) q2[n++] = p_v->JC()->P().Abs2();
  if (n<2) return Complex(1.,0.);
  return m_prop(q2[0]) * m_prop(q2[1]);
}

// ----------------------------------------------------------------------
// Getter declarations
// ----------------------------------------------------------------------

DECLARE_GETTER(FFTFF,"TFF",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFTFF>::
operator()(const Vertex_Key &args) const { return new FFTFF(args); }
void ATOOLS::Getter<Form_Factor,Vertex_Key,FFTFF>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"gamma* gamma* -> P transition form factor"; }

DECLARE_GETTER(FFVMD,"VMD",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFVMD>::
operator()(const Vertex_Key &args) const {
  msg_Debugging()<<METHOD<<"(size = "<<args.m_j.size()<<")\n";
  return new FFVMD(args);
}
void ATOOLS::Getter<Form_Factor,Vertex_Key,FFVMD>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"VMD"; }