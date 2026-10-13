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
    // hep-ph/0512180 Eq. (10) terms E and F: gamma -> rho(1450)/rho(1700) pi,
    // each carrying ONLY its partner isoscalar pole (phi(1020) for E,
    // omega(1631) for F) rather than the full isoscalar sum V_pi uses - see
    // FFVMD::ConstructVectorPionFormFactorRhoE/F and the production vertices
    // in LowEnergy_Model::InitThreePionVertices.
    V_pi_rho1450    = 302,
    V_pi_rho1700    = 303,
    // gamma* -> pi+ pi- pi0 as one contact current (VPPP_LC.C), carrying the
    // isoscalar tower in q^2 times the rho tower in the two-pion sub-masses.
    ee3pi           = 401,
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
    case V_pi_rho1450: return str<<"V_pi_rho1450 (hep-ph/0512180 Eq.10 term E)";
    case V_pi_rho1700: return str<<"V_pi_rho1700 (hep-ph/0512180 Eq.10 term F)";
    case ee3pi    : return str<<"ee3pi (gamma* -> pi+ pi- pi0 contact current)";
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
    // V_pi_rho1450/V_pi_rho1700 (hep-ph/0512180 Eq.10 terms E/F) use this
    // instead of m_props: their HepPh0512180_BreitWigner is NOT normalised
    // to 1 at s=0 and must be evaluated raw, which Summed_Propagator's
    // divide-by-weight-sum default cannot do without corrupting the value -
    // see the class comment in Propagator.H. Owned; NULL for every other
    // mode.
    METOOLS::Propagator_Base * p_raw_prop;
    Complex      m_raw_weight;
    FF_0_PP_mode m_mode;
    int          m_pos;
    void FixMode(const Vertex_Key &key);
    void Construct();
    void ConstructPionFormFactor();
    void ConstructKaonFormFactor();
    void ConstructThreePionFormFactor();
    void ConstructVectorPionFormFactor();
    void ConstructVectorPionFormFactorSubtracted();
    void ConstructEE3PiFormFactor();
    // ee3pi only: the isoscalar a(q^2) and one rho propagator per two-pion
    // sub-mass.  Owned; NULL in every other mode.
    METOOLS::Propagator_Base * p_iso3pi;
    METOOLS::Propagator_Base * p_rho3pi;
    // hep-ph/0512180 eq. (10) terms E and F.  Each pairs ONE isoscalar pole
    // with ONE excited-rho sub-current, so they cannot be folded into the
    // A-D isoscalar or into the shared rho tower - eq. (10) is a sum of
    // products, not a product of sums.  NULL unless the eq-10 isoscalar is
    // active.
    METOOLS::Propagator_Base * p_isoE, * p_rhoE;   // phi(1020) x rho(1450)
    METOOLS::Propagator_Base * p_isoF, * p_rhoF;   // omega(1631) x rho(1700)
    Complex m_cE, m_cF;
    double m_alphaA;
    Complex EE3PiFF();
    METOOLS::Propagator_Base * BuildEE3PiRhoTower();
    void ConstructVectorPionFormFactorRhoE();
    void ConstructVectorPionFormFactorRhoF();
  public:
    FFVMD(const Vertex_Key &key);
    ~FFVMD();
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
  m_props(), p_raw_prop(NULL), m_raw_weight(1.,0.),
  p_iso3pi(NULL), p_rho3pi(NULL),
  p_isoE(NULL), p_rhoE(NULL), p_isoF(NULL), p_rhoF(NULL),
  m_cE(0.,0.), m_cF(0.,0.), m_alphaA(1.),
  m_mode(FF_0_PP_mode::unknown),
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
  const Complex f0(m_mode==FF_0_PP_mode::ee3pi ? Complex(m_alphaA,0.) :
                   p_raw_prop!=NULL ? m_raw_weight*(*p_raw_prop)(0.) : m_props(0.));
  if (report) FFTableFoot("F(0)", f0);
  // p_raw_prop modes (Eq.10 E/F) are not normalised to 1 at s=0 by design -
  // see the HepPh0512180_BreitWigner class comment in Propagator.H - so the
  // usual F(0)=1 check does not apply to them.
  // ee3pi's F(0) is the anomaly alpha_A, a physical GeV^-3 scale, not 1.
  if (p_raw_prop==NULL && m_mode!=FF_0_PP_mode::ee3pi)
    CheckFFNormalisation(tag.str(), f0, IsChargeFormFactor(m_mode));
  msg_Debugging()<<METHOD<<":\n";
}

FFVMD::~FFVMD() {
  delete p_raw_prop; delete p_iso3pi; delete p_rho3pi;
  delete p_isoE; delete p_rhoE; delete p_isoF; delete p_rhoF;
}

void FFVMD::FixMode(const Vertex_Key &key) {
  multiset<kf_code> kfs;
  for (size_t i(0); i<key.p_mv->in.size(); ++i)
    if (key.p_mv->in[i].IsHadron()) kfs.insert(key.p_mv->in[i].Kfcode());

  if (kfs == multiset<kf_code>{kf_pi, kf_pi_plus, kf_pi_plus}) {
    // Both the old three-leg contact vertex and the new four-leg
    // gamma-pi-pi-pi current present the same hadron multiset, so the leg
    // count is what separates them: 4 legs means the photon is part of the
    // vertex and this is the full ee -> 3pi current.
    m_mode = (key.p_mv->in.size()==4) ? FF_0_PP_mode::ee3pi
                                      : FF_0_PP_mode::pipipi_0;
    return;
  }
  if (kfs == multiset<kf_code>{kf_rho_770, kf_pi} ||
      kfs == multiset<kf_code>{kf_rho_770_plus, kf_pi_plus}) {
    m_mode = FF_0_PP_mode::V_pi;
    return;
  }
  if (kfs == multiset<kf_code>{kf_rho_1450, kf_pi} ||
      kfs == multiset<kf_code>{kf_rho_1450_plus, kf_pi_plus}) {
    m_mode = FF_0_PP_mode::V_pi_rho1450;
    return;
  }
  if (kfs == multiset<kf_code>{kf_rho_1700, kf_pi} ||
      kfs == multiset<kf_code>{kf_rho_1700_plus, kf_pi_plus}) {
    m_mode = FF_0_PP_mode::V_pi_rho1700;
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
  case int(FF_0_PP_mode::ee3pi):      ConstructEE3PiFormFactor(); break;
  case int(FF_0_PP_mode::V_pi_rho1450): ConstructVectorPionFormFactorRhoE(); break;
  case int(FF_0_PP_mode::V_pi_rho1700): ConstructVectorPionFormFactorRhoF(); break;
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

// The once-subtracted isoscalar of arXiv:1907.01556 eq. (2.14),
//
//     a(q^2) = alpha_A [ 1 + Sum_V d_V ( BW_V(q^2) - BW_V(0) ) ]
//
// ported from tau-decays-in-sherpa (METOOLS-hadron-currents),
// FF_0_Isoscalar3Pi::Construct().  OPT-IN, default off: the summed form below
// is what the existing THREE_PION_FORM_FACTOR amplitudes were tuned against,
// so switching the functional form silently would invalidate that tune.
// Enable with THREE_PION_FORM_FACTOR: {Subtracted: true}.
//
// Why this is not the same as the default path with different weights:
//   Summed_Propagator gives Sum w_V BW_V / Sum w_V, which is 1 at q^2 = 0 and
//   falls as 1/q^2;
//   this gives alpha_A at q^2 = 0 and tends to the CONSTANT
//   alpha_A (1 - Sum_V d_V BW_V(0)), which is what holds the form factor up
//   between resonances.
// The 1/q^2 asymptote of the default path is the likely origin of the
// long-standing high-mass tail problem the omega(1375)/omega(1631) amplitudes
// were added to patch - see the note by them below.
//
// TWO DEVIATIONS from the reference, both deliberate and both meaning the
// d_V below may need refitting rather than being trusted as-is:
//  1. The reference takes all four isoscalars as RUNNING-width registry line
//     shapes.  Our Omega_1420_Lineshape/Omega_1600_Lineshape exist but have
//     their entire constructor bodies commented out (Omega_Decays.C:72,92),
//     so they would register zero channels, give a zero running width and
//     leave an undamped pole.  Fixed-width Breit-Wigners with the masses and
//     widths from our own parameter table are used for all four instead.
//  2. alpha_A is NOT applied here.  It is an overall constant, and FF() is a
//     shape that the AVVP vertex coupling g_gammarhopi multiplies; putting a
//     ~11.4 GeV^-3 constant in the form factor while the coupling still
//     carries its own tuned normalisation would double count.  The subtracted
//     SHAPE is what this provides; the scale stays in the coupling.
void FFVMD::ConstructVectorPionFormFactorSubtracted() {
  msg_Debugging()<<METHOD<<": subtracted three-pion isoscalar\n";
  ThreePionFF::ReadAll();

  // c0 = 1: the shape is normalised so a(0) = 1, with the physical anomaly
  // scale left in g_gammarhopi (see 2. above).
  METOOLS::Subtracted_Propagator * iso =
    new METOOLS::Subtracted_Propagator(Complex(1.,0.));
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Omega782().Mass(),
             ThreePionFF::Omega782().Width()),
           ThreePionFF::Omega782().Weight());
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Phi1020().Mass(),
             ThreePionFF::Phi1020().Width()),
           ThreePionFF::Phi1020().Weight());
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Omega1375().Mass(),
             ThreePionFF::Omega1375().Width()),
           ThreePionFF::Omega1375().Weight());
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Omega1631().Mass(),
             ThreePionFF::Omega1631().Width()),
           ThreePionFF::Omega1631().Weight());
  // Evaluated directly rather than through m_props: m_props normalises by the
  // weight sum, which is exactly what this form exists to avoid.
  p_raw_prop = iso;
  m_raw_weight = Complex(1.,0.);
}

// gamma* -> pi+ pi- pi0 as one current, following arXiv:1907.01556 as
// implemented in tau-decays-in-sherpa's FF_0_Isoscalar3Pi:
//
//     F(q^2; s_+-, s_+0, s_-0) = a(q^2) * [ BW(s_+-) + BW(s_+0) + BW(s_-0) ]
//     a(q^2) = alpha_A [ 1 + Sum_V d_V ( BW_V(q^2) - BW_V(0) ) ]
//
// The rho sum is ALL-PLUS and that is forced, not chosen: an I = 0 three-pion
// state is totally antisymmetric in isospin and the momentum epsilon is
// totally antisymmetric, so |A|^2 must be totally symmetric under permuting
// the pions.  Writing the three terms here, rather than letting three
// separate Comix diagrams supply them, is the whole point of the contact
// construction - in the factorised gamma->rho pi chain the relative weights
// come out (+2,-1,+1) instead, and have to be corrected by hand.
//
// alpha_A is the Wess-Zumino-Witten anomaly and is PREDICTED, not fitted.  It
// is applied here rather than in the vertex coupling, because in this
// construction FF() is the entire dynamical factor and the vertex carries
// only the epsilon structure.
void FFVMD::ConstructEE3PiFormFactor() {
  msg_Debugging()<<METHOD<<": ee -> 3pi contact current\n";
  ThreePionFF::ReadAll();
  PionFF::ReadAll();

  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"] };
  const double fpi(s["F_pi"].SetDefault(0.1307).Get<double>()/sqrt(2.));
  const double F3pi(1./(4.*sqr(M_PI)*pow(fpi,3.)));
  m_alphaA = F3pi/3.*1.066;

  // a(q^2)/alpha_A, subtracted so that it is 1 at the origin and tends to a
  // constant rather than falling as 1/q^2.
  //
  // The couplings d_V are RELATIVE TO THE ANOMALY, not to omega(782): the
  // q^2 -> 0 value is fixed at 1 by the low-energy theorem, so d_V measures
  // how much resonance strength sits on top of it.  They are therefore NOT
  // the same quantity as the omega-relative THREE_PION_FORM_FACTOR
  // amplitudes the non-contact path uses, and the two must not be mixed.
  //
  // Defaults are the reference implementation's own fit to BaBar 2021
  // (arXiv:2110.00520) over 0.62-1.80 GeV - tau-decays-in-sherpa,
  // FF_0_Isoscalar3Pi::FixParameters.  Carried over because they belong to
  // exactly this functional form; they are a starting point for a refit
  // here, not a guarantee, since that fit used running-width registry line
  // shapes where these are fixed-width (see below).
  //
  // TWO DISTINCT SETS, not one - do not mix them.  The reference's "fitted"
  // mode uses the complex set below (amplitude and phase). Its "dispersive"
  // mode, the one the sum rule belongs to, uses a SEPARATE real-only set
  // instead - RealCoupling(), not Coupling(), in FF_0_Isoscalar3Pi.C - and
  // that is deliberate, not an oversight: a nonzero phase puts
  // Im a(infinity) != 0, which the sum rule's own derivation requires to
  // vanish (eq. 2.14's dispersion integral is otherwise logarithmically
  // divergent). Confusing the two - which the first version of this function
  // did, feeding the complex fitted set into the sum-rule machinery - breaks
  // the sum rule's premise before it even solves anything, and was the
  // likely cause of a 3-7x inflation measured over 1.3-3.5 GeV against
  // BELLE_2024_I2775022.
  const bool sumrule(s["SumRule"].SetDefault(false).Get<bool>());
  Complex w[4];
  if (sumrule) {
    const double rv[4] = { 0.409504, -0.032589, -0.003088, -0.047799 };
    const char * rvn[4] = { "d_omega782","d_phi1020","d_omega1420","d_omega1650" };
    for (size_t i(0);i<4;++i)
      w[i] = Complex(s[rvn[i]].SetDefault(rv[i]).Get<double>(),0.);
  }
  else {
    const double dv[4][2] = { {1.0160,5.7518}, {0.0628,2.0436},
                              {0.0203,2.3891}, {0.0661,1.3573} };
    const char * dvn[4] = { "d_omega782","d_phi1020","d_omega1420","d_omega1650" };
    for (size_t i(0);i<4;++i) {
      Scoped_Settings ds{ s[dvn[i]] };
      const double a(ds["Amplitude"].SetDefault(dv[i][0]).Get<double>());
      const double p(ds["Phase"].SetDefault(dv[i][1]).Get<double>());
      w[i] = a*Complex(cos(p),sin(p));
    }
  }

  // Summed (default) vs subtracted isoscalar.  They agree wherever the
  // resonances dominate and differ only in the asymptotics:
  //   summed      a = Sum_V w_V BW_V(q^2) / Sum_V w_V, 1 at q^2 = 0 and
  //               falling as 1/q^2;
  //   subtracted  a = 1 + Sum_V d_V ( BW_V(q^2) - BW_V(0) ), also 1 at the
  //               origin but tending to the CONSTANT 1 - Sum_V d_V BW_V(0).
  //
  // SUMMED IS THE DEFAULT, and that is a measured choice.  Compared point by
  // point against hep-ph/0512180 eq. (10) - the variant demonstrated to track
  // BELLE_2024_I2775022 over the full 1.05-3.5 GeV range - using
  // doc/examples/Soft_QCD/LowEnergy/ff3pi_check.C:
  //
  //     q (GeV)      1.32   1.62   1.92   2.22   2.52   2.82   3.12   3.42
  //     summed/ref   0.82   0.90   1.30   1.33   1.39   1.31   1.40   1.99
  //     subtr/ref    0.96   0.94   1.57   3.19   5.01   6.37   8.68  15.34
  //
  // Below ~1.7 GeV either will do.  Above it the subtracted form's constant
  // is the entire disagreement: eq. (10) has NO subtraction constant, so it
  // dies as 1/q^2, while ours plateaus.  The subtracted construction is not
  // wrong physics - it is arXiv:1907.01556 eq. (2.14), where the anomaly sum
  // rule is what cancels that constant - but our sum rule's c_i are not
  // fitted for this model, so the constant currently survives.  Until they
  // are, the summed form is the one that reproduces the reference.
  //
  // NOTE the weights mean different things in the two cases - d_V are
  // relative to the anomaly, the summed weights are relative to each other -
  // so a fit done for one does not carry over to the other.
  // ---- Eq. (10): the reference's own isoscalar, DEFAULT ----
  //
  // hep-ph/0512180 eq. (10)'s A-D parenthesis, on UNNORMALISED Breit-Wigners
  // 1/(s/m^2-1+i Gamma/m), with that paper's own real coefficients.  This is
  // the variant demonstrated to track BELLE_2024_I2775022 across the whole
  // 1.05-3.5 GeV range, and it is the default because measuring it against
  // the alternatives left no contest.
  //
  // What it does that neither of ours did, measured against our summed form:
  //     q (GeV)     1.1   1.2   1.3   1.4   1.5   1.6   1.7   1.8   2.0
  //     ref/ours   1.21  1.23  1.26  1.20  1.01  0.84  0.64  0.44  0.22
  // i.e. MORE strength through 1.1-1.4 and a far steeper falloff above 1.6 -
  // exactly the two-sided correction our scans asked for, and exactly what no
  // single normalisation could supply.
  //
  // The mechanism is the RELATIVE SIGNS.  All of B, C, D are negative against
  // A = +18.20, so the omega(1375)/omega(1631) terms interfere destructively
  // with the omega(782) tail: that flattens the region between them into a
  // plateau and then cancels hard at large q^2.  Our fitted-mode d_V are
  // complex with phases fitted for the SUBTRACTED form against BaBar, a
  // different interference pattern entirely, which is why reusing them here
  // reproduced neither feature.
  //
  // Absolute scale: eq. (10) carries its own (GeV^-3), so this REPLACES
  // alpha_A * a(q^2) rather than multiplying it - see EE3PiFF().
  // DEFAULT OFF.  Measured against BELLE_2024_I2775022 the summed isoscalar
  // below is the better of the two overall: its tail agrees (median
  // model/data 0.88 above 1.9 GeV) and its only defect is a smooth
  // undershoot over 1.06-1.5 GeV.  Eq. (10) fixes that mid-range but its
  // A-D combination has a destructive-interference trough near 2.1 GeV
  // which shows up as a DIP AT 2 GeV; terms E and F fill most of it but not
  // all, and eq. (10) also pairs A-D with H_rho(770) alone whereas our A-D
  // multiplies the full GS tower, so some excited-rho strength is counted
  // twice.  Enable with THREE_PION_FORM_FACTOR: {EqTen: true} once that
  // pairing is sorted out.
  const bool eqten(s["EqTen"].SetDefault(false).Get<bool>());
  if (eqten) {
    const double cA(s["C_A"].SetDefault( 18.20).Get<double>());
    const double cB(s["C_B"].SetDefault( -0.87).Get<double>());
    const double cC(s["C_C"].SetDefault( -0.77).Get<double>());
    const double cD(s["C_D"].SetDefault( -1.12).Get<double>());
    // Table-1 masses and widths of that paper, NOT our particle table: the
    // coefficients only mean anything with the poles they were fitted to.
    METOOLS::Summed_Propagator * iso = new METOOLS::Summed_Propagator(false);
    iso->Add(new METOOLS::HepPh0512180_BreitWigner(0.7824 ,0.00869),cA);
    iso->Add(new METOOLS::HepPh0512180_BreitWigner(1.01924,0.00414),cB);
    iso->Add(new METOOLS::HepPh0512180_BreitWigner(1.3750 ,0.2500 ),cC);
    iso->Add(new METOOLS::HepPh0512180_BreitWigner(1.6310 ,0.2450 ),cD);
    p_iso3pi = iso;
    p_rho3pi = BuildEE3PiRhoTower();
    // alpha_A is NOT applied on top of this one.
    m_alphaA = 1.;

    // ---- terms E and F ----
    // Without them the A-D combination has a destructive-interference trough:
    // |A-D| falls to a minimum of ~0.65 around 2.10-2.15 GeV before
    // flattening, which shows up as a DIP AT 2 GeV in the cross section.
    // E and F attach phi(1020) to the rho(1450) sub-current and omega(1631)
    // to the rho(1700) one; both excited rhos peak exactly where A-D is
    // collapsing, which is what fills the trough.  That pairing is why
    // eq. (10) is a sum of products - folding these into the shared rho
    // tower would spread each across every isoscalar instead.
    const double cE(s["C_E"].SetDefault(-0.72).Get<double>());
    const double cF(s["C_F"].SetDefault(-0.59).Get<double>());
    m_cE = Complex(cE,0.); m_cF = Complex(cF,0.);
    p_isoE = new METOOLS::HepPh0512180_BreitWigner(1.01924,0.00414);
    p_isoF = new METOOLS::HepPh0512180_BreitWigner(1.6310 ,0.2450 );
    // eq. (10)'s own excited-rho poles, Table 1 - not our particle table,
    // since the coefficients only mean anything with the poles they were
    // fitted against.
    p_rhoE = new METOOLS::HepPh0512180_BreitWigner(1.465,0.310);
    p_rhoF = new METOOLS::HepPh0512180_BreitWigner(1.700,0.235);

    msg_Out()<<METHOD<<": isoscalar = hep-ph/0512180 eq.(10), "
	     <<"A="<<cA<<" B="<<cB<<" C="<<cC<<" D="<<cD
	     <<" E="<<cE<<" F="<<cF<<"\n";
    return;
  }

  const bool subtracted(s["Subtracted"].SetDefault(false).Get<bool>());
  // The sum rule only means anything for the subtracted form - it exists to
  // cancel that form's constant asymptote.  Asking for it with the summed
  // isoscalar would silently do nothing except swap in the real coupling
  // set (selected above), which is a quiet way to get wrong numbers, so
  // refuse the combination instead.
  if (sumrule && !subtracted)
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: SumRule needs Subtracted: true "
	  "- the sum rule cancels the subtracted form's constant asymptote and "
	  "has no meaning for the summed form.");
  METOOLS::Propagator_Base * isoprop(NULL);
  if (!subtracted) {
    METOOLS::Summed_Propagator * sum = new METOOLS::Summed_Propagator(true);
    sum->Add(new METOOLS::FixedBreitWigner(NULL,
               ThreePionFF::Omega782().Mass(),
               ThreePionFF::Omega782().Width()), w[0]);
    sum->Add(new METOOLS::FixedBreitWigner(NULL,
               ThreePionFF::Phi1020().Mass(),
               ThreePionFF::Phi1020().Width()), w[1]);
    sum->Add(new METOOLS::FixedBreitWigner(NULL,
               ThreePionFF::Omega1375().Mass(),
               ThreePionFF::Omega1375().Width()), w[2]);
    sum->Add(new METOOLS::FixedBreitWigner(NULL,
               ThreePionFF::Omega1631().Mass(),
               ThreePionFF::Omega1631().Width()), w[3]);
    p_iso3pi = sum;
    p_rho3pi = BuildEE3PiRhoTower();
    msg_Out()<<METHOD<<": isoscalar = SUMMED (falls as 1/q^2), alpha_A = "
	     <<m_alphaA<<" GeV^-3\n";
    return;
  }
  msg_Out()<<METHOD<<": isoscalar = SUBTRACTED (tends to a constant)\n";

  METOOLS::Subtracted_Propagator * iso =
    new METOOLS::Subtracted_Propagator(Complex(1.,0.));
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Omega782().Mass(),
             ThreePionFF::Omega782().Width()), w[0]);
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Phi1020().Mass(),
             ThreePionFF::Phi1020().Width()), w[1]);
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Omega1375().Mass(),
             ThreePionFF::Omega1375().Width()), w[2]);
  iso->Add(new METOOLS::FixedBreitWigner(NULL,
             ThreePionFF::Omega1631().Mass(),
             ThreePionFF::Omega1631().Width()), w[3]);
  // ---- the anomaly sum rule, arXiv:1907.01556 eq. (2.16) ----
  //
  // Without this a(q^2) tends to the CONSTANT 1 - Sum_V d_V P_V(0) and the
  // three-pion cross section flattens out above ~2 GeV, where the data falls
  // steeply - measured against BELLE_2024_I2775022 d02, the plateau sits at
  // ~1.7-2.0 while the data drops to 0.02-0.08, i.e. wrong by up to ~90x at
  // 3.3 GeV.  The sum rule drives a(q^2) -> 0 instead, as an unsubtracted
  // dispersion relation demands.
  //
  // Two constraints, both linear in the coefficients, solved for the last
  // two of them (so at least three are needed):
  //   (I)  a(infinity) = 0, i.e. C_p(infinity) = Sum_V d_V P_V(0) - 1;
  //   (II) Sum_i i c_i = 0, removing the square-root cusp at s_inel.
  //
  // NOTE the units differ from the reference: its Subtracted_Propagator
  // carries c0 = alpha_A, so its target is (Sum_V alpha_A d_V P_V(0) -
  // alpha_A) in GeV^-3.  Ours carries c0 = 1 with alpha_A applied in
  // EE3PiFF(), so the same condition divided through by alpha_A - i.e. the
  // target below is dimensionless.  Getting this wrong silently rescales the
  // constraint.
  //
  // OFF by default: the reference flags this construction as needing a refit
  // before it is believed, and its c_i were fitted with running-width line
  // shapes where these are fixed-width.  THREE_PION_FORM_FACTOR: {SumRule: true}.
  if (sumrule) {
    std::vector<double> cp;
    cp.push_back(s["cpoly_1"].SetDefault(-0.12).Get<double>());
    cp.push_back(s["cpoly_2"].SetDefault(-1.36).Get<double>());
    cp.push_back(s["cpoly_3"].SetDefault(-0.47).Get<double>());
    cp.push_back(0.);   // determined
    cp.push_back(0.);   // determined
    METOOLS::Conformal_Polynomial * cpoly =
      new METOOLS::Conformal_Polynomial(cp,1.,-1.);
    const size_t p(cp.size());
    const double z0(cpoly->Z0());

    // P_V(0) is EVALUATED, not assumed to be 1: true for a Breit-Wigner,
    // false for Gounaris-Sakurai, and Subtracted_Propagator measures it the
    // same way, so the two agree by construction.
    const double mm[4] = { ThreePionFF::Omega782().Mass(),
                           ThreePionFF::Phi1020().Mass(),
                           ThreePionFF::Omega1375().Mass(),
                           ThreePionFF::Omega1631().Mass() };
    const double gg[4] = { ThreePionFF::Omega782().Width(),
                           ThreePionFF::Phi1020().Width(),
                           ThreePionFF::Omega1375().Width(),
                           ThreePionFF::Omega1631().Width() };
    double sub(0.);
    for (size_t i(0);i<4;++i) {
      METOOLS::FixedBreitWigner bw(NULL,mm[i],gg[i]);
      sub += (w[i]*bw(0.)).real();
    }
    const double want(sub-1.);            // required C_p(infinity)

    const size_t i1(p-1), i2(p);          // the two 1-based powers solved
    double r1(0.), r2(0.);
    for (size_t k(0);k+2<p;++k) {
      const double i = double(k+1);
      r1 += i*cp[k];
      r2 += cp[k]*(((k+1)%2) ? -1.-pow(z0,i) : 1.-pow(z0,i));
    }
    const double f1 = ((i1%2) ? -1.-pow(z0,double(i1)) : 1.-pow(z0,double(i1)));
    const double f2 = ((i2%2) ? -1.-pow(z0,double(i2)) : 1.-pow(z0,double(i2)));
    const double det = double(i1)*f2-double(i2)*f1;
    if (dabs(det)<1.e-12)
      THROW(fatal_error,"ee3pi: conformal-polynomial constraints degenerate.");
    const double b1(-r1), b2(want-r2);
    const double c1(( b1*f2-double(i2)*b2)/det);
    const double c2((double(i1)*b2-b1*f1)/det);
    cpoly->SetCoefficient(p-2,c1);
    cpoly->SetCoefficient(p-1,c2);
    iso->Add(cpoly,Complex(1.,0.));
    msg_Out()<<METHOD<<": sum rule wants C_p(inf) = "<<want
	     <<"; solved c_"<<i1<<" = "<<c1<<", c_"<<i2<<" = "<<c2
	     <<"; C_p(inf) = "<<cpoly->AtInfinity()<<"\n";
  }

  p_iso3pi = iso;

  // The rho tower in the two-pion sub-mass - see BuildEE3PiRhoTower().
  p_rho3pi = BuildEE3PiRhoTower();

  msg_Out()<<METHOD<<": alpha_A = "<<m_alphaA<<" GeV^-3 (WZW anomaly, F_pi = "
	   <<fpi<<")\n";
}

// F = alpha_A * (a(q^2)/alpha_A) * Sum_ij BW_rho(s_ij), evaluated from the
// three pion momenta.  Unlike every other mode here, this needs the two-pion
// sub-masses and not just the photon virtuality.
// The SAME rho tower as ConstructPionFormFactor's rho(770)/rho(1450)/
// rho(1700) mixture (weight 1 for rho(770), PionFF::Rho1450()/Rho1700() for
// the excited states) - this is what the reference's LineShapes->GetChannel
// (hadronic_channel::pipi_Pwave) actually returns (tau-decays-in-sherpa,
// Line_Shapes.C: a Summed_Propagator over exactly those three, with exactly
// PION_FORM_FACTOR's own amplitudes), so this is not a new fit, it is
// reusing the one already in this model.  A bare GS(770), which is what
// EE3PiFF used before this, has no rho(1450)/rho(1700) contribution at all
// and would undershoot in the 1.3-2 GeV region where they matter most.
METOOLS::Propagator_Base * FFVMD::BuildEE3PiRhoTower() {
  METOOLS::Summed_Propagator * rho = new METOOLS::Summed_Propagator(true);
  rho->Add(new METOOLS::GounarisSakurai(
             LineShapes->Get(Flavour(kf_rho_770)),
             resonance_type::GS,
             PionFF::Rho().Mass(),
             PionFF::Rho().Width()),
           Complex(1.,0.));
  rho->Add(new METOOLS::GounarisSakurai(
             LineShapes->Get(Flavour(kf_rho_1450)),
             resonance_type::GS,
             PionFF::Rho1450().Mass(),
             PionFF::Rho1450().Width()),
           PionFF::Rho1450().Weight());
  rho->Add(new METOOLS::GounarisSakurai(
             LineShapes->Get(Flavour(kf_rho_1700)),
             resonance_type::GS,
             PionFF::Rho1700().Mass(),
             PionFF::Rho1700().Width()),
           PionFF::Rho1700().Weight());
  return rho;
}

Complex FFVMD::EE3PiFF() {
  Vec4D pp(0.,0.,0.,0.), pm(0.,0.,0.,0.), p0(0.,0.,0.,0.);
  bool hp(false), hm(false), h0(false);
  for (size_t i(0); i<p_v->J().size(); ++i) {
    const Flavour & fl(p_v->J(i)->Flav());
    if (!fl.IsHadron()) continue;
    const Vec4D p(p_v->J(i)->P());
    if      (fl.Kfcode()==kf_pi)                      { p0=p; h0=true; }
    else if (fl.Kfcode()==kf_pi_plus && fl.IntCharge()>0) { pp=p; hp=true; }
    else if (fl.Kfcode()==kf_pi_plus && fl.IntCharge()<0) { pm=p; hm=true; }
  }
  // The outgoing leg is not among p_v->J(); recover it from the others.
  if (!(hp&&hm&&h0)) {
    const Flavour & fc(p_v->JC()->Flav());
    const Vec4D pc(p_v->JC()->P());
    if      (!h0 && fc.Kfcode()==kf_pi)      { p0=pc; h0=true; }
    else if (!hp && fc.IntCharge()>0)        { pp=pc; hp=true; }
    else if (!hm && fc.IntCharge()<0)        { pm=pc; hm=true; }
  }
  // The vertex is declared {gamma, pi+, pi-, pi0}, so all three pions are
  // always present and this never fires - but returning a bare 1 for a
  // dimensionful (GeV^-3) form factor would be a silently-wrong number that
  // does not crash, so fail loudly instead of guessing.
  if (!(hp&&hm&&h0))
    THROW(fatal_error,"EE3PiFF: could not identify all three pions on the "
	  "gamma-pi+-pi--pi0 vertex.");

  const double q2((pp+pm+p0).Abs2());
  const double s12((pp+pm).Abs2()), s13((pp+p0).Abs2()), s23((pm+p0).Abs2());
  const Complex rho = (*p_rho3pi)(s12)+(*p_rho3pi)(s13)+(*p_rho3pi)(s23);
  Complex F(m_alphaA*(*p_iso3pi)(q2)*rho);
  // eq. (10) terms E and F: each is ONE isoscalar pole times its OWN
  // excited-rho sub-current, summed over the three pion pairings.  Present
  // only for the eq-10 isoscalar; they fill its 2 GeV interference trough.
  if (p_isoE!=NULL)
    F += m_cE*(*p_isoE)(q2)*((*p_rhoE)(s12)+(*p_rhoE)(s13)+(*p_rhoE)(s23));
  if (p_isoF!=NULL)
    F += m_cF*(*p_isoF)(q2)*((*p_rhoF)(s12)+(*p_rhoF)(s13)+(*p_rhoF)(s23));
  return F;
}

void FFVMD::ConstructVectorPionFormFactor() {
  msg_Debugging()<<METHOD<<": three-pion isoscalar parameters\n";
  ThreePionFF::ReadAll();   // omega(782), phi(1020), J/psi

  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"] };
  if (s["Subtracted"].SetDefault(false).Get<bool>()) {
    ConstructVectorPionFormFactorSubtracted();
    return;
  }

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
  // The two higher isoscalars of hep-ph/0512180 Eq. (10) - see the note by
  // their parameters in FormFactorParameters.C. They enter exactly like the
  // omega and phi above, so this stays a sum of isoscalars in q^2 multiplying
  // the rho chain built in ConstructThreePionFormFactor.
  //
  // NOTE what this does NOT add. Eq. (10) also has two terms in which an
  // EXCITED rho chain is paired with one specific isoscalar - rho(1450) with
  // phi(1020), rho(1700) with omega(1631), coefficients E and F. That is a sum
  // of products, whereas this class computes a product of two sums, so those
  // two terms cannot be expressed here by adding parameters: a rho added to
  // the chain would multiply every isoscalar rather than just its partner.
  // Reproducing the full Eq. (10) needs FF() restructured into a sum over
  // (isoscalar x its own rho chain) terms.
  //
  // Measured against Belle 2024 (BELLE_2024_I2775022, chi2/ndf on d02-x01-y03,
  // 4000 events per energy point):
  //
  //                       without   with   x0.5
  //     1.85-3 GeV tail     419.2    5.2   81.4
  //     1.1-2.0 GeV          68.1  616.0  150.8
  //     1-3 GeV             216.6  398.4  135.4
  //
  // So these two do fix what they were added for - the tail runs away to
  // 3-10 x data without them and lands on 0.9-1.6 x with them, which says the
  // masses, widths and shape are right.  The 1.1-2.0 GeV damage is the missing
  // E and F: C and D were fitted TOGETHER with them and lose a cancellation on
  // their own.  Halving both amplitudes is the best single number over
  // 1-3 GeV but trades the tail straight back, so it is a compromise and not a
  // fix - the fix is structural.
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,
              ThreePionFF::Omega1375().Mass(),
              ThreePionFF::Omega1375().Width()),
            ThreePionFF::Omega1375().Weight());
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,
              ThreePionFF::Omega1631().Mass(),
              ThreePionFF::Omega1631().Width()),
            ThreePionFF::Omega1631().Weight());
}

// hep-ph/0512180 Eq. (10) term E: phi(1020) x H_rho(1450). This is the
// gamma -> rho(1450) pi production vertex's own form factor - a SINGLE
// isoscalar pole (phi(1020)'s own mass/width), not the full isoscalar sum
// V_pi carries, because E multiplies ONLY the rho(1450) sub-current, not
// every isoscalar. See the note in ConstructVectorPionFormFactor above and
// the ThreePionFF::Phi1020RhoE() note in FormFactorParameters.C.
//
// Uses HepPh0512180_BreitWigner (Propagator.H), NOT FixedBreitWigner inside
// m_props: FixedBreitWigner(0)=1 by construction and Summed_Propagator then
// divides by the weight sum, forcing F(0)=1 - but eq. (10)'s own BW is NOT
// normalised to 1 at s=0 (about -1 there for a narrow resonance) and must be
// used raw. Confirmed by a direct 3-way comparison against the
// tau-decays-in-sherpa fork's own eq-10 output: the FixedBreitWigner/
// Summed_Propagator version inflated the pi+pi-pi0 cross section 57x at
// 3 GeV against BELLE data, while the fork's raw formula tracks BELLE's
// 2-3.5 GeV tail well. p_raw_prop is evaluated directly in FF(), bypassing
// m_props entirely for this mode.
void FFVMD::ConstructVectorPionFormFactorRhoE() {
  msg_Debugging()<<METHOD<<": hep-ph/0512180 Eq.10 term E parameters\n";
  ThreePionFF::ReadAll();
  p_raw_prop = new METOOLS::HepPh0512180_BreitWigner(
              ThreePionFF::Phi1020RhoE().Mass(),
              ThreePionFF::Phi1020RhoE().Width());
  m_raw_weight = ThreePionFF::Phi1020RhoE().Weight();
}

// hep-ph/0512180 Eq. (10) term F: omega(1631) x H_rho(1700). Same structure
// as term E above, but for the rho(1700) sub-current with its own partner
// isoscalar, omega(1631).
void FFVMD::ConstructVectorPionFormFactorRhoF() {
  msg_Debugging()<<METHOD<<": hep-ph/0512180 Eq.10 term F parameters\n";
  ThreePionFF::ReadAll();
  p_raw_prop = new METOOLS::HepPh0512180_BreitWigner(
              ThreePionFF::Omega1631RhoF().Mass(),
              ThreePionFF::Omega1631RhoF().Width());
  m_raw_weight = ThreePionFF::Omega1631RhoF().Weight();
}

Complex FFVMD::FF() {
  if (m_mode==FF_0_PP_mode::ee3pi) return EE3PiFF();
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
  if (p_raw_prop!=NULL) return m_raw_weight*(*p_raw_prop)(Q2);
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