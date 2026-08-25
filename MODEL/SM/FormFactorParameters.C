#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "MODEL/SM/FormFactorParameters.H"
#include <iomanip>
#include <sstream>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

namespace METOOLS {

  const int s_ffwidth = 80;

  // ------------------------------------------------------------------
  // ResonanceParameters implementation
  // ------------------------------------------------------------------
  ResonanceParameters::ResonanceParameters(const string& block,
                                           const string& tag, double m,
                                           double Gamma, double c, double phi,
                                           ResonanceType type, kf_code flav)
      : m_block(block), m_sub(""), m_tag(tag), m_m0(m), m_G0(Gamma), m_c0(c),
        m_phi0(phi), m_m(m), m_G(Gamma), m_c(c), m_phi(phi), m_type(type),
        m_flav(flav), m_shown(false)
  {
  }

  ResonanceParameters::ResonanceParameters(const string& block,
                                           const string& sub, const string& tag,
                                           double m, double Gamma, double c,
                                           double phi, ResonanceType type,
                                           kf_code flav)
      : m_block(block), m_sub(sub), m_tag(tag), m_m0(m), m_G0(Gamma), m_c0(c),
        m_phi0(phi), m_m(m), m_G(Gamma), m_c(c), m_phi(phi), m_type(type),
        m_flav(flav), m_shown(false)
  {
  }

  void ResonanceParameters::Read()
  {
    Settings& main = Settings::GetMainSettings();
    Scoped_Settings s =
        m_sub.empty() ? main[m_block][m_tag] : main[m_block][m_sub][m_tag];
    m_m = s["Mass"].SetDefault(m_m0).Get<double>();
    m_G = s["Width"].SetDefault(m_G0).Get<double>();
    m_c = s["Amplitude"].SetDefault(m_c0).Get<double>();
    m_phi = s["Phase"].SetDefault(m_phi0).Get<double>();
    if (m_shown) return;
    m_shown = true;
    PrintLine();
  }

  void ResonanceParameters::PrintLine() const
  {
    ostringstream line;
    line << setw(16) << m_tag << setw(14) << m_m << setw(14) << m_G << setw(16)
         << m_c << setw(14) << m_phi;
    msg_Out() << Frame_Line(line.str(), s_ffwidth);
  }

  // ------------------------------------------------------------------
  // ParameterRegistry
  // ------------------------------------------------------------------
  ParameterRegistry& ParameterRegistry::Instance()
  {
    static ParameterRegistry inst;
    return inst;
  }

  
void ParameterRegistry::Add(const std::string& sector,
                            ResonanceParameters* params)
{
  if (params != nullptr)
    m_data[sector].push_back(params);
}

const std::vector<ResonanceParameters*>&
ParameterRegistry::Get(const std::string& sector) const
{
  static const std::vector<ResonanceParameters*> empty;

  auto it = m_data.find(sector);

  if (it == m_data.end())
    return empty;

  return it->second;
}

#define DEFINE_PARAM(name, block, sub, tag, m, G, c, phi, type, flav) \
    static ResonanceParameters s_##name(block, sub, tag, m, G, c, phi, type, flav)
  // ---------- Pion form factor ----------
DEFINE_PARAM(rho,        "PION_FORM_FACTOR", "",               "rho(770)",   0.77456, 0.14832, 1.0,     0.0,     ResonanceType::GounarisSakurai, kf_rho_770);
DEFINE_PARAM(rho1450,    "PION_FORM_FACTOR", "",               "rho(1450)",  1.4859,  0.37360, 0.14104, 3.7797,  ResonanceType::GounarisSakurai, kf_rho_1450);
DEFINE_PARAM(rho1700,    "PION_FORM_FACTOR", "",               "rho(1700)",  1.8668,  0.30334, 0.0614,  1.429,   ResonanceType::GounarisSakurai, kf_rho_1700);
DEFINE_PARAM(rho2150,    "PION_FORM_FACTOR", "",               "rho(2150)",  2.2645,  0.11327, 0.0047,  0.921,   ResonanceType::GounarisSakurai, kf_rho_2150);
DEFINE_PARAM(omega782,   "PION_FORM_FACTOR", "",               "omega(782)", 0.78248, 0.00855, 0.00158, 0.075,   ResonanceType::FixedBreitWigner, kf_omega_782);
DEFINE_PARAM(phi1020,    "PION_FORM_FACTOR", "",               "phi(1020)",  1.01947, 0.00425, 0.00045, 2.888,   ResonanceType::FixedBreitWigner, kf_phi_1020);

// ---------- Kaon form factor ----------
DEFINE_PARAM(K_rho,       "KAON_FORM_FACTOR", "rho(770)",      "rho(770)",   0.77456, 0.14832,  1.195,     0.0,     ResonanceType::GounarisSakurai, kf_rho_770);
DEFINE_PARAM(K_rho1450,   "KAON_FORM_FACTOR", "rho(1450)",     "rho(1450)",  1.4859,  0.37360, -0.112,    M_PI,    ResonanceType::GounarisSakurai, kf_rho_1450);
DEFINE_PARAM(K_rho1700,   "KAON_FORM_FACTOR", "rho(1700)",     "rho(1700)",  1.8668,  0.30334, -0.083,    M_PI,    ResonanceType::GounarisSakurai, kf_rho_1700);
DEFINE_PARAM(K_rho2150,   "KAON_FORM_FACTOR", "rho(2150)",     "rho(2150)",  2.2645,  0.11327, -4.22045e-3,0.0,     ResonanceType::GounarisSakurai, kf_rho_2150);
DEFINE_PARAM(K_omega,     "KAON_FORM_FACTOR", "omega(782)",    "omega(782)", 0.78248, 0.00855,  1.195,     0.0,     ResonanceType::FixedBreitWigner, kf_omega_782);
DEFINE_PARAM(K_omega1420, "KAON_FORM_FACTOR", "omega(1420)",   "omega(1420)",1.410,   0.290,   -0.112,    M_PI,    ResonanceType::FixedBreitWigner, kf_omega_1420);
// DEFINE_PARAM(K_omega1650, "KAON_FORM_FACTOR", "omega(1650)",   "omega(1650)",1.67,    0.315,    1-1.195+0.112, M_PI, ResonanceType::FixedBreitWigner, kf_omega_1650);
DEFINE_PARAM(K_phi,       "KAON_FORM_FACTOR", "phi(1020)",     "phi(1020)",  1.01947, 0.00434,  1.018,     0.0,     ResonanceType::FixedBreitWigner, kf_phi_1020);
DEFINE_PARAM(K_phi1680,   "KAON_FORM_FACTOR", "phi(1680)",     "phi(1680)",  1.680,   0.150,    1 - 1.018, M_PI,    ResonanceType::FixedBreitWigner, kf_phi_1680);

// ---------- Three-pion form factor ----------
DEFINE_PARAM(v_omega782,  "THREE_PION_FORM_FACTOR", "omega(782)", "omega(782)", 0.78266,  0.00868, 1.0,    0.0,   ResonanceType::FixedBreitWigner, kf_omega_782);
DEFINE_PARAM(v_phi1020,   "THREE_PION_FORM_FACTOR", "phi(1020)",  "phi(1020)",  1.019461, 0.004249,0.042,  M_PI,  ResonanceType::FixedBreitWigner, kf_phi_1020);
DEFINE_PARAM(v_jpsi,      "THREE_PION_FORM_FACTOR", "J/psi",      "J/psi",      3.0969,   0.0000926,0.001, 0.0,   ResonanceType::FixedBreitWigner, kf_J_psi_1S);
// The two higher isoscalars of hep-ph/0512180 Eq. (10), which carry the
// structure between roughly 1.1 and 2 GeV - without them the three-pion form
// factor has nothing above the phi until the J/psi.
//
// Masses and widths are that paper's Table 1. Its omega(1375) and omega(1631)
// are PDG's omega(1420) and omega(1650), whose kf codes are kf_omega_1420 and
// kf_omega_1600 (30223 - registered in Hadron_Init.C under the name
// "omega(1650)"; the symbol is the older omega(1600) spelling). The kaon
// sector's own omega(1650) line above is commented out because it reaches for
// a kf_omega_1650 that does not exist - kf_omega_1600 is the one to use.
//
// The flavour is only a label here: these are built as
// FixedBreitWigner(NULL, m, Gamma), so the mass and width come from this table
// and not from the particle entry, which differs slightly (1.670/0.31 there
// against this fit's 1.631/0.245) and must not be substituted - the amplitudes
// only mean anything with the masses they were fitted against.
//
// The amplitudes are ITS coefficients C and D divided by its A, because Eq.
// (10) normalises to an absolute A = 18.20 whereas this block is relative to
// omega(782) = 1:
//     C/A = -0.77/18.20 = -0.0423     D/A = -1.12/18.20 = -0.0615
// and, following phi(1020) above, a negative coefficient is written as a
// positive amplitude with phase pi.
//
// Sanity check on the convention: the same reduction applied to that paper's
// B gives B/A = -0.87/18.20 = -0.0478 for phi(1020), against the 0.042 tuned
// here - so the two normalisations do line up.
DEFINE_PARAM(v_omega1375, "THREE_PION_FORM_FACTOR", "omega(1375)","omega(1375)",1.3750,   0.2500,  0.042308, M_PI, ResonanceType::FixedBreitWigner, kf_omega_1420);
DEFINE_PARAM(v_omega1631, "THREE_PION_FORM_FACTOR", "omega(1631)","omega(1631)",1.6310,   0.2450,  0.061538, M_PI, ResonanceType::FixedBreitWigner, kf_omega_1600);

// ---------- Two-photon transition ----------
// pi0
DEFINE_PARAM(T_pi0_rho,   "TWO_PHOTON_FORM_FACTOR", "pi0",        "rho(770)",   0.77526,  0.1474,  1.0,     0.0,     ResonanceType::FixedBreitWigner, kf_rho_770);
DEFINE_PARAM(T_pi0_omega, "TWO_PHOTON_FORM_FACTOR", "pi0",        "omega(782)", 0.78266,  0.00868, 1.0109,  0.3972,  ResonanceType::FixedBreitWigner, kf_omega_782);
DEFINE_PARAM(T_pi0_phi,   "TWO_PHOTON_FORM_FACTOR", "pi0",        "phi(1020)",  1.019461, 0.004249,0.0676, -2.9058, ResonanceType::FixedBreitWigner, kf_phi_1020);
// eta
DEFINE_PARAM(T_eta_rho,   "TWO_PHOTON_FORM_FACTOR", "eta",        "rho(770)",   0.77526,  0.1474,  1.0,     0.0,     ResonanceType::FixedBreitWigner, kf_rho_770);
DEFINE_PARAM(T_eta_omega, "TWO_PHOTON_FORM_FACTOR", "eta",        "omega(782)", 0.78266,  0.00868, 0.0350,  0.5872,  ResonanceType::FixedBreitWigner, kf_omega_782);
DEFINE_PARAM(T_eta_phi,   "TWO_PHOTON_FORM_FACTOR", "eta",        "phi(1020)",  1.019461, 0.004249,0.1135,  2.3182,  ResonanceType::FixedBreitWigner, kf_phi_1020);
// eta'
DEFINE_PARAM(T_etap_rho,  "TWO_PHOTON_FORM_FACTOR", "eta'",       "rho(770)",   0.77526,  0.1474,  1.0,     0.0,     ResonanceType::FixedBreitWigner, kf_rho_770);
DEFINE_PARAM(T_etap_omega,"TWO_PHOTON_FORM_FACTOR", "eta'",       "omega(782)", 0.78266,  0.00868, 1.0,     0.0,     ResonanceType::FixedBreitWigner, kf_omega_782);
DEFINE_PARAM(T_etap_phi,  "TWO_PHOTON_FORM_FACTOR", "eta'",       "phi(1020)",  1.019461, 0.004249,1.0,     0.0,     ResonanceType::FixedBreitWigner, kf_phi_1020);

// ---------- Baryon form factor ----------
DEFINE_PARAM(N1440,    "BARYON_FORM_FACTOR", "N(1440)",      "N(1440)",   1.44,  0.35, 1.0, 0.0, ResonanceType::FixedBreitWigner, kf_none);
DEFINE_PARAM(N1520,    "BARYON_FORM_FACTOR", "N(1520)",      "N(1520)",   1.52,  0.12, 1.0, 0.0, ResonanceType::FixedBreitWigner, kf_none);
DEFINE_PARAM(N1535,    "BARYON_FORM_FACTOR", "N(1535)",      "N(1535)",   1.535, 0.15, 1.0, 0.0, ResonanceType::FixedBreitWigner, kf_none);
DEFINE_PARAM(Delta1232,"BARYON_FORM_FACTOR", "Δ(1232)",      "Δ(1232)",   1.232, 0.118,1.0, 0.0, ResonanceType::FixedBreitWigner, kf_none);
DEFINE_PARAM(Delta1600,"BARYON_FORM_FACTOR", "Δ(1600)",      "Δ(1600)",   1.60,  0.35, 1.0, 0.0, ResonanceType::FixedBreitWigner, kf_none);
// ------------------------------------------------------------------
// Sector accessor implementations
// ------------------------------------------------------------------
#define DEFINE_SECTOR_ACCESSORS(Sector, prefix)                                \
  namespace Sector {                                                           \
    ResonanceParameters& prefix##Rho() { return s_##prefix##rho; }             \
    ResonanceParameters& prefix##Rho1450() { return s_##prefix##rho1450; }     \
    ResonanceParameters& prefix##Rho1700() { return s_##prefix##rho1700; }     \
    ResonanceParameters& prefix##Rho2150() { return s_##prefix##rho2150; }     \
    ResonanceParameters& prefix##Omega782() { return s_##prefix##omega782; }   \
    ResonanceParameters& prefix##Phi1020() { return s_##prefix##phi1020; }     \
    void ReadAll()                                                             \
    {                                                                          \
      prefix##Rho().Read();                                                    \
      prefix##Rho1450().Read();                                                \
      prefix##Rho1700().Read();                                                \
      prefix##Rho2150().Read();                                                \
      prefix##Omega782().Read();                                               \
      prefix##Phi1020().Read();                                                \
    }                                                                          \
  }

  // Pion sector (already has specific names)
  namespace PionFF {
    ResonanceParameters& Rho() { return s_rho; }
    ResonanceParameters& Rho1450() { return s_rho1450; }
    ResonanceParameters& Rho1700() { return s_rho1700; }
    ResonanceParameters& Rho2150() { return s_rho2150; }
    ResonanceParameters& Omega782() { return s_omega782; }
    ResonanceParameters& Phi1020() { return s_phi1020; }
    void ReadAll()
    {
      Rho().Read();
      Rho1450().Read();
      Rho1700().Read();
      Rho2150().Read();
      Omega782().Read();
      Phi1020().Read();
    }
  } // namespace PionFF

  // Kaon sector
  namespace KaonFF {
    ResonanceParameters& Rho() { return s_K_rho; }
    ResonanceParameters& Rho1450() { return s_K_rho1450; }
    ResonanceParameters& Rho1700() { return s_K_rho1700; }
    ResonanceParameters& Rho2150() { return s_K_rho2150; }
    ResonanceParameters& Omega() { return s_K_omega; }
    ResonanceParameters& Omega1420() { return s_K_omega1420; }
    // ResonanceParameters& Omega1650() { return s_K_omega1650; }
    ResonanceParameters& Phi() { return s_K_phi; }
    ResonanceParameters& Phi1680() { return s_K_phi1680; }
    void ReadAll()
    {
      Rho().Read();
      Rho1450().Read();
      Rho1700().Read();
      Rho2150().Read();
      Omega().Read();
      Omega1420().Read();
      Omega1650().Read();
      Phi().Read();
      Phi1680().Read();
    }
  } // namespace KaonFF

  // Three-pion sector
  namespace ThreePionFF {
    ResonanceParameters& Omega782() { return s_v_omega782; }
    ResonanceParameters& Phi1020() { return s_v_phi1020; }
    ResonanceParameters& Jpsi() { return s_v_jpsi; }
    ResonanceParameters& Omega1375() { return s_v_omega1375; }
    ResonanceParameters& Omega1631() { return s_v_omega1631; }
    void ReadAll()
    {
      Omega782().Read();
      Phi1020().Read();
      Jpsi().Read();
      Omega1375().Read();
      Omega1631().Read();
    }
  } // namespace ThreePionFF

  // Two-photon sector
  namespace TwoPhotonFF {
    // pi0
    ResonanceParameters& Pi0Rho() { return s_T_pi0_rho; }
    ResonanceParameters& Pi0Omega() { return s_T_pi0_omega; }
    ResonanceParameters& Pi0Phi() { return s_T_pi0_phi; }
    // eta
    ResonanceParameters& EtaRho() { return s_T_eta_rho; }
    ResonanceParameters& EtaOmega() { return s_T_eta_omega; }
    ResonanceParameters& EtaPhi() { return s_T_eta_phi; }
    // eta'
    ResonanceParameters& EtapRho() { return s_T_etap_rho; }
    ResonanceParameters& EtapOmega() { return s_T_etap_omega; }
    ResonanceParameters& EtapPhi() { return s_T_etap_phi; }
    void ReadAll()
    {
      Pi0Rho().Read();
      Pi0Omega().Read();
      Pi0Phi().Read();
      EtaRho().Read();
      EtaOmega().Read();
      EtaPhi().Read();
      EtapRho().Read();
      EtapOmega().Read();
      EtapPhi().Read();
    }
  } // namespace TwoPhotonFF

  // Baryon sector
  namespace BaryonFF {
    ResonanceParameters& N1440() { return s_N1440; }
    ResonanceParameters& N1520() { return s_N1520; }
    ResonanceParameters& N1535() { return s_N1535; }
    ResonanceParameters& Delta1232() { return s_Delta1232; }
    ResonanceParameters& Delta1600() { return s_Delta1600; }
    void ReadAll()
    {
      N1440().Read();
      N1520().Read();
      N1535().Read();
      Delta1232().Read();
      Delta1600().Read();
    }
  } // namespace BaryonFF

  // ------------------------------------------------------------------
  // Utility functions (identical to before, but using s_ffwidth)
  // ------------------------------------------------------------------
  void FFTableHead(const string& title)
  {
    msg_Out() << Frame_Header(s_ffwidth) << Frame_Line(title, s_ffwidth)
              << Frame_Separator(s_ffwidth);
    ostringstream line;
    line << setw(16) << "Resonance" << setw(14) << "Mass [GeV]" << setw(14)
         << "Width [GeV]" << setw(16) << "Amplitude" << setw(14) << "Phase";
    msg_Out() << Frame_Line(line.str(), s_ffwidth)
              << Frame_Separator(s_ffwidth);
  }

  void FFTableFoot(const string& label, const Complex& f0)
  {
    ostringstream line;
    line << label << " = " << f0;
    msg_Out() << Frame_Separator(s_ffwidth) << Frame_Line(line.str(), s_ffwidth)
              << Frame_Footer(s_ffwidth);
  }

  void CheckFFNormalisation(const string& tag, const Complex& f0,
                            bool requireunity)
  {
    if (IsNan(f0)) {
      msg_Error() << om::red << "Warning in " << tag << ": F(0) = " << f0
                  << " is not a number." << om::reset << "\n";
      return;
    }
    if (!requireunity) {
      msg_Debugging() << "  " << tag << ": F(0) = " << f0
                      << " (not a charge form factor, F(0)=1 not required)\n";
      return;
    }
    const double dev = std::abs(f0 - Complex(1., 0.));
    if (dev > 1e-9)
      msg_Error() << om::red << "Warning in " << tag << ": F(0) = " << f0
                  << ", expected 1 (off by " << dev << ").\n"
                  << "  Charge normalisation is violated, which breaks gauge\n"
                  << "  invariance in processes with a radiated photon."
                  << om::reset << "\n";
  }

} // namespace METOOLS