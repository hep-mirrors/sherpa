#include "MODEL/SM/ThreePion_FittedCurrent.H"
#include "METOOLS/FormFactors/Line_Shapes.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include <cmath>
#include <set>
#include <sstream>

using namespace METOOLS;
using namespace ATOOLS;

struct Fitted_ThreePion_Current::Impl {
  /*
    M^2/(M^2 - s - i M Gamma) at every s, s = 0 included, as FixedBreitWigner
    of the fitting branch. (METOOLS::FixedBreitWigner returns 1 at s = 0, which
    differs by 1/(1 - i Gamma/M) - 15% for the omega(1631).)
  */
  class Constant_Width_BreitWigner: public Propagator_Base {
  public:
    Constant_Width_BreitWigner(const double & M,const double & Gamma) :
      Propagator_Base(NULL,resonance_type::fixed,M,Gamma) {}
    const Complex operator()(const double & s)
    { return m_M2/Complex(m_M2-s,-m_M*m_Gamma); }
    const Complex Normalised(const double & s)  { return (*this)(s); }
    const double  Normalised2(const double & s) { return norm((*this)(s)); }
  };

  /*
    Gounaris-Sakurai (Phys. Rev. Lett. 21 (1968) 244) with the imaginary part
    M Gamma(s), GounarisSakuraiM of the fitting branch
    (METOOLS/HadronCurrents/FormFactors/Propagator.C there):
      m_num/(M^2 - s + f(s) - i M Gamma(s)),
      Gamma(s) = Gamma (M/sqrt(s)) (k/k0)^3,  k = sqrt(s/4 - m_pi^2),
    real-valued below threshold (k = 0), m_num = M^2 (1 + d Gamma/M).
  */
  class Gounaris_Sakurai_MGamma {
    double m_M, m_M2, m_G0, m_mpi, m_mpi2;
    double m_k0, m_h0, m_dh0, m_d, m_num;
    double K(const double & s) const {
      if (s<=4.*m_mpi2) return 0.;
      return 0.5*sqrt(s-4.*m_mpi2);
    }
    double H(const double & s) const {
      const double k(K(s));
      if (k<=0. || s<=0.) return 0.;
      const double rs(sqrt(s));
      return (2./M_PI)*(k/rs)*log((rs+2.*k)/(2.*m_mpi));
    }
    double F(const double & s) const {
      const double k(K(s));
      return m_G0*m_M2/pow(m_k0,3.)*(sqr(k)*(H(s)-m_h0)+(m_M2-s)*sqr(m_k0)*m_dh0);
    }
    double Width(const double & s) const {
      const double k(K(s));
      if (k<=0. || s<=0.) return 0.;
      return m_G0*(m_M/sqrt(s))*pow(k/m_k0,3.);
    }
  public:
    Gounaris_Sakurai_MGamma(const double & M,const double & G0,const double & mpi) :
      m_M(M), m_M2(sqr(M)), m_G0(G0), m_mpi(mpi), m_mpi2(sqr(mpi))
    {
      m_k0  = K(m_M2);
      m_h0  = H(m_M2);
      m_dh0 = m_h0*(1./(8.*sqr(m_k0))-1./(2.*m_M2))+1./(2.*M_PI*m_M2);
      m_d   = (3./M_PI)*(m_mpi2/sqr(m_k0))*log((m_M+2.*m_k0)/(2.*m_mpi))
            + m_M/(2.*M_PI*m_k0) - (m_mpi2*m_M)/(M_PI*pow(m_k0,3.));
      m_num = m_M2*(1.+m_d*m_G0/m_M);
    }
    Complex operator()(const double & s) const {
      return m_num/Complex(m_M2-s+F(s),-m_M*Width(s));
    }
  };

  // the couplings A - F of hep-ph/0512180 eq. (10), and the isoscalar poles
  enum class Coupling { A = 0, B, C, D, E, F };
  enum class Isoscalar { omega782 = 0, phi1020, omega1375, omega1631 };
  struct Pole { std::string name; double mass, width; };

  double coupling[6];
  Pole   pole[4];
  double rho770_mass, rho770_width, pion_mass;
  // rho-omega mixing on the omega term (Fit: RhoOmegaMixing), off for epsilon = 0
  double mix_eps, mix_kappa, mix_rho_mass, mix_rho_width;
  std::unique_ptr<Propagator_Base> iso[3];
  std::unique_ptr<Gounaris_Sakurai_MGamma> rho770;
  // rho(1450), rho(1700): [neutral, charged]
  std::unique_ptr<BreitWigner> rho1450[2], rho1700[2];
  std::string description;

  double C(const Coupling & c) const { return coupling[size_t(c)]; }
  const Pole & P(const Isoscalar & v) const { return pole[size_t(v)]; }

  //! a THREE_PION_FORM_FACTOR key that configures the LowEnergy model's own
  //! current and is set to something the fitted current cannot honour
  static bool ConflictsWithFittedCurrent(Scoped_Settings & s,const std::string & key);
  static void CheckSettings();
};

bool Fitted_ThreePion_Current::Impl::
ConflictsWithFittedCurrent(Scoped_Settings & s,const std::string & key) {
  // keys the fitted current reads, or that do not touch its shape (the chain
  // couplings cancel in it; Contact is refused by the model; RhoOmegaMixing
  // adds its own term)
  static const std::set<std::string> compatible
    { "Current", "Fit", "Contact", "RhoOmegaMixing", "F_pi",
      "g_gammarhopi", "g_rhopipi", "g_rho1450pipi", "g_rho1700pipi" };
  if (compatible.count(key)) return false;
  // switches that are harmless when off
  static const std::set<std::string> off_switches
    { "Subtracted", "SumRule", "EqTen" };
  if (off_switches.count(key)) return s[key].SetDefault(false).Get<bool>();
  if (key=="RhoWidth") return s[key].SetDefault("Fixed").Get<std::string>()!="Fixed";
  if (key=="rho0_isospin") return s[key].SetDefault(1.).Get<double>()!=1.;
  // every other key is a parameter of the LowEnergy current (its isoscalar
  // poles and amplitudes, the sum-rule polynomial, the contact isoscalar)
  return true;
}

void Fitted_ThreePion_Current::Impl::CheckSettings() {
  Settings & main(Settings::GetMainSettings());
  if (main["THREE_PION_CONTACT"].SetDefault(false).Get<bool>())
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: Current: Fitted needs the "
          "gamma* -> 3 pi current with its epsilon structure; THREE_PION_CONTACT "
          "builds the old scalar vertex instead.");
  Scoped_Settings s{ main["THREE_PION_FORM_FACTOR"] };
  std::string bad;
  for (const std::string & key : s.GetKeys())
    if (ConflictsWithFittedCurrent(s,key)) bad += (bad.empty()?"":", ")+key;
  if (!bad.empty())
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: Current: Fitted brings its own "
          "isoscalar poles, couplings and rho line shapes (parameters under Fit); "
          "these keys configure the LowEnergy current and would be ignored or "
          "applied twice: "+bad+". Remove them.");
}

bool Fitted_ThreePion_Current::Selected() {
  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"] };
  const std::string name(s["Current"].SetDefault("LowEnergy").Get<std::string>());
  if (name!="LowEnergy" && name!="Fitted")
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: Current must be LowEnergy or "
          "Fitted, not '"+name+"'.");
  return name=="Fitted";
}

Fitted_ThreePion_Current & Fitted_ThreePion_Current::Instance() {
  static Fitted_ThreePion_Current current;
  return current;
}

Fitted_ThreePion_Current::Fitted_ThreePion_Current() : p_impl(new Impl) {
  Impl::CheckSettings();
  Impl & d(*p_impl);
  // the fit of the update_decay_currents branch: couplings A - F refitted
  // there against Belle II 2024 below 2.5 GeV with exactly the line shapes of
  // this file (commit 16462ab6e, FF_0_Isoscalar3Pi::FixParameters) [GeV^-3].
  // The published values of C, D, F (-0.77, -1.12, -0.59) belong to the
  // paper's own propagators.
  const char * const coupling_names[6] = { "A", "B", "C", "D", "E", "F" };
  const double fitted_couplings[6] = { 18.20, -0.87, -0.5068, -1.1262, -0.1954, -0.6368 };
  // the isoscalar poles, constant-width Breit-Wigners: omega(782) and
  // phi(1020) at the particle-table masses of the fitting branch with the
  // widths it gave them, omega' and omega'' at hep-ph/0512180's own states
  const Impl::Pole fitted_poles[4] = { { "omega(782)",  0.78194, 0.00868  },
                                       { "phi(1020)",   1.01941, 0.004249 },
                                       { "omega(1375)", 1.375,   0.250    },
                                       { "omega(1631)", 1.631,   0.245    } };
  Scoped_Settings fit{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"]["Fit"] };
  for (size_t i(0);i<6;++i)
    d.coupling[i] = fit[coupling_names[i]].SetDefault(fitted_couplings[i]).Get<double>();
  for (size_t i(0);i<4;++i) {
    d.pole[i] = fitted_poles[i];
    Scoped_Settings p{ fit[fitted_poles[i].name] };
    d.pole[i].mass  = p["Mass"].SetDefault(fitted_poles[i].mass).Get<double>();
    d.pole[i].width = p["Width"].SetDefault(fitted_poles[i].width).Get<double>();
  }
  // rho(770): the fitting branch's particle-table rho(770) and pi+, for both
  // charges (the two have the same entry there)
  Scoped_Settings r{ fit["rho(770)"] };
  d.rho770_mass  = r["Mass"].SetDefault(0.77).Get<double>();
  d.rho770_width = r["Width"].SetDefault(0.1507).Get<double>();
  d.pion_mass    = 0.13957;
  // arXiv:2307.02546 eq. (31) in narrow-width form: omega -> omega [1 + eps
  // kappa s/(M_rho^2 - s - i M_rho Gamma_rho)], kappa = (g_omegagamma/
  // g_rhogamma)^2 = (16.2/4.9)^2 (text after their eq. 30)
  Scoped_Settings mix{ fit["RhoOmegaMixing"] };
  d.mix_eps       = mix["Epsilon"].SetDefault(0.).Get<double>();
  d.mix_kappa     = mix["Kappa"].SetDefault(sqr(16.2/4.9)).Get<double>();
  d.mix_rho_mass  = mix["Mass"].SetDefault(0.77526).Get<double>();
  d.mix_rho_width = mix["Width"].SetDefault(0.1491).Get<double>();
  // the two rho-omega mixing terms describe the same physics
  const bool chainmix(Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"]
                      ["RhoOmegaMixing"]["Active"].SetDefault(false).Get<bool>());
  if (d.mix_eps!=0. && chainmix)
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: RhoOmegaMixing: {Active: true} and "
          "Fit: RhoOmegaMixing: {Epsilon != 0} both add rho-omega mixing; switch "
          "one of them off.");

  for (size_t k(0);k<3;++k)
    d.iso[k].reset(NewIsoscalar(Tower(k),Complex(1.,0.)));
  d.rho770.reset(new Impl::Gounaris_Sakurai_MGamma(d.rho770_mass,d.rho770_width,
                                                   d.pion_mass));

  if (LineShapes==NULL)
    THROW(fatal_error,"Fitted_ThreePion_Current: no line-shape registry; the "
          "fitted current needs the LowEnergy model.");
  const kf_code rho_kf[2][2] = { { kf_rho_1450, kf_rho_1450_plus },
                                 { kf_rho_1700, kf_rho_1700_plus } };
  for (size_t c(0);c<2;++c) {
    Total_Width_Base * w1450(LineShapes->Get(Flavour(rho_kf[0][c])));
    Total_Width_Base * w1700(LineShapes->Get(Flavour(rho_kf[1][c])));
    if (w1450==NULL || w1700==NULL)
      THROW(fatal_error,"Fitted_ThreePion_Current: the line-shape registry has no "
            "running width for the rho(1450) or rho(1700).");
    d.rho1450[c].reset(new BreitWigner(w1450,resonance_type::running));
    d.rho1700[c].reset(new BreitWigner(w1700,resonance_type::running));
  }

  std::ostringstream str;
  str<<"fitted current (hep-ph/0512180 eq. (10), update_decay_currents fit): A-F =";
  for (size_t i(0);i<6;++i) str<<" "<<d.coupling[i];
  str<<"; constant-width isoscalars";
  for (size_t i(0);i<4;++i)
    str<<" "<<d.pole[i].name<<" ["<<d.pole[i].mass<<", "<<d.pole[i].width<<"]";
  str<<"; rho(770) Gounaris-Sakurai (M Gamma(s)) ["<<d.rho770_mass<<", "<<d.rho770_width
     <<"], rho(1450)/rho(1700) running registry widths";
  if (d.mix_eps!=0.)
    str<<"; rho-omega mixing eps = "<<d.mix_eps<<", kappa = "<<d.mix_kappa
       <<", rho ["<<d.mix_rho_mass<<", "<<d.mix_rho_width<<"]";
  d.description = str.str();
  msg_Info()<<"Fitted_ThreePion_Current: "<<d.description<<"\n";
}

Fitted_ThreePion_Current::~Fitted_ThreePion_Current() {}

Propagator_Base * Fitted_ThreePion_Current::NewIsoscalar(const Tower & k,
                                                         const Complex & scale) const {
  typedef Impl::Isoscalar Iso;
  typedef Impl::Coupling  Cpl;
  const Impl & d(*p_impl);
  Summed_Propagator * iso(new Summed_Propagator(false));
  auto add = [&](const Iso & v,const Cpl & c) {
    iso->Add(new Impl::Constant_Width_BreitWigner(d.P(v).mass,d.P(v).width),scale*d.C(c));
  };
  switch (k) {
  case Tower::rho770:
    if (d.mix_eps==0.) add(Iso::omega782,Cpl::A);
    else {
      /*
        rho-omega mixing with a constant-width rho keeps the factor a sum of
        Breit-Wigners: with L = M^2 - i M Gamma,
          BW_om [1 + eps kappa s/(L_rho - s)]
            = (1 - eps kappa L_om/(L_om - L_rho)) BW_om
              + eps kappa M_om^2 L_rho/(M_rho^2 (L_om - L_rho)) BW_rho
      */
      const Impl::Pole & om(d.P(Iso::omega782));
      const double Mr(d.mix_rho_mass), Gr(d.mix_rho_width);
      const Complex Lom(sqr(om.mass),-om.mass*om.width), Lrho(sqr(Mr),-Mr*Gr);
      const Complex ek(d.mix_eps*d.mix_kappa);
      iso->Add(new Impl::Constant_Width_BreitWigner(om.mass,om.width),
               scale*d.C(Cpl::A)*(1.-ek*Lom/(Lom-Lrho)));
      iso->Add(new Impl::Constant_Width_BreitWigner(Mr,Gr),
               scale*d.C(Cpl::A)*ek*sqr(om.mass)*Lrho/(sqr(Mr)*(Lom-Lrho)));
    }
    add(Iso::phi1020,Cpl::B);
    add(Iso::omega1375,Cpl::C);
    add(Iso::omega1631,Cpl::D);
    break;
  case Tower::rho1450:
    add(Iso::phi1020,Cpl::E);
    break;
  case Tower::rho1700:
    add(Iso::omega1631,Cpl::F);
    break;
  }
  return iso;
}

Complex Fitted_ThreePion_Current::RhoLineShape(const Tower & k,const bool & charged,
                                               const double & s) {
  Impl & d(*p_impl);
  const size_t c(charged ? 1 : 0);
  switch (k) {
  case Tower::rho770:  return (*d.rho770)(s);
  case Tower::rho1450: return (*d.rho1450[c])(s);
  case Tower::rho1700: return (*d.rho1700[c])(s);
  }
  return Complex(0.,0.);
}

Complex Fitted_ThreePion_Current::operator()(const double & q2,const double & spm,
                                             const double & sp0,const double & sm0) {
  // sum over the towers of the isoscalar factor times H_k, the rho line shape
  // over the three pion pairs (neutral rho at s+-, charged at s+0, s-0)
  Complex res(0.,0.);
  for (size_t k(0);k<3;++k)
    res += (*p_impl->iso[k])(q2)*
      (RhoLineShape(Tower(k),false,spm)+RhoLineShape(Tower(k),true,sp0)+
       RhoLineShape(Tower(k),true,sm0));
  return res;
}

std::string Fitted_ThreePion_Current::Describe() const { return p_impl->description; }

// ----------------------------------------------------------------------
// FFThreePionFittedRho: the fitted current's rho line shapes in the rho pi chain
// ----------------------------------------------------------------------
// On the rho -> pi pi vertices (SSV) of the rho(770), rho(1450), rho(1700),
// with THREE_PION_FORM_FACTOR: {Current: Fitted}. Comix propagates the rho with
// a fixed width, i/(m^2 - s - i m Gamma); this multiplies it by
//   R(s) = LS_k(s) (m^2 - s - i m Gamma)/M_k^2,
// LS_k the fitted current's line shape and M_k the neutral rho_k mass of the
// production normalisation N_k (LowEnergy_Model::InitThreePionVertices), so that
// every chain becomes (i/M_k^2) LS_k(s) with the production form factor
// carrying the isoscalar factor of its tower (FFVMD, V_pi modes).
namespace METOOLS {
  class FFThreePionFittedRho: public Form_Factor {
  private:
    kf_code m_rho;
    Fitted_ThreePion_Current::Tower m_tower;
    bool    m_charged;
    double  m_m2, m_mgamma, m_norm2;
  public:
    FFThreePionFittedRho(const Vertex_Key &key);
    Complex FF();
  };
}

FFThreePionFittedRho::FFThreePionFittedRho(const Vertex_Key &key):
  Form_Factor("ThreePionFittedRho",key), m_rho(0),
  m_tower(Fitted_ThreePion_Current::Tower::rho770),
  m_charged(false), m_m2(0.), m_mgamma(0.), m_norm2(1.)
{
  for (size_t i(0); i<key.p_mv->in.size(); ++i) {
    const Flavour &fl(key.p_mv->in[i]);
    if (fl.IntSpin()!=2) continue;
    m_rho = fl.Kfcode();
    m_charged = fl.IntCharge()!=0;
    m_m2 = sqr(fl.Mass());
    m_mgamma = fl.Mass()*fl.Width();
  }
  kf_code neutral(0);
  switch (m_rho) {
  case kf_rho_770:  case kf_rho_770_plus:
    m_tower = Fitted_ThreePion_Current::Tower::rho770;  neutral = kf_rho_770;  break;
  case kf_rho_1450: case kf_rho_1450_plus:
    m_tower = Fitted_ThreePion_Current::Tower::rho1450; neutral = kf_rho_1450; break;
  case kf_rho_1700: case kf_rho_1700_plus:
    m_tower = Fitted_ThreePion_Current::Tower::rho1700; neutral = kf_rho_1700; break;
  default:
    THROW(fatal_error,"ThreePionFittedRho form factor on a vertex without a "
          "rho(770), rho(1450) or rho(1700).");
  }
  m_norm2 = sqr(Flavour(neutral).Mass());
}

Complex FFThreePionFittedRho::FF() {
  const Current *rho(NULL);
  for (size_t i(0); i<p_v->J().size() && rho==NULL; ++i)
    if (p_v->J(i)->Flav().Kfcode()==m_rho) rho = p_v->J(i);
  if (rho==NULL) rho = p_v->JC();
  const double s(rho->P().Abs2());
  return Fitted_ThreePion_Current::Instance().RhoLineShape(m_tower,m_charged,s)*
    Complex(m_m2-s,-m_mgamma)/m_norm2;
}

DECLARE_GETTER(FFThreePionFittedRho,"ThreePionFittedRho",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFThreePionFittedRho>::
operator()(const Vertex_Key &args) const { return new FFThreePionFittedRho(args); }
void ATOOLS::Getter<Form_Factor,Vertex_Key,FFThreePionFittedRho>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"rho -> pi pi: the fitted three-pion current's rho line shape"; }
