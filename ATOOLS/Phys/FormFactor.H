#ifndef ATOOLS_Phys_FormFactor_H
#define ATOOLS_Phys_FormFactor_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

struct ff {
  enum mode { off = 0, pion = 1, kplus = 2, kls = 3 };
  enum type { factored = 0, vertex = 1 };
};

std::istream &operator>>(std::istream &str, ff::mode &sm);
std::ostream &operator<<(std::ostream &str, const ff::mode &sm);
std::istream &operator>>(std::istream &str, ff::type &sm);
std::ostream &operator<<(std::ostream &str, const ff::type &sm);

class FormFactor {
public:
  FormFactor();
  ~FormFactor();

  void RegisterDefaultsPion();
  void RegisterDefaultsKaon();
  // void Initalize();
  double Eval(const double &q2);
  void Calculate(const double &q2);

  Complex ppi(const double &q2);
  Complex d(const double &m);
  Complex h(const double &q2);
  Complex dh(const double &q2);
  inline Complex Form() { return m_form; }

  Complex f(const double &q2, const double &m, const double &g);
  Complex gamma(const double &q2, const double &m, const double &g);

  inline bool On() { return m_form_mode != ff::off; }
  inline ff::type Type() { return m_form_type ; }

  Complex b(const double &q2, const Flavour &fl);
  Complex BW_GS(const double &q2, const Flavour &fl);
  Complex Pion(const double &q2);
  Complex KPlus(const double &q2);
  void Init();
  double m_crho, m_crhop, m_crhopp, m_crhoppp;
  double m_comega, m_comegap, m_comegapp;
  double m_cphi, m_cphip, m_cphipp;
  double m_prho, m_prhop, m_prhopp, m_prhoppp, m_pomega, m_pphi;
  double m_omega_mass, m_omega_g, m_omegap_mass, m_omegap_g, m_omegapp_mass,
      m_omegapp_g;
  double m_phi_mass, m_phi_g, m_phip_mass, m_phip_g;
  double m_rho_mass, m_rho_g, m_rhop_mass, m_rhop_g, m_rhopp_mass, m_rhopp_g,
      m_rhoppp_mass, m_rhoppp_g;
  bool m_debug, m_initialize;
  double m_xmin, m_xmax;
  Flavour m_flv, m_rho, m_rhop, m_rhopp, m_rhoppp;
  Flavour m_omega, m_omegap, m_omegapp;
  Flavour m_phi, m_phip;

  Particle_Info m_rhoi, m_rhopi, m_rhoppi, m_rhopppi;
  Particle_Info m_omegai, m_omegapi, m_omegappi;
  Particle_Info m_phii, m_phipi;

  Complex m_I, m_form;
  ff::mode m_form_mode;
  ff::type m_form_type;
};
extern FormFactor *pionform;
inline std::string to_string(ff::mode fs) {
  std::ostringstream os;
  os << fs;
  return os.str();
}
inline std::string to_string(ff::type fs) {
  std::ostringstream os;
  os << fs;
  return os.str();
}

} // namespace ATOOLS

#endif
