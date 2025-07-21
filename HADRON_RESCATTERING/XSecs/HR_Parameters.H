#ifndef HADRON_RESCATTERING_Tools_HR_Parameters_H
#define HADRON_RESCATTERING_Tools_HR_Parameters_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <list>
#include <string>

namespace HADRON_RESCATTERING {
  enum class HR_Res {
    none    = 0,
    pip_pip = 1,
    pip_pim = 2,
    pip_pi0 = 3,
    pi0_pi0 = 4,
    Kp_pip  = 11,
    Kp_pim  = 12,
    Kp_pi0  = 13,
    K0_pip  = 14,
    K0_pi0  = 15
  };
    
  class HR_Resonance {
  private:
    ATOOLS::Flavour m_flav;
    double m_BR, m_m2, m_symmfac;
    double m_mass, m_mass2, m_width, m_width2, m_prest, m_spinfac;
    int    m_spin, m_spin1, m_spin2;
  public:
    HR_Resonance(ATOOLS::Flavour flav,
		 ATOOLS::Flavour fl1,ATOOLS::Flavour fl2,
		 const double & symmfac=1., const double & BR=1.,
		 const double & mass=-1., const double & width=-1.);
    const double Width2(const double & s) const;
    const double XStot_R(const double & s) const;
  };
  typedef std::list<HR_Resonance> HR_resonances;

  class HR_Parameters {
  private:
    std::map<std::string,double> m_parameters;
  public:
    HR_Parameters();
    ~HR_Parameters();

    double operator()(const std::string& keyword) const;
  };

  extern const HR_Parameters * hrpars;
}
#endif
