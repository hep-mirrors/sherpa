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

    double  m_fpi_pom, m_fpi_f2, m_fpi_rho, m_fpi_I2;
    double  m_fK_pom, m_fK_f2, m_fK_rho;
    double  m_s0, m_b;
    std::vector<double> m_alpha_pom, m_alpha_f2, m_alpha_rho;
    Complex m_xi_pom, m_xi_f2, m_xi_rho, m_xi_I2;
  public:
    HR_Parameters();
    ~HR_Parameters();




    inline const double s0() const { return m_s0; }
    inline const double alpha_pom(const double & t) const {
      return m_alpha_pom[0] + m_alpha_pom[1]*t;
    }
    inline const double alpha_f2(const double & t) const {
      return m_alpha_f2[0] + m_alpha_f2[1]*t + m_alpha_f2[2]*t*t/2.;
    }
    inline const double alpha_rho(const double & t) const {
      return m_alpha_rho[0] + m_alpha_rho[1]*t + m_alpha_rho[2]*t*t/2.;
    }
    inline const double alpha_I2(const double & t) const {
      return alpha_rho(t)+m_alpha_rho[0]-1.;
    }
    inline const double beta_pi_pom(const double & t) const {
      double alpha_pom_t = alpha_pom(t);
      return ( m_fpi_pom * sqrt(alpha_pom_t*(1.+alpha_pom_t)/2.) *
	       exp(m_b*t/2.) );
    }
    inline const double beta_pi_f2(const double & t) const {
      double alpha_f2_t = alpha_f2(t);
      return ( m_fpi_f2 *
	       sqrt(alpha_f2_t/m_alpha_f2[0]*
		    (1.+alpha_f2_t)/(1.+m_alpha_f2[0])) *
	       exp(m_b*t/2.) );
    }
    inline const double beta_pi_rho(const double & t) const {
      return ( m_fpi_rho *
	       sqrt( ((1.5+1.)*exp(m_b*t)-1.5) *
		     (1.+alpha_rho(t))/(1.+m_alpha_rho[0]) ) );
    }
    inline const double beta_pi_I2(const double & t) const {
      return m_fpi_I2 * exp(m_b*t/2.);
    }
    inline const double beta_K_pom(const double & t) const {
      return m_fK_pom/m_fpi_pom * beta_pi_pom(t);
    }
    inline const double beta_K_f2(const double & t) const {
      return m_fK_f2/m_fpi_f2 * beta_pi_f2(t);
    }
    inline const double beta_K_rho(const double & t) const {
      return m_fK_rho/m_fpi_rho * beta_pi_rho(t);
    }
    inline const Complex & Xi_pom() const { return m_xi_pom; } 
    inline const Complex & Xi_f2()  const { return m_xi_f2; } 
    inline const Complex & Xi_rho() const { return m_xi_rho; } 
    inline const Complex & Xi_I2()  const { return m_xi_I2; } 
    
    double operator()(const std::string& keyword) const;
  };

  extern const HR_Parameters * hrpars;
}
#endif
