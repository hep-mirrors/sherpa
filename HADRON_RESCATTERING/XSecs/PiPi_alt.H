#ifndef HADRON_RESCATTERING_XSecs_PiPi_H
#define HADRON_RESCATTERING_XSecs_PiPi_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"

//////////////////////////////////////////////////////////////////////////////
//
// All total cross sections in mb.
// For details regarding parametrizations etc., cf the C-file.
//
//////////////////////////////////////////////////////////////////////////////

namespace HADRON_RESCATTERING {  
  class PiPi {
  public:
    static double s_B_S[6], s_z02_S, s_K_S[4], s_d_S[3], s_eps_S[3];
  private:
    double  m_mpi, m_mpi2, m_mK, m_mK2, m_alpha, m_s0;
    Complex m_sp;
    
    
    inline const double sigma(const double & s) const {
      return sqrt(1.-4.*m_mpi2/s);
    }
    inline const double q(const double & s) const {
      return sqrt(s-4.*m_mpi2)/2.;
    }
    inline const Complex omega(const double & s) const {
      Complex sqrt_s   = Complex(sqrt(s),0.);
      Complex sqrt_s0s = ATOOLS::csqrt(m_s0-s);
      return (sqrt_s-m_alpha*sqrt_s0s)/(sqrt_s+m_alpha*sqrt_s0s);
    }
    inline const Complex omega1(const Complex & s) const {
      return 2.*(sqrt(s)-2.*m_mK)/(1.5-2.*m_mK)-1.;
    }
    const Complex Phi00(const double & s) const;
    const Complex f(const Complex & s) const;
    const Complex JBar(const Complex & s,const double & mi2) const;

    const Complex t00(const double & s) const;
    const Complex t00_conf(const double & s) const;
    const Complex t0_f0(const double & s) const;

    void Test();
  public:
    PiPi();
    ~PiPi() {}
    double XStot(const double & s);
    double XSel(const double & s);
  };
}

#endif
