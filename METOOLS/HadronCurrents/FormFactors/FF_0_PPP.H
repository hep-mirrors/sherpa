#ifndef METOOLS_Current_Library_FF_0_PPP_H
#define METOOLS_Current_Library_FF_0_PPP_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  enum FF_0_PPP_mode {
    piP_pi0_pi0       = 1,
    piM_piP_piP       = 2,
    K_pi_K            = 10,
    K0_pi_K0b         = 11,
    KS_pi_KS          = 12,
    KS_pi_KL          = 13,
    KL_pi_KL          = 14,
    K_pi0_K0          = 15,
    pi0_pi0_K         = 16,
    K_pi_pi           = 17,
    pi_K0b_pi0        = 18,
    unknown           = 999
  };
  
  class FF_0_PPP_Base : public FormFactor_Base {
  protected:
    FF_0_PPP_mode       m_mode;
    Complex             m_norm;

    void    FixMode();
    virtual Complex FF_KS(const double & s123,const double & s12,
                          const double & s13,const double & s23);
    virtual Complex FF_RChiPT(const double & s123,const double & s12,
                              const double & s13,const double & s23);
  public:
    FF_0_PPP_Base(const FF_Parameters & params) :
      FormFactor_Base(params),
      m_mode(FF_0_PPP_mode::unknown),
      m_norm(Complex(0.,0.)) {
      FixMode();
    }
    ~FF_0_PPP_Base() {}

    Complex operator()(const ATOOLS::Vec4D_Vector& moms);
    inline const Complex & GetNorm() const { return m_norm; }
  };
}

#endif
