#ifndef METOOLS_Current_Library_FF_0_PPP_H
#define METOOLS_Current_Library_FF_0_PPP_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  enum FF_0_PPP_mode {
    piP_pi0_pi0       = 1,
    piM_piP_piP       = 2,
    unknown           = 999
  };
  
  class FF_0_PPP_Base : public FormFactor_Base {
  protected:
    FF_0_PPP_mode       m_mode;
    Complex             m_norm;

    void    FixMode();
    virtual Complex FF_KS(const double & s123,const double & s1,const double & s2);
    virtual Complex FF_RChiPT(const double & s123,const double & s1,const double & s2);
  public:
    FF_0_PPP_Base(const FF_Parameters & params) :
      FormFactor_Base(params),
      m_mode(FF_0_PPP_mode::unknown),
      m_norm(Complex(0.,0.)) {
      m_ffmodel = ff_model::KS;
      FixMode();
    }
    ~FF_0_PPP_Base() {}

    Complex operator()(const ATOOLS::Vec4D_Vector& moms);
    inline const Complex & GetNorm() const { return m_norm; }
  };
}

#endif
