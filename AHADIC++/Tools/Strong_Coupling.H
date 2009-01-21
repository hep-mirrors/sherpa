#ifndef Strong_Coupling_H
#define Strong_Coupling_H

#include "Function_Base.H"
#include "Return_Value.H"

namespace AHADIC {

  struct asform {
    enum code {
      constant           = 0,
      fall_off           = 1,
      fixed              = 2,
      fall_off_times_fix = 3,
      GDH_inspired       = 10,
      IR_regulated_GDH   = 11,
      unknown            = 99
    };
  };

  class Strong_Coupling : public ATOOLS::Function_Base {
  private:
    asform::code m_form;
    double       m_beta0, m_Lambda, m_Lambda2, m_kappa2, m_pt02, m_lastpt2, m_asmax;

    double m_gamma, m_a, m_b, m_c, m_d, m_m2;
    const double mg2(const double) const;
    const double n(const double)   const;
  public:
    Strong_Coupling(const asform::code,const double);

    const double operator()(double q2) const;     
    const double Inverse(const double as) const;
    const double MaxValue() const { return m_asmax; }
    const double SelectPT(const double pt2max,const bool pole=true);
    const double GetPTArgument() const { return m_lastpt2; }
  };
}

#endif
