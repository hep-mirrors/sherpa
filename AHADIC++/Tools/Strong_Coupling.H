#ifndef AHADIC_Tools_Strong_Coupling_H
#define AHADIC_Tools_Strong_Coupling_H

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/Return_Value.H"

namespace AHADIC {

  struct asform {
    enum code {
      constant           = 0,
      IRregularised      = 1,
      IRregularised_IR0  = 2,
      GDH_inspired       = 8
    };
  };

  class Strong_Coupling : public ATOOLS::Function_Base {
  private:
    asform::code m_form;
    double       m_beta0, m_Lambda, m_Lambda2, m_kappa2, m_pt02, m_pt2max;
    double       m_lastpt2, m_asmax, m_eta;
    double       m_gamma, m_a, m_b, m_c, m_d, m_m2;
    double mg2(const double) const;
    double n(const double)   const;
  public:
    Strong_Coupling(const asform::code);

    double operator()(double q2,bool reweight=false) const;     
    double MaxValue() const { return m_asmax; }
    double SelectPT(const double & pt2max,const double & pt2min=0.);
    double GetPTArgument() const { return m_lastpt2; }
  };
}

#endif
