#ifndef Strong_Coupling_H
#define Strong_Coupling_H

#include "Function_Base.H"
#include "Return_Value.H"

namespace AHADIC {
  class Hadronisation_Parameters;

  struct asform {
    enum code {
      constant = 0,
      fall_off = 1,
      unkonwn  = 99
    };
  };

  class Strong_Coupling : public ATOOLS::Function_Base {
  private:
    asform::code m_form;
    double       m_beta0, m_Lambda2, m_kappa2, m_pt02, m_asmax;
  public:
    Strong_Coupling(const asform::code,const double);

    const double operator()(double q2) const;     
    const double MaxValue() const { return m_asmax; }
  };
}

#endif
