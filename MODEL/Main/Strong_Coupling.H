#ifndef MODEL_Main_Strong_Coupling_H
#define MODEL_Main_Strong_Coupling_H

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/Return_Value.H"

namespace MODEL {

  struct asform {
    enum code {
      constant     = 0,
      frozen       = 1,
      smooth       = 2,
      IR0          = 3,
      GDH_inspired = 10
    };
  };

  class Strong_Coupling : public ATOOLS::Function_Base {
  private:
    Running_AlphaS * p_as; 
    asform::code     m_form;
    double           m_pt02;
    double           m_beta0, m_Lambda, m_Lambda2, m_kappa2, m_pt2max;
    //m_Qas2,m_maxalphaS
    double       m_lastpt2, m_asmax, m_eta;
    double       m_gamma, m_a, m_b, m_c, m_d, m_m2;
    double mg2(const double) const;
    double n(const double)   const;
  public:
    Strong_Coupling(Running_AlphaS * as,const asform::code & asf,
		    const double & pt02);
    double operator()(double q2,bool reweight=false) const;     
    double Weight(const double & pt2,const bool & noless=false) const { 
      return (*this)(noless?ATOOLS::Max(m_pt02,pt2):pt2)/m_asmax; 
    }
    const double & MaxValue() const { return m_asmax; }
    const double & PT02()     const { return m_pt02; }
    const asform::code Form() const { return m_form; }
  };
}

#endif
