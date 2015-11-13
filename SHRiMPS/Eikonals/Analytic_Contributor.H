#ifndef SHRIMPS_Eikonals_Analytic_Contributor_H
#define SHRIMPS_Eikonals_Analytic_Contributor_H

#include "SHRiMPS/Eikonals/Form_Factors.H"

namespace SHRIMPS {
  class Analytic_Contributor : public ATOOLS::Function_Base {
  private:
    Form_Factor * p_ff;
    double        m_Delta, m_Y;
    int           m_sign;
  public:
    Analytic_Contributor(Form_Factor * ff,const double & Delta,const double & Y,
			 const int & sign) :
      p_ff(ff), m_Delta(Delta), m_Y(Y), m_sign(sign) {}
    double operator()(const double & b,const double & y) const;
  };
} 

#endif
