#ifndef XS_EW_H
#define XS_EW_H

#include "Single_XS.H"

namespace EXTRAXS {

  class XS_q1q2b_q3q4b: public Single_XS {  
  private:
    
    double m_ckm2[2], m_mw2, m_ww2, m_aqed, m_sin2tw, m_kappa, m_colfac;
    bool   m_barred[2];

  public:
 
    // constructor
    XS_q1q2b_q3q4b(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

    double KFactor(double scale);

  };// end of class XS_q1q2b_llb
 
}// end of namespace EXTRAXS

#endif
