#ifndef XS_EW_H
#define XS_EW_H

#include "Single_XS.H"

namespace EXTRAXS {

  class XS_f1f2b_f3f4b: public Single_XS {  
  private:
    double m_ckm2[2], m_mw2, m_ww2, m_aqed, m_sin2tw, m_kappa, m_colfac;
    bool   m_barred[2], m_W_on;
  public:
    // constructor
    XS_f1f2b_f3f4b(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

    double KFactor(double scale);
  };

  class XS_f1f1_f1f1: public Single_XS {  
  private:
    double m_aqed, m_eq, m_mz2, m_wz2, m_sin2tw, m_cos2tw, m_pref_qed, m_pref_Z;
    double m_kappa, m_colfac;
    double M_t, M_u, M_mix;
    bool   m_Z_on, m_P_on;
    int    m_anti, m_y3f;
  public:
    // constructor
    XS_f1f1_f1f1(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(double scale) { return 1.; };
  };
 


}// end of namespace EXTRAXS

#endif
