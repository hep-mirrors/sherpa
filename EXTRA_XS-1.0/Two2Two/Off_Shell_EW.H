#ifndef Off_Shell_EW_H
#define Off_Shell_EW_H

#include "Single_XS.H"

namespace EXTRAXS {

  class Off_Shell_qqb_llb: public Single_XS {  
  private:
    
    bool   barred;
    double qe, qf, ae, af, ve, vf, mass;
    double kappa, sin2tw, MZ2, GZ2, alpha;
    double chi1, chi2, term1, term2, colfac;
    int    kswitch;
    double fac,fin;

  public:

    // constructor
    Off_Shell_qqb_llb(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

    double KFactor(double scale);

  };// end of class Off_Shell_qqb_llb
 
  class Off_Shell_q1q2b_lnulb: public Single_XS {  
  private:
    
    double m_ckm2, m_mw2, m_ww2, m_aqed, m_sin2tw, m_kappa, m_colfac;
    bool   m_barred;

  public:

    // constructor
    Off_Shell_q1q2b_lnulb(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

    double KFactor(double scale);

  };// end of class Off_Shell_q1q2b_llb
 
}// end of namespace EXTRAXS

#endif






