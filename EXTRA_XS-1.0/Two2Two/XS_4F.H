#ifndef XS_EW_H
#define XS_EW_H

#include "Single_XS.H"

namespace EXTRAXS {

  class XS_f1f1_f1f1: public Single_XS {  
  private:
    bool   m_Z_on, m_P_on;
    int    m_anti;
    double m_mz2, m_wz2;
    double m_sin2tw, m_cos2tw, m_eq, m_y3f, m_v, m_a;
    double m_aqed, m_pref_qed, m_pref_Z, m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f1_f1f1(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };
 
  class XS_f1f1b_f1f1b: public Single_XS {  
  private:
    bool   m_Z_on, m_P_on;
    int    m_anti1, m_anti2;
    double m_mz2, m_wz2;
    double m_sin2tw, m_cos2tw, m_eq, m_y3f, m_v, m_a;
    double m_aqed, m_pref_qed, m_pref_Z, m_kappa, m_colfac;
    double M_t, M_s, M_mix;
  public:
    // constructor
    XS_f1f1b_f1f1b(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };

  class XS_f1f1b_f2f2b: public Single_XS {  
  private:
    bool   m_Z_on, m_P_on, m_W_on;
    int    m_anti1, m_anti2;
    double m_mz2, m_wz2, m_mw2, m_ww2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    Complex m_ckm;
    double m_aqed, m_pref_qed, m_pref_Z, m_pref_W, m_kappa, m_colfac;
    double M_t, M_s, M_mix;
  public:
    // constructor
    XS_f1f1b_f2f2b(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };

  class XS_f1f2_f1f2: public Single_XS {  
  private:
    bool   m_Z_on, m_P_on, m_W_on;
    int    m_anti, m_rev;
    double m_mz2, m_wz2, m_mw2, m_ww2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    Complex m_ckm;
    double m_aqed, m_pref_qed, m_pref_Z, m_pref_W, m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f2_f1f2(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };

  class XS_f1f2b_f1f2b: public Single_XS {  
  private:
    bool   m_Z_on, m_P_on, m_W_on;
    int    m_anti, m_rev;
    double m_mz2, m_wz2, m_mw2, m_ww2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    Complex m_ckm;
    double m_aqed, m_pref_qed, m_pref_Z, m_pref_W, m_kappa, m_colfac;
    double M_t, M_s, M_mix;
  public:
    // constructor
    XS_f1f2b_f1f2b(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };

  class XS_f1f2_f3f4: public Single_XS {  
  private:
    int    m_anti, m_rev;
    double m_mw2, m_ww2, m_pref_W;
    Complex m_ckm1, m_ckm2;
    double m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f2_f3f4(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };

  class XS_f1f2b_f3f4b: public Single_XS {  
  private:
    int    m_anti, m_rev, m_schannel;
    double m_mw2, m_ww2, m_pref_W;
    Complex m_ckm1, m_ckm2;
    double m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f2b_f3f4b(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);
    bool   SetColours();

    double KFactor(const double scale);
  };


}// end of namespace EXTRAXS

#endif
