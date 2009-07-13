#ifndef MODEL_Main_Running_Alpha_HV_H
#define MODEL_Main_Running_Alpha_HV_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"

namespace MODEL {

  //! Contains data for alpha_S running up to 2nd order
  struct AsDataSet {
    double low_scale,high_scale;
    double as_low,as_high;
    int nf;
    double lambda2;
    double b[2];
    double beta0;
  };
  //! coefficients of the Beta functions for QCD
  // cf. T. van Ritbergen, J. Vermaseren, S. Larin PLB 400 (1997) 379

  class Running_Alpha_HV : public ATOOLS::Function_Base {
  protected:

    int         m_order;
    int         m_nth,m_mzset;
    double      m_TR,m_CF,m_CA;
    double      m_as_MZ,m_m2_MZ;
    double      m_fac;
    std::string m_group;
    double      m_Nc;
    double      m_cutq2;

    AsDataSet * p_thresh;

    double Beta0(const int);
    double Beta1(const int);
    double Lambda2(const int);
    double ZetaOS2(const double,const double,const double,const int);
    double InvZetaOS2(const double,const double,const double,const int);
    double Alpha_HV_Lam(const double,const int);
    void   ContinueAlpha_HV(int &);

  public:

    // constructor
    Running_Alpha_HV(const double as_MZ,const double m2_MZ,const int order=0, const double fac=1.,
		     const std::string group="SU", const double Nc=3.);
    
    // destructor
    ~Running_Alpha_HV();

    // member functions
    double operator()(double q2); 
    double Alpha_HV(const double q2);
    void   SelfTest();

    // inline functions
    
    inline int    Order()                { return m_order;      }
    inline double Beta0(const double q2) { return Beta0(Nf(q2)); }
    inline double Beta1(const double q2) { return Beta1(Nf(q2)); }
    inline double AsMZ()                 { return m_as_MZ;      }
    inline double AsFixed()              { return m_defval;     }
    int    Nf(const double q2);

    inline double CutQ2() { return m_cutq2; }    
    inline double ScaleFactor() { return m_fac; }
  };

  extern Running_Alpha_HV * as_HV;
}



#endif
