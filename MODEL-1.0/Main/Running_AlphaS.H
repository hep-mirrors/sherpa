#ifndef running_alphas_h
#define running_alphas_h

#include "Flavour.H"
#include "Function_Base.H"

namespace MODEL {

  //! Contains data for alpha_S running up to 3rd order
  struct AsDataSet {
    double low_scale,high_scale;
    double as_low,as_high;
    int nf;
    double lambda2;
    // b_i = beta_i / beta_0;   i>=1
    double b[4];
    // beta_0 = 1/12 * ( 33 - 2*nf)
    double beta0;
    // alpha_s at high and low scale 
    // (note there might exist discontinuties at 2nd order)
  };


  //! QCD SM Running of AlphaS
  /*!
    K.~G.~Chetyrkin, B.~A.~Kniehl and M.~Steinhauser,
   ``Strong coupling constant with flavour thresholds at four loops in the  MS-bar scheme,''
    Phys.\ Rev.\ Lett.\  {\bf 79} (1997) 2184
    [arXiv:hep-ph/9706430].

    Note, the following restriction apply.

    - Working is one, two and three loop running, i.e. four loop is 
      still experimentally. The order is fixed in the constructor where
      order=0 corresponds to lowest order (one loop). Default is order=1, ie
      two loop.
    - The thresholds are fixed to the onshell masses of the correspondig quark.
      At three loop (and any higer order) the running will not be continues any 
      longer.
    - In order to avoid the Landau Pole, the running is stopped as soon as 
      alpha_S becomes one. (Below that scale alpha_s is set linearly to zero.)
  */

  class Running_AlphaS : public AMATOOLS::Function_Base {
  protected:
    int         m_order;
    int         m_nth,m_mzset;
    double      m_CF,m_CA;
    double      m_as_MZ,m_m2_MZ;

    AsDataSet * p_thresh;

    double Beta0(const int);
    double Beta1(const int);
    double Beta2(const int);
    double Beta3(const int);
    double Lambda2(const int);
    double ZetaOS2(const double,const double,const double,const int);
    double InvZetaOS2(const double,const double,const double,const int);
    double AlphaSLam(const double,const int);
    void   ContinueAlphaS(int &);
  public:
    Running_AlphaS(double,double);
    ~Running_AlphaS();
    double operator()(double); 
    double AsMZ()       { return m_as_MZ; }
    double AsFixed()    { return m_defval; }

    double AlphaS(const double);
    int    Nf(const double);
    double Beta0(const double);

    void   SelfTest();
  };

  inline int Running_AlphaS::Nf(const double t) {
    for (int i=0;i<=m_nth;++i) {
      if (t<=p_thresh[i].high_scale && t>p_thresh[i].low_scale )
	return p_thresh[i].nf;
    }
    return m_nth;
  };

  inline double Running_AlphaS::Beta0(const double t){
    return Beta0(Nf(t)); 
  };



  extern Running_AlphaS * as;     // * alphaS


  // 0 - no runing, 1 - running
  //! order 0==one loop running, 1==two loop running, 2==three loop running
  //! main data 
  //! defines thresholds
  //! coefficients of the Beta functions for QCD
  // cf. T. van Ritbergen, J. Varmaseren, S. Larin PLB 400 (1997) 379
  // up to factor 4^(i+1)
  // @{
  // @}
  
  /*
  //! number of active quark flavours (QCD)
  int Nf(double t);
  */
  
  //! calculates lambda^2  and all beta coefficients of a given domain
  
  
  
  //! calculates the downward--decoupling constant
  //! calculates the upward--decoupling constant
  //! calculates alpha_S with lambda of the "nr"-dataset
  //! initialises continuation over landau pole
  // running coupling (independent of switches)
  //! alpha at MZ
  //! alpha at the beam energy (calculated during Init() or read from Parameterfile)
  //! returns alpha_s at a given scale Q^2  (i.e. running coupling)
  /*
    behavior should be influenced be rpa-switches
  */
  // performs a number of selftests
}



#endif
