#ifndef running_alphas_h
#define running_alphas_h

#include "Flavour.H"
#include "Function_Base.H"
//#include "MathTools.H"

namespace APHYTOOLS {

  //! contains date for alpha_S running up to 3rd order
  struct AsDataSet {
    double low_scale;
    double high_scale;
    int nf;
    double lambda2;
    // b_i = beta_i / beta_0;   i>=1
    double b[4];
    // beta_0 = 1/12 * ( 33 - 2*nf)
    double beta0;
    // alpha_s at high and low scale 
    // (note there might exist discontinuties at 2nd order)
    double as_low;
    double as_high;
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
    // 0 - no runing, 1 - running
    int mode;

    //! order 0==one loop running, 1==two loop running, 2==three loop running
    int order;

    int     nth;
    double  CF,CA;
    double  as_MZ;
    double  m2_MZ;
    int     mzset;

    double  as_eff;

    //! main data 
    AsDataSet * thresh;

    //! defines thresholds
    void    Init();

    //! coefficients of the Beta functions for QCD
    // cf. T. van Ritbergen, J. Varmaseren, S. Larin PLB 400 (1997) 379
    // up to factor 4^(i+1)
    // @{
    double Beta0(int nf);
    double Beta1(int nf);
    double Beta2(int nf);
    double Beta3(int nf);
    // @}

    /*
    //! number of active quark flavours (QCD)
    int Nf(double t);
    */

    //! calculates lambda^2  and all beta coefficients of a given domain
    double Lambda2(int nr);
    //! calculates the downward--decoupling constant
    double ZetaOS2(double as, double mass2_os, double mu2, int nl);
    //! calculates the upward--decoupling constant
    double InvZetaOS2(double as, double mass2_os, double mu2, int nl);
    //! calculates alpha_S with lambda of the "nr"-dataset
    double AlphaSLam(double Q2, int nr);
    //! initialises continuation over landau pole
    void ContinueAlphaS(int & nr);
  public:
    Running_AlphaS();
    ~Running_AlphaS();

    // running coupling (independent of switches)
    double operator()(double t); 

    //! alpha at MZ
    double AsMZ()       { return as_MZ; }

    //! alpha at the beam energy (calculated during Init() or read from Parameterfile)
    double AsFixed()    { return as_eff; }

    //! returns alpha_s at a given scale Q^2  (i.e. running coupling)
    /*
      behavior should be influenced be rpa-switches
    */
    double AlphaS(double);

    // performs a number of selftests
    void SelfTest();


    int    Nf(double);
    double Beta0(double);
  };

  inline int Running_AlphaS::Nf(double t){

    for (int i=0; i<=nth; ++i) {
      if (t<= thresh[i].high_scale  && t > thresh[i].low_scale )
	return thresh[i].nf;
    }
    return nth;
  };

  inline double Running_AlphaS::Beta0(double t){
    return Beta0(Nf(t)); 
  };

  extern Running_AlphaS * as;     // * alphaS
}



#endif
