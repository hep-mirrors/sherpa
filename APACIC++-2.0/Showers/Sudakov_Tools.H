#ifndef Sudakov_Tools_H
#define Sudakov_Tools_H

#include "Model_Base.H"
#include "Function_Base.H"
#include "MathTools.H"


namespace APACIC {
  /*!
    This is the class containing all generic helpers needed for the
    Sudakov form factors and the shower evolution. Inside this
    class, there're methods for (running) coupling constants.
  */
  class Sudakov_Tools {
  private :
    ATOOLS::Function_Base * p_as, * p_aqed;
    int    scheme;
    double lambda2, beta0, scalefac;
    double alphaSmax;
    double alphaQEDmax;
  private :
    //! Determines lambda_QCD in LO from alphaS
    void FixLambda2(double);
    //! Set a scale factor to play it safe.
    void Setscalefac(double);
    //! Ouputs some data.
    void Output();
  public :
    /*!
      Standard Constructor, initializes the running alpha_s
    */
    Sudakov_Tools(MODEL::Model_Base *);
    /*!
      Constructor specifying the scheme for alpha_s (fixed or running)
      and the boundaries for the splittings. Gives maximum values of 
      alpha_s at upper and lwoer boundaries.
    */
    Sudakov_Tools(int,MODEL::Model_Base *,double,double);
    ~Sudakov_Tools() {};

    //! Estimate AlphaS from lambda
    double  CrudeAlphaS(double);
    //! Exact alpha_S
    double  AlphaS(double);
    //! Exact alpha_QED
    double  Alpha(double);

    const double GetLambda2() const { return lambda2; }
    const double GetBnorm()   const { return 2.*M_PI*beta0/scalefac; }
    const double GetASmax()   const { return alphaSmax; }
    const double GetAQEDmax() const { return alphaQEDmax; }
  };

}


#endif



