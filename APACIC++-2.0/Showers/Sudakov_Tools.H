#ifndef Sudakov_Tools_H
#define Sudakov_Tools_H

#include "Model_Base.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"
#include "Run_Parameter.H"


namespace APACIC {

  class Sudakov_Tools {
  private :

    ATOOLS::Function_Base *p_as, *p_aqed;

    int    m_scheme;
    double m_lambda2, m_beta0, m_scalefac;
    double m_alphaSmax, m_alphaQEDmax;

    void FixLambda2(double);
    void Setscalefac(double);

    void Output();

  public :

    Sudakov_Tools(MODEL::Model_Base *);

    ~Sudakov_Tools() {};    

    void CalculateMaxCouplings(const int,const double,const double);

    double CrudeAlphaS(double) const;
    double AlphaS(double) const;
    double Nf(double) const;
    double Alpha(double) const;

    inline double GetLambda2() const { return m_lambda2;                  }
    inline double GetBnorm() const   { return 2.*M_PI*m_beta0/m_scalefac; }
    inline double GetASmax() const   { return m_alphaSmax;                }
    inline double GetAQEDmax() const { return m_alphaQEDmax;              }
  };

  /*!
    \file
    \brief contains the class APACIC::Sudakov_Tools
  */

  /*!
    \class Sudakov_Tools 
    \brief Helper class for the shower evolution, esp. the 
          evaluation of (running) coupling constants.

    This is the class containing all generic helpers needed for the
    Sudakov form factors and the shower evolution. Inside this
    class, there're methods for (running) coupling constants.
  */

}

#endif



