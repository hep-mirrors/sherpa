#ifndef Sudakov_Tools_H
#define Sudakov_Tools_H

#include "MathTools.H"


namespace APACIC {

  class Sudakov_Tools {
    int m_scheme;
    double m_lambda2, m_beta0, m_scalefac;
    double m_alpha_s_max;
    double m_alpha_qed_max;
  private :
    void FixLambda2(double);
    void Setscalefac(double);
    void Output();
  public :
    Sudakov_Tools();
    Sudakov_Tools(int,double,double);
    ~Sudakov_Tools() {};

    double  CrudeAlphaS(double);
    double  AlphaS(double);
    double  Alpha(double);

    const double GetLambda2() const { return m_lambda2; }
    const double GetBnorm()   const { return 2.*M_PI*m_beta0/m_scalefac; }
    const double GetASmax()   const { return m_alpha_s_max; }
    const double GetAQEDmax() const { return m_alpha_qed_max; }
  };

  /*!
    \file 
    \brief contains the class APACIC::Sudakov_Tools
  */

  /*!
    \class Sudakov_Tools
    \brief Provides all information for an approximated running AlphaS 
    
    This is the class containing all generic helpers needed for the
    Sudakov form factors and the shower evolution. Inside this
    class, there're methods for (running) coupling constants.
  */

  /*!
    \fn void Sudakov_Tools::FixLambda2(double);
    \brief Determines lambda_QCD in LO from alphaS
  */

  /*!
    \fn void Sudakov_Tools::Setscalefac(double);
    \brief Set a scale factor to play it safe.
  */

  /*!
    \fn void Sudakov_Tools::Output();
    \brief Ouputs some data.
  */

  /*!
    \fn Sudakov_Tools::Sudakov_Tools();
    \brief Standard Constructor
  */

  /*!
    \fn Sudakov_Tools::Sudakov_Tools(int,double,double);
    \brief Constructor specifying the scheme for alpha_s (fixed or running)
    and the boundaries for the splittings. Gives maximum values of 
    alpha_s at upper and lower boundaries.
  */

  /*!
    \fn Sudakov_Tools::~Sudakov_Tools();
    \brief Destructor
  */

  /*!
    \fn double  Sudakov_Tools::CrudeAlphaS(double);
    \brief upper estimate AlphaS from lambda
  */

  /*!
    \fn double  Sudakov_Tools::AlphaS(double);
    \brief returns exact alpha_s, normaly only an interface to ATOOLS::Running_AlphaS  
  */

  /*!
    \fn double  Sudakov_Tools::Alpha(double);
    \brief returns exact alpha_QED, normaly only an interface to ATOOLS::Running_AlphaQED
  */

  /*!
    \fn const double Sudakov_Tools::GetLambda2() const { return lambda2; }
    \brief returns lambda_qcd value usd in the approximative alphaS routine
  */

  /*!
    \fn const double Sudakov_Tools::GetBnorm()   const { return 2.*M_PI*beta0/scalefac; }
    \brief accesses betafactor used in the appoximative alphaS routine
  */


  /*
    \fn const double Sudakov_Tools::GetASmax()   const { return alphaSmax; }
    \brief returns maximum possible value of running strong coupling \f$\alpha_S\f$
  */

  /*!
    \fn const double Sudakov_Tools::GetAQEDmax() const { return alphaQEDmax; }
    \brief returns maximum possible value of running electromagnetic
    coupling \f$\alpha_{\rm QED}\f$
  */

}


#endif



