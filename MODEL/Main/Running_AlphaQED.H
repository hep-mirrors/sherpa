#ifndef MODEL_Main_Running_AlphaQED_H
#define MODEL_Main_Running_AlphaQED_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"

namespace MODEL {
  class Running_AlphaQED : public ATOOLS::Function_Base {
    const static double m_A[4],m_B[4],m_C[4];
    double m_alpha0,m_alpha_MZ,m_alpha_eff;
    double m_MZ2,m_split_scale;
  

    double  PiGamma(const ATOOLS::Flavour &, double);
  public:
    Running_AlphaQED(const double,const double,const double=0.);

    double operator()(double); 
    double AqedThomson()  { return m_alpha0; }
    double AqedMZ()       { return m_alpha_MZ; }
    double AqedFixed()    { return m_defval; }
    double Aqed(double);
  };

  extern Running_AlphaQED * aqed;

  /*!
    \class Running_AlphaQED
    \brief The class for the (running) electromagnetic coupling constant.

    This is an implementation of the
    <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?key=2184940">
    electromagnetic coupling constant </A> by R. Kleiss et al. with the 
    <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?key=2076233">
    hadronic component</A> by H. Burkhardt et al..
  */
  /*!
    \var const static double Running_AlphaQED::m_A[4]
    The \f$A_i\f$ parameters needed to calculate the hadronic component of \f$\alpha_{QED}\f$
    \f[A_0 = 0.0,\; A_1 = 0.0,\; A_2 = 0.00165,\; A_3 = 0.00221\,.\f]
  */
  /*!
    \var const static double Running_AlphaQED::m_B[4]
    The \f$A_i\f$ parameters needed to calculate the hadronic component of \f$\alpha_{QED}\f$
    \f[B_0 = 0.00835,\; B_1 = 0.00238,\; B_2 = 0.00299,\; B_3 = 0.00293\,.\f]
  */
  /*!
    \var const static double Running_AlphaQED::m_C[4]
    The \f$A_i\f$ parameters needed to calculate the hadronic component of \f$\alpha_{QED}\f$
    \f[C_0 = 1.0,\; C_1 = 3.927,\; C_2 = 1.0,\; C_3 = 1.0\f]
  */
  /*!
    \var double Running_AlphaQED::m_alpha0
    \f$\alpha_{QED}\f$ in the Thomson limit, 
    \f[\alpha^{(0)}_{QED} = 1/137.03599976\,.\f]
  */
  /*!
    \var double Running_AlphaQED::m_alpha_MZ
    \f$\alpha_{QED}\f$ at the Z-pole.    
  */
  /*!
    \var double Running_AlphaQED::m_alpha_eff
    \todo Not used at the moment.
  */
  /*!
    \var double Running_AlphaQED::m_MZ2
    The square of the Z-mass to calculate \f$\alpha_{QED}\f$ at the Z-pole. 
  */
  /*!
    \var double Running_AlphaQED::m_ecms2
    The square of the c.m. energy as given by the beams.
    \todo Not used at the moment.
  */
  /*!
    \var double Running_AlphaQED::m_split_scale
    \todo Not used at the moment.
  */
  /*!
    \fn double Running_AlphaQED::PiGamma(const ATOOLS::Flavour &, double)
    With the mass of the particle and the scale the vacuum polarization is given by
    \f[\Pi_\gamma(m^2,s)] = \left\{ \begin{array}{lcl}
                    -5/3-\log(m^2/s)\;& \;\mbox{\rm if}\;& \; 4m^2/s<10^{-3}\\
                    1/3-(1+2m^2/s)
                    \left[2+\sqrt{1-4m^2/s}\cdot
		          \log\frac{1-\sqrt{1-4m^2/s}}{1+\sqrt{1-4m^2/s}}\right]
                    \;&\;\mbox{\rm if}\;&\; 4m^2/s<1\\
                    0\;& \;\mbox{\rm if}\;& \; 4m^2/s\ge 1
		    \end{array} \right.\f]
  */
  /*!
    \fn Running_AlphaQED::Running_AlphaQED(const double,const double,const double=0.);
    Initialises \f$\alpha_{QED}\f$, mainly the scales and default values.
  */
  /*!
    \fn double Running_AlphaQED::operator()(double); 
    Returns the value for running \f$\alpha_{QED}\f$.
    \f[\alpha_{QED}(Q^2) = \frac{\alpha^{(0)}_{QED}}{1-\sigma}\,,\f]
    where
    \f[\sigma = \sigma_{\rm lepton} + \sigma_{\rm hadron} + \sigma_{\rm top}\f]
    The leptonic component is given by
    \f[\sigma_{\rm lepton}(Q^2) = \sum\limits_{l=e^-,\mu^-,\tau^-} 
      \frac{\alpha^{(0)}_{QED}}{3\pi}\cdot\Pi_\gamma(l,Q^2)\,.\f]
    \f[\sigma_{\rm hadron}(Q^2) = A[i] + B[i] \log\left(1+C[i]\cdot Q^2\right)\,.\f]
    \f[\sigma_{\rm top}(Q^2) = \frac{\alpha^{(0)}_{QED}}{3\pi}\cdot\Pi_\gamma(t,Q^2)\,.\f]
  */
  /*!
    \fn double Running_AlphaQED::AqedThomson()
    Returns \f$\alpha_{QED}\f$ in the Thomson limit.
  */
  /*!
    \fn double Running_AlphaQED::AqedMZ()
    Returns \f$\alpha_{QED}\f$ at the Z-pole.
  */
  /*!
    \fn double Running_AlphaQED::AqedFixed()
    Returns the default value for \f$\alpha_{QED}\f$.
  */
  /*!
    \fn double Running_AlphaQED::Aqed(double);
    Returns the operator, Running_AlphaQED::operator()(double).
  */
}
#endif
