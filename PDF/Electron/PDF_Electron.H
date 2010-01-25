#ifndef PDF_Main_PDF_Electron_H
#define PDF_Main_PDF_Electron_H

#include "PDF/Main/PDF_Base.H"

namespace PDF {
  class PDF_Electron : public PDF_Base {
    double m_mass;
    double m_alpha,m_beta;
    double m_xpdf;
    int    m_izetta,m_order,m_init;
  public:
    PDF_Electron(const ATOOLS::Flavour,const int,const int);
    ~PDF_Electron() {}
    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double GetXPDF(const ATOOLS::Flavour );
    
    bool EWOn() { return true; }

  };
  /*!
    \class PDF_Electron
    \brief This is a pure QED structure function for an electron or, more general, for a lepton. 

    This class houses the QED structure function for an electron, or, more general, charged
    leptons in the parametrization of 
    <A HREF=""></A>. 
    It is based on the leading log
    approximation leading to exponentiation including further higher order terms up to
    \f${\cal O}(\alpha^3)\f$ in the electromagnetic coupling constant. Furthermore, different
    exponentiation schemes are available, see also 
    <A HREF=""></A>.
  */
  /*!
    \var double PDF_Electron::m_alpha
    The electromagnetic coupling constant, taken at the scale at which the structure 
    function is to be evaluated.
  */
  /*!
    \var double PDF_Electron::m_beta
    The characteristic exponent of the lepton PDF. It is given by
    \f[
    \beta = \alpha(m_l^2)/\pi (\log(E^2/m_l^2)-1)
    \f]
  */
  /*!
    \var double PDF_Electron::m_mass
    The mass of the lepton.
  */
  /*!
    \var int PDF_Electron::m_order
    The order in alpha for the caluclation of the structure function.
  */
  /*!
    \var int PDF_Electron::m_izetta
    The \f$\zeta\f$-scheme for the definition of how the logarithms enter the exponentiation.
  */
  /*!
    \fn PDF_Electron::PDF_Electron(const ATOOLS::Flavour,const int,const int)
    The constructor, initializes all constant parameters for evaluation in Calculate.
  */
  /*!
    \fn PDF_Base * PDF_Electron::GetCopy()
    A method to initialize another electron PDF as exact copy of the current one.
    This is needed for the initial state shower of APACIC.
  */
  /*!
    \fn void PDF_Electron::Calculate(const double, const double);
    Here, the following expression is evaluated:
    \f[
    \begin{array}{l}
    f(x,Q^2) = \\ \\ \\ \\ \\
    \end{array}
    \begin{array}{l}
    (1-x)^{\frac{\beta}{2}-1}\cdot
    \frac{\beta\exp\left(\frac12\Gamma_E*\beta+\frac38\beta_S\right)}{2\gamma}\\
    - \frac{\beta_H}{4}(1+x)
    - \frac{\beta_H^2}{32}
    \left[\frac{1+3x^2}{1-x}\log(x) + 4(1+x)\log(1-x)+5+x\right]\\
    - \frac{\beta_H^3}{384}
	    \left[\vphantom{\frac32}
	          (1+x)\left(6\mbox{\rm Li}_2(x)+12\log^2(1-x)-3\pi^2\right) \right.\\
	    \;\;\;\;\;\;\;\; 
	          + \frac{1}{1-1x}\left(\frac{3(1+8x+3x^2)}{2}\log(x)
		                        + 6(x+5)(1-x)\log(1-x)
					+12(1+x^2)\log(x)\log(1-x) \right.\\
	     \;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\left.\left.				       
					+12(1+x^2)\log(x)\log(1-x)
					-\frac{1+7x^2}{2}\log^2(x)
					+\frac{39-24x-15x^2}{4}\right)\right]\,,
    \end{array}
    \f]
    where a number of choices are possible for the various \f$\beta\f$.
    Defining 
    \f[
    L = 2\log\frac{Q}{m_e}\;,\;\;
    \beta_e = \frac{2\alpha (L-1)}{\pi}\;,\;\; 
    \eta = \frac{2\alpha L}{\pi} 
    \f]
    for various values of \f$\zeta\f$ (m_izetta) the \f$\beta\f$ are given by:
    \f[
    \zeta = \left\{
            \begin{array}{l} 0 \\ 1 \\ \mbox{\rm else}\end{array}
            \begin{array}{l}  \beta = \beta_e\,,\;\; \beta_H = \beta_S = \eta\,,\\
	                      \beta = \beta_S = \beta_e\,,\;\; \beta_H = \eta\,,\\
	                      \beta = \beta_S = \beta_H = \beta_e\,.
	    \end{array}\right.
    \f]
    The above expression for \f$f(x,Q^2)\f$ is valid for
    \f[ 
    x\in [0,0.9999]\,,
    \f]
    and the pdf weight yields
    \f[
    {\cal W}(x,Q^2) = x f(x,Q^2)\,.
    \f]
    For
    \f[ 
    x\in [0.9999,0.999999]
    \f]
    the pdf is replaced by
    \f[
    {\cal W}(x,Q^2) = x f(x,Q^2)\cdot \frac{100^{\beta/2}}{100^{\beta/2}-1}\,.
    \f]
    For higher $x$ values a zero is returned. Hence, the modification for the high \f$x\f$
    range basically amounts to moving a good portion of the contributions from the potentially 
    numerically instable region close to 1 to a lower range.
  */
} 
#endif // PDF_Electron_H
