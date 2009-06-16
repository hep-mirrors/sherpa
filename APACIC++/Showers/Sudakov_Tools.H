#ifndef APACIC_Showers_Sudakov_Tools_H
#define APACIC_Showers_Sudakov_Tools_H

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"


namespace APACIC {

  class Sudakov_Tools {
  private :

    ATOOLS::Function_Base *p_as, *p_aqed;

    int    m_scheme;
    double m_lambda2, m_beta0;
    double m_alphaSmax, m_alphaQEDmax;

    void FixLambda2(double);

    void Output();

  public :

    Sudakov_Tools(MODEL::Model_Base *);

    ~Sudakov_Tools() {};    

    void CalculateMaxCouplings(const int,const double,const double);

    double AlphaS(double) const;
    double Nf(double) const;
    double Alpha(double) const;

    inline double GetLambda2() const { return m_lambda2;                  }
    inline double GetBnorm() const   { return 2.*M_PI*m_beta0;            }
    inline double GetASmax() const   { return m_alphaSmax;                }
    inline double GetAQEDmax() const { return m_alphaQEDmax;              }
  };

  /*!
    \file
    \brief contains the class APACIC::Sudakov_Tools

  */
  /*!

    \class Sudakov_Tools 
    \brief Helper class for the shower evolution, especially the evaluation of 
           (running) coupling constants.

    This is the class containing all generic helpers needed for the Sudakov form factors 
    and the shower evolution, in particular the running coupling constants.  

    The strong coupling constant is (over-)estimated at one loop through
    \f[
    \frac{1}{\alpha_s(t)} \approx \frac{1}{\kappa}\beta_{0}\log\left(\frac{\mu_{ren} t}{\Lambda^2}\right),,
    \f]
    where
    \f$\kappa\f$ is an extra factor to ensure that the actual \f$\alpha_s\f$ is always smaller than the 
    estimate, \f$\mu_{ren}\f$ is a global scale factor applied to all strong coupling constants throughout 
    the full program, \f$\beta_0\f$ is the corresponding \f$\beta\f$ factor at LO, which depends on the number 
    of flavours and thus the actual scale \f$t\f$, and \f$\Lambda\f$ is the corresponding QCD scale.

    The separation of coupling constants may become important because various ways of how they are taken 
    in the parton shower could be encoded here; also here scale factors modifying the scales for the 
    coupling constants inside the shower only should be taken care of.  At the moment, however, the scale 
    factors are common to both the matrix element and the parton shower and they are transmitted
    by ATOOLS::Run_Parameter::gen::RenormalizationScaleFactor() .

  */
  /*!

    \var ATOOLS::Function_Base *Sudakov_Tools::p_as
    \brief A pointer to the running strong coupling constant \f$\alpha_s\f$ provided in the MODEL 
           phasespace, MODEL::Running_AlphaS .

  */
  /*!

    \var ATOOLS::Function_Base *Sudakov_Tools::p_aqed
    \brief A pointer to the running electromagnetic coupling constant \f$\alpha_{QED}\f$ provided in 
           the MODEL phasespace, MODEL::Running_AlphaQED .

  */
  /*!

    \var int Sudakov_Tools::m_scheme
    \brief This could be used to have different ways of running the shower, for instance with
           fixed vs. running coupling constants.  But it has not yet been fully implemented.
    \todo  Implement different schemes for the evaluation of coupling constants inside the shower.

  */
  /*!

    \var double Sudakov_Tools::m_lambda2
    \brief The QCD scale \f$\Lambda^2\f$ used in the crude \f$\alpha_s\f$ estimate for the shower, see above.  

  */
  /*!

    \var double Sudakov_Tools::m_beta0
    \brief The parameter govering the LO running of the strong coupling constant, \f$\beta_0\f$.  It's used 
    in the crude \f$\alpha_s\f$ estimate for the shower, see above.  

  */
  /*!

    \var double Sudakov_Tools::m_scalefac
    \brief The extra factor \f$\kappa\f$ ensuring that the crude estimate for \f$\alpha_s\f$ is larger
           than the actual value, see above.

  */
  /*!

    \var double Sudakov_Tools::m_alphaSmax
    \brief The maximal value of \f$\alpha_s\f$ in the scale region accessible to the shower.

  */
  /*!

    \var double Sudakov_Tools::m_alphaQEDmax
    \brief The maximal value of \f$\alpha_{QED}\f$ in the scale region accessible to the shower.

  */
  /*!

    \fn Sudakov_Tools::Sudakov_Tools(MODEL::Model_Base *)
    \brief Initialisation of the class, basically filling the pointers to the coupling constants
           and setting the correction factor \f$\kappa=\f$Sudakov_Tools::m_scalefac to one.  

  */
  /*!

    \fn Sudakov_Tools::~Sudakov_Tools()
    \brief Declared but with no function.  

  */
  /*!

    \fn void Sudakov_Tools::FixLambda2(double)    
    \brief With \f$\alpha_s(t)\f$ known, \f$\Lambda^2\f$ follows from the equation above.  The only
           difference here is that the correction factor \f$\kappa\f$ is not taken into account.

  */
  /*!

    \fn void Sudakov_Tools::Setscalefac(double)
    \brief Calculates the correction factor \f$\kappa\f$ from above as the ratio of actual and crude
           \f$\alpha_s\f$ at the lowest possible scale.

  */
  /*!

    \fn void Sudakov_Tools::Output()
    \brief Some output to check what's going on, mainly used for debugging purposes. 

  */
  /*!

    \fn void Sudakov_Tools::CalculateMaxCouplings(const int,const double,const double)
    \brief Determines the maximal values for the running coupling constants.  In general, the maximal value 
           of \f$\alpha_s\f$ is given at the lowest accessible scale, in this case this is the minimum of the 
           minimal scale in the parton shower times the renormalisation scale factor and the lowest scale 
           permitted in the definition of \f$\alpha_s\f$ in MODEL::Running_AlphaS .  For the electromagnetic 
           coupling constant the highest scale is relevant, given by the collider c.m. energy.  

  */
  /*!

    \fn double Sudakov_Tools::CrudeAlphaS(double t) const
    \brief Returns the estimated \f$\alpha_s(t)\f$ from above, i.e. including the
           renormalisation scale factor.

  */
  /*!

    \fn double Sudakov_Tools::Nf(double t) const
    \brief Returns the number of active flavours at the given scale from MODEL::Running_AlphaS::Nf(double) .  

  */
  /*!

    \fn double Sudakov_Tools::AlphaS(double t) const
    \brief Returns the actual \f$\alpha_s(\mu_{ren} t)\f$ from MODEL::Running_AlphaS::operator(double).

  */
  /*!

    \fn double Sudakov_Tools::Alpha(double t) const
    \brief Returns the actual \f$\alpha_{QED}(\mu_{ren} t)\f$ from MODEL::Running_AlphaQED::operator(double).

  */
  /*!

    \fn inline double Sudakov_Tools::GetLambda2() const
    \brief Returns the QCD scale \f$\Lambda^2\f$.

  */
  /*!

    \fn inline double Sudakov_Tools::GetBnorm() const  
    \brief Returns a suitably normalised \f$\beta_0\f$, i.e. \f$2\pi\beta_0/\kappa\f$.

  */
  /*!

    \fn inline double Sudakov_Tools::GetASmax() const  
    \brief Returns the maximal QCD coupling constant over the shower regime, Sudakov_Tools::m_alphaSmax.

  */
  /*!

    \fn inline double Sudakov_Tools::GetAQEDmax() const
    \brief Returns the maximal QED coupling constant over the shower regime, Sudakov_Tools::m_alphaQEDmax.
  */
}

#endif



