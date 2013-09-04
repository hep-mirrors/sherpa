#ifndef MODEL_Main_Running_Fermion_Mass_H
#define MODEL_Main_Running_Fermion_Mass_H

#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Function_Base.H"

namespace MODEL {
  class Running_Fermion_Mass : public ATOOLS::Function_Base {
  protected:
    int              m_order, m_runbelowpole;
    double           m_polemass, m_a;
    Running_AlphaS * p_as;
    double Series(const double &a,const int nf) const;
    double Beta0(const double &nf) const;
    double Beta1(const double &nf) const;
    double Beta2(const double &nf) const;
    double Gamma0(const double &nf) const;
    double Gamma1(const double &nf) const;
    double Gamma2(const double &nf) const;
  public:
    Running_Fermion_Mass(ATOOLS::Flavour,double,Running_AlphaS *);
    ~Running_Fermion_Mass() {}
    double operator()(double); 
    void SelfTest();
  };

  /*!
    \class Running_Fermion_Mass
    \brief A class for running quark masses 

    This class is for the runnnig of fermion masses, in its current implementation,
    however, only quark masses are dealt with.
    \todo Include leptons as well and repeat the exercise for other types of particles.
  */
  /*!
    \var int Running_Fermion_Mass::m_order
    The order of the running masses, at the moment fixed to one.
    \todo More orders ?
  */
  /*!
    \var double Running_Fermion_Mass::m_gamma
    The anomalous dimension, \f$\gamma_m = 1/\pi \f$. 
  */
  /*!
    \var double Running_Fermion_Mass::m_beta
    The beta function to lowest order at the scale of the quark, 
    \f$\beta_0 = (33-2n_f)/(12\pi)\f$.
  */
  /*!
    \var double Running_Fermion_Mass::m_polemass
    The pole mass of the quark.
  */
  /*!
    \var double Running_Fermion_Mass::m_alphahat
    \f$\alpha_s\f$ at the scale of the polemass.
  */
  /*!
    \var Running_AlphaS * Running_Fermion_Mass::Running_Fermion_Mass::p_as
    A pointer to the running \f$\alpha_s\f$ used in the code.
  */
  /*!
    \fn Running_Fermion_Mass::Running_Fermion_Mass(ATOOLS::Flavour,double,Running_AlphaS *)
    In the constructor of this class, the quark in question, its polemass and a poniter to
    the running \f$\alpha_s\f$ are needed.
  */
  /*!
    \fn double Running_Fermion_Mass::operator()(double) 
    Returns he vaule of the running mass:
    \f[ m(t) = m_{\rm pole}\;\;\; \mbox{\rm if}\;\; |t|<m_{\rm pole}\f]
    \f[ m(t) = m_{\rm pole}\cdot
               \left[\frac{\alpha_s(t)}{\alpha_s(m_{\rm pole})}\right]^{\gamma_m/\beta_0}
               \;\;\; \mbox{\rm if}\;\; |t|\ge m_{\rm pole}\,.\f]
  */
  /*!
    \fn void Running_Fermion_Mass::SelfTest()
    Prints out a range for the running mass. 
  */
}


#endif
