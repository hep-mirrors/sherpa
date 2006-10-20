#ifndef Pair_Popper_H
#define Pair_Popper_H

#include "Constituents.H"
#include "Flavour.H"
#include "Return_Value.H"

namespace AHADIC {
  class Hadronisation_Parameters;

  class Pair_Popper {
  private:
    Constituents * p_constituents;
    int            m_ptmode, m_flavmode;
    double         m_tension, m_ptexp;

    double         m_totweight;    
    double         m_maxmass,m_minmass;
  public:
    Pair_Popper();
    ~Pair_Popper() {}
    ATOOLS::Flavour SelectFlavour(const double=100.,const bool=true);
    double          SelectPT(const double=100.);
    double          PopWeight(ATOOLS::Flavour & flav);
  };
  /*!
    \file Pair_Popper.H
    \brief Contains the class AHADIC::Pair_Popper. 
  */

  /*! 
    \class Pair_Popper
    \brief Is used to "pop" quark-antiquark or diquark-antidiquark pairs out of the vacuum,
    subject to kinematical conditions.

    
    The Pair_Popper can run in different modes (Default: m_flavmode=1, m_ptmode=2):
    - m_flavmode describes the way, flavours are selected:  
      - m_flavmode=1 only uses the a-priori weights given by the AHADIC::Constituents.
      - m_flavmode=2 attachs a flavour-mass specific weight to the total a priori weights given 
        in the constituents (accessed through Constituents::TotWeight(ATOOLS::Flavour & flav)).  
	This extra mass-dependent weight is given by
        \f[
        wt = \exp\left[-\frac{\pi m^2}{\sigma}\right]\,,
        \f]
        where \f$\sigma\f$ is one of the parameters of the model ('COLOUR_TENSION').
    - m_ptmode defines, how \f$p_\perp\f$ is selected:
      - m_ptmode=0 corresponds to an exponential distribution
                   \f$p_\perp^2 \propto \exp\left[-\frac{\pi p_\perp^2}{\sigma}\right]\f$
      - m_ptmode=1 corresponds to a simple Gaussian
                   \f$p_\perp \propto \exp\left[-\frac{\pi p_\perp^2}{\sigma}\right]\f$
      - m_ptmode=2 corresponds to a shifted Gaussian
                   \f$p_\perp \propto 
                   \exp\left[-\frac{\pi (p_\perp-\langle p_\perp\rangle)^2}{\sigma}\right]\f$,
        where \f$\langle p_\perp\rangle\f$ is given by the parameter '<PT_SHIFT>'.

    AHADIC::Cluster_Part and AHADIC::Soft_Decay_Handler employ this class.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Pair_Popper
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Constituents * Pair_Popper::p_constituents
    \brief A pointer to the constituents that may pop out of the vacuum.
  */
  /*!
    \var double Pair_Popper::m_totweight
    \brief The static normalisation, calculated in the constructor.
  */
  /*!
    \var int Pair_Popper::m_ptmode
    \brief Defines, how \f$p_\perp\f$ is selected, see in the class description.
    The default value is m_ptmode=2.
  */
  /*!
    \var int Pair_Popper::m_flavmode
    \brief Defines, how the flavour of a popped pair is selected, see in the 
    class description.  The default value is m_flavmode=1.
  */
  /*!
    \var double Pair_Popper::m_tension
    \brief This is \f$\sigma\f$, given by 'COLOUR_TENSION'
  */ 
  /*!
    \var double Pair_Popper::m_ptexp
    \brief This is \f$\langle p_\perp\rangle\f$, relevant if m_ptmode=2.  Given by '<PT_SHIFT>'.
  */ 
  /*!
    \var double Pair_Popper::m_totweight
    \brief The sum of all flavour popping weights if there's no upper mass for the popped pair.
  */
  /*!
    \var double Pair_Popper::m_maxmass
    \brief The mass of the heaviest constituent that can be popped.
  */
  /*!
    \var double Pair_Popper::m_minmass
    \brief The mass of the lightest constituent that can be popped.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Pair_Popper
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Pair_Popper::Pair_Popper()
    \brief Constructor, fills in the p_constituents and calculates the normalisation 
    for popping a pair.
  */
  /*!
    \fn Pair_Popper::~Pair_Popper()
    \brief Destructor, no task so far.
  */
  /*!
    \fn ATOOLS::Flavour Pair_Popper::SelectFlavour(const double=100.,const bool=true)
    \brief Selects a flavour.

    This method selects a flavour that is lighter than the upper limit given by the 
    double.  In case the bool is set to true the flavours are restricted to quarks only.
  */
  /*!
    \fn double Pair_Popper::SelectPT(const double=100.)
    \brief Returns a \f$p_\perp\f$ smaller than an upper limit given as argument.
  */
  /*!
    \fn double Pair_Popper::PopWeight(ATOOLS::Flavour & flav)
    \brief Returns the mass-dependent part of the weight related to a pair-popping, i.e.
    \f$wt = \exp\left[-\frac{\pi m^2}{\sigma}\right]\f$.
  */
}
#endif
