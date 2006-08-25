#ifndef Pair_Popper_H
#define Pair_Popper_H

#include "Constituents.H"
#include "Flavour.H"
#include "Return_Value.H"

namespace AHADIC {
  class Pair_Popper {
  private:
    Constituents * p_constituents;
    double         m_totweight;
    double         m_tension;
    double         m_maxmass,m_minmass;
  public:
    Pair_Popper();
    ~Pair_Popper() {}
    ATOOLS::Flavour SelectFlavour(const double=100.);
    double          SelectPT(const double=100.);
    double          PopWeight(ATOOLS::Flavour & flav);
  };
  /*!
    \file Pair_Popper.H
    \brief Contains the class AHADIC::Pair_Popper. 
  */

  /*! 
    \class Pair_Popper
    \brief Is used to "pop" quark-antiquark or diquark-antidiquark pairs out of the vacuum.
    The AHADIC::Cluster_Part employs this method.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Pair_Popper
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Constituents * Pair_Popper::p_constituents
    \brief The constituents that may pop out of the vacuum.
  */
  /*!
    \var double Pair_Popper::m_totweight
    \brief The static normalisation, calculated in the constructor.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Pair_Popper
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Pair_Popper::Pair_Popper()
    \brief Constructor, fills in the p_constituents and calcualtes the normalisation for popping a pair.
  */
  /*!
    \fn Pair_Popper::~Pair_Popper()
    \brief Destructor, no task so far.
  */
  /*!
    \fn bool Pair_Popper::Pop(ATOOLS::Flavour &)
    \brief Gives the flavour for a pair of particles popped out of the vacuum.

    The flavour is determined according to the relative probabilities given by the static, flavour-dependent
    weight provided by AHADIC::Constituents. The popped out flavour msut have a mass smaller than the 
    mass of the argument. If no such flavour can be found, "false" is returned. In all other cases,
    "true" is returned and the argument is changed to yield the selected flavour.
  */
}
#endif
