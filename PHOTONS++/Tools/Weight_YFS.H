#ifndef PHOTONS_Tools_Weight_YFS_H
#define PHOTONS_Tools_Weight_YFS_H

#include "PHOTONS++/Tools/Weight_Base.H"
#include "PHOTONS++/Tools/YFS_Form_Factor.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

namespace PHOTONS {
  class Weight_YFS: public Weight_Base {
    private:
      double          m_nbar;
      double          m_Y;
      double          m_Ymax;

      virtual void CalculateWeight();
      virtual void CalculateMax();

    public:
      Weight_YFS(const Particle_Vector, const Particle_Vector, const double, const double);
      virtual ~Weight_YFS();
  };


  

  /*!
    \file Weight_YFS.H
    \brief contains the class Weight_YFS
  */

  /*!
    \class Weight_YFS
    \brief calculates \f$ W_{YFS} \f$ of the the event given
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variable of Weight_YFS
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Weight_YFS::m_nbar
    \brief avarage photon number of the configuration
  */

  /*!
    \var double Weight_YFS::m_Y
    \brief YFS_Form_Factor of the event
  */

  /*!
    \var double Weight_YFS::m_Ymax
    \brief YFS_Form_Factor of an event with \f$ K=0 \f$
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Weight_YFS
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_YFS::Weight_YFS(const Particle_Vector, const Particle_Vector, const double, const double)
    \brief initialises weight variables, calculates both form factors and calculates the weight and the maximum weight

    The values to be passed are the multipole after corrections, the multipole before corrections, the infrared cut-off used in the photon generation and the avarage photon number of the configuration.

    The maximum weight actually contains the weight of a event with \f$ K=0 \f$ which is allways smaller than the actual weight of any event with finite photon momentum, but this contribution is needed for the maximum weight of the overall event.

    Uses the <tt>YFS_Form_Factor</tt> class to calculate the YFS-Form-Factor.
  */

  /*!
    \fn void Weight_YFS::CalculateWeight()
    \brief calculates the weight of the event

    The weight is \f$ W_{YFS} = \exp(Y+\bar{n}) \f$.
  */

  /*!
    \fn void Weight_YFS::CalculateMax()
    \brief calculates the weight of an event with \f$ K=0 \f$

    The "maximum weight" is actually a minimum weight if one looks only at the YFS_Form_Factor contribution. But in terms of the overall weight, the value at \f$ K=0 \f$ is needed.
  */
}

#endif
