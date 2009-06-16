#ifndef PHOTONS_Tools_Weight_Base_H
#define PHOTONS_Tools_Weight_Base_H

#include "PHOTONS++/Main/Dipole_Type.H"

namespace PHOTONS {
  class Weight_Base {
    protected:
      double            m_weight;
      double            m_maxweight;
      Dipole_Type::code m_dtype;

      virtual void CalculateWeight() = 0;
      virtual void CalculateMax() = 0;

    public:
      Weight_Base();
      virtual ~Weight_Base();

      inline double Get()       { return m_weight; }
      inline double GetMax()    { return m_maxweight; }
  };

  

  /*!
    \file Weight_Base.H
    \brief contains the class Weight_Base
  */

  /*!
    \class Weight_Base
    \brief provides basic variables and methods for all weight classes
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables of Weight_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Weight_Base::m_weight
    \brief weight of the event given
  */

  /*!
    \var double Weight_Base::m_maxweight
    \brief maximum weight of the configuration given
  */

  /*!
    \var Dipole_Type::code Weight_Base::m_dtype
    \brief Dipole_Type of the configuration
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member methods of Weight_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_Base::Weight_Base
    \brief sets all variables to default values
  */

  /*!
    \fn void Weight_Base::CalculateWeight()
    \brief purely virtual in Weight_Base
  */

  /*!
    \fn void Weight_Base::CalculateMax()
    \brief purely virtual in Weight_Base
  */

  /*!
    \fn double Weight_Base::Get()
    \brief returns the weight of the event given
  */

  /*!
    \fn double Weight_Base::GetMax()
    \brief returns the maximum weight of the configuration
  */
}

#endif
