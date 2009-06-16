#ifndef PHOTONS_Main_Dipole_Type_H
#define PHOTONS_Main_Dipole_Type_H

// #define PHOTONS_DEBUG

#include "ATOOLS/Phys/Particle.H"

using namespace ATOOLS;

namespace PHOTONS {

  typedef std::vector<Particle_Vector>    Particle_Vector_Vector;
  typedef std::vector<Particle_Int_Map>   Particle_Int_Map_Vector;

  class Dipole_Type {
    public:
      enum code {
        ii     =   1,     // initial-initial
        fi     =   2,     // final-initial
        ff     =   3,     // final-final

        unknown =   0
      };
  };


  

  /*!
    \file Dipole_Type.H
    \brief contains the class Dipole_Type and additional typedefs for Particle handling
  */

  /*!
    \class Dipole_Type
    \brief defines a numbering scheme for the different dipole types
  */

  /*!
    \enum Dipole_Type::code
    \brief enumeration of the Dipole_Type

    - ii  -- initial-initial; all charged particles are initial state
    - fi  -- final-initial; charged particles are both initial and final state
    - ff  -- final-final; all charged particles are final state

    - unknown
    .
  */

  /*!
    \typedef std::vector<Particle_Vector>    Particle_Vector_Vector
    \brief builds a vector out of Particle_Vector's
  */

  /*!
    \typedef std::vector<Particle_Int_Map>   Particle_Int_Map_Vector
    \brief builds a vector out of Particle_Int_Map's
  */

}

#endif
