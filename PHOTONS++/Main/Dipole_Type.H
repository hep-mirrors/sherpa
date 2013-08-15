#ifndef PHOTONS_Main_Dipole_Type_H
#define PHOTONS_Main_Dipole_Type_H

#include <ostream>
#include <cstdlib>
#include <vector>

// #define PHOTONS_DEBUG

namespace ATOOLS {
  class Particle;
  typedef std::vector<Particle* > Particle_Vector;
}

namespace PHOTONS {

  typedef std::vector<ATOOLS::Particle_Vector> Particle_Vector_Vector;

  struct IdPair {
    public:
      size_t i,j;
      IdPair(size_t _i, size_t _j): i(_i), j(_j) {}
  };

  std::ostream& operator<<(std::ostream&,const IdPair&);

  struct IdPairNbar {
    public:
      IdPair ij;
      double nbar;
      IdPairNbar(IdPair _ij, double _nbar) : 
        ij(_ij), nbar(_nbar) {}
      IdPairNbar(size_t _i, size_t _j, double _nbar) : 
        ij(IdPair(_i,_j)), nbar(_nbar) {}
  };

  std::ostream& operator<<(std::ostream&,const IdPairNbar&);

  typedef std::vector<IdPairNbar> IdPairNbarVector;

  struct Dipole_Type {
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
