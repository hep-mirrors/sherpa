#ifndef Kinematics_h
#define Kinematics_h

#include "Particle.H"
#include "Flavour.H"
#include "Vector.H"
#include <string>


namespace SHERPA {
  class Kinematics {
  public:
    void ShuffleMomenta(std::vector<ATOOLS::Particle *>);
    void ShuffleMomenta(std::vector<ATOOLS::Vec4D>,std::vector<double>);
  };
}

#endif
