#ifndef METOOLS_SpinCorrelations_Decay_Matrix_H
#define METOOLS_SpinCorrelations_Decay_Matrix_H

#include "METOOLS/SpinCorrelations/Amplitude2_Matrix.H"

namespace METOOLS {
  class Amplitude2_Tensor;

  class Decay_Matrix : public Amplitude2_Matrix {
    inline Decay_Matrix(const Amplitude2_Matrix& D) : Amplitude2_Matrix(D) {}
  public:
    Decay_Matrix(ATOOLS::Particle* p);
    Decay_Matrix(ATOOLS::Particle* p, Amplitude2_Tensor* amps);
  };

}


#endif
