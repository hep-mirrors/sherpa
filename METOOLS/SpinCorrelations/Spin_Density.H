#ifndef METOOLS_SpinCorrelations_Spin_Density_H
#define METOOLS_SpinCorrelations_Spin_Density_H

#include <vector>
#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Matrix.H"

namespace ATOOLS {
  class Particle;
}

namespace METOOLS {
  class Amplitude2_Tensor;

  class Spin_Density : public Amplitude2_Matrix {
  public:
    Spin_Density(ATOOLS::Particle* p);
    Spin_Density(ATOOLS::Particle* p, const Amplitude2_Tensor* amps);
    Spin_Density(ATOOLS::Particle* p, const Spin_Density* s,
                 const Amplitude2_Tensor* amps);
    ~Spin_Density();
  };

}


#endif
