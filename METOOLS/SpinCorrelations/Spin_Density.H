#ifndef METOOLS_SpinCorrelations_Spin_Density_H
#define METOOLS_SpinCorrelations_Spin_Density_H

#include <vector>
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Vector.H"
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
    Spin_Density(const Spin_Density& s);
    ~Spin_Density();

    Spin_Density* CopyAll();
  };

}

typedef std::vector<std::pair<std::pair<ATOOLS::Flavour,ATOOLS::Vec4D>, METOOLS::Spin_Density*> > SpinDensityMap;
namespace ATOOLS {
  template <> Blob_Data<SpinDensityMap*>::~Blob_Data();

  template <> Blob_Data_Base* Blob_Data<SpinDensityMap*>::ClonePtr();
}

#endif
