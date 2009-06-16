#ifndef Analysis_Observables_Momentum_Tensor_H
#define Analysis_Observables_Momentum_Tensor_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Phys/Particle_List.H"

namespace ANALYSIS {
  class Momentum_Tensor {
  protected:
    double m_power;
    double m_eigenvalues[3];
    ATOOLS::Vec3D m_eigenvectors[3];
  public:
    Momentum_Tensor(double power);
    void Calculate(const ATOOLS::Particle_List & pl);
    void Evaluate(const ATOOLS::Matrix<3> & m);
    double EigenValue(const int i){ return m_eigenvalues[i]; }
    const ATOOLS::Vec3D & EigenVector(const int i) { return m_eigenvectors[i]; }
  };
}


#endif
