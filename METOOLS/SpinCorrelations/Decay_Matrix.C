#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace METOOLS;
using namespace ATOOLS;

Decay_Matrix::Decay_Matrix(ATOOLS::Particle* p) :
  Amplitude2_Matrix(p)
{
  // create diagonal normalised matrix
  Complex OneOverN=Complex(1.0/double(m_nhel), 0.0);
  for (size_t i(0); i<m_nhel; ++i) (*this)[(m_nhel+1)*i]=OneOverN;
}

Decay_Matrix::Decay_Matrix(ATOOLS::Particle* p, Amplitude2_Tensor* amps) :
  Amplitude2_Matrix(p)
{
  DEBUG_FUNC(p->RefFlav());
  (*this)=amps->ReduceToMatrix(p);
  Normalise();
  DEBUG_VAR(*this);
}
