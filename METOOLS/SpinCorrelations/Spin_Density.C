#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Particle.H"

using namespace METOOLS;
using namespace ATOOLS;

Spin_Density::Spin_Density(ATOOLS::Particle* p) :
  Amplitude2_Matrix(p)
{
  // create diagonal normalised matrix
  Complex OneOverN=Complex(1.0/double(m_nhel), 0.0);
  for (size_t i(0); i<m_nhel; ++i) (*this)[(m_nhel+1)*i]=OneOverN;
}

Spin_Density::Spin_Density(ATOOLS::Particle* p, const Amplitude2_Tensor* amps) :
  Amplitude2_Matrix(amps->ReduceToMatrix(p))
{
  Normalise();
}

Spin_Density::Spin_Density(ATOOLS::Particle* p, const Spin_Density* sigma0,
                           const Amplitude2_Tensor* amps) :
  Amplitude2_Matrix(p)
{
  if (amps->Next().size()!=sigma0->size()) THROW(fatal_error, "Internal1.");
  for (size_t i(0); i<sigma0->size(); ++i) {
    this->Add(amps->Next()[i]->ReduceToMatrix(p), (*sigma0)[i]);
  }
  Normalise();
}


Spin_Density::~Spin_Density()
{
}
