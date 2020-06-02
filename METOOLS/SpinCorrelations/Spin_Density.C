#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"

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

Spin_Density::Spin_Density(const Spin_Density& s) : Amplitude2_Matrix(s)
{
}


Spin_Density::~Spin_Density()
{
}

typedef std::vector<std::pair<std::pair<ATOOLS::Flavour,ATOOLS::Vec4D>, Spin_Density*> > SpinDensityMap;
namespace ATOOLS {
  template <> Blob_Data<SpinDensityMap*>::~Blob_Data()
  {
    for (SpinDensityMap::iterator it=m_data->begin(); it!=m_data->end(); it++) {
      delete it->second;
    }
    delete m_data; m_data=NULL;
  }

  template <> Blob_Data_Base* Blob_Data<SpinDensityMap*>::ClonePtr()
  {
    SpinDensityMap* newdata = new SpinDensityMap();
    for (SpinDensityMap::iterator it = m_data->begin(); it!=m_data->end(); ++it) {
      std::pair<ATOOLS::Flavour,ATOOLS::Vec4D> first = std::make_pair(it->first.first, it->first.second);
      newdata->push_back(make_pair(first, new Spin_Density(*it->second)));
    }
    return new Blob_Data(newdata);
  }

  template class Blob_Data<SpinDensityMap*>;
  template SpinDensityMap*&Blob_Data_Base::Get<SpinDensityMap*>();
}
