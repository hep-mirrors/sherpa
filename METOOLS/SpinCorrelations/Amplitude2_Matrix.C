#include "METOOLS/SpinCorrelations/Amplitude2_Matrix.H"
#include "METOOLS/Main/Polarization_Tools.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"

using namespace METOOLS;
using namespace ATOOLS;

Amplitude2_Matrix::Amplitude2_Matrix(const ATOOLS::Particle* p) :
  p_part(p)
{
  m_nhel=p_part->RefFlav().IntSpin()+1;
  if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

  resize(m_nhel*m_nhel, Complex(0.0,0.0));
}

Amplitude2_Matrix::~Amplitude2_Matrix()
{
}

Complex Amplitude2_Matrix::Trace() const
{
  Complex sum(0.,0.);
  for (size_t i(0); i<m_nhel; i++) 
    sum+=(*this)[(m_nhel+1)*i];
  return sum;
}

void Amplitude2_Matrix::Add(const Amplitude2_Matrix& sigma, const Complex& fac)
{
  if (size()!=sigma.size()) THROW(fatal_error, "Internal error.");
  for (size_t i(0); i<size(); ++i) {
    (*this)[i]+=fac*sigma[i];
  }
}

void Amplitude2_Matrix::Normalise()
{
  Complex factor(Complex(1.0,0.0)/Trace());
  for (size_t i(0); i<size(); ++i) (*this)[i]*=factor;
}

void Amplitude2_Matrix::PolBasisTrafo(vector<vector<Complex>> coeff, vector<vector<Complex>> conj_coeff) {
  METOOLS::Amplitude2_Matrix old_matrix(*this);
  for (size_t a(0); a < m_nhel*m_nhel; ++a) {
    Complex value(0, 0);
    // calculate one new element of amplitude matrix in new spin basis
    for (size_t c(0); c < m_nhel * m_nhel; ++c) {
      value += coeff[a - (a / m_nhel) * m_nhel][c - (c / m_nhel) * m_nhel] * conj_coeff[a / m_nhel][c / m_nhel] *
               old_matrix[c];
    }
    (*this)[a] = value;
  }
}

Complex Amplitude2_Matrix::operator*(const Amplitude2_Matrix& sigma) const
{
  Complex result(0.0,0.0);
  for (size_t i=0; i<m_nhel; ++i) {
    for (size_t j=0; j<m_nhel; ++j) {
      result+=(*this)(i,j)*sigma(i,j);
    }
  }
  return result;
}

namespace METOOLS {
  std::ostream& operator<<(std::ostream& ostr, const Amplitude2_Matrix& m) {
    ostr<<"   Matrix with "<<m.m_nhel<<" spin combinations for "
        <<(m.Particle()?m.Particle()->RefFlav():Flavour(kf_none))<<":"<<std::endl;
    for(size_t i=0;i<m.m_nhel;i++) {
      for(size_t j=0;j<m.m_nhel;j++) {
        ostr<<m(i,j)<<", ";
      }
      ostr<<std::endl;
    }
    return ostr;
  }
}
