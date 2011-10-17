#ifndef METOOLS_SpinCorrelations_Amplitude2_Matrix_H
#define METOOLS_SpinCorrelations_Amplitude2_Matrix_H

#include <vector>
#include "ATOOLS/Math/MyComplex.H"

namespace ATOOLS {
  class Particle;
}

namespace METOOLS {

  class Amplitude2_Matrix : public std::vector<Complex> {
  protected:
    const ATOOLS::Particle* p_part;
    size_t m_nhel;

  public:
    Amplitude2_Matrix(const ATOOLS::Particle* p);
    ~Amplitude2_Matrix();

    Complex Trace() const;
    void Normalise();
    void Add(const Amplitude2_Matrix& sigma, const Complex& factor);
    inline const ATOOLS::Particle* Particle() const { return p_part; }
    inline void SetParticle(const ATOOLS::Particle* p) { p_part=p; }

    Complex operator*(const Amplitude2_Matrix& sigma) const;
    inline const Complex& operator()(size_t lambda, size_t lambdap) const {
      return (*this)[lambda+m_nhel*lambdap];
    }

    friend std::ostream& operator<<(std::ostream&, const Amplitude2_Matrix&);
  };

}

#endif
