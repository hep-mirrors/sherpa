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
    Amplitude2_Matrix(const ATOOLS::Particle* p, Complex factor);
    Amplitude2_Matrix(const Amplitude2_Matrix& amp) : std::vector<Complex>(amp), p_part(amp.p_part), m_nhel(amp.m_nhel) {}
    ~Amplitude2_Matrix();

    Complex Trace() const;
    void Normalise();
    void Add(const Amplitude2_Matrix& sigma, const Complex& factor);
    /**
    * @param coeff vector containing the transformation coefficients, the vector contains one vector for each final
         * helicity which includes all related transformation coefficients to describe the final polarization object in
         * terms of the polarization objects in an old polarization definition (typically the polarization definition
         * used for matrix element calculation), the entries should be ordered according to the polarization ordering
         * in the Amplitude2_Matrix which should be transformed
    * @param conj_coeff vector containing the coefficients for the transformation of the complex conjugate
     * matrix elements, form analog to coeff
    *
    * method for transforming the Amplitude2_Matrix to another basis for polarization definition
    */
    void PolBasisTrafo(vector<vector<Complex>> coeff, vector<vector<Complex>> conj_coeff);
    inline const ATOOLS::Particle* Particle() const { return p_part; }
    inline void SetParticle(const ATOOLS::Particle* p) { p_part=p; }
    inline size_t GetSpin() { return m_nhel; }

    Complex operator*(const Amplitude2_Matrix& sigma) const;
    inline const Complex& operator()(size_t lambda, size_t lambdap) const {
      return (*this)[lambda+m_nhel*lambdap];
    }

    friend std::ostream& operator<<(std::ostream&, const Amplitude2_Matrix&);
  };

}

#endif
