#include <vector>
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/Matrix.H"

#ifndef METOOLS_Main_SpinFuncs_H
#define METOOLS_Main_SpinFuncs_H

// Implement gamma-matrices corresponding to COMIX_DEFAULT_GAUGE: 0 for the calculation of the Weyl spinors
// For the other gauges different bases for the gamma matrices need to be added here
namespace METOOLS{
  template <class Scalar>
  class PauliVector{
    ATOOLS::TCMatrix<Scalar> sigma0;
    ATOOLS::TCMatrix<Scalar> sigma1;
    ATOOLS::TCMatrix<Scalar> sigma2;
    ATOOLS::TCMatrix<Scalar> sigma3;
  public:
    typedef std::complex<Scalar> SComplex;
    PauliVector();
    ATOOLS::TCMatrix<Scalar> operator[](int i) const;
  };

  template <class Scalar>
  class Gamma{
    ATOOLS::TCMatrix<Scalar> gamma0;
    ATOOLS::TCMatrix<Scalar> gamma1;
    ATOOLS::TCMatrix<Scalar> gamma2;
    ATOOLS::TCMatrix<Scalar> gamma3;
  public:
    typedef std::complex<Scalar> SComplex;
    Gamma();
    ATOOLS::TCMatrix<Scalar> operator[](int i) const;
    ATOOLS::TCMatrix<Scalar> operator*(const ATOOLS::Vec4<double> &p);
  };
}
#endif