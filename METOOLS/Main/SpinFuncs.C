#include "METOOLS/Main/SpinFuncs.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

template<class Scalar>
METOOLS::PauliVector<Scalar>::PauliVector() : sigma0(ATOOLS::TCMatrix(2, SComplex(0.0))),
sigma1(ATOOLS::TCMatrix(2, SComplex(0.0))), sigma2(ATOOLS::TCMatrix(2, SComplex(0.0))),
sigma3(ATOOLS::TCMatrix(2, SComplex(0.0)))  {
  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  int gauge = s["COMIX_DEFAULT_GAUGE"].Get<int>();
  if ( gauge != 0)
    THROW(not_implemented, "Basis for gamma matrices for the chosen COMIX_DEFAULT_GAUGE is not implemented")
  sigma0[0][0] = sigma0[1][1] = sigma1[0][1] = sigma1[1][0] = sigma3[0][0] = SComplex(1);
  sigma3[1][1]=SComplex(-1);
  sigma2[0][1]=SComplex(0, -1);
  sigma2[1][0]=SComplex(0, 1);
}

template<class Scalar>
ATOOLS::TCMatrix<Scalar> METOOLS::PauliVector<Scalar>::operator[](int i) const {
  if (i>3) THROW(fatal_error, "There are only four Pauli matrices!")
  if (i==0) return sigma0;
  if (i==1) return sigma1;
  if (i==2) return sigma2;
  if (i==3) return sigma3;
}

template<class Scalar>
METOOLS::Gamma<Scalar>::Gamma() : gamma0(ATOOLS::TCMatrix(4, SComplex(0.0))),
gamma1(ATOOLS::TCMatrix(4, SComplex(0.0))), gamma2(ATOOLS::TCMatrix(4, SComplex(0.0))),
gamma3(ATOOLS::TCMatrix(4, SComplex(0.0))) {

  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  int gauge = s["COMIX_DEFAULT_GAUGE"].Get<int>();
  if ( gauge != 0)
    THROW(not_implemented, "Basis for gamma matrices for the chosen COMIX_DEFAULT_GAUGE is not implemented")

  // Gamma vector in Weyl basis
  METOOLS::PauliVector<Scalar> paulivector = PauliVector<Scalar>();
  for (int i(0); i<2; ++i){
    for (int j(0); j<2; ++j){
        gamma0[i][j+2] = paulivector[0][i][j];
        gamma0[i+2][j] = paulivector[0][i][j];
        gamma1[i][j+2] = paulivector[1][i][j];
        gamma1[i+2][j] = -paulivector[1][i][j];
        gamma2[i][j+2] = paulivector[2][i][j];
        gamma2[i+2][j] = -paulivector[2][i][j];
        gamma3[i][j+2] = paulivector[3][i][j];
        gamma3[i+2][j] = -paulivector[3][i][j];
    }
  }
}

template<class Scalar>
ATOOLS::TCMatrix<Scalar> METOOLS::Gamma<Scalar>::operator[](int i) const{
  if (i>3) THROW(fatal_error, "There are only four Pauli matrices!")
  if (i==0) return gamma0;
  if (i==1) return gamma1;
  if (i==2) return gamma2;
  if (i==3) return gamma3;
}

// Feynman slash
template <class Scalar>
ATOOLS::TCMatrix<Scalar> METOOLS::Gamma<Scalar>::operator*(const ATOOLS::Vec4<double> &p) {
  return (*this)[0]*SComplex(p[0])+(*this)[1]*SComplex(-p[1])+(*this)[2]*SComplex(-p[2])+(*this)[3]*SComplex(-p[3]);
}

//=============================
//  Explicit instantiations.
//=============================
template class METOOLS::PauliVector<double>;
template class METOOLS::PauliVector<long double>;
template class METOOLS::Gamma<double>;
template class METOOLS::Gamma<long double>;

