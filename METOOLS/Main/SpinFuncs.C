#include "METOOLS/Main/SpinFuncs.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Spinor.H"

template<class Scalar>
METOOLS::PauliVector<Scalar>::PauliVector() : sigma0(ATOOLS::TCMatrix(2, SComplex(0.0))),
sigma1(ATOOLS::TCMatrix(2, SComplex(0.0))), sigma2(ATOOLS::TCMatrix(2, SComplex(0.0))),
sigma3(ATOOLS::TCMatrix(2, SComplex(0.0)))  {
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

  std::vector<unsigned int> gauge_vec = GetGauge();
  // Gamma vector in Weyl basis
  METOOLS::PauliVector<Scalar> paulivector = PauliVector<Scalar>();
  for (int i(0); i<2; ++i){
    for (int j(0); j<2; ++j){
        gamma0[i][j+2] = paulivector[0][i][j];
        gamma0[i+2][j] = paulivector[0][i][j];
        gamma1[i][j+2] = paulivector[gauge_vec[1]][i][j];
        gamma1[i+2][j] = -paulivector[gauge_vec[1]][i][j];
        gamma2[i][j+2] = paulivector[gauge_vec[2]][i][j];
        gamma2[i+2][j] = -paulivector[gauge_vec[2]][i][j];
        gamma3[i][j+2] = paulivector[gauge_vec[3]][i][j];
        gamma3[i+2][j] = -paulivector[gauge_vec[3]][i][j];
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
  return (*this)[0] * SComplex(p[0]) + (*this)[1] * SComplex(-p[1])
         + (*this)[2] * SComplex(-p[2]) + (*this)[3] * SComplex(-p[3]);
}

template <class Scalar>
std::vector<unsigned int> METOOLS::Gamma<Scalar>::GetGauge() {
  // The gauge of the Weyl spinors is implemented such that $p^\m$$\gamma_\mu$ always has the following form (with
  // variables corresponding to methods in Spinor.C) when the momentum entries are shuffled according to s_ri:
  // 0        0       PPMinus  -PTC    `
  // 0        0       -PT      PPPlus   `
  // PPPlus  PTC      0          0      ´
  // PT      PPMinus  0          0    ´
  // One does not receive the same result if instead the $\gamma_\mu$s are shuffled with the same s_ri. Therefore,
  // the s_ri's are changed here to achieve that since the shuffled gamma's are necessary for multiplying with
  // wave functions
  // The gauge of the Weyl spinors in Sherpa can be understood as a rotation of the coordinate system for the
  // gamma matrices.
  std::vector<unsigned int> gauge_vec{0, ATOOLS::Spinor<Scalar>::R1(), ATOOLS::Spinor<Scalar>::R2(), ATOOLS::Spinor<Scalar>::R3()};
  if (gauge_vec[1]==2 && gauge_vec[2]==3 && gauge_vec[3]==1){
    gauge_vec[1] = 3; gauge_vec[2] = 1; gauge_vec[3] = 2;
  }
  else if (gauge_vec[1]==3 && gauge_vec[2]==1 && gauge_vec[3]==2){
    gauge_vec[1] = 2; gauge_vec[2] = 3; gauge_vec[3] = 1;
  }
  else if (gauge_vec[1]==1 && gauge_vec[2]==2 && gauge_vec[3]==3) {}
  else
    THROW(fatal_error,"Gauge choice not implemented");
  return gauge_vec;
}

//=============================
//  Explicit instantiations.
//=============================
template class METOOLS::PauliVector<double>;
template class METOOLS::PauliVector<long double>;
template class METOOLS::Gamma<double>;
template class METOOLS::Gamma<long double>;

