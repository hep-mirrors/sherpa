#include "METOOLS/Main/SpinFuncs.H"

METOOLS::PauliVector::PauliVector(){
  // sigma_0
  (*this)[0] = ATOOLS::CMatrix(2);
  (*this)[0][0][0]=Complex(1);
  (*this)[0][1][1]=Complex(1);

  // sigma_1
  (*this)[1] = ATOOLS::CMatrix(2);
  (*this)[1][1][0]=Complex(1);
  (*this)[1][0][1]=Complex(1);

  // sigma_2
  (*this)[2] = ATOOLS::CMatrix(2);
  (*this)[2][0][1]=Complex(0, -1);
  (*this)[2][1][0]=Complex(0, 1);

  // sigma_3
  (*this)[3] = ATOOLS::CMatrix(2);
  (*this)[3][0][0]=Complex(1);
  (*this)[3][1][1]=Complex(-1);
}

METOOLS::Gamma::Gamma() {
  // Gamma vector in Weyl basis
  METOOLS::PauliVector paulivector = PauliVector();
  for (int i(0); i<4; ++i){
    (*this)[i] = ATOOLS::CMatrix(4);
    double prefactor(-1);
    if (i==0) prefactor=1;
    for (int j(0); j<2; ++j){
      for (size_t k(0); k<2; ++k){
        (*this)[i][j][k+2] = paulivector[i][j][k];
        (*this)[i][j+2][k] = prefactor*paulivector[i][j][k];
      }
    }
  }
}

// Feynman slash
ATOOLS::CMatrix METOOLS::Gamma::operator*(const ATOOLS::Vec4<double> &p) {
  return (*this)[0]*Complex(p[0])+(*this)[1]*Complex(-p[1])+(*this)[2]*Complex(-p[2])+(*this)[3]*Complex(-p[3]);
}
