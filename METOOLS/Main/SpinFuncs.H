#include <vector>
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/Matrix.H"

#ifndef METOOLS_Main_SpinFuncs_H
#define METOOLS_Main_SpinFuncs_H

namespace METOOLS{
  class PauliVector: public ATOOLS::Vec4<ATOOLS::CMatrix> {
  public:
    PauliVector();
  };


  class Gamma: public ATOOLS::Vec4<ATOOLS::CMatrix>{
  public:
    Gamma();
    ATOOLS::CMatrix operator*(const ATOOLS::Vec4<double> &p);
  };
}
#endif