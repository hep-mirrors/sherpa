#ifndef HELICITIES_Loops_BOX_4_H
#define HELICITIES_Loops_BOX_4_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__GOLEM95
#include "HELICITIES/Loops/Divergence_Array.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Math/MyComplex.H"

namespace HELICITIES {

class Box_4 : public Divergence_Array<ATOOLS::Lorentz_Ten4<Complex> > {
  const std::vector<std::vector<double> >& m_smat;
  const std::vector<ATOOLS::Vec4D>& m_delta_a1;
  const std::vector<ATOOLS::Vec4D>& m_delta_a2;
  const std::vector<ATOOLS::Vec4D>& m_delta_a3;
  const std::vector<ATOOLS::Vec4D>& m_delta_a4;
  
public:
  Box_4(const std::vector<std::vector<double> >& smat,
        const std::vector<ATOOLS::Vec4D>& delta_a1,
        const std::vector<ATOOLS::Vec4D>& delta_a2,
        const std::vector<ATOOLS::Vec4D>& delta_a3,
        const std::vector<ATOOLS::Vec4D>& delta_a4);
  ~Box_4();
  
  void Calc();
};

}

#endif
#endif
