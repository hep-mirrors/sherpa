#ifndef TRIANGLE_2_H
#define TRIANGLE_2_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__GOLEM95
#include "HELICITIES/Loops/Divergence_Array.H"
#include "ATOOLS/Math/Vector.H"

namespace HELICITIES {

class Triangle_2 {
  const std::vector<std::vector<double> >& m_smat;
  const std::vector<ATOOLS::Vec4D>& m_delta_a0;
  const std::vector<ATOOLS::Vec4D>& m_delta_a1;
  
  std::vector<std::vector<DivArrC> > m_A32;
  DivArrC m_B32;
public:
  Triangle_2(const std::vector<std::vector<double> >& smat,
             const std::vector<ATOOLS::Vec4D>& delta_a0,
             const std::vector<ATOOLS::Vec4D>& delta_a1);
  ~Triangle_2();
  
  void CalcFFs();
  DivArrC Contract(const ATOOLS::Vec4C& v0, const ATOOLS::Vec4C& v1);
};

}

#endif
#endif
