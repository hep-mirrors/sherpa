#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__GOLEM95

#include "HELICITIES/Loops/Triangle_2.H"
#include "HELICITIES/Loops/Golem95_Wrapper.H"

namespace HELICITIES {

Triangle_2::Triangle_2(const std::vector<std::vector<double> >& smat,
                       const std::vector<ATOOLS::Vec4D>& delta_a0,
                       const std::vector<ATOOLS::Vec4D>& delta_a1) :
  m_smat(smat), m_delta_a0(delta_a0), m_delta_a1(delta_a1),
  m_A32(3, std::vector<DivArrC>(3, DivArrC()))
{
}

Triangle_2::~Triangle_2()
{
}

void Triangle_2::CalcFFs()
{
  DEBUG_FUNC("");
  golem95::init(3);
  golem95::set_s_mat(m_smat);
  for (size_t l1=1; l1<4; ++l1) {
    for (size_t l2=1; l2<4; ++l2) {
      golem95::a32(l1, l2, m_A32[l1-1][l2-1]);
      DEBUG_INFO("m_A32["<<l1-1<<"]["<<l2-1<<"]="<<m_A32[l1-1][l2-1].Finite());
    }
  }
  golem95::b32(m_B32);
  golem95::finish();
}

DivArrC Triangle_2::Contract(const ATOOLS::Vec4C& v0, const ATOOLS::Vec4C& v1)
{
  DivArrC res(vector<Complex>(6,Complex(0.0,0.0)));
  for (size_t l0=0; l0<4; ++l0) {
    for (size_t l1=0; l1<4; ++l1) {
      res+=m_A32[l0][l1]*((m_delta_a0[l0]*v0)*(m_delta_a1[l1]*v1));
    }
  }
  res+=m_B32*(v0*v1);
  return res;
}

}

#endif
