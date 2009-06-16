#include "HELICITIES/Loops/Box_4.H"
#include "HELICITIES/Loops/Golem95_Wrapper.H"
#include "ATOOLS/Org/Exception.H"

#ifdef USING__GOLEM95
using namespace std;
using namespace ATOOLS;

namespace HELICITIES {

Box_4::Box_4(const vector<vector<double> >& smat,
             const vector<ATOOLS::Vec4D>& delta_a1,
             const vector<ATOOLS::Vec4D>& delta_a2,
             const vector<ATOOLS::Vec4D>& delta_a3,
             const vector<ATOOLS::Vec4D>& delta_a4) :
  Divergence_Array<Lorentz_Ten4<Complex> >(Lorentz_Ten4<Complex>(Complex(0.0,0.0))),
  m_smat(smat),
  m_delta_a1(delta_a1), m_delta_a2(delta_a2),
  m_delta_a3(delta_a3), m_delta_a4(delta_a4)
{
}

Box_4::~Box_4()
{
}

void Box_4::Calc()
{
  // eq. (3.20) in JHEP10(2005)015
  DEBUG_FUNC("");
  golem95::init(4);
  golem95::set_s_mat(m_smat);
  for (size_t l1=0; l1<4; ++l1) {
    for (size_t l2=0; l2<4; ++l2) {
      for (size_t l3=0; l3<4; ++l3) {
        for (size_t l4=0; l4<4; ++l4) {
          DivArrC A44;
          golem95::a44(l1+1, l2+1, l3+1, l4+1, A44);
          Lorentz_Ten4<double> ten=BuildTensor(m_delta_a1[l1],
                                               m_delta_a2[l2],
                                               m_delta_a3[l3],
                                               m_delta_a4[l4]);
          for (size_t i=0; i<this->GetResult().size(); ++i)
            (*this)[i]+=ten*A44[i];
        }
      }
      DivArrC B44;
      golem95::b44(l1+1, l2+1, B44);
      Lorentz_Ten4<double> ten=
           BuildTensor(MetricTensor(),1,2, BuildTensor(m_delta_a3[l1],m_delta_a4[l2]),3,4);
      ten+=BuildTensor(MetricTensor(),1,3, BuildTensor(m_delta_a2[l1],m_delta_a4[l2]),2,4);
      ten+=BuildTensor(MetricTensor(),1,4, BuildTensor(m_delta_a2[l1],m_delta_a3[l2]),2,3);
      ten+=BuildTensor(MetricTensor(),2,3, BuildTensor(m_delta_a1[l1],m_delta_a4[l2]),1,4);
      ten+=BuildTensor(MetricTensor(),2,4, BuildTensor(m_delta_a1[l1],m_delta_a3[l2]),1,3);
      ten+=BuildTensor(MetricTensor(),3,4, BuildTensor(m_delta_a1[l1],m_delta_a2[l2]),1,2);
      for (size_t i=0; i<this->GetResult().size(); ++i) (*this)[i]+=ten*B44[i];
    }
  }
  DivArrC C44;
  Lorentz_Ten4<double> ten=
       BuildTensor(MetricTensor(),1,2, MetricTensor(),3,4);
  ten+=BuildTensor(MetricTensor(),1,3, MetricTensor(),2,4);
  ten+=BuildTensor(MetricTensor(),2,3, MetricTensor(),1,4);
  golem95::c44(C44);
  for (size_t i=0; i<this->GetResult().size(); ++i) (*this)[i]+=ten*C44[i];
  golem95::finish();
}


}

#endif
