#include "EXTRA_XS/NLO/Loop_ME_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;

namespace EXTRAXS {
  class DY_QCD_Loop : public Loop_ME_Base {
  public:
    DY_QCD_Loop(const Process_Info& pi, const Flavour_Vector& flavs) :
      Loop_ME_Base(pi, flavs)
    {
    }

    ~DY_QCD_Loop() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
  };
}


void DY_QCD_Loop::Calc(const Vec4D_Vector& momenta) {
  vector<int> spins(4);
  for (spins[0]=0; spins[0]<2; ++spins[0]) {
    for (spins[1]=0; spins[1]<2; ++spins[1]) {
      for (spins[2]=0; spins[2]<2; ++spins[2]) {
        for (spins[3]=0; spins[3]<2; ++spins[3]) {
          DivArrC amp;
          m_res.Insert(amp, spins);
        }
      }
    }
  }
}


//DECLARE_LOOPME_GETTER(DY_QCD_Loop_Getter,"DY_QCD_Loop")
//Loop_ME_Base *DY_QCD_Loop_Getter::operator()(const Process_Info &pi) const
//{
//  if (pi.m_nloqcdtype==nlo_type::loop && pi.m_nloewtype==nlo_type::lo) {
//    Flavour_Vector fl=pi.ExtractFlavours();
//    if (fl.size()!=4) return NULL;
//    if ((fl[2].IsLepton() && fl[3]==fl[2].Bar() &&
//         fl[0].IsQuark()  && fl[1]==fl[0].Bar()) ||   
//        (fl[0].IsLepton() && fl[1]==fl[0].Bar() &&
//         fl[2].IsQuark()  && fl[3]==fl[2].Bar())) {
//      if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
//        return new DY_QCD_Loop(pi, fl);
//      }
//    }
//  }
//  return NULL;
//}
