#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "EXTRA_XS/NLO/Loop_ME_Base.H"
#include "ATOOLS/Org/Exception.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace HELICITIES;

namespace EXTRAXS {
  class PPPP_QED_Virtual : public Virtual_ME2_Base {
    Loop_ME_Base* p_loop_me;
  public:
    PPPP_QED_Virtual(Loop_ME_Base* loop_me,
                     const Process_Info& pi, const Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs), p_loop_me(loop_me)
    {
    }

    ~PPPP_QED_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    bool SetColours(const ATOOLS::Vec4D_Vector& momenta) { return false; }
  };
}


void PPPP_QED_Virtual::Calc(const Vec4D_Vector& momenta) {
  p_loop_me->Calc(momenta);
  
  DivArrC sum(vector<Complex>(5, Complex(0.0, 0.0)));
  for (size_t ihel=0; ihel<p_loop_me->Result().size(); ++ihel) {
    sum+=conj(p_loop_me->Result()[ihel])*p_loop_me->Result()[ihel];
  }
  if (!IsZero(imag(sum))) THROW(fatal_error, "imag(sum tree*loop)!=0");
  m_res=2.0*real(sum);
}


DECLARE_VIRTUALME2_GETTER(PPPP_QED_Virtual_Getter,"PPPP_QED_Virtual")
Virtual_ME2_Base *PPPP_QED_Virtual_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloewtype==nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    for (size_t i=0; i<4; ++i) if (!fl[i].IsPhoton()) return NULL;

    Loop_ME_Base* loop_me=Loop_ME_Base::GetME(pi);
    if (!loop_me) return NULL;
    return new PPPP_QED_Virtual(loop_me, pi, fl);
  }
  return NULL;
}
