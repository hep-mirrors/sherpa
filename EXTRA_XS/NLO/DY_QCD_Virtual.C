#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#define CF 1.33333333333333333
#define CA 3.
#define TR 0.5

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class DY_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_fac;
  public:
    DY_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs)
    {
      m_fac = CF;
    }

    ~DY_QCD_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);

  };
}

using namespace EXTRAXS;

void DY_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  // 1/epsIR
  m_res.IR()=-3.*m_fac;
  // 1/epsIR2
  m_res.IR2()=-2.*m_fac;
  // finite
  m_res.Finite()=(-8.+sqr(M_PI))*m_fac;
}

bool DY_QCD_Virtual::SetColours(const ATOOLS::Vec4D_Vector& momenta) {
  return true;
}

DECLARE_VIRTUALME2_GETTER(DY_QCD_Virtual_Getter,"DY_QCD_Virtual")
Virtual_ME2_Base *DY_QCD_Virtual_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    if ((fl[2].IsLepton() && fl[3]==fl[2].Bar() &&
         fl[0].IsQuark()  && fl[1]==fl[0].Bar()) ||   
        (fl[0].IsLepton() && fl[1]==fl[0].Bar() &&
         fl[2].IsQuark()  && fl[3]==fl[2].Bar())) {
      if ((pi.m_oqcd==1 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        return new DY_QCD_Virtual(pi, fl);
      }
    }
    if ((fl[2].IsLepton() && fl[3].LeptonFamily()==fl[2].LeptonFamily() &&
         fl[0].IsQuark()  && fl[1].QuarkFamily()==fl[0].QuarkFamily())) {
      if ((pi.m_oqcd==1 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        return new DY_QCD_Virtual(pi, fl);
      }
    }
  }
  return NULL;
}
