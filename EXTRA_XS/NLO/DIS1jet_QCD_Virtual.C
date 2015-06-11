#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class DIS1jet_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_fac;
  public:
    DIS1jet_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs) {
      m_fac=4.0/3.0;
    }

    ~DIS1jet_QCD_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& mom);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 2.*M_PI*m_mur2/(mom[0]*mom[2]);
    }
  };
}

using namespace EXTRAXS;

void DIS1jet_QCD_Virtual::Calc(const Vec4D_Vector& mom) {
  m_res.IR()=-3.*m_fac;
  m_res.IR2()=-2.*m_fac;
  m_res.Finite()=(-8.)*m_fac;
}


DECLARE_VIRTUALME2_GETTER(DIS1jet_QCD_Virtual,"DIS1jet_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,DIS1jet_QCD_Virtual>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    if (fl[0].IsLepton() && fl[1].IsQuark() && fl[2]==fl[0]  && fl[3]==fl[1]) {
      if (pi.m_maxcpl[0]==1 && pi.m_maxcpl[1]==2 &&
	  pi.m_mincpl[0]==1 && pi.m_mincpl[1]==2) {
        return new DIS1jet_QCD_Virtual(pi, fl);
      }
    }
  }
  return NULL;
}
