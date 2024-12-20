#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class PionPionVirtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin;
  public:
    PionPionVirtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep) :
      Virtual_ME2_Base(pi, flavs),
      m_eps2(ep2), m_eps(ep)
    {
    }

    ~PionPionVirtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };
}

using namespace EXTRAXS;

void PionPionVirtual::Calc(const Vec4D_Vector& momenta) {
  double factor(1.);
  if      (m_stype&sbt::qcd) factor=2*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2*M_PI/AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");
  // 1/epsIR
  // m_res.IR()=m_eps*factor;
  // // 1/epsIR2
  // m_res.IR2()=m_eps2*factor;
  // // finite
  m_res.Finite()=0;
}

DECLARE_VIRTUALME2_GETTER(EXTRAXS::PionPionVirtual,"PionPionVirtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::PionPionVirtual>::
operator()(const Process_Info &pi) const
{
  PRINT_INFO("HERE");
  if (pi.m_loopgenerator.find("Internal")!=0) return NULL;
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  Flavour_Vector fl(pi.ExtractFlavours());
  if(fl.size()!=4) return NULL;
  // if(fl[2])
  return new PionPionVirtual(pi,fl,0,0);
  // }
  return NULL;
}
