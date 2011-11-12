#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class Dummy_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
  public:
    Dummy_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs)
    {
    }

    ~Dummy_QCD_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };
}

using namespace EXTRAXS;

void Dummy_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  // 1/epsIR
  m_res.IR()=0.3;
  // 1/epsIR2
  m_res.IR2()=0.3;
  // finite
  m_res.Finite()=0.3;
}

DECLARE_VIRTUALME2_GETTER(Dummy_QCD_Virtual_Getter,"Dummy_QCD_Virtual")
Virtual_ME2_Base *Dummy_QCD_Virtual_Getter::operator()(const Process_Info &pi) const
{
  Data_Reader read(" ",";","!","=");
  if (read.GetValue<int>("USE_DUMMY_VIRTUAL",0)==0) return NULL;
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    PRINT_INFO("Using Dummy_QCD_Virtual for "<<fl);
    return new Dummy_QCD_Virtual(pi, fl);
  }
  return NULL;
}
