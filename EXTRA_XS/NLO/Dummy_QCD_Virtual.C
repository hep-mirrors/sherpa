#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class Dummy_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin;
  public:
    Dummy_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                      const double& ep2, const double& ep, const double& fn) :
      Virtual_ME2_Base(pi, flavs), m_eps2(ep2), m_eps(ep), m_fin(fn)
    {
    }

    ~Dummy_QCD_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };
}

using namespace EXTRAXS;

void Dummy_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  double factor=2*M_PI/AlphaQCD();
  // 1/epsIR
  m_res.IR()=m_eps*factor;
  // 1/epsIR2
  m_res.IR2()=m_eps2*factor;
  // finite
  m_res.Finite()=m_fin*factor;
}

DECLARE_VIRTUALME2_GETTER(Dummy_QCD_Virtual,"Dummy_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,Dummy_QCD_Virtual>::
operator()(const Process_Info &pi) const
{
  Data_Reader read(" ",";","!","=");
  std::vector<double> helpvd;
  if (!read.VectorFromFile(helpvd,"USE_DUMMY_VIRTUAL")) return NULL;
  if (!(helpvd.size()>0 && helpvd[0]==1)) return NULL;
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    double eps2(0.3), eps(0.3), fin(0.3);
    if (helpvd.size()>1) fin=helpvd[1];
    if (helpvd.size()>2) eps=helpvd[2];
    if (helpvd.size()>3) eps2=helpvd[3];
    Flavour_Vector fl(pi.ExtractFlavours());
    msg_Info()<<om::bold<<om::green<<"Caution: Using Dummy_QCD_Virtual for "
              <<fl<<om::reset<<std::endl;
    return new Dummy_QCD_Virtual(pi,fl,eps2,eps,fin);
  }
  return NULL;
}
