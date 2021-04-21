#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class Dummy_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin;
  public:
    Dummy_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep, const double& fn) :
      Virtual_ME2_Base(pi, flavs),
      m_eps2(ep2), m_eps(ep), m_fin(fn)
    {
    }

    ~Dummy_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };
}

using namespace EXTRAXS;

void Dummy_Virtual::Calc(const Vec4D_Vector& momenta) {
  double factor(1.);
  if      (m_stype&sbt::qcd) factor=2*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2*M_PI/AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");
  // 1/epsIR
  m_res.IR()=m_eps*factor;
  // 1/epsIR2
  m_res.IR2()=m_eps2*factor;
  // finite
  m_res.Finite()=m_fin*factor;
}

DECLARE_VIRTUALME2_GETTER(Dummy_Virtual,"Dummy_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,Dummy_Virtual>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator.find("Dummy")!=0) return NULL;
  double def(pi.m_fi.m_nlocpl[0]==1.?0.3:0.03);
  double eps2(def), eps(def), fin(def);
  std::string args(pi.m_loopgenerator);
  size_t pos(args.find('['));
  if (pos!=std::string::npos) {
    args=args.substr(pos+1);
    pos=args.rfind(']');
    if (pos!=std::string::npos) args=args.substr(0,pos);
    else THROW(fatal_error,"No closing ] in Dummy parameters.");
    Data_Reader read(",",";","!","=");
    read.SetString(args);
    std::vector<double> helpvd;
    if (!read.VectorFromString(helpvd,""))
      THROW(fatal_error,"Error in readin.");
    if (helpvd.size()>0) fin=helpvd[0];
    if (helpvd.size()>1) eps=helpvd[1];
    if (helpvd.size()>2) eps2=helpvd[2];
  }
  if (pi.m_fi.m_nlotype==nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    msg_Info()<<om::bold<<om::green<<"Caution: Using dummy virtual ( e2 = "
              <<eps2<<", e = "<<eps<<", f = "<<fin<<" ) for "
              <<fl<<" @ O("<<pi.m_maxcpl[0]<<","<<pi.m_maxcpl[1]<<")"
              <<om::reset<<std::endl;
    return new Dummy_Virtual(pi,fl,eps2,eps,fin);
  }
  return NULL;
}
