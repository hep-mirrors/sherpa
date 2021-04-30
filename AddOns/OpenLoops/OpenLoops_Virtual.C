#include "OpenLoops_Virtual.H"

#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"


using namespace OpenLoops;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

OpenLoops_Virtual::OpenLoops_Virtual(const Process_Info& pi,
                                     const Flavour_Vector& flavs,
                                     int ol_id) :
  Virtual_ME2_Base(pi, flavs), m_ol_id(ol_id), m_ismapped(false),
  m_modebackup(m_mode)
{
  DEBUG_FUNC("");
  msg_Debugging()<<PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi)
                 <<" -> "<<OpenLoops_Interface::s_procmap[m_ol_id]
                 <<" ("<<m_ol_id<<")"<<std::endl;
  m_modebackup=m_mode=OpenLoops_Interface::s_vmode;
}

void OpenLoops_Virtual::SwitchMode(const int mode)
{
  OpenLoops_Interface::SwitchMode(mode);
}

void OpenLoops_Virtual::Calc(const Vec4D_Vector& momenta)
{
  m_mode=m_modebackup;

  OpenLoops_Interface::SetParameter("alpha", AlphaQED());
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());
  OpenLoops_Interface::SetParameter("mu", sqrt(m_mur2));

  bool shouldprinttime(msg_LevelIsDebugging());
  MyTiming* timing;
  if (shouldprinttime) {
    timing = new MyTiming();
    timing->Start();
  }

  OpenLoops_Interface::EvaluateLoop(m_ol_id, momenta, m_born, m_res);

  // factor which by Sherpa convention has to be divided out at this stage
  // if both types of subtraction, still take alphas out
  double coupling(1.);
  if      (m_stype&sbt::qcd) coupling=AlphaQCD();
  else if (m_stype&sbt::qed) coupling=AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");

  if (shouldprinttime) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
               <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    double factor=coupling/2.0/M_PI;
    msg_Out()<<"B     = "<<m_born<<std::endl;
    msg_Out()<<"V_fin = "<<m_res.Finite()<<" -> "<<m_res.Finite()/m_born/factor<<std::endl;
    msg_Out()<<"V_e1  = "<<m_res.IR()<<" -> "<<m_res.IR()/m_born/factor<<std::endl;
    msg_Out()<<"V_e2  = "<<m_res.IR2()<<" -> "<<m_res.IR2()/m_born/factor<<std::endl;
    delete timing;
  }

  // if Born vanishes, do not divide by it, reset mode for this event
  if (!(m_mode&1) && m_born==0.) {
    m_mode|=1;
    msg_Tracking()<<METHOD<<"(): switch to mode 1, Born vanishes"<<std::endl;
  }
  double factor=((m_mode&1)?1.:m_born)*coupling/2.0/M_PI;
  m_res.Finite()/=factor;
  m_res.IR()/=factor;
  m_res.IR2()/=factor;
}

bool OpenLoops_Virtual::IsMappableTo(const PHASIC::Process_Info& pi)
{
  Process_Info looppi(pi);
  if (looppi.m_fi.m_nlotype!=nlo_type::lo) looppi.m_fi.m_nlotype=nlo_type::loop;
  std::string name(PHASIC::Process_Base::GenerateName(looppi.m_ii,looppi.m_fi));
  DEBUG_FUNC(name);
  if (OpenLoops_Interface::s_procmap[m_ol_id]==name) m_ismapped=true;
  else                                               m_ismapped=false;
  msg_Debugging()<<(m_ismapped?"yes":"no")<<std::endl;
  return m_ismapped;
}

DECLARE_VIRTUALME2_GETTER(OpenLoops_Virtual,"OpenLoops_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,OpenLoops_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="OpenLoops") return NULL;
  if (!(pi.m_fi.m_nlotype==nlo_type::loop)) return NULL;

  DEBUG_VAR(pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0]);
  DEBUG_VAR(pi.m_fi.m_nlocpl[0]);
  DEBUG_VAR(pi.m_maxcpl[1]-pi.m_fi.m_nlocpl[1]);
  DEBUG_VAR(pi.m_fi.m_nlocpl[1]);
  OpenLoops_Interface::SetParameter
    ("coupling_qcd_0", (int) pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0]);
  OpenLoops_Interface::SetParameter
    ("coupling_qcd_1", (int) pi.m_fi.m_nlocpl[0]);
  OpenLoops_Interface::SetParameter
    ("coupling_ew_0", (int) pi.m_maxcpl[1]-pi.m_fi.m_nlocpl[1]);
  OpenLoops_Interface::SetParameter
    ("coupling_ew_1", (int) pi.m_fi.m_nlocpl[1]);

  int id = OpenLoops_Interface::RegisterProcess(pi.m_ii, pi.m_fi, 11);
  if (id>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new OpenLoops_Virtual(pi, flavs, id);
  }
  else {
    return NULL;
  }
}

