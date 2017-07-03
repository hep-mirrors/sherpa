#include "GoSam_Virtual.H"

#include "AddOns/GoSam/GoSam_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"


using namespace GoSam;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

GoSam_Virtual::GoSam_Virtual(const Process_Info& pi,
                             const Flavour_Vector& flavs,
                             int gs_id) :
  Virtual_ME2_Base(pi, flavs), m_gs_id(gs_id), m_ismapped(false),
  m_modebackup(m_mode)
{
  DEBUG_FUNC("");
  msg_Debugging()<<PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi)
                 <<" -> "<<GoSam_Interface::s_procmap[m_gs_id]
                 <<" ("<<m_gs_id<<")"<<std::endl;
  m_modebackup=m_mode=GoSam_Interface::s_vmode;
}

void GoSam_Virtual::Calc(const Vec4D_Vector& momenta)
{
  DEBUG_FUNC(m_stype);
  m_mode=m_modebackup;

  GoSam_Interface::SetParameter("alpha", AlphaQED());
  GoSam_Interface::SetParameter("alphas", AlphaQCD());

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
  }

  GoSam_Interface::EvaluateLoop(m_gs_id, momenta, m_mur2, m_born, m_res, m_accu);
  if (msg_LevelIsDebugging()) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
               <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    msg_Out()<<"B     = "<<m_born<<std::endl;
    msg_Out()<<"V_fin = "<<m_res.Finite()<<std::endl;
    msg_Out()<<"V_e1  = "<<m_res.IR()<<std::endl;
    msg_Out()<<"V_e2  = "<<m_res.IR2()<<std::endl;
  }

  // factor which by Sherpa convention has to be divided out at this stage
  // if both types of subtraction, still take alphas out
  double coupling(1.);
  if      (m_stype&sbt::qcd) coupling=AlphaQCD();
  else if (m_stype&sbt::qed) coupling=AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");

  // if Born vanishes, do not divide by it, reset mode for this event
  if (!(m_mode&1) && m_born==0.) {
    m_mode|=1;
    msg_Tracking()<<METHOD<<"(): switch to mode 1, Born vanishes."<<std::endl;
  }
  double factor=((m_mode&1)?1.:m_born)*coupling/2.0/M_PI;
  msg_Debugging()<<"cpl="<<coupling/2.0/M_PI<<std::endl;
  m_res.Finite()/=factor;
  m_res.IR()/=factor;
  m_res.IR2()/=factor;
}

bool GoSam_Virtual::IsMappableTo(const PHASIC::Process_Info& pi)
{
  Process_Info looppi(pi);
  if (looppi.m_fi.m_nlotype!=nlo_type::lo) looppi.m_fi.m_nlotype=nlo_type::loop;
  std::string name(PHASIC::Process_Base::GenerateName(looppi.m_ii,looppi.m_fi));
  DEBUG_FUNC(name);
  if (GoSam_Interface::s_procmap[m_gs_id]==name) m_ismapped=true;
  else                                           m_ismapped=false;
  msg_Debugging()<<(m_ismapped?"yes":"no")<<std::endl;
  return m_ismapped;
}

DECLARE_VIRTUALME2_GETTER(GoSam_Virtual,"GoSam_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,GoSam_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="GoSam") return NULL;
  if (!(pi.m_fi.m_nlotype==nlo_type::loop)) return NULL;

  DEBUG_VAR(pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0]);
  DEBUG_VAR(pi.m_fi.m_nlocpl[0]);
  DEBUG_VAR(pi.m_maxcpl[1]-pi.m_fi.m_nlocpl[1]);
  DEBUG_VAR(pi.m_fi.m_nlocpl[1]);
//  GoSam_Interface::SetParameter
//    ("CouplingPower QCD Born", (int) pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0]);
//  GoSam_Interface::SetParameter
//    ("CouplingPower QCD Loop", (int) pi.m_fi.m_nlocpl[0]);
//  GoSam_Interface::SetParameter
//    ("CouplingPower QED Born", (int) pi.m_maxcpl[1]-pi.m_fi.m_nlocpl[1]);
//  GoSam_Interface::SetParameter
//    ("CouplingPower QED Loop", (int) pi.m_fi.m_nlocpl[1]);

  int id = GoSam_Interface::RegisterProcess(pi.m_ii, pi.m_fi);
  if (id>=0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new GoSam_Virtual(pi, flavs, id);
  }
  else {
    return NULL;
  }
}

