#include "OpenLoops_Virtual.H"

#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"


using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace OpenLoops {

OpenLoops_Virtual::OpenLoops_Virtual(const Process_Info& pi,
                                     const Flavour_Vector& flavs,
                                     int ol_id) :
  Virtual_ME2_Base(pi, flavs), m_ol_id(ol_id)
{
}

void OpenLoops_Virtual::Calc(const Vec4D_Vector& momenta) {

  OpenLoops_Interface::SetParameter("alpha", AlphaQED());
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());
  OpenLoops_Interface::SetParameter("mu", sqrt(m_mur2));
  
  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
  }
  OpenLoops_Interface::EvaluateLoop(m_ol_id, momenta, m_born, m_res);
  if (msg_LevelIsDebugging()) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
               <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
  }

  // factor which by Sherpa convention has to be divided out at this stage
  double factor=m_born*AlphaQCD()/2.0/M_PI;
  m_res.Finite()/=factor;
  m_res.IR()/=factor;
  m_res.IR2()/=factor;
}

}


using namespace OpenLoops;

DECLARE_VIRTUALME2_GETTER(OpenLoops_Virtual,"OpenLoops_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,OpenLoops_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="OpenLoops") return NULL;
  if ((pi.m_fi.m_nloewtype==nlo_type::loop) != (pi.m_fi.m_nloqcdtype!=nlo_type::loop)) return NULL;

  OpenLoops_Interface::SetParameter
    ("coupling_qcd_0", (int) pi.m_maxcpl[0]-(pi.m_fi.m_nloqcdtype==nlo_type::loop));
  OpenLoops_Interface::SetParameter
    ("coupling_qcd_1", (int) pi.m_fi.m_nloqcdtype==nlo_type::loop);
  OpenLoops_Interface::SetParameter
    ("coupling_ew_0", (int) pi.m_maxcpl[1]-(pi.m_fi.m_nloewtype==nlo_type::loop));
  OpenLoops_Interface::SetParameter
    ("coupling_ew_1", (int) pi.m_fi.m_nloewtype==nlo_type::loop);

  int id = OpenLoops_Interface::RegisterProcess(pi.m_ii, pi.m_fi, 11);
  if (id>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new OpenLoops_Virtual(pi, flavs, id);
  }
  else {
    return NULL;
  }
}
