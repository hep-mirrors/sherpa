#include "Griffin_Virtual.H"

#include "AddOns/Griffin/Griffin_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;
using namespace Griffin;


Griffin::Griffin_Virtual::Griffin_Virtual(const Process_Info& pi,
        const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi, flavs)
{
  m_mode=1;
}
  
void Griffin::Griffin_Virtual::Calc(const Vec4D_Vector& momenta) {
 Griffin_Interface::EvaluateLoop(momenta,m_res);
 double s = (momenta[2]+momenta[3]).Abs2();
 m_res.Finite() *= 2*s*s;
}



DECLARE_VIRTUALME2_GETTER(Griffin::Griffin_Virtual,"Griffin_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,Griffin::Griffin_Virtual>::
operator()(const PHASIC::Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Griffin") return NULL;

  if (pi.m_fi.m_nlotype!=nlo_type::loop) return NULL;
  Griffin_Interface::RegisterProcess(pi);
  Flavour_Vector flavs = pi.ExtractFlavours();
  if(flavs.size()==4) return new Griffin_Virtual(pi, flavs);
  else return NULL;
}
