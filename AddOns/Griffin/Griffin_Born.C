#include "AddOns/Griffin/Griffin_Born.H"
#include "AddOns/Griffin/Griffin_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"


using namespace PHASIC;

namespace Griffin {


  Griffin_Born::Griffin_Born(const External_ME_Args& args, const GriffinProcess& proc) :
    Tree_ME2_Base(args), m_process(proc)
  {
  }


  double Griffin_Born::Calc(const Vec4D_Vector& momenta)
  {
    double res = 0.;
    Griffin_Interface::Instance().EvaluateBorn(momenta, m_process, res);
    return res;
  }

}

using namespace Griffin;

DECLARE_TREEME2_GETTER(Griffin::Griffin_Born,
		       "Griffin_Born")

Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,
			      Griffin::Griffin_Born>::
operator()(const External_ME_Args& args) const
{
  if(args.m_source.length() &&
     args.m_source != "Griffin") return nullptr;
  if(args.Flavours().size()!=4) return nullptr;
  return new Griffin_Born(args, Griffin_Interface::MakeProcess(args));
}
