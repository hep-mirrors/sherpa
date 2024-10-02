#include "AddOns/Griffin/Griffin_Born.H"
#include "AddOns/Griffin/Griffin_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"


using namespace PHASIC;

namespace Griffin {


  Griffin_Born::Griffin_Born(const External_ME_Args& args) :
    Tree_ME2_Base(args)
  {
    // m_symfac =Flavour::FSSymmetryFactor(args.m_outflavs);
    // m_symfac*=Flavour::ISSymmetryFactor(args.m_inflavs);
  }


  double Griffin_Born::Calc(const Vec4D_Vector& momenta) 
  {
    double res;
    Griffin_Interface::EvaluateBorn(momenta, res);
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
     args.m_source != "Griffin") return NULL;
  Griffin_Interface::RegisterProcess(args); 
  return new Griffin_Born(args);
}

