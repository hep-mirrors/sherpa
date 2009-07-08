#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class No_KFactor_Setter: public KFactor_Setter_Base {
  public:

    No_KFactor_Setter(Process_Base *const proc);

    double KFactor();

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(No_KFactor_Setter_Getter,"NO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *No_KFactor_Setter_Getter::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new No_KFactor_Setter(args.p_proc);
}

void No_KFactor_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme\n";
}

No_KFactor_Setter::No_KFactor_Setter
(Process_Base *const proc): KFactor_Setter_Base(proc) {}

double No_KFactor_Setter::KFactor() 
{
  if (!m_kfkey.Assigned()) {
     std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),2,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  m_kfkey<<1.0;
  return 1.0;
}

