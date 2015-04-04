#include "AMISIC++/Tools/MPI_KFactor_Setter.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace PHASIC;
using namespace ATOOLS;

double MPI_KFactor_Setter::s_pt0=-1.0;

DECLARE_GETTER(MPI_KFactor_Setter,"MPI",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,MPI_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new MPI_KFactor_Setter(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    MPI_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MPI kfactor scheme\n";
}

MPI_KFactor_Setter::MPI_KFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args)
{
  msg_Debugging()<<METHOD<<"(): p_{T,0} = "<<s_pt0<<".\n";
}

double MPI_KFactor_Setter::KFactor() 
{
  double pt2=p_proc->ScaleSetter()->Momenta()[2].PPerp2(), mt2=pt2+sqr(s_pt0);
  return m_weight=sqr(pt2/mt2*(*MODEL::as)(mt2)/(*MODEL::as)(pt2));
}
