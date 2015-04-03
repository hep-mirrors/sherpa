#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace AMISIC {

  class MPI_KFactor_Setter: public PHASIC::KFactor_Setter_Base {
  private:

    double m_pt0;

  public:

    MPI_KFactor_Setter(const PHASIC::KFactor_Setter_Arguments &args);

    double KFactor();

  };// end of class MPI_KFactor_Setter

}// end of namespace AMISIC

using namespace AMISIC;
using namespace PHASIC;
using namespace ATOOLS;

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
  size_t pos(args.m_kfac.find('{'));
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid p_{T,0} '"+args.m_kfac+"'");
  std::string kftag=args.m_kfac.substr(pos+1);
  pos=kftag.rfind('}');
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid p_{T,0} '"+args.m_kfac+"'");
  m_pt0=ToType<double>(kftag.substr(0,pos));
  msg_Debugging()<<METHOD<<"(): p_{T,0} = "<<m_pt0<<".\n";
}

double MPI_KFactor_Setter::KFactor() 
{
  double pt2=p_proc->ScaleSetter()->Momenta()[2].PPerp2(), mt2=pt2+sqr(m_pt0);
  return m_weight=sqr(pt2/mt2*(*MODEL::as)(mt2)/(*MODEL::as)(pt2));
}
