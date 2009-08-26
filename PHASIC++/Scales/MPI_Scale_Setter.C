#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class MPI_Scale_Setter: public Scale_Setter_Base {
  public:

    MPI_Scale_Setter(Process_Base *const proc,
		     const std::string &scale);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(MPI_Scale_Setter_Getter,"MPI",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *MPI_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new MPI_Scale_Setter(args.p_proc,args.m_scale);
}

void MPI_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mpi scale scheme\n";
}

MPI_Scale_Setter::MPI_Scale_Setter
(Process_Base *const proc,const std::string &scale): 
  Scale_Setter_Base(proc) {}

double MPI_Scale_Setter::CalculateScale(const std::vector<ATOOLS::Vec4D> &momenta) 
{
  if (!m_kfkey.Assigned()) {
    std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),3,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  double s((momenta[0]+momenta[1]).Abs2());
  double t((momenta[0]-momenta[2]).Abs2());
  double u((momenta[0]-momenta[3]).Abs2());
  m_scale[stp::fac]=m_scale[stp::ren]=-1.0/(1.0/s+1.0/t+1.0/u);
  msg_Debugging()<<METHOD<<"(): Set \\mu_r = "
		 <<sqrt(m_scale[stp::ren])<<", \\mu_f = "
		 <<sqrt(m_scale[stp::fac])<<".\n";
  m_kfkey[0]=m_scale[stp::ren];
  m_kfkey[2]=m_kfkey[1]=m_scale[stp::fac];
  return m_scale[stp::fac];
}

