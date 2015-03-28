#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class MPI_Scale_Setter: public Scale_Setter_Base {
  public:

    MPI_Scale_Setter(const Scale_Setter_Arguments &args);

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(MPI_Scale_Setter,"MPI",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,MPI_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new MPI_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    MPI_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mpi scale scheme";
}

MPI_Scale_Setter::MPI_Scale_Setter(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args)
{
  SetCouplings();
}

double MPI_Scale_Setter::Calculate
(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  double s((momenta[0]+momenta[1]).Abs2());
  double t((momenta[0]-momenta[2]).Abs2());
  double u((momenta[0]-momenta[3]).Abs2());
  m_scale[stp::fac]=m_scale[stp::ren]=-1.0/(1.0/s+1.0/t+1.0/u);
  m_scale[stp::res]=m_scale[stp::fac];
  msg_Debugging()<<METHOD<<"(): Set \\mu_r = "
		 <<sqrt(m_scale[stp::ren])<<", \\mu_f = "
		 <<sqrt(m_scale[stp::fac])<<".\n";
  return m_scale[stp::fac];
}

