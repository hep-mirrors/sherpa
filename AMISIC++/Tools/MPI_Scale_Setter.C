#include "AMISIC++/Tools/MPI_Scale_Setter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
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
  str<<"MPI scale scheme: -(1/s + 1/t + 1/u)^{-1}";
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

