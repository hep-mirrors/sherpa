#include "PHASIC++/Scales/Scale_Setter_Base.H"

namespace PHASIC {

  class Pilot_Scale_Setter: public Scale_Setter_Base {

    public:

    Pilot_Scale_Setter(const Scale_Setter_Arguments &args);

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode);

  };

}

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Pilot_Scale_Setter,"PILOT",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,Pilot_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Pilot_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    Pilot_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"pilot scale scheme";
}

Pilot_Scale_Setter::Pilot_Scale_Setter
(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args)
{
}

double Pilot_Scale_Setter::Calculate
(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  return m_scale[stp::fac];
}
