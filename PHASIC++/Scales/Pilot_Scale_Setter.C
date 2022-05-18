#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace PHASIC {

  class Pilot_Scale_Setter: public Scale_Setter_Base {

    public:

    Pilot_Scale_Setter(const Scale_Setter_Arguments &args);

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode);

    private:

    Scale_Setter_Base* p_var_scale_setter;
    Scale_Setter_Base* p_mets_scale_setter;

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
  // Construct VAR scale setter from string argument passed to the PILOT scale
  // setter.
  Scale_Setter_Arguments varargs = args;
  varargs.m_scale = StringReplace(args.m_scale, "PILOT", "VAR");
  p_var_scale_setter = Scale_Setter_Base::Scale_Getter_Function::GetObject(
      varargs.m_scale, varargs);
  if (p_var_scale_setter == NULL)
    THROW(fatal_error,
          "Could not construct VAR scale setter within PILOT one.");

  // Construct STRICT_METS scale setter.
  Scale_Setter_Arguments metsargs = args;
  metsargs.m_scale = "STRICT_METS";
  p_mets_scale_setter = Scale_Setter_Base::Scale_Getter_Function::GetObject(
      metsargs.m_scale, metsargs);
  if (p_mets_scale_setter == NULL)
    THROW(fatal_error,
          "Could not construct STRICT_METS scale setter within PILOT one.");
}

double Pilot_Scale_Setter::Calculate
(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  return m_scale[stp::fac];
}
