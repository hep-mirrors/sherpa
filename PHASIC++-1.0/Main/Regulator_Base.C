#include "Regulator_Base.H"

#include "Regulator.H"
#include "Message.H"

using namespace PHASIC;

std::ostream &PHASIC::operator<<(std::ostream &ostr,const rf::code code)
{
  switch (code) {
  case rf::none:        return ostr<<"None";
  case rf::identity:    return ostr<<"Identity";
  case rf::qcd_trivial: return ostr<<"QCD Trivial";
  }
  return ostr;
}

Regulator_Base::Regulator_Base(Integrable_Base *const xs,
			       const std::vector<double> &parameters,
			       const rf::code type):
  m_type(type),
  m_parameters(parameters),
  p_xs(xs) {}

Regulator_Base::~Regulator_Base() {}

double Regulator_Base::operator()(const double dsigma) const
{
  ATOOLS::msg.Error()<<"Regulator_Base::operator()(..): "
		     <<"Virtual method called."<<std::endl;
  return 0.0;
}

double Regulator_Base::operator[](const double scale) const
{
  ATOOLS::msg.Error()<<"Regulator_Base::operator[](..): "
		     <<"Virtual method called."<<std::endl;
  return scale;
}

Regulator_Base *Regulator_Base::GetRegulator(Integrable_Base *const xs,
					     const std::string &regulator,
					     const std::vector<double> &parameters)
{
  Regulator_Base *function=NULL;
  if ((function=SelectRegulator<LO_QCD_Regulator>(xs,regulator,parameters))!=NULL);
  else function=SelectRegulator<Identity>(xs,regulator,parameters);
  return function;
}
