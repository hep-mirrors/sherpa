#include "Regulator_Base.H"

#include "Regulator.H"
#include "Message.H"

using namespace EXTRAXS;

std::ostream &EXTRAXS::operator<<(std::ostream &ostr,const rf::code code)
{
  switch (code) {
  case rf::none:               return ostr<<"None";
  case rf::identity:           return ostr<<"Identity";
  case rf::massive_propagator: return ostr<<"Massive Propagator";
  }
  return ostr;
}

Regulator_Base::Regulator_Base(XS_Base *const xs,const std::vector<double> &parameters,
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

Regulator_Base *Regulator_Base::GetRegulator(XS_Base *const xs,const std::string &regulator,
					     const std::vector<double> &parameters)
{
  Regulator_Base *function=NULL;
  if ((function=SelectRegulator<Massive_Propagator>(xs,regulator,parameters))!=NULL);
  else function=SelectRegulator<Identity>(xs,regulator,parameters);
  return function;
}
