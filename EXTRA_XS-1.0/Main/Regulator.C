#include "Regulator.H"

using namespace EXTRAXS;

template <> Regulator_Base*
Regulator_Base::SelectRegulator<Identity>(Single_XS *const xs,const std::string &regulator,
					  const std::vector<double> &parameters)
{
  if (regulator==std::string("Identity")) {
    return new Identity(xs,parameters);
  }
  return NULL;
}

Identity::Identity(Single_XS *const xs,const std::vector<double> &parameters):
  Regulator_Base(xs,parameters) {}

double Identity::operator()(const double dsigma) const
{
  return dsigma;
}

template <> Regulator_Base*
Regulator_Base::SelectRegulator<Massive_Propagator>(Single_XS *const xs,const std::string &regulator,
						    const std::vector<double> &parameters)
{
  if (regulator==std::string("Massive_Propagator") && parameters.size()>0) {
    return new Massive_Propagator(xs,parameters);
  }
  return NULL;
}

Massive_Propagator::Massive_Propagator(Single_XS *const xs,
				       const std::vector<double> &parameters):
  Regulator_Base(xs,parameters) {}

double Massive_Propagator::operator()(const double dsigma) const
{
  double pperp2=p_xs->Momenta()[2].PPerp2();
  return dsigma*ATOOLS::sqr(pperp2)/
    (pperp2+ATOOLS::sqr(m_parameters[0]));
}
