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

double Identity::operator[](const double scale) const
{
  return scale;
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

double Massive_Propagator::operator[](const double scale) const
{
  const ATOOLS::Vec4D *p=p_xs->Momenta();
  double pperp2=p[2].PPerp2(), oldscale=pperp2, s=0.0, t=0.0, u=0.0;
  switch (p_xs->ScaleScheme()) {
  case 11:
    s=(p[0]+p[1]).Abs2();
    s=(p[0]-p[2]).Abs2();
    s=(p[0]-p[3]).Abs2();
    oldscale=2.*s*t*u/(s*s+t*t+u*u);
    break;
  default:
    oldscale=pperp2;
    break;
  }
  return scale*(pperp2+ATOOLS::sqr(m_parameters[0]))/oldscale;
}
