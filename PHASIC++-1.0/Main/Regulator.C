#include "Regulator.H"

using namespace PHASIC;

template <> Regulator_Base*
Regulator_Base::SelectRegulator<Identity>(Integrable_Base *const xs,
					  const std::string &regulator,
					  const std::vector<double> &parameters)
{
  if (regulator==std::string("Identity")) {
    return new Identity(xs,parameters);
  }
  return NULL;
}

Identity::Identity(Integrable_Base *const xs,const std::vector<double> &parameters):
  Regulator_Base(xs,parameters,rf::identity) {}

double Identity::operator()(const double dsigma) const
{
  return dsigma;
}

double Identity::operator[](const double scale) const
{
  return scale;
}

template <> Regulator_Base*
Regulator_Base::SelectRegulator<LO_QCD_Regulator>(Integrable_Base *const xs,
						  const std::string &regulator,
						  const std::vector<double> &parameters)
{
  if (regulator==std::string("QCD_Trivial") && parameters.size()>0) {
    return new LO_QCD_Regulator(xs,parameters);
  }
  return NULL;
}

LO_QCD_Regulator::LO_QCD_Regulator(Integrable_Base *const xs,
				   const std::vector<double> &parameters):
  Regulator_Base(xs,parameters,rf::qcd_trivial) {}

double LO_QCD_Regulator::operator()(const double dsigma) const
{
  double pperp2=p_xs->Momenta()[2].PPerp2();
  return dsigma*ATOOLS::sqr(pperp2)/ATOOLS::sqr(pperp2+ATOOLS::sqr(m_parameters[0]));
}

double LO_QCD_Regulator::operator[](const double scale) const
{
  p_xs->SetMomenta();
  const ATOOLS::Vec4D *p=p_xs->Momenta();
  double oldscale=0.0;
  switch (p_xs->ScaleScheme()) {
  case 21: {
    double s=(p[0]+p[1]).Abs2();
    double t=(p[0]-p[2]).Abs2();
    double u=(p[0]-p[3]).Abs2();
    oldscale=2.*s*t*u/(s*s+t*t+u*u);
    break;
  }
  default:
    oldscale=p[2].PPerp2();
    break;
  }
  return scale*(oldscale+ATOOLS::sqr(m_parameters[0]))/oldscale;
}
