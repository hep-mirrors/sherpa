#include "Effective_Higgs_Coupling.H"
#include "Message.H"
#include "MathTools.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Effective_Higgs_Coupling::Effective_Higgs_Coupling(const double mhiggs) : m_mhiggs(mhiggs)
{}

Complex Effective_Higgs_Coupling::f(double tau)
{
  if (tau<=0.) return Complex(0.,0.);
  if (tau>=1.) return Complex(ATOOLS::sqr(::asin(sqrt(1./tau))),0.);

  double eta=sqrt(1.-tau);
  Complex a(log((1.+eta)/(1.-eta)),-M_PI);
  return -.25*a*a;
}

double Effective_Higgs_Coupling::GetFermionContribution(double mass)
{
  if (mass<=0.) return 2./3.;
  double tau=ATOOLS::sqr(2.*mass/m_mhiggs);
  return tau*(1.+(1.-tau)*real(f(tau)));
}
