#include "MODEL/Main/Effective_Higgs_Coupling.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

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

Complex Effective_Higgs_Coupling::GetFermionContribution(double mass,bool pseudoscalar)
{
  if (!pseudoscalar) {
    if (mass<=0.) return Complex(2./3.,0.);
    double tau=ATOOLS::sqr(2.*mass/m_mhiggs);
    return tau*(Complex(1.,0.)+(1.-tau)*f(tau));
  }
  if (mass<=0.) return Complex(-1.,0.);
  double tau=ATOOLS::sqr(2.*mass/m_mhiggs);
  return tau*f(tau);  
}

Complex Effective_Higgs_Coupling::GetVectorContribution(double mass)
{
  if (mass<=0.) return -3.5;
  double tau=ATOOLS::sqr(2.*mass/m_mhiggs);
  return Complex(-1.-1.5*tau,0.)-1.5*tau*(2.-tau)*f(tau);
}

Complex Effective_Higgs_Coupling::GetScalarContribution(double mass)
{
  if (mass<=0.) return 1./6.;
  double tau=ATOOLS::sqr(2.*mass/m_mhiggs);
  return 0.5*tau*(tau*f(tau)-Complex(1.,0.));
}
