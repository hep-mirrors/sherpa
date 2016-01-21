#include "SHRiMPS/Cross_Sections/Sigma_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_Base::Calculate(Omega_ik * eikonal) {
  SetEikonal(eikonal);
  ATOOLS::Gauss_Integrator integrator(this);
  double bmax(MBpars.GetEikonalParameters().bmax);
  double accu(MBpars.GetEikonalParameters().accu);
  return m_sigma = integrator.Integrate(0.,bmax,accu,1)*rpa->Picobarn();
}

double Sigma_Base::operator()(double B) { 
  if (p_eikonal!=NULL) return 2.*M_PI*B*GetValue(B); 
  return 2.*M_PI*B*GetCombinedValue(B); 
}

