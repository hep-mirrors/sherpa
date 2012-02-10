#include "SINIC++/Cross_Sections/Sigma_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

using namespace SINIC;
using namespace ATOOLS;

double Sigma_Base::Calculate(const double & Bmin, const double & Bmax,
			     Omega_ik * eikonal) {
  SetEikonal(eikonal);
  ATOOLS::Gauss_Integrator integrator(this);
  return m_sigma = integrator.Integrate(Bmin,Bmax,m_accu,1)*rpa->Picobarn();
}

double Sigma_Base::operator()(double B) { 
  if (p_eikonal!=NULL) return 2.*M_PI*B*GetValue(B); 
  return 2.*M_PI*B*GetCombinedValue(B); 
}

