#include "SHRiMPS/Cross_Sections/Sigma_Total.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_Tot::GetValue(const double & B) { 
  return 2.*p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.)); 
}

double Sigma_Tot::GetCombinedValue(const double & B) { 
  double value(0.); //, pref, eik;
  for (size_t i=0;i<p_eikonals->size();i++) {
    for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
      Omega_ik * eikonal = (*p_eikonals)[i][j];
      value += 2.*eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)/2.)); 
    }
  }
  return value;
}

double Sigma_Tot::Test() {
  const Eikonal_Parameters & eikparams(MBpars.GetEikonalParameters());
  const FormFactor_Parameters & ffparams(MBpars.GetFFParameters());
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double a(ffparams.Lambda2/(8.*(1.+ffparams.kappa)));
  double c(eikparams.beta02*ffparams.Lambda2*(1.+ffparams.kappa)*
	   exp(2.*eikparams.Delta*eikparams.Ymax)/(8.*M_PI));
  double alpha(4.*M_PI*ffparams.norm);
  ExpInt expint;
  double ei(expint.GetExpInt(-c/2.));
  return alpha*(EulerGamma-ei+log(c/2.))/(2.*a)*rpa->Picobarn();
}




double Elastic_Slope::GetValue(const double & B) { 
  return B*B*p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.))/m_stot; 
}

double Elastic_Slope::GetCombinedValue(const double & B) { 
  double value(0.); //, pref, eik;
  for (size_t i=0;i<p_eikonals->size();i++) {
    for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
      Omega_ik * eikonal = (*p_eikonals)[i][j];
      value += eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)/2.)); 
    }
  }
  return B*B*value/m_stot;
}

