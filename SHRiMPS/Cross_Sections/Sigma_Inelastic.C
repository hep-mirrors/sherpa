#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_Inelastic::GetValue(const double & B) { 
  return p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B))); 
}

double Sigma_Inelastic::GetCombinedValue(const double & B) { 
  double value(0.);
  for (size_t i=0;i<p_eikonals->size();i++) {
    for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
      Omega_ik * eikonal = (*p_eikonals)[i][j];
      value += eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)));
    }
  }
  return value;
}

double Sigma_Inelastic::GetValuePerChannel(const int k, const int l, const double &B) {
    if (k >= int(p_eikonals->size()) || l >= int(p_eikonals->size())) {
        msg_Error()<<"Error in "<<METHOD<<": requested GW state does not exist, will return zero for cross section.\n";
        return 0.;
    }
    double value(0.);
    if (k >= 0 && l >= 0) {
        Omega_ik * eikonal = (*p_eikonals)[k][l];
        value = eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)));
        value *= sqr(p_eikonals->size());
    }
    if (k >= 0 && l < 0) {
        for (size_t j=0;j<(*p_eikonals)[k].size();j++) {
            Omega_ik * eikonal = (*p_eikonals)[k][j];
            value += eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)));
        }
        value *= p_eikonals->size();
    }
    if (k < 0 && l >= 0) {
        for (size_t i=0;i<(*p_eikonals)[l].size();i++) {
            Omega_ik * eikonal = (*p_eikonals)[i][l];
            value += eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)));
        }
        value *= p_eikonals->size();
    }
    if (k < 0 && l < 0) value = GetCombinedValue(B);
    return value;
}

std::vector<double> * Sigma_Inelastic::FillBGrid(Omega_ik * eikonal) {
  p_eikonal = eikonal;
  std::vector<double> * grid = new std::vector<double>;
  double deltaB(eikonal->DeltaB()), B(0.), sigma(0.), val1(0.), val2(0.);
  grid->push_back(0.);
  do {
    B     += deltaB;
    val2   = 2.*M_PI*B*GetValue(B);
    sigma += deltaB*(val1+val2)/2.;
    grid->push_back(sigma);
    val1   = val2;
  } while (B<MBpars.GetEikonalParameters().bmax);
  return grid;
}  
  
double Sigma_Inelastic::Test() {
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double a(MBpars.GetFFParameters().Lambda2/
	   (8.*(1.+MBpars.GetFFParameters().kappa)));
  double c(MBpars.GetEikonalParameters().beta02*
	   MBpars.GetFFParameters().Lambda2*
	   (1.+MBpars.GetFFParameters().kappa)*
	   exp(2.*MBpars.GetEikonalParameters().Delta*
	       MBpars.GetEikonalParameters().Ymax)/(8.*M_PI));
  double alpha(2.*M_PI*MBpars.GetFFParameters().norm);
  return alpha*(EulerGamma+log(c))/(2.*a)*rpa->Picobarn();
}
