#include "SHRiMPS/Eikonals/Analytic_Eikonal.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Analytic_Eikonal::Analytic_Eikonal() :
  m_norm(1./(2.*(1.+MBpars.GetFFParameters().kappa))),
  m_prefactor(MBpars.GetEikonalParameters().beta02*
	      MBpars.GetFFParameters().Lambda2/(4.*M_PI)*
	      exp(2.*MBpars.GetEikonalParameters().Delta*
		  MBpars.GetEikonalParameters().Ymax)*
	      m_norm*sqr(1.+MBpars.GetFFParameters().kappa)),
  m_expnorm(1./4.*MBpars.GetFFParameters().Lambda2*m_norm)
{}

double Analytic_Eikonal::operator()(const double & B) const {
  if (B<0.) return 0.;
  return m_prefactor * exp(-B*B*m_expnorm);
}
