#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"

#include <iostream>


using namespace MODEL;
using namespace ATOOLS;


namespace MODEL {
  Running_AlphaQED * aqed  = 0 ;

  const double Running_AlphaQED::m_A[4]={0.0,0.0,0.00165,0.00221};
  const double Running_AlphaQED::m_B[4]={0.00835,0.00238,0.00299,0.00293};
  const double Running_AlphaQED::m_C[4]={1.0,3.927,1.0,1.0};
}


Running_AlphaQED::Running_AlphaQED(const double _alpha0) :
  m_alpha0(_alpha0)
{
  m_type = std::string("Running Coupling");
  m_name  = "Alpha_QED";
  m_defval = _alpha0;
}


double Running_AlphaQED::operator()(double t)
{
  double Q2    = t;
  if (t<0.) Q2 = -t; 

  int i = 3;
  if (Q2<0.3)        i=0;
  else if (Q2<3.0)   i=1;
  else if (Q2<100.0) i=2;
    
  double sig_lep_gg = m_alpha0/(3.*M_PI) * 
    (PiGamma(Flavour(kf_e),Q2)+PiGamma(Flavour(kf_mu),Q2)+PiGamma(Flavour(kf_tau),Q2));
  double sig_ha_gg  = m_A[i] + m_B[i]*log(1+m_C[i]*Q2);
  double sig_top_gg = m_alpha0/(3.*M_PI) * 3. * (PiGamma(Flavour(kf_t),Q2));
  double sigma_gg   = sig_lep_gg+sig_ha_gg+sig_top_gg;

  return m_alpha0/(1.-sigma_gg);
}  

double Running_AlphaQED::PiGamma(const Flavour & fl,double scale) {
  double mass2  = sqr(fl.Mass(true)); // onshell mass
  if(mass2==0.) THROW(fatal_error, "Cannot evolve QED coupling with zero fermion masses");    
  double mqs    = mass2/scale;
  if (scale==0.) return 0.;
  if (4.*mqs<1.e-3) return (-5./3.-log(mqs));
  else if (4.*mqs<=1.) {
    double beta = sqrt(1.-4.*mqs);
    return (1./3.-(1.+2.*mqs)*(2.+beta*log((1.-beta)/(1.+beta))));
  } 
  else { 
    return 0.; 
  }
}
