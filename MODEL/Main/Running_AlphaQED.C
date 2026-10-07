#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"


#include <iostream>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;



namespace MODEL {
  Running_AlphaQED * aqed  = 0 ;

  const double Running_AlphaQED::m_A[4]={0.0,0.0,0.00165,0.00221};
  const double Running_AlphaQED::m_B[4]={0.00835,0.00238,0.00299,0.00293};
  const double Running_AlphaQED::m_C[4]={1.0,3.927,1.0,1.0};
 #ifdef USING__HADALPHAQED
    extern "C"{
      void hadr5x_(double *e, double *st2,
                  double *der, double *errdersta,
                  double *errdersys, double *deg,
                  double *errdegsta, double *errdegsys);
    }
    #endif

}



Running_AlphaQED::Running_AlphaQED(const double _alpha0) :
  m_alpha0(_alpha0)
{
  m_type = std::string("Running Coupling");
  m_name  = "Alpha_QED";
  m_defval = _alpha0;
  Scoped_Settings s{ Settings::GetMainSettings()[m_name] };
  m_mode     = s["VPMODE"].SetDefault("Default").Get<vpmode::code>();
#ifndef USING__HADALPHAQED
  if (m_mode==vpmode::full || m_mode==vpmode::hp || m_mode==vpmode::lp)
    THROW(fatal_error, "Alpha_QED: VPMODE: Full, HP and LP need Jegerlehner's "
                       "hadr5x, i.e. Sherpa configured with "
                       "-DSHERPA_ENABLE_ALPHAQEDHAD=ON.");
#endif
}


double Running_AlphaQED::operator()(double t)
{
  double Q2    = t;
  if (t<0.) Q2 = -t; 

  int i = 3;
  if (Q2<0.3)        i=0;
  else if (Q2<3.0)   i=1;
  else if (Q2<100.0) i=2;

  if (IsZero(t) || m_mode==vpmode::off) return m_alpha0;

  double sig_lep_gg = m_alpha0/(3.*M_PI) * 
    (PiGamma(Flavour(kf_e),Q2)+PiGamma(Flavour(kf_mu),Q2)+PiGamma(Flavour(kf_tau),Q2));
  double sig_ha_gg  = m_A[i] + m_B[i]*log(1+m_C[i]*Q2);
  double sig_top_gg = m_alpha0/(3.*M_PI) * 3. * (PiGamma(Flavour(kf_t),Q2));
  double sigma_gg   = sig_lep_gg+sig_ha_gg+sig_top_gg;

  #ifdef USING__HADALPHAQED
    double delta_r(0.),errdersta,errdersys,deg,errdegsta,errdegsys;
    // Reference value of sin^2(theta_W) for hadr5x (hadr5x19.f: "st2=0.23153 is
    // the reference value")
    double sin2 = 0.23153;
    if (m_mode==vpmode::full || m_mode==vpmode::hp) {
      // The routine takes E = signed ENERGY in GeV:
      //   "E energy ( momentum transfer ): E>0 timelike , E<0 spacelike"
      //   (hadr5x19.f:40)
      double e = (t<0. ? -sqrt(Q2) : sqrt(Q2));
      hadr5x_(&e, &sin2, &delta_r, &errdersta,
              &errdersys, &deg, &errdegsta, &errdegsys);
    }
    switch(m_mode){
      case vpmode::full:
        sigma_gg = sig_lep_gg + delta_r + sig_top_gg;
        break;
      case vpmode::lp:
        sigma_gg = sig_lep_gg;
        break;
      case vpmode::hp:
        sigma_gg = delta_r;
        break;
      case vpmode::off:
      case vpmode::builtin:
        break;   // off returned above; builtin keeps the built-in sigma_gg
      }
  #endif

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

void Running_AlphaQED::PrintSummary()
{
  msg_Info()<<"Set \\alpha according to EW scheme"
            <<"\n  1/\\alpha(0)   = "<<1./m_alpha0
            <<"\n  1/\\alpha(def) = "<<1./m_defval<<"\n";
}

std::istream &MODEL::operator>>(std::istream &str,vpmode::code &vp)
{
  std::string tag;
  str>>tag;
  if      (tag=="None")    vp=vpmode::off;
  else if (tag=="Full")    vp=vpmode::full;
  else if (tag=="HP")      vp=vpmode::hp;
  else if (tag=="LP")      vp=vpmode::lp;
  else if (tag=="Default") vp=vpmode::builtin;
  else
    THROW(fatal_error, "Unknown Alpha_QED: VPMODE '"+tag+
                       "'. Use None, Full, HP or LP, or leave it unset.");
  return str;
}
