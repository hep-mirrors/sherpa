#include "Running_AlphaQED.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MathTools.H"

#include <iostream>


using namespace MODEL;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;



namespace MODEL {
  Running_AlphaQED * aqed  = 0 ;

  const double Running_AlphaQED::m_A[4]={0.0,0.0,0.00165,0.00221};
  const double Running_AlphaQED::m_B[4]={0.00835,0.00238,0.00299,0.00293};
  const double Running_AlphaQED::m_C[4]={1.0,3.927,1.0,1.0};
}


Running_AlphaQED::Running_AlphaQED(const double _alpha0,const double _MZ2,const double _split_scale) :
  m_alpha0(_alpha0),m_MZ2(_MZ2),m_split_scale(_split_scale)
{
  m_type = std::string("Running Coupling");

  if (m_alpha0!=1./137.03599976) { // Thomson limit (basis for running coupling)
    msg.Error()<<" WARNING: using value for alpha_QED(0)="<<m_alpha0<<std::endl;
  } 

  m_split_scale = sqr(m_split_scale);
  m_ecms2       = sqr(rpa.gen.Ecms());
  m_alpha_MZ    = Aqed(m_MZ2);
  m_defval      = Aqed(m_MZ2);

  // default mode =0 and E_splitt=0
  // if splitted:
  // (a) E_splitt < 0  => always use alpha_eff
  // (b) E_splitt about 1 GeV
  //                   => use alpha0 for small scales (E_splitt ==0 => alpha0 for Q^2=0 only)
  //                   => use alpah_eff for large scales
  // (d) E_splitt >> E_cms
  //                   => always use alpha0


  // status
  if (rpa.gen.Tracking()) {
    msg.Out()<<" alpha_QED(0)    = 1./"<<1./m_alpha0<<std::endl;
    msg.Out()<<" alpha_QED(MZ)   = 1./"<<1./m_alpha_MZ<<std::endl;
    msg.Out()<<" alpha_QED_fixed = 1./"<<1./m_defval<<std::endl;
  }
}


// possibly running coupling (if switched on)
double  Running_AlphaQED::Aqed(double t){
  return operator()(t);

  if (t>m_split_scale) return m_defval;
  return m_alpha0;
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
    (PiGamma(Flavour(kf::e),Q2)+PiGamma(Flavour(kf::mu),Q2)+PiGamma(Flavour(kf::tau),Q2));
  double sig_ha_gg  = m_A[i] + m_B[i]*log(1+m_C[i]*Q2);
  double sig_top_gg = m_alpha0/(3.*M_PI) * 3. * (PiGamma(Flavour(kf::t),Q2));
  double sigma_gg   = sig_lep_gg+sig_ha_gg+sig_top_gg;

  return m_alpha0/(1.-sigma_gg);
}  

double Running_AlphaQED::PiGamma(const Flavour & fl,double scale) {
  double mass2  = sqr(fl.PSMass()); // onshell mass
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
