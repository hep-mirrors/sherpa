#include "Running_AlphaQED.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MathTools.H"

#include <iostream>


using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;



namespace APHYTOOLS {
  Running_AlphaQED * aqed  = 0 ;

  const double Running_AlphaQED::A[4]={0.0,0.0,0.00165,0.00221};
  const double Running_AlphaQED::B[4]={0.00835,0.00238,0.00299,0.00293};
  const double Running_AlphaQED::C[4]={1.0,3.927,1.0,1.0};
}


Running_AlphaQED::Running_AlphaQED()
{
  initialised=0;
  Init();
}


void Running_AlphaQED::Init()
{
  if (rpa.consts.IsAQEDRunning()==Switch::Off)
    mode = 0;
  // (mode==0)
  // fixed (splitted) reading parameter "alpha0" and "alpha_eff"
  // (mode==2)
  // fixed (splitted) reading parameter "alpha0" and calculating "alpha_eff"
  else 
    mode = 1;
  // (mode==1)
  // running, reading parameter "alpha0" and and calculating all alpha(Q^2)
  //             (no splitting)


  ecms2         = sqr(rpa.gen.Ecms());

  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());
  // Note: all alphas become inverted (see below) !!!
  alpha_MZ      = dr.GetValue<double>("alpha_QED(MZ)");
  alpha0        = dr.GetValue<double>("alpha_QED(0)");
  alpha_eff     = dr.GetValue<double>("alpha_QED_fixed");
  splitt_scale  = dr.GetValue<double>("splitt_scale");

  // default values
  if (alpha0==NotDefined<double>()) { // Thomson limit (basis for running coupling)
    alpha0=1./137.03599976;  // cf. PDG 2002
    if (mode==1) 
      msg.Error()<<" ERROR: alpha_QED(0) must be defined in oder to use running coupling!"<<std::endl;
    msg.Out()<<" WARNING: using default value for alpha_QED(0)="<<alpha0<<std::endl;
  } 
  else 
    alpha0= 1./alpha0;

  // if scheme default:

  if (alpha_MZ==NotDefined<double>()) {
    alpha_MZ=Aqed(sqr(Flavour(kf::Z).PSmass()));
    if (alpha_eff==NotDefined<double>()) { // alpha_QED fixed
      alpha_eff=Aqed(ecms2);
      msg.Out()<<" WARNING: using default value for alpha_eff="<<alpha_eff;
      msg.Out()<<"          at scale = "<<ecms2<<std::endl;
    }
    else alpha_eff=1./alpha_eff;
  }
  else {
    alpha_MZ=1./alpha_MZ;
 
    if (alpha_eff==NotDefined<double>()) { // alpha_QED fixed
      // use alpha_MZ:
      alpha_eff=alpha_MZ;
      msg.Out()<<" WARNING: using alpha_QED(MZ) value for alpha_QED_fixed="<<alpha_eff<<std::endl;
    }
    else alpha_eff=1./alpha_eff;
  }

  // if scheme not default (e.g. G_mu)
  // calculate alpha_mz (and alpha_eff)


  if (splitt_scale==NotDefined<double>())
    splitt_scale=0.;
  splitt_scale=sqr(splitt_scale);

  // default mode =0 and E_splitt=0
  // if splitted:
  // (a) E_splitt < 0  => always use alpha_eff
  // (b) E_splitt about 1 GeV
  //                   => use alpha0 for small scales (E_splitt ==0 => alpha0 for Q^2=0 only)
  //                   => use alpah_eff for large scales
  // (d) E_splitt >> E_cms
  //                   => always use alpha0

  initialised=1;

  // status
  if (rpa.gen.Tracking()) {
    msg.Out()<<" alpha_QED(0)    = 1./"<<1./alpha0<<std::endl;
    msg.Out()<<" alpha_QED(MZ)   = 1./"<<1./alpha_MZ<<std::endl;
    msg.Out()<<" alpha_QED_fixed = 1./"<<1./alpha_eff<<std::endl;
    if (mode==0) {
      if (splitt_scale>ecms2) 
	msg.Out()<<" always using  alpha_QED(0) "<<std::endl;
      else if (splitt_scale>=0) 
	msg.Out()<<" using  alpha_QED(0) up to Q^2="<<splitt_scale<<std::endl;
      else
	msg.Out()<<" always using  alpha_QED_fixed "<<std::endl;
    }
    else
      msg.Out()<<" using running alpha_QED(Q^2) "<<std::endl;

  }

}


// possibly running coupling (if switched on)
double  Running_AlphaQED::Aqed(double t){
  if (mode) {
    // fixed aQED!
    if (t>splitt_scale)
      return alpha_eff;
    
    return alpha0;
  }
  // running coupling
  return operator()(t);
}


// running Alpha_QED (independent of switches)
double Running_AlphaQED::operator()(double t)
{
  // according to 
  //   R. Kleiss et al.: CERN yellow report 89-08, vol.3 p.129
  //   Hadronic component from: H. Burkhardt et al.: Z. Phys C43 (89) 497
  // cf. HERWIG / PYTHIA
    
  double Q2;
  if (t<0.) Q2=-t; else Q2=t;
  int i;
  if (Q2<0.3) i=0;
  else if (Q2<3.0) i=1;
  else if (Q2<100.0) i=2;
  else i=3;
    
  // leptonic components:
  double sig_lep_gg=alpha0/(3.*M_PI) * 
    (PiGamma(Flavour(kf::e),Q2)+PiGamma(Flavour(kf::mu),Q2)+PiGamma(Flavour(kf::tau),Q2));
  // hadronic component (and light quarks):
  double sig_ha_gg=A[i]+B[i]*log(1+ C[i]*Q2);
  // top quark
  double sig_top_gg=alpha0/(3.*M_PI) * 3. * 
    (PiGamma(Flavour(kf::t),Q2));
  // sum:
  double sigma_gg=sig_lep_gg+sig_ha_gg+sig_top_gg;
  return alpha0/(1.-sigma_gg);
}  

double Running_AlphaQED::PiGamma(const Flavour & fl,double scale) {
  double mass2=sqr(fl.PSmass()); // onshell mass
  double mqs=mass2/scale;
  if (scale==0.) return 0.;
  if (4.*mqs<1.e-3) {
    double test=-5./3.-log(mqs); 
    return test; 
  }
  else if (4.*mqs<=1.) {
    double beta=sqrt(1.-4.*mqs);
    double test=1./3.-(1.+2.*mqs)*(2.+beta*log((1.-beta)/(1.+beta)));
    return test;
  } 
  else { // 4.*mqs>1. 
    return 0.; 
    /*
    double beta=sqrt(4.*mqs-1.);
    return 1./3.-(1.+2.*mqs)*(2.-2.*beta*atan(1./beta));
    */
  }

};
