#include "PDF_Electron.H"
#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Message.H"
#include "MathTools.H"
//#include <iostream>


using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace PDF;

PDF_Electron::PDF_Electron(int mode) {
  if (mode==0) beam = Flavour(kf::e);
          else beam = Flavour(kf::e).bar();
  partons.push_back(beam);
  
  mass     = beam.PSmass();
  alpha    = (*aqed)(sqr(rpa.gen.Ecms()));

  double L = log(sqr(rpa.gen.Ecms()/beam.PSmass()));
  beta     = (*aqed)(sqr(beam.PSmass()))/M_PI*(L-1.);
}


void PDF_Electron::Calculate(const double x, const double Q2) 
{
  xpdf  = 0.;
  alpha = (*aqed)(Q2);
  if (x>=0.999999) return;

  double L       = 2.*log(sqrt(Q2)/mass);
  double beta_e  = 2.*alpha/M_PI*(L-1.);
  double eta     = 2.*alpha/M_PI*L;

  // parameters :
  int izetta=1;
  int order=1;

  double betaS,betaH;
  switch(izetta) {
  case 0:
    beta = beta_e;
    betaS = betaH = eta;
    break;
  case 1:
    beta = betaS = beta_e;
    betaH = eta;
    break;
  default:
    beta = betaS = betaH = beta_e; 
  }
  double gamma   = ::exp(Gammln(1.+beta/2.)); 
  // beta,betaS,betaH


  // Produces collinear bremsstrahlung in exponentiated LLA,
  double S=0.,h0=0.,h1=0.,h2=0.;
  S  = exp(-.5*GAMMA_E*beta+.375*betaS)/gamma*beta/2.; // *pow(1.-x,beta/2.-1.);
  h0 = -.25*(1.+x)*betaH;

  if (order>=2) {  
    h1   = -1./32.*betaH*betaH*
      ((1.+3.*x*x)/(1.-x)*log(x)+4.*(1.+x)*log(1.-x)+5.+x); 
  }
  if (order==3) {
    int i = 1;
    double Li = 0.;
    double del= 1.;
    double del2=1.;
    for (i=1;del2>rpa.gen.Accu();++i) {  
      del *=x;
      del2 = del/double(i*i);
      Li  += del2;
    }
    h2 = -1./384.*betaH*betaH*betaH*
      ((1.+x)*(6.*Li+12.*sqr(log(1.-x))-3.*M_PI*M_PI)+
       1./(1.-x)*(1.5*(1.+8.*x+3.*x*x)*log(x)+6.*(x+5.)*(1.-x)*log(1.-x)+
       12.*(1.+x*x)*log(x)*log(1.-x)-(.5+3.5*x*x)*sqr(log(x))+
		  .25*(39.-24.*x-15.*x*x)));                
  } 

  xpdf = x * (S*pow(1.-x,beta/2.-1.)+(h0+h1+h2));  
  if (x>0.9999) xpdf *= pow(100.,beta/2)/(pow(100.,beta/2)-1.);
}









