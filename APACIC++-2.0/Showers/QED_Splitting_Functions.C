#include "QED_Splitting_Functions.H"

using namespace APACIC;

// --------- class f_fp -----------------------------
// * fermion to fermion + photon  splitting function
// --------------------------------------------------
f_fp::f_fp(APHYTOOLS::Flavour fermionflavour) {
  flavs[0] = fermionflavour; // a
  flavs[1] = fermionflavour; // b
  flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // c
  qsqr     = AMATOOLS::sqr(fermionflavour.Charge());
  alpha    = 1.;
}

f_fp::f_fp(APHYTOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
  tools (_tools) {
  flavs[0] = fermionflavour; // a
  flavs[1] = fermionflavour; // b
  flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // c
  qsqr     = AMATOOLS::sqr(fermionflavour.Charge());
  alpha    = tools->GetAQEDmax();
}

double f_fp::operator()(double z) {return qsqr*(1.+z*z)/(1.-z);}; 
double f_fp::GetZ()      {
  return 1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),AMATOOLS::ran.Get());
}

double f_fp::GetCoupling()         { return alpha;}

double f_fp::GetCoupling(double t) { return tools->Alpha(t);}


double f_fp::GetWeight(double z,double pt2,bool massterm) { 
  if (!massterm) return 0.5*(1.+z*z);
  return ( (1.+z*z)/2. - 
	   z*AMATOOLS::sqr((1.-z)*flavs[0].PSMass())/
	   (pt2+AMATOOLS::sqr((1.-z)*flavs[0].PSMass())) );
}

double f_fp::CrudeInt(double _zmin, double _zmax) {
  zmin = _zmin;
  zmax = _zmax;
  return 2.*qsqr*alpha*log((1.-zmin)/(1.-zmax));
}



// --------- class f_pf -----------------------------
// * fermion to photon + fermion splitting function 
//   (only used in Initial State Shower)
// --------------------------------------------------

f_pf::f_pf(APHYTOOLS::Flavour fermionflavour) {
  flavs[0] = fermionflavour; // a
  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // b
  flavs[2] = fermionflavour; // c
  qsqr     = AMATOOLS::sqr(fermionflavour.Charge());
  alpha    = 1.;
}

f_pf::f_pf(APHYTOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
  tools (_tools) {
  flavs[0] = fermionflavour; // a
  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // b
  flavs[2] = fermionflavour; // c
  qsqr     = AMATOOLS::sqr(fermionflavour.Charge());
  alpha    = tools->GetAQEDmax();
}

double f_pf::operator()(double z) {return qsqr*(1.+(1.-z)*(1.-z))/z;} 

double f_pf::GetZ()      
{
  return zmin*pow(zmax/zmin,AMATOOLS::ran.Get());
}

double f_pf::GetCoupling()         { return alpha;}
double f_pf::GetCoupling(double t) { return tools->Alpha(t);}
double f_pf::GetWeight(double z,double pt2,bool massterm) { 
  if (!massterm) return 0.5*(1.+AMATOOLS::sqr(1.-z));
  return ( (1.+AMATOOLS::sqr(1.-z))/2. - 
	   (1.-z)*AMATOOLS::sqr(z*flavs[0].PSMass())/
	   (pt2+AMATOOLS::sqr(z*flavs[0].PSMass())) );
}

double f_pf::CrudeInt(double _zmin, double _zmax) {
  zmin = _zmin;
  zmax = _zmax;
  return 2.*qsqr*alpha*log(zmax/zmin);
}


// --------- class f_pf -----------------------------
// * photon to fermion + anti-fermion splitting function
// --------------------------------------------------

p_ff::p_ff(APHYTOOLS::Flavour fermionflavour) {
  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // a
  flavs[1] = fermionflavour; // b
  flavs[2] = fermionflavour.Bar(); // c
  qsqr     = AMATOOLS::sqr(fermionflavour.Charge());
  alpha    = 1.;
}

p_ff::p_ff(APHYTOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
  tools (_tools) {
  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // a
  flavs[1] = fermionflavour; // b
  flavs[2] = fermionflavour.Bar(); // c
  qsqr     = AMATOOLS::sqr(fermionflavour.Charge());
  alpha    = tools->GetAQEDmax();
}

double p_ff::operator()(double z) 
{
  return qsqr*(z*z+ AMATOOLS::sqr(1-z));
}

double p_ff::GetZ()      
{
  return zmin+(zmax-zmin)*AMATOOLS::ran.Get();
}

double p_ff::GetCoupling()         { return alpha;}
double p_ff::GetCoupling(double t) { return tools->Alpha(t);}

double p_ff::GetWeight(double z,double pt2,bool masses) 
{ 
  if (masses) return (*this)(z)/qsqr;
  return (1. - 2.*z*(1.-z)*pt2/(pt2+AMATOOLS::sqr(flavs[1].PSMass())));
}                 

double p_ff::CrudeInt(double _zmin, double _zmax) 
{
  zmin = _zmin;
  zmax = _zmax;
  return (zmax-zmin)*alpha*qsqr;
}
