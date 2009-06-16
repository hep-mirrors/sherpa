#include "APACIC++/Showers/QED_Splitting_Functions.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "APACIC++/Showers/Sudakov_Tools.H"

using namespace APACIC;

// --------- class f_fp -----------------------------
// * fermion to fermion + photon  splitting function
// --------------------------------------------------
f_fp::f_fp(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour fermionflavour):
  Splitting_Function(ms)
{
  m_flavs[0] = fermionflavour; // a
  m_flavs[1] = fermionflavour; // b
  m_flavs[2] = ATOOLS::Flavour(kf_photon); // c
  m_qsqr     = ATOOLS::sqr(fermionflavour.Charge());
  m_alpha    = 1.;
}

f_fp::f_fp(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
  Splitting_Function(ms), p_tools (_tools) {
  m_flavs[0] = fermionflavour; // a
  m_flavs[1] = fermionflavour; // b
  m_flavs[2] = ATOOLS::Flavour(kf_photon); // c
  m_qsqr     = ATOOLS::sqr(fermionflavour.Charge());
  m_alpha    = p_tools->GetAQEDmax();
}

double f_fp::operator()(double z) {return m_qsqr*(1.+z*z)/(1.-z);} 
double f_fp::GetZ()      
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());
}

double f_fp::GetCoupling()         { return m_alpha;}

double f_fp::GetCoupling(double t) { return p_tools->Alpha(t);}


double f_fp::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+z*z);
  return ( (1.+z*z)/2. - 
	   z*ATOOLS::sqr((1.-z)*p_ms->Mass(m_flavs[0]))/
	   (pt2+ATOOLS::sqr((1.-z)*p_ms->Mass(m_flavs[0]))) );
}

double f_fp::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*m_qsqr*m_alpha*log((1.-m_zmin)/(1.-m_zmax));
}



// --------- class f_pf -----------------------------
// * fermion to photon + fermion splitting function 
//   (only used in Initial State Shower)
// --------------------------------------------------

f_pf::f_pf(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour fermionflavour) :
  Splitting_Function(ms)
{
  m_flavs[0] = fermionflavour; // a
  m_flavs[1] = ATOOLS::Flavour(kf_photon); // b
  m_flavs[2] = fermionflavour; // c
  m_qsqr     = ATOOLS::sqr(fermionflavour.Charge());
  m_alpha    = 1.;
}

f_pf::f_pf(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
  Splitting_Function(ms), p_tools (_tools) 
{
  m_flavs[0] = fermionflavour; // a
  m_flavs[1] = ATOOLS::Flavour(kf_photon); // b
  m_flavs[2] = fermionflavour; // c
  m_qsqr     = ATOOLS::sqr(fermionflavour.Charge());
  m_alpha    = p_tools->GetAQEDmax();
}

double f_pf::operator()(double z) {return m_qsqr*(1.+(1.-z)*(1.-z))/z;} 

double f_pf::GetZ()      
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double f_pf::GetCoupling()         { return m_alpha;}
double f_pf::GetCoupling(double t) { return p_tools->Alpha(t);}

double f_pf::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+ATOOLS::sqr(1.-z));
  return ( (1.+ATOOLS::sqr(1.-z))/2. - 
	   (1.-z)*ATOOLS::sqr(z*p_ms->Mass(m_flavs[0]))/
	   (pt2+ATOOLS::sqr(z*p_ms->Mass(m_flavs[0]))) );
}

double f_pf::CrudeInt(double _zmin, double _zmax) {
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*m_qsqr*m_alpha*log(m_zmax/m_zmin);
}


// --------- class f_pf -----------------------------
// * photon to fermion + anti-fermion splitting function
// --------------------------------------------------

p_ff::p_ff(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour fermionflavour) :
  Splitting_Function(ms)
{
  m_flavs[0] = ATOOLS::Flavour(kf_photon); // a
  m_flavs[1] = fermionflavour; // b
  m_flavs[2] = fermionflavour.Bar(); // c
  m_qsqr     = ATOOLS::sqr(fermionflavour.Charge());
  m_alpha    = 1.;
}

p_ff::p_ff(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
  Splitting_Function(ms), p_tools (_tools) 
{
  m_flavs[0] = ATOOLS::Flavour(kf_photon); // a
  m_flavs[1] = fermionflavour; // b
  m_flavs[2] = fermionflavour.Bar(); // c
  m_qsqr     = ATOOLS::sqr(fermionflavour.Charge());
  m_alpha    = p_tools->GetAQEDmax();
}

double p_ff::operator()(double z) 
{
  return m_qsqr*(z*z+ ATOOLS::sqr(1-z));
}

double p_ff::GetZ()      
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran.Get();
}

double p_ff::GetCoupling()         { return m_alpha;}
double p_ff::GetCoupling(double t) { return p_tools->Alpha(t);}

double p_ff::GetWeight(double z,double pt2,bool masses) 
{ 
  if (masses) return (*this)(z)/m_qsqr;
  return (1. - 2.*z*(1.-z)*pt2/(pt2+ATOOLS::sqr(p_ms->Mass(m_flavs[1]))));
}                 

double p_ff::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return (m_zmax-m_zmin)*m_alpha*m_qsqr;
}
