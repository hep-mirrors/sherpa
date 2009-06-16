#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

namespace CSSHOWER {
  
  class LF_VVV_FF: public SF_Lorentz {
  public:

    inline LF_VVV_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,int mode=0);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

    double J(const double);

  };

  class LF_VVV_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVV_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,int mode=0);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

    double J(const double,const double,const double);

  };

  class LF_VVV_IF: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVV_IF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,int mode=0);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

    double J(const double,const double,const double);

  };

  class LF_VVV_II: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVV_II(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,int mode=0);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

    double J(const double,const double,const double);

  };

}

#include "ATOOLS/Math/Random.H"
#include "PDF/Main/PDF_Base.H"

using namespace CSSHOWER;
using namespace PDF;
using namespace ATOOLS;

double LF_VVV_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,int mode)
{
  
  double muk2  = sqr(p_ms->Mass(m_flspec))/Q2;
  //the massless case
  double massless = 2. * ( 1./(1.-z+z*y) + 1./(z+y-z*y) -2. + z*(1.-z) );
  if (muk2==0.) {
    double value = p_cf->Coupling(scale,0) * massless;
    if (mode&1) return value/2.0;
    return value * J(y);
  }
  else {
    //the massive case
    double vijk  = sqrt(sqr(2.*muk2+(1.-muk2)*(1.-y))-4.*muk2)/((1.-muk2)*(1.-y));
    double zm = 0.5*(1.- vijk);  
    double zp = 0.5*(1.+ vijk);
    double massive = 2. * ( 1./(1.-z+z*y) + 1./(z+y-z*y) + 
				   (z*(1.-z) - zp*zm - 2.)/vijk );
    if (massive < 0.) {
      //std::cout<<" g -> gg FF mass correction : "<<massive/massless<<"\n"; 
      return 0.;
    }
    massive *= (1.-muk2)/sqrt(Lambda(1.,0.,muk2));
    double value = p_cf->Coupling(scale,0) * massive;
    if (mode&1) return value/2.0;
    return value * J(y);
  }
}

double LF_VVV_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  return 2.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax)));
}

double LF_VVV_FF::OverEstimated(const double z,const double y)
{
  return 2.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) );
}

double LF_VVV_FF::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmax) *
	     pow(m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax),ATOOLS::ran.Get()));
}

double LF_VVV_FF::J(const double y)
{ 
  return (1.-y);
}

double LF_VVV_FI::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,int mode)
{
  double value = 2.0*p_cf->Coupling(scale,0) * ( 1./(1.-z+y) + 1./(z+y) -2. + z*(1.-z) );
  if (mode&1) return value/2.0;
  return value * J(y,eta,scale);
}

double LF_VVV_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax=2.;
  return 2.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
}

double LF_VVV_FI::OverEstimated(const double z,const double y)
{
  return 2.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) ) * m_Jmax;
}

double LF_VVV_FI::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double LF_VVV_FI::J(const double y,const double eta,const double scale)
{ 
  if (scale<sqr(p_ms->Mass(m_flspec)) ||
      scale<p_pdf[m_beam]->Q2Min() || eta/(1.-y)>1.)   return 0.;
  p_pdf[m_beam]->Calculate(eta/(1.-y),scale);
  double fresh = p_pdf[m_beam]->GetXPDF(m_flspec);
  p_pdf[m_beam]->Calculate(eta,scale);
  double old = p_pdf[m_beam]->GetXPDF(m_flspec);
  if (fresh<0.0 || old<0.0 || IsZero(old) || IsZero(fresh)) return 0.; 
  return (1.-y) * fresh/old;
}

double LF_VVV_IF::operator() 
  (const double z,const double y,const double eta,
   const double scale,const double Q2,int mode)
{
  double muk2 = sqr(p_ms->Mass(m_flspec))/Q2;
  double massless = 2. * ( 1./(1.-z+y) + 1./z - 2. +z*(1.-z));
  if (muk2==0.) {
    //the massless case
    double value = p_cf->Coupling(scale,0) * massless;
    if (mode&1) return value/2.0;
    return value * J(z,eta,scale);
  }
  else {
    //the massive case
    double massive = massless - 2. * muk2*y/(z*(1.-y));
    if (massive < 0.) {
      return 0.;
  }
    double value = p_cf->Coupling(scale,0) * massive;
    if (mode&1) return value/2.0;
    return value * J(z,eta,scale);
  }
}

double LF_VVV_IF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 5.; 
  return 2.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
}

double LF_VVV_IF::OverEstimated(const double z,const double y)
{
  return 2.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) ) * m_Jmax;
}

double LF_VVV_IF::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double LF_VVV_IF::J(const double z,const double eta,const double scale) { 
  if (scale<p_pdf[m_beam]->Q2Min() || eta/z>1.)   return 0.;
  p_pdf[m_beam]->Calculate(eta/z,scale);
  double fresh = p_pdf[m_beam]->GetXPDF(m_flavs[0]);
  p_pdf[m_beam]->Calculate(eta,scale);
  double old = p_pdf[m_beam]->GetXPDF(m_flavs[0]);
  if (fresh<0.0 || old<0.0 || IsZero(old) || IsZero(fresh)) return 0.; 
  return fresh/old;
}

double LF_VVV_II::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,int mode)
{
  double value = 2.0 * p_cf->Coupling(scale,0) * ( 1./(1.-z) + 1./z -2. +z*(1.-z));
  if (mode&1) return value/2.0;
  return value * J(z,eta,scale);
}

double LF_VVV_II::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 5.;
  return 2.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
}

double LF_VVV_II::OverEstimated(const double z,const double y)
{
  return 2.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) ) * m_Jmax;
}

double LF_VVV_II::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double LF_VVV_II::J(const double z,const double eta,const double scale)
{ 
  if (scale<p_pdf[m_beam]->Q2Min() || eta/z>1.)   return 0.;
  p_pdf[m_beam]->Calculate(eta/z,scale);
  double fresh = p_pdf[m_beam]->GetXPDF(m_flavs[0]);
  p_pdf[m_beam]->Calculate(eta,scale);
  double old = p_pdf[m_beam]->GetXPDF(m_flavs[0]);
  if (fresh<0.0 || old<0.0 || IsZero(old) || IsZero(fresh)) return 0.; 
  return fresh/old;
}

DECLARE_GETTER(LF_VVV_Getter,"Gauge3",SF_Lorentz,SF_Key);

SF_Lorentz *LF_VVV_Getter::operator()
  (const Parameter_Type &args) const
{
  switch (args.m_type) {
  case cstp::FF: return new LF_VVV_FF(args);
  case cstp::FI: return new LF_VVV_FI(args);
  case cstp::IF: return new LF_VVV_IF(args);
  case cstp::II: return new LF_VVV_II(args);
  case cstp::none: break;
  }
  return NULL;
}

void LF_VVV_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"vvv lorentz functions";
}
