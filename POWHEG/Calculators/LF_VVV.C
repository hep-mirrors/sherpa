#include "POWHEG/Showers/Splitting_Function_Base.H"

namespace POWHEG {
  
  class LF_VVV_FF: public SF_Lorentz {
  public:

    inline LF_VVV_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,ATOOLS::Cluster_Amplitude *const sub);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VVV_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVV_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,ATOOLS::Cluster_Amplitude *const sub);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VVV_IF: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVV_IF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,ATOOLS::Cluster_Amplitude *const sub);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VVV_II: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VVV_II(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double,ATOOLS::Cluster_Amplitude *const sub);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

}

#include "ATOOLS/Math/Random.H"

using namespace POWHEG;
using namespace PDF;
using namespace ATOOLS;

double LF_VVV_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,Cluster_Amplitude *const sub)
{
  
  double muk2  = sqr(p_ms->Mass(m_flspec))/Q2;
  //the massless case
  double massless = 2. * ( 1./(1.-z+z*y) + 1./(z+y-z*y) -2. + z*(1.-z) );
  if (muk2==0.) {
    double value = 2.0 * p_cf->Coupling(scale,0,sub) * massless;
    return value * JFF(y);
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
    double value = 2.0 * p_cf->Coupling(scale,0,sub) * massive;
    return value * JFF(y);
  }
}

double LF_VVV_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  return 4.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax)));
}

double LF_VVV_FF::OverEstimated(const double z,const double y)
{
  return 4.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) );
}

double LF_VVV_FF::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow(m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax),ATOOLS::ran.Get()));
}

double LF_VVV_FI::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,Cluster_Amplitude *const sub)
{
  double value = 4.0*p_cf->Coupling(scale,0,sub) * ( 1./(1.-z+y) + 1./(z+y) -2. + z*(1.-z) );
  return value * JFI(y,eta,scale,sub);
}

double LF_VVV_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax=m_flspec.Kfcode()<3?5.:1.;
  return 4.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
}

double LF_VVV_FI::OverEstimated(const double z,const double y)
{
  return 4.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) ) * m_Jmax;
}

double LF_VVV_FI::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double LF_VVV_IF::operator() 
  (const double z,const double y,const double eta,
   const double scale,const double Q2,Cluster_Amplitude *const sub)
{
  double muk2 = sqr(p_ms->Mass(m_flspec))/Q2;
  double massless = 2. * ( 1./(1.-z+y) + 1./z - 2. +z*(1.-z));
  if (muk2==0.) {
    //the massless case
    double value = 2.0 * p_cf->Coupling(scale,0,sub) * massless;
    return value * JIF(z,y,eta,scale,sub);
  }
  else {
    //the massive case
    double massive = massless - 2. * muk2*y/(z*(1.-y));
    if (massive < 0.) {
      return 0.;
  }
    double value = 2.0 * p_cf->Coupling(scale,0,sub) * massive;
    return value * JIF(z,y,eta,scale,sub);
  }
}

double LF_VVV_IF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1.; 
  return 4.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
}

double LF_VVV_IF::OverEstimated(const double z,const double y)
{
  return 4.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) ) * m_Jmax;
}

double LF_VVV_IF::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double LF_VVV_II::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,Cluster_Amplitude *const sub)
{
  double value = 4.0 * p_cf->Coupling(scale,0,sub) * ( 1./(1.-z) + 1./z -2. +z*(1.-z));
  return value * JII(z,y,eta,scale,sub);
}

double LF_VVV_II::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1.;
  return 4.*p_cf->MaxCoupling(0) * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
}

double LF_VVV_II::OverEstimated(const double z,const double y)
{
  return 4.*p_cf->MaxCoupling(0) * ( 1./(z*(1.-z)) ) * m_Jmax;
}

double LF_VVV_II::Z()
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

namespace POWHEG {

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

}
