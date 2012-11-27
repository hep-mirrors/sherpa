#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

// #define USING_DIS_MEC 

namespace CSSHOWER {
  
  class LF_FFV_FF: public SF_Lorentz {
  public:

    inline LF_FFV_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FFV_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

    double CDIS(const double z,const double y);
    double CDISMax();

  public:

    inline LF_FFV_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FFV_IF: public SF_Lorentz {
  protected:

    double m_Jmax;

    double CDIS(const double z,const double y);
    double CDISMax();

  public:

    inline LF_FFV_IF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FFV_II: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_FFV_II(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FVF_FF: public SF_Lorentz {
  public:

    inline LF_FVF_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FVF_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

    double CDIS(const double z,const double y);
    double CDISMax();

  public:

    inline LF_FVF_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FVF_IF: public SF_Lorentz {
  protected:

    double m_Jmax;

    double CDIS(const double z,const double y);
    double CDISMax();

  public:

    inline LF_FVF_IF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_FVF_II: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_FVF_II(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VFF_FF: public SF_Lorentz {
  public:

    inline LF_VFF_FF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VFF_FI: public SF_Lorentz {
  protected:

    double m_Jmax;

    double CDIS(const double z,const double y);
    double CDISMax();

  public:

    inline LF_VFF_FI(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VFF_IF: public SF_Lorentz {
  protected:

    double m_Jmax;

    double CDIS(const double z,const double y);
    double CDISMax();

  public:

    inline LF_VFF_IF(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

  class LF_VFF_II: public SF_Lorentz {
  protected:

    double m_Jmax;

  public:

    inline LF_VFF_II(const SF_Key &key): SF_Lorentz(key) {}

    double operator()(const double,const double,const double,
		      const double,const double);
    double OverIntegrated(const double,const double,
			  const double,const double);
    double OverEstimated(const double,const double);
    double Z();

  };

}

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "ATOOLS/Math/Random.H"

using namespace CSSHOWER;
using namespace ATOOLS;

double LF_FFV_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double muij2 = sqr(p_ms->Mass(m_flavs[0]))/Q2;
  double mi2   = sqr(p_ms->Mass(m_flavs[1]));
  double mui2  = mi2/Q2;
  double muk2  = sqr(p_ms->Mass(m_flspec))/Q2;
  //the massless case
  double massless = ( 2./(1.-z+z*y) - (1.+z) );
  if (muij2==0. && mui2==0. && muk2==0.) {
    double longpol = 0.5 * ( 1. - z );
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JFF(y);
  }
  else {
    //the massive case
    double vtijk = sqrt(Lambda(1.,muij2,muk2))/(1.-muij2-muk2);
    double vijk  = sqrt(sqr(2.*muk2+(1.-mui2-muk2)*(1.-y))-4.*muk2)/((1.-mui2-muk2)*(1.-y));
    double pipj  = Q2*(1.0-mui2-muk2)*y/2.0;
    double massive = ( 2./(1.-z+z*y) - vtijk/vijk * (1.+z + mi2/pipj) );
    if (massive < 0.) {
      //std::cout<<" Q -> Qg FF mass correction : "<<massive/massless<<" )\n"; 
      return 0.;
    }
    massive *= 1./((1.-mui2-muk2)+1./y*(mui2-muij2))*sqr(1.-mui2-muk2)/sqrt(Lambda(1.,muij2,muk2));
    double longpol = 0.5 * ( 1. - z );
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JFF(y);
  } 
}

double LF_FFV_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) *log((1.-zmin)/(1.-zmax));
}

double LF_FFV_FF::OverEstimated(const double z,const double y)
{
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))/(1.-z);
}

double LF_FFV_FF::Z()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran->Get());
}

double LF_FFV_FI::CDIS(const double z,const double y)
{
#ifdef USING_DIS_MEC 
  return y*(1.+3.*z*(1.-y));
#else
  return 0.0;
#endif
}

double LF_FFV_FI::CDISMax()
{
  return 0.;
}

double LF_FFV_FI::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{  
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  //the massless case
  double massless = ( 2./(1.-z+y) - 1.-z + CDIS(z,y) );
  if (mi2==0.) {
    double longpol = 0.5 * ( 1. - z );
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JFI(y,eta,scale);
  }
  else {
    //the massive case
    double yt = y/(1.0-y), pipj = yt*(Q2+mi2)/2.0;
    double massive = massless - mi2/pipj;
    if (massive < 0.) {
      //std::cout<<" Q -> Qg FI mass correction : "<<massive<<" / "<<massless<<"\n";
      return 0.;
    }
    double longpol = 0.5 * ( 1. - z );
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JFI(y,eta,scale);
  }
}

double LF_FFV_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax=m_flspec.Kfcode()<3?5.:1.;
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1))*log((1.-zmin)/(1.-zmax)) * m_Jmax;
}

double LF_FFV_FI::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1))/(1.-z) * m_Jmax;
}

double LF_FFV_FI::Z()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran->Get());
}

double LF_FFV_IF::CDIS(const double z,const double y)
{
#ifdef USING_DIS_MEC 
  return y*(1.+3.*z*(1.-y));
#else
  return 0.0;
#endif
}

double LF_FFV_IF::CDISMax()
{
  return 0.;
}

double LF_FFV_IF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double value = 2.0 * p_cf->Coupling(scale,0) * ( 2./(1.-z+y) - (1.+z) + CDIS(z,y) )
    + p_cf->Coupling(scale,1) * 0.5 * ( 1. - z );
  return value * JIF(z,y,eta,scale);
}

double LF_FFV_IF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode()<3?5.:1.; 
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1)) * log((1.-zmin)/(1.-zmax)) * m_Jmax;
}

double LF_FFV_IF::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1))/(1.-z) * m_Jmax;
}

double LF_FFV_IF::Z()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran->Get());
}

double LF_FFV_II::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double value = 2.0 * p_cf->Coupling(scale,0) * ( 2./(1.-z) - (1.+z) )
    + p_cf->Coupling(scale,1) * 0.5 * ( 1. - z );
  return value * JII(z,y,eta,scale);
}

double LF_FFV_II::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax; 
  m_Jmax = m_flavs[0].Kfcode()<3?5.:1.;
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) * log((1.-zmin)/(1.-zmax)) * m_Jmax;
}

double LF_FFV_II::OverEstimated(const double z,const double y)
{
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))/(1.-z) * m_Jmax;
}

double LF_FFV_II::Z()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran->Get());
}

double LF_FVF_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double muij2 = sqr(p_ms->Mass(m_flavs[0]))/Q2;
  double mj2   = sqr(p_ms->Mass(m_flavs[2]));
  double muj2  = mj2/Q2;
  double muk2  = sqr(p_ms->Mass(m_flspec))/Q2;
  //the massless case
  double massless = ( 2./(z+y-z*y) - 2. + z );
  if (muij2==0. && muj2==0. && muk2==0.) {
    double longpol = 0.5 * z;
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JFF(y);
  }
  else {
    //the massive case
    double vtijk = sqrt(Lambda(1.,muij2,muk2))/(1.-muij2-muk2);
    double vijk  = sqrt(sqr(2.*muk2+(1.-muj2-muk2)*(1.-y))-4.*muk2)/((1.-muj2-muk2)*(1.-y));
    double pipj  = Q2*(1.0-muj2-muk2)*y/2.0;
    double massive = ( 2./(z+y-z*y) - vtijk/vijk * (2.+z + mj2/pipj) );
    if (massive < 0.) {
      //std::cout<<" Q -> gQ FF mass correction : "<<massive/massless<<"\n"; 
      return 0.;
    }
    massive *= 1./((1.-muj2-muk2)+1./y*(muj2-muij2))*sqr(1.-muj2-muk2)/sqrt(Lambda(1.,muij2,muk2));
    double longpol = 0.5 * z;
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JFF(y);
  }
}

double LF_FVF_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))*log(zmax/zmin);
}

double LF_FVF_FF::OverEstimated(const double z,const double y)
{
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))/z;
}

double LF_FVF_FF::Z()
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran->Get());
}

double LF_FVF_FI::CDIS(const double z,const double y)
{
#ifdef USING_DIS_MEC 
  return y*(1.+3.*z*(1.-y));
#else
  return 0.0;
#endif
}

double LF_FVF_FI::CDISMax()
{
  return 0.;
}

double LF_FVF_FI::operator() (const double z,const double y,
			    const double eta, const double scale,const double Q2) {
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  //the massless case
  double massless = (2./(z+y) - 2.+z + CDIS(1.-z,y));
  if (mj2==0.) {
    double longpol = 0.5 * z;
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JFI(y,eta,scale);
  }
  else {
    //the massive case
    double yt = y/(1.0-y), pipj = yt*(Q2+mj2)/2.0;
    double massive = massless - mj2/pipj;
    if (massive < 0.) {
      //std::cout<<" Q -> gQ FI mass correction : "<<massive/massless<<"\n";
      return 0.;
    }
    double longpol = 0.5 * z;
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JFI(y,eta,scale);
  }
}
double LF_FVF_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax=m_flspec.Kfcode()<3?5.:1.;
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1)) * log(zmax/zmin) * m_Jmax;
}

double LF_FVF_FI::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0)*(2.0+CDISMax()) + 0.5*p_cf->MaxCoupling(1))/z * m_Jmax;
}

double LF_FVF_FI::Z()
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran->Get());
}

double LF_FVF_IF::CDIS(const double z,const double y)
{
#ifdef USING_DIS_MEC 
  return y*(1.+3.*z*(1.-y));
#else
  return 0.0;
#endif
}

double LF_FVF_IF::CDISMax()
{
  return 0.;
}

double LF_FVF_IF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double mk2  = sqr(p_ms->Mass(m_flspec));
  double muk2 = mk2*z/(Q2+mk2); 
  double massless = ( 2./z - 2. +z + CDIS(1.-z,y));
  if (muk2==0.) {
    //the massless case
    double longpol = 0.5 * z;
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JIF(z,y,eta,scale);
  }
  else {
    //the massive case
    double massive = massless - 2.*muk2*y/(z*(1.-y));
    if (massive < 0.) {
      //std::cout<<" q -> gq IF mass correction : "<<massive/massless<<" ( massive : "<<massive<<" massless : "<<massless<< " )\n"; 
      return 0.;
    }
    double longpol = 0.5 * z;
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JIF(z,y,eta,scale);
  }
}

double LF_FVF_IF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  double fresh = p_sf->GetXPDF(scale,xbj,m_flavs[0],m_beam);
  double old   = p_sf->GetXPDF(scale,xbj,m_flavs[1],m_beam,1);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.;
  m_Jmax = 5.*fresh/old;
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1)) * log(zmax/zmin) * m_Jmax;
}

double LF_FVF_IF::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0)*(2.+CDISMax()) + 0.5*p_cf->MaxCoupling(1))/z * m_Jmax;
}

double LF_FVF_IF::Z()
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran->Get());
}

double LF_FVF_II::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double value = 2.0 * p_cf->Coupling(scale,0) * ( 2./z - 2. +z )
    + p_cf->Coupling(scale,1) * 0.5 * z;
  return value * JII(z,y,eta,scale);
}

double LF_FVF_II::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
 {
  m_zmin = zmin; m_zmax = zmax;
  double fresh = p_sf->GetXPDF(scale,xbj,m_flavs[0],m_beam);
  double old   = p_sf->GetXPDF(scale,xbj,m_flavs[1],m_beam,1);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  m_Jmax = 5.*fresh/old;
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))* log(zmax/zmin) * m_Jmax;
}

double LF_FVF_II::OverEstimated(const double z,const double y)
{
  return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1))/z * m_Jmax;
}

double LF_FVF_II::Z()
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran->Get());
}

double LF_VFF_FF::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double mui2  = sqr(p_ms->Mass(m_flavs[1]))/Q2;
  double muj2  = sqr(p_ms->Mass(m_flavs[2]))/Q2;
  double muk2  = sqr(p_ms->Mass(m_flspec))/Q2;
  //the massless case 
  double massless = (1.-2.*z*(1.-z));
  double longpol = 0.5;
  if (mui2==0. && muj2==0. && muk2==0.) {
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JFF(y);
  }
  else {
    //the massive case
    double fac  = 1.-mui2-muj2-muk2;
    double viji = sqrt(sqr(fac*y)-4.*mui2*muj2)/(fac*y+2.*mui2);
    double vijk = sqrt(sqr(2.*muk2+fac*(1.-y))-4.*muk2)/(fac*(1.-y));
    double frac = (2.*mui2+fac*y)/(2.*(mui2+muj2+fac*y));
    double zm = frac*(1.- viji*vijk);  
    double zp = frac*(1.+ viji*vijk);
    double massive = vijk * (1.- 2.*(z*(1.-z) - 2.* zp*zm));
    massive *= 1./((1.-mui2-muj2-muk2)+1./y*(mui2+muj2))*sqr(1.-mui2-muj2-muk2)/sqrt(Lambda(1.,0.,muk2));
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JFF(y);
  }
}
  
double LF_VFF_FF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  return (2.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) * (m_zmax-m_zmin);
}

double LF_VFF_FF::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1));
}

double LF_VFF_FF::Z() {
  return m_zmin + (m_zmax-m_zmin)*ATOOLS::ran->Get();
}

double LF_VFF_FI::CDIS(const double z,const double y)
{
#ifdef USING_DIS_MEC 
  return 4.*y*z*(1.-z);
#else
  return 0.0;
#endif
}

double LF_VFF_FI::CDISMax()
{
#ifdef USING_DIS_MEC 
  return 0.5;
#else
  return 0.0;
#endif
}

double LF_VFF_FI::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double muQ2 = sqr(p_ms->Mass(m_flavs[1]))*(1.-y)/Q2;
  //the massless case 
  double massless = ( (1.-2.*z*(1.-z))*(1.-0.5/z*CDIS(y,z)) + CDIS(z,y) );
  double longpol = 0.5;
  if (muQ2==0.) {
    double value = 2.0 * p_cf->Coupling(scale,0) * massless + p_cf->Coupling(scale,1) * longpol;
    return value * JFI(y,eta,scale);
  }
  else {
    //the massive case
    double delta   = sqrt(sqr(y-2.*muQ2)-4.*muQ2*muQ2)/y; 
    double zp      = 0.5 * (1. + delta);
    double zm      = 0.5 * (1  - delta); 
    double massive = (1.-2.*(zp-z)*(z-zm));
    if (massive < 0.) std::cout<<" massive V_FF FI < 0. "<<massive<<std::endl; 
    double value = 2.0 * p_cf->Coupling(scale,0) * massive + p_cf->Coupling(scale,1) * longpol;
    return value * JFI(y,eta,scale);
  }
}

double LF_VFF_FI::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax=m_flspec.Kfcode()<3?5.:1.;
  return (2.0*p_cf->MaxCoupling(0)*(1.+CDISMax()) + 0.5*p_cf->MaxCoupling(1))* (m_zmax-m_zmin) * m_Jmax;
}

double LF_VFF_FI::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0)*(1.+CDISMax()) + 0.5*p_cf->MaxCoupling(1))* m_Jmax;
}

double LF_VFF_FI::Z()
{
  return m_zmin + (m_zmax-m_zmin)*ATOOLS::ran->Get();
}

double LF_VFF_IF::CDIS(const double z,const double y)
{
#ifdef USING_DIS_MEC 
  return 4.*y*z*(1.-z);
#else
  return 0.0;
#endif
}

double LF_VFF_IF::CDISMax()
{
#ifdef USING_DIS_MEC 
  return 0.5;
#else
  return 0.0;
#endif
}

double LF_VFF_IF::operator() 
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double value = 2.0 * p_cf->Coupling(scale,0) * ( (1.-2.*z*(1.-z))*(1.0-0.5/z*CDIS(y,z)) + CDIS(z,y) )
    + p_cf->Coupling(scale,1) * 0.5;
  return value * JIF(z,y,eta,scale);
}

double LF_VFF_IF::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  double fresh = p_sf->GetXPDF(scale,xbj,m_flavs[0],m_beam);
  double old   = p_sf->GetXPDF(scale,xbj,m_flavs[1],m_beam,1);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  m_Jmax = 5.*fresh/old; 
  return (2.0*p_cf->MaxCoupling(0)*(1.+CDISMax()) + 0.5*p_cf->MaxCoupling(1)) * (m_zmax-m_zmin) * m_Jmax;
}

double LF_VFF_IF::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0)*(1.+CDISMax()) + 0.5*p_cf->MaxCoupling(1)) * m_Jmax;
}

double LF_VFF_IF::Z()
{
  return m_zmin + (m_zmax-m_zmin)*ATOOLS::ran->Get();
}

double LF_VFF_II::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2)
{
  double value = 2.0 * p_cf->Coupling(scale,0) * (1.-2.*z*(1.-z))
    + p_cf->Coupling(scale,1) * 0.5;
  return value * JII(z,y,eta,scale);
}

double LF_VFF_II::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  m_zmin = zmin; m_zmax = zmax;
  double fresh = p_sf->GetXPDF(scale,xbj,m_flavs[0],m_beam);
  double old   = p_sf->GetXPDF(scale,xbj,m_flavs[1],m_beam,1);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  m_Jmax = 5.*fresh/old;
  return (2.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) * (m_zmax-m_zmin) * m_Jmax;
}

double LF_VFF_II::OverEstimated(const double z,const double y)
{
  return (2.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) * m_Jmax;
}

double LF_VFF_II::Z()
{
  return m_zmin + (m_zmax-m_zmin)*ATOOLS::ran->Get();
}

DECLARE_GETTER(LF_FFV_Getter,"Gamma",SF_Lorentz,SF_Key);

SF_Lorentz *LF_FFV_Getter::operator()
  (const Parameter_Type &args) const
{
  if (args.m_col<0) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_FFV_FF(args);
    case cstp::FI: return new LF_FFV_FI(args);
    case cstp::IF: return new LF_FFV_IF(args);
    case cstp::II: return new LF_FFV_II(args);
    case cstp::none: break;
    }
  }
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==1) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_FVF_FF(args);
    case cstp::FI: return new LF_FVF_FI(args);
    case cstp::IF: return new LF_FVF_IF(args);
    case cstp::II: return new LF_FVF_II(args);
    case cstp::none: break;
    }
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    switch (args.m_type) {
    case cstp::FF: return new LF_VFF_FF(args);
    case cstp::FI: return new LF_VFF_FI(args);
    case cstp::IF: return new LF_VFF_IF(args);
    case cstp::II: return new LF_VFF_II(args);
    case cstp::none: break;
    }
  }
  return NULL;
}

void LF_FFV_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ffv lorentz functions";
}
