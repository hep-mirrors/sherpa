#include "METOOLS/Explicit/Dipole_Terms.H"

#include "METOOLS/Explicit/Dipole_Info.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

double Lam(const double &a,const double &b,const double &c)
{
  return sqr(a-b-c)-4.0*b*c;
}

I_Args::I_Args(const ATOOLS::Vec4D &pij,const ATOOLS::Vec4D &pk,
	       const double &_mij,const double &_mk):
  mij(_mij), mij2(mij*mij), mk(_mk), mk2(mk*mk),
  r(0.0), rj2(0.0), rk2(0.0)
{
  is=pij[0]<0.0||pk[0]<0.0;
  s=dabs(2.0*pij*pk);
  Q=sqrt(Q2=(is?s:s+mij2+mk2));
  v=sqrt(1.0-sqr(2.0*mij*mk)/(s*s));
  if (mij && mk) {
    r=sqrt(r2=(1.0-v)/(1.0+v));
    rj2=(1.0-v+2.0*mij2/s)/(1.0+v+2.0*mij2/s);
    rk2=(1.0-v+2.0*mk2/s)/(1.0+v+2.0*mk2/s);
  }
  Qa2=s;
}

void I_Args::Swap()
{
  std::swap<double>(mij,mk);
  std::swap<double>(mij2,mk2);
  std::swap<double>(rj2,rk2);
}

NLO_Value METOOLS::FFVS(const I_Args &a,const Dipole_Info *info)
{
  if (a.mij==0.0 && a.mk==0.0) return NLO_Value(1.0,0.0,0.0);
  if (a.mij==0 ^ a.mk==0.0) {
    double m(a.mij?a.mij:a.mk), m2(m*m);
    return NLO_Value
      (0.5,0.5*log(m2/a.s),
       -0.25*sqr(log(m2/a.s))-sqr(M_PI)/12.0
       -0.5*log(m2/a.s)*log(a.s/a.Q2)
       -0.5*log(m2/a.Q2)*log(a.s/a.Q2));
  }
  return NLO_Value
    (0.0,log(a.r)/a.v,1.0/a.v*
     (-0.25*sqr(log(a.rj2))-0.25*sqr(log(a.rk2))
      -sqr(M_PI)/6.0+log(a.r)*log(a.Q2/a.s)));
}

double METOOLS::FFAE(const I_Args &ia,const Dipole_Info *info)
{
  if ((ia.is&&info->Massive()) || info->AMax()==1.0) return 0.0;
  if (ia.mij==0.0 && ia.mk==0.0) return -0.5*sqr(log(info->AMax()));
  if (ia.mij==0.0) {
    double al(info->AMax()), yp((ia.Q-ia.mk)/(ia.Q+ia.mk));
    double xp(yp*(1.0-al)+sqrt((1.0-al)*(1.0-al*yp*yp)));
    return 0.5*sqr(log((1.0-yp*yp+2.0*xp*yp)/(1.0-sqr(yp-xp))))
      -sqr(log((1.0+yp-xp)/(1.0+yp)))
      +2.0*(log((1.0+yp)/2.0)*log((1.0-yp+xp)/(1.0-yp))
	    +log((1.0+yp)/(2.0*yp))*log((1.0-yp*yp+2.0*xp*yp)/(1.0-yp*yp))
	    +DiLog((1.0-yp)/(1.0+yp))-DiLog((1.0-yp*yp+2.0*xp*yp)/sqr(1.0+yp))
	    +DiLog((1.0-yp+xp)/2.0)-DiLog((1.0-yp)/2.0));
  }
  if (ia.mk==0.0) {
    double al(info->AMax());
    return -log(al)*log(ia.mij2/ia.Q2)
      -DiLog(1.0-ia.Q2/ia.mij2)+DiLog(al*(1.0-ia.Q2/ia.mij2));
  }
  double al(info->AMax()), d(0.5*ia.s/ia.Q2), a(ia.mk/ia.Q/d);
  double b((1.0-ia.mk/ia.Q)/d), c(ia.mk*(ia.Q-ia.mk)/ia.Q2/d), yp(1.0-c);
  double x(yp*(1.0-al)+sqrt
	   (yp*(1.0-al)*(1.0/yp-al*yp+
			 4.0*ia.mij2*ia.mk2/(ia.mij2-sqr(ia.Q-ia.mk))/ia.s)));
  double xm((sqr(ia.Q-ia.mk)-ia.mij2-sqrt(Lam(ia.Q2,ia.mij2,ia.mk2)))/ia.s);
  double xp((sqr(ia.Q-ia.mk)-ia.mij2+sqrt(Lam(ia.Q2,ia.mij2,ia.mk2)))/ia.s);
  return (-DiLog((a+x)/(a+xp))+DiLog(a/(a+xp))+DiLog((xp-x)/(xp-b))-DiLog(xp/(xp-b))
	  +DiLog((c+x)/(c+xp))-DiLog(c/(c+xp))+DiLog((xm-x)/(xm+a))-DiLog(xm/(xm+a))
	  -DiLog((b-x)/(b-xm))+DiLog(b/(b-xm))-DiLog((xm-x)/(xm+c))+DiLog(xm/(xm+c))
	  +DiLog((b-x)/(b+a))-DiLog(b/(b+a))-DiLog((c+x)/(c-a))+DiLog(c/(c-a))
	  +log(c+x)*log((a-c)*(xp-x)/(a+x)/(c+xp))-log(c)*log((a-c)*xp/a/(c+xp))
	  +log(b-x)*log((a+x)*(xm-b)/(a+b)/(xm-x))-log(b)*log(a*(xm-b)/(a+b)/xm)
	  -log((a+x)*(b-xp))*log(xp-x)+log(a*(b-xp))*log(xp)
	  +log(d)*log((a+x)*xp*xm/a/(xp-x)/(xm-x))+log((xm-x)/xm)*log((c+xm)/(a+xm))
	  +0.5*log((a+x)/a)*log(a*(a+x)*sqr(a+xp)))/ia.v;
}

NLO_Value METOOLS::FFQQ(const I_Args &a,const Dipole_Info *info)
{
  NLO_Value i(FFVS(a,info)+FFGQQ(a,info));
  i.m_f+=FFVNSQQ(a,info)-sqr(M_PI)/3.0;
  i.m_f+=7.0/2.0-sqr(M_PI)/6.0+2.0*FFAE(a,info)+FFACQQ(a,info);
  return i;
}

NLO_Value METOOLS::FFGQQ(const I_Args &a,const Dipole_Info *info)
{
  if (a.mij==0.0) return NLO_Value(0.0,3.0/2.0,3.0/2.0);
  double lmu2(log(info->Mu2()/a.s));
  return NLO_Value(0.0,1.0,0.5*log(a.mij2/info->Mu2())-2.0
		   -lmu2+3.0/2.0*(1.0+lmu2));
}

double METOOLS::FFVNSQQ(const I_Args &a,const Dipole_Info *info)
{
  if (a.mij==0.0 && a.mk==0.0) return 0.0;
  if (a.mk==0.0) {
    return 3.0/2.0*log(a.s/a.Q2)+sqr(M_PI)/6-DiLog(a.s/a.Q2)
      -2.0*log(a.s/a.Q2)-a.mij2/a.s*log(a.mij2/a.Q2);
  }
  if (a.mij==0.0) {
    return 3.0/2.0
      *(log(1.0-a.mk2/a.Q2)-2.0*log(1.0-a.mk/a.Q)-2.0*a.mk/(a.Q+a.mk))
       +sqr(M_PI)/6-DiLog(1.0-a.mk2/a.Q2);
  }
  return 3.0/2.0*log(a.s/a.Q2)
    +1.0/a.v*(log(a.r2)*log(1.0+a.r2)+2.0*DiLog(a.r2)
	    -DiLog(1.0-a.rj2)-DiLog(1.0-a.rk2)-sqr(M_PI)/6.0)
    +log(1.0-a.mk/a.Q)-2.0*log((sqr(a.Q-a.mk)-a.mij2)/a.Q2)
    -2.0*a.mij2/a.s*log(a.mij/(a.Q-a.mk))
    -a.mk/(a.Q-a.mk)+2.0*a.mk*(2.0*a.mk-a.Q)/a.s+sqr(M_PI)/2.0;
}

double METOOLS::FFACQQ(const I_Args &a,const Dipole_Info *info)
{
  if ((a.is&&info->Massive()) || info->AMax()==1.0) return 0.0;
  if (a.mij==0.0 && a.mk==0.0)
    return 3.0/2.0*(info->AMax()-1.0-log(info->AMax()));
  if (a.mij==0.0) {
    double yp((a.Q-a.mk)/(a.Q+a.mk));
    return 3.0/2.0*(yp*(info->AMax()-1.0)-log(info->AMax()));
  }
  if (a.mk==0.0) {
    double al(info->AMax()), muq2(a.mij2/a.Q2);
    return 0.5*(3.0*al-2.0
		-(3.0-muq2)/(1.0-muq2)*log(al+(1.0-al)*muq2)
		-al/(al+(1.0-al)*muq2))
      -2.0*log(al)+2.0*log(al+(1.0-al)*muq2)/(1.0-muq2);
  }
  double al(info->AMax()*(1.0-2.0*a.mk*(a.Q-a.mk)/a.s));
  return 3.0/2.0*(1.0+al)+a.Q/(a.Q-a.mk)-2.0*(2.0*a.Q2-2.0*a.mij2-a.mk*a.Q)/a.s
    +0.5*(1.0-al)*a.mij2/(a.mij2+al*a.s)-2.0*log(al*a.s/(sqr(a.Q-a.mk)-a.mij2))
    +0.5*(a.Q2+a.mij2-a.mk2)/a.s*log((a.mij2+al*a.s)/sqr(a.Q-a.mk));
}

NLO_Value METOOLS::FFGQ(const I_Args &a,const Dipole_Info *info,const double &m)
{
  NLO_Value i(FFGGQ(a,info,m));
  i.m_f+=FFVNSGQ(a,info,m);
  i.m_f+=(m?0.0:-10.0/9.0)+FFACGQ(a,info,m);
  return i;
}

NLO_Value METOOLS::FFGGQ(const I_Args &a,const Dipole_Info *info,const double &m)
{
  if (m==0.0) return NLO_Value(0.0,-2.0/3.0,-2.0/3.0);
  return NLO_Value(0.0,0.0,-2.0/3.0*log(m*m/a.Qa2));
}

double METOOLS::FFVNSGQ(const I_Args &a,const Dipole_Info *info,const double &m)
{
  if (m==0.0) {
    if (a.mk==0.0) return 0.0;
    return -2.0/3.0*(log(a.s/a.Q2)-2.0*log((a.Q-a.mk)/a.Q)-2.0*a.mk/(a.Q+a.mk))
      +(info->Kappa()-2.0/3.0)*a.mk2/a.s*2.0*log(2.0*a.mk/(a.Q+a.mk));
  }
  double v=2.0/3.0*log(m*m/a.Qa2);
  if (a.s<4.0*m*(m+a.mk)) return v;
  double r1(sqrt(1.0-4.0*m*m/sqr(a.Q-a.mk))), r2(sqrt(1.0-4.0*m*m/(a.Q2-a.mk2)));
  v+=4.0/3.0*(log((a.Q-a.mk)/a.Q)+a.mk*r1*r1*r1/(a.Q+a.mk)+log(0.5*(1.0+r1))
	      -r1*(1.0+r1*r1/3.0)-0.5*log(m*m/a.Q2));
  if (a.mk==0.0) return v;
  return v+(info->Kappa()-2.0/3.0)*a.mk2/a.s*2.0*
    (r2*r2*r2*log((r2-r1)/(r2+r1))-log((1.0-r1)/(1.0+r1))-8.0*r1*m*m/a.s);
}

double METOOLS::FFACGQ(const I_Args &ia,const Dipole_Info *info,const double &m)
{
  if ((ia.is&&info->Massive()) || info->AMax()==1.0) return 0.0;
  if (m==0.0) {
    if (ia.mk==0.0)
      return -2.0/3.0*(info->AMax()-1.0-log(info->AMax()));
    double muk(ia.mk/ia.Q), al(info->AMax()*(1.0-muk)/(1.0+muk));
    return 2.0/3.0*((1.0-muk-al*(1.0+muk))/(1.0+muk)
		    +log(al*(1.0+muk)/(1.0-muk)))
      +2.0*(info->Kappa()-2.0/3.0)*muk*muk/(1.0-muk*muk)
      *log((1.0-al)*(1.0+muk)/(2.0*muk));
  }
  if (ia.s<4.0*m*(m+ia.mk)) return 0.0;
  double yp(1.0-2.0*ia.mk*(ia.Q-ia.mk)/(ia.Q2-2.0*m*m-ia.mk2));
  double al(info->AMax()), muj2(m*m/ia.Q2), muk2(ia.mk2/ia.Q2);
  if (ia.mk==0.0) {
    double f(sqrt(1.0-4.0*muj2));
    double e(sqrt(sqr(al*(1.0-2.0*muj2))-4.0*muj2*muj2));
    return -2.0/3.0*
      (2.0*e/(2.0*(al-1.0)*muj2-al)+e+(2.0*muj2-1.0)
       *(-log((e+al*(2.0*muj2-1.0))/(f+2.0*muj2-1.0))
	 +2.0*atan(2.0*muj2/e)-2.0*atan(2.0*muj2/f))+f);
  }
  double c(-1.0+2.0*muj2+muk2), d(sqrt(sqr(al*c*yp)-4.0*muj2*muj2));
  double b(sqrt(sqr(c*yp)-4.0*muj2*muj2)), a(sqrt(1.0-muk2));
  return -((-8.0*a*muj2*muj2+2.0*a*c*(c+1.0)+4.0*a*muj2)
	   *log((al*c*c*yp-d*sqrt(c*c-4.0*muj2*muj2)-4.0*muj2*muj2)/
		(c*c*yp-b*sqrt(c*c-4.0*muj2*muj2)-4.0*muj2*muj2))
	   -2.0*a*(c*c+c-4.0*muj2*muj2+2.0*muj2)*log((1.0-al*yp)/(1.0-yp))
	   +(-3.0*c*c-2.0*c+4.0*c*muj2)
	   *sqrt(2.0*muj2-c)*log((al*c*yp+d)/(c*yp+b))
	   +2.0*sqrt(2.0*muj2-c)*(c*c-2.0*(c+1.0)*muj2+4.0*muj2*muj2)
	   *(atan(2.0*muj2/d)-atan(2.0*muj2/b))
	   +(c*sqrt(2.0*muj2-c)*c*c*yp*(al*al*b*yp-2.0*al*b-d*(yp-2.0))
	     +4.0*c*muj2*(b*(al*yp-1.0)-d*(yp-1.0))+4.0*muj2*muj2*(b-d))/(b*d)
	   )/(3.0*c*pow(2.0*muj2-c,3.0/2.0));
}

NLO_Value METOOLS::FFGG(const I_Args &a,const Dipole_Info *info)
{
  NLO_Value i(FFVS(a,info)+FFGGG(a,info));
  i.m_f+=FFVNSGG(a,info)-sqr(M_PI)/3.0;
  i.m_f+=67.0/18.0-sqr(M_PI)/6.0+2.0*FFAE(a,info)+FFACGG(a,info);
  return i;
}

NLO_Value METOOLS::FFGGG(const I_Args &a,const Dipole_Info *info)
{
  return NLO_Value(0.0,11.0/6.0,11.0/6.0);
}

double METOOLS::FFVNSGG(const I_Args &a,const Dipole_Info *info)
{
  if (a.mk==0.0) return 0.0;
  return 11.0/6.0*(log(a.s/a.Q2)-2.0*log((a.Q-a.mk)/a.Q)-2.0*a.mk/(a.Q+a.mk))
    +sqr(M_PI)/6.0-DiLog(a.s/a.Q2)
    -(info->Kappa()-2.0/3.0)*a.mk2/a.s*log(2.0*a.mk/(a.Q+a.mk));
}

double METOOLS::FFACGG(const I_Args &a,const Dipole_Info *info)
{
  if ((a.is&&info->Massive()) || info->AMax()==1.0) return 0.0;
  if (a.mk==0.0)
    return 11.0/6.0*(info->AMax()-1.0-log(info->AMax()));
  double muk(a.mk/a.Q), al(info->AMax()*(1.0-muk)/(1.0+muk));
  return -11.0/6.0*((1.0-muk-al*(1.0+muk))/(1.0+muk)
		    +log(al*(1.0+muk)/(1.0-muk)))
    -(info->Kappa()-2.0/3.0)*muk*muk/(1.0-muk*muk)
    *log((1.0-al)*(1.0+muk)/(2.0*muk));
}

namespace METOOLS {

  std::ostream &operator<<(std::ostream &str,const NLO_Value &v)
  {
    return str<<'{'<<v.m_e2<<','<<v.m_e1<<','<<v.m_f<<'}';
  }

}
