#include "Off_Shell_QCD.H"

#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Random.H"
#include "Flow.H"

#define NC 3

using namespace EXTRAXS;

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_gg_qqb>(const size_t nin,const size_t nout,
						   const ATOOLS::Flavour *flavours)
{
  if (flavours[2].IsQuark() && flavours[3]==flavours[2].Bar() && 
      flavours[0].IsGluon() && flavours[1].IsGluon()){ 
    return new Off_Shell_gg_qqb(nin,nout,flavours); 
  }
  return NULL;
}

Off_Shell_gg_qqb::Off_Shell_gg_qqb(const size_t nin,const size_t nout,
				   const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours) 
{
  m_alphas=(*MODEL::as)(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms()));
  m_nvector=m_nvector+2;
  CreateMomenta(m_nvector);
}

double Off_Shell_gg_qqb::operator()(double s,double t,double u) 
{
  ATOOLS::Vec4D *p=p_momenta;
  double S=p[4]*p[5], M2=p[2].Abs2();
  double a3=p[4]*p[2], a4=p[4]*p[3];
  double b3=p[5]*p[2], b4=p[5]*p[3];
  double k12=-p[0].Abs2(), k22=-p[1].Abs2();
  double delta=2.*S*(2.*(a4*b3-a3*b4)/S-k12*a3/(p[4]*p[1])+k22*b3/(p[5]*p[0])+p[2]*(p[0]-p[1]));
  double abelian=4.*(S*S/((t-M2)*(u-M2))-1./(k12*k22)*ATOOLS::sqr(S+2.*(a3*b4/(t-M2)+b3*a4/(u-M2))));
  double nonabelian=4.*S*S*(S/(s*(p[4]*p[1])*(p[5]*p[0]))-1./((t-M2)*(u-M2))
			    -(1/(t-M2)-1/(u-M2))/s*(a4/(p[4]*p[1])-b4/(p[5]*p[0])))
    +2./(k12*k22)*(S+4.*a3*b4/(t-M2)-delta/s)*(S+4.*a4*b3/(u-M2)+delta/s);
  return ATOOLS::sqr(4.*M_PI*m_alphas*(p[0]*p[5])*(p[1]*p[4])/(S*S))*
    (abelian/(2.*NC)+nonabelian*NC/(2.*(NC*NC-1.)));
}

bool Off_Shell_gg_qqb::SetColours(double s,double t,double u) 
{ 
  RestoreInOrder();
  int r=(int)p_flavours[2].IsAnti();
  double Mt=u/t;
  double Mu=t/u;
  m_scale[PHASIC::stp::fac]=(2.*s*t*u)/(s*s+t*t+u*u);
  p_colours[0][0]=ATOOLS::Flow::Counter();
  p_colours[0][1]=ATOOLS::Flow::Counter();
  if (Mt*(1-r)+Mu*r>(Mt+Mu)*ATOOLS::ran.Get()) {
    p_colours[2+r][0]=p_colours[0][0];
    p_colours[3-r][1]=p_colours[1][1]=ATOOLS::Flow::Counter();
    p_colours[1][0]=p_colours[0][1];
  }
  else {
    p_colours[2+r][0]=p_colours[1][0]=ATOOLS::Flow::Counter();
    p_colours[3-r][1]=p_colours[0][1];
    p_colours[1][1]=p_colours[0][0];
  }
  return true;
}

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_gg_gg>(const size_t nin,const size_t nout,
						  const ATOOLS::Flavour *flavours)
{
  if (flavours[2].IsGluon() && flavours[3].IsGluon() && 
      flavours[0].IsGluon() && flavours[1].IsGluon()){ 
    return new Off_Shell_gg_gg(nin,nout,flavours); 
  }
  return NULL;
}

Off_Shell_gg_gg::Off_Shell_gg_gg(const size_t nin,const size_t nout,
				   const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours) 
{
  m_alphas=(*MODEL::as)(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms()));
  m_nvector=m_nvector+2;
  CreateMomenta(m_nvector);
}

static double res, num;

double Off_Shell_gg_gg::operator()(double s,double t,double u) 
{
  const ATOOLS::Vec4D *p=p_momenta;
  const ATOOLS::Vec4D &p1=p[4], &p2=p[5];
  const ATOOLS::Vec4D &k1=p[0], &k2=p[1], &p3=p[2], &p4=p[3];
  ATOOLS::Vec4D k=k1-p3;
  double S=2.*p1*p2, S2=2./S;
  double a=p2*k*S2, b=p1*k*S2;
  double a1=p2*k1*S2, b2=p1*k2*S2;
  ATOOLS::Vec4D A1=
    (k1.Perp()*k.Perp())*(2.*k1-p3)
    +((2.*p3-k1)*k1.Perp())*k.Perp()
    -((k1+p3)*k.Perp())*k1.Perp();
  ATOOLS::Vec4D A2=
    (k2.Perp()*k.Perp())*(p4-2.*k2)
    +((k2+p4)*k.Perp())*k2.Perp()
    -((k2-2.*p4)*k2.Perp())*k.Perp();
  double M11=2*(A1*p3)*(A1*p4)/(p3*p4)-A1*A1;
  double M22=2*(A2*p4)*(A2*p3)/(p4*p3)-A2*A2;
  double result=M11*M22*
    ATOOLS::sqr(4./(S*S*k.PPerp2()*a1*a*b*b2));
  return ATOOLS::sqr(4.*M_PI*m_alphas)*(NC*NC)/(NC*NC-1)/(2.*4.)*result;
}

bool Off_Shell_gg_gg::SetColours(double s,double t,double u) 
{ 
  return true; 
}

