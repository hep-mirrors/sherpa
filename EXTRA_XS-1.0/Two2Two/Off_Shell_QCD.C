#include "Off_Shell_QCD.H"

#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Running_AlphaS.H"

#define NC 3

using namespace EXTRAXS;

Off_Shell_gg_qqb::Off_Shell_gg_qqb(const size_t nin,const size_t nout,
				   const ATOOLS::Flavour *flavours,const int scalescheme,
				   const int kfactorscheme,const double scalefactor):
  Single_XS(nin,nout,flavours) 
{
  m_scalescheme=scalescheme;
  m_kfactorscheme=kfactorscheme;
  m_scalefactor=scalefactor;
  m_alphas=(*MODEL::as)(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms()));
  m_nvector=m_nvector+2;
  delete [] p_momenta;
  p_momenta = new ATOOLS::Vec4D[m_nvector];
}

double Off_Shell_gg_qqb::operator()(double s,double t,double u) 
{
  ATOOLS::Vec4D *p=p_momenta;
  double S=p[4]*p[5], M2=p[2].Abs2();
  double a3=p[4]*p[2], a4=p[4]*p[3];
  double b3=p[5]*p[2], b4=p[5]*p[3];
  double k12=-p[0].Abs2(), k22=-p[1].Abs2();
  double delta=2.*S*(2.*(a3*b4-a4*b3)/S-k12*a4/(p[4]*p[1])+k22*b4/(p[5]*p[0])+p[3]*(p[0]-p[1]));
  double abelian=4.*(S*S/((t-M2)*(u-M2))-1./(k12*k22)*ATOOLS::sqr(S+2.*(a4*b3/(t-M2)+b4*a3/(u-M2))));
  double nonabelian=4.*S*S*(S/(s*(p[4]*p[1])*(p[5]*p[0]))-1./((t-M2)*(u-M2))
			    -(1/(t-M2)-1/(u-M2))/s*(a3/(p[4]*p[1])-b3/(p[5]*p[0])))
    +2./(k12*k22)*(S+4.*a4*b3/(t-M2)-delta/s)*(S+4.*a3*b4/(u-M2)+delta/s);
  return ATOOLS::sqr(4.*M_PI*m_alphas*(p[0]*p[5])*(p[1]*p[4])/(S*S))*
    (abelian/(2.*NC)+nonabelian*NC/(2.*(NC*NC-1.)));
}

bool Off_Shell_gg_qqb::SetColours(double s,double t,double u) 
{ 
  return true; 
}

