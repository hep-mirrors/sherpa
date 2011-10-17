#include "PHASIC++/Process/Flavour_Kernels.H"

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

Flavour_Kernels::Flavour_Kernels()
{
  m_nf = Flavour(kf_quark).Size()/2;
  SetNC(3.0);
}

void Flavour_Kernels::SetNC(const double &nc)
{
  m_TR=0.5;
  m_CA=m_NC=nc;
  m_CF=m_TR*(m_NC-1.0/m_NC);
  m_g1=1.5*m_CF;
  m_g2=11.0/6.0*m_CA-2.0/3.0*m_TR*m_nf;
}

void Flavour_Kernels::SetAlpha(const double &a)
{
  m_alpha=a; m_loga=log(a);
}

double Flavour_Kernels::Kb1(int type,double x)
{
  switch(type) {
  case 1:
    return 2.*m_CF/(1.-x)*log((1.-x)/x);
  case 4:
    return 2.*m_CA/(1.-x)*log((1.-x)/x);
  }
  return 0.;
}

double Flavour_Kernels::Kb2(int type)
{
  if (type==2||type==3) return 0.;
  switch(type) {
  case 1:
    return m_CF*(sqr(M_PI)-5.+2.*sqr(m_loga))-m_g1*(m_alpha-1-m_loga);
  case 4:
    return m_CA*(sqr(M_PI)-50./9.+2.*sqr(m_loga))-m_g2*(m_alpha-1-m_loga)+16./9.*m_TR*m_nf;
  }
  return 0.;
}

double Flavour_Kernels::Kb3(int type,double x)
{
  double at=0.;
  if (m_alpha<1.&& (type==1||type==4)) {
    if (x<1.-m_alpha) at=log((1.-x)/(2.-x));
    at+=log(m_alpha*(2.-x)/(1.+m_alpha-x));
    at/=(1.-x);
  }
  switch(type) {
  case 1:
    return m_CF*(1.-x-(1.+x)*log(m_alpha*(1.-x)/x)+2.*at);
  case 2:
    return m_CF*((1.+sqr(1.-x))/x*log(m_alpha*(1.-x)/x)+x);
  case 3:
    return m_TR*((x*x+sqr(1.-x))*log(m_alpha*(1.-x)/x)+2.*x*(1.-x));
  case 4:
    return 2.*m_CA*(((1.-x)/x-1.+x*(1.-x))*log(m_alpha*(1.-x)/x)+at);
  }
  return 0.;
}

double Flavour_Kernels::hfnc1(double x)
{
  double l1x=log(1.-x);
  return -0.5*sqr(l1x)+l1x*log(x)+DiLog(x);
}

double Flavour_Kernels::Kb4(int type,double x)
{
  switch(type) {
  case 1:
    return 2.*m_CF*hfnc1(x);
  case 4:
    return 2.*m_CA*hfnc1(x);
  }
  return 0.;
}


double Flavour_Kernels::t1(double x)
{
  if (x<1.-m_alpha) return 0.;
  return 1./(1.-x);
}

double Flavour_Kernels::t2()
{
  return m_alpha;
}

double Flavour_Kernels::t4(double x)
{
  if (x<1.-m_alpha) return 0.;
  return -log(1.-x)+m_loga;
}

double Flavour_Kernels::ft(int type)
{
  switch(type) {
  case 1:
    return m_g1/m_CF;
  case 2:
    return m_g2/m_CA;
  }
  return 0.;
}


double Flavour_Kernels::Kt1(int type,double x)
{
  if (x<1.-m_alpha) return 0.;
  switch(type) {
  case 1:
    return 2./(1.-x)*log(1.-x);
  case 4:
    return 2./(1.-x)*log(1.-x);
  }
  return 0.;
}

double Flavour_Kernels::Kt2(int type)
{
  if (type==1||type==4) return -sqr(M_PI)/3.;
  return 0.;
}

double Flavour_Kernels::Kt3(int type,double x)
{
  double at=0.,ax=0.;
  if (m_alpha<1.) {
    if (type==1||type==4) {
      if (x>1.-m_alpha) at=-log(2.-x);
      at+=log(1.+m_alpha-x)-m_loga;
      at*=2./(1.-x);
    }
    if (m_alpha<(1.-x)) ax=log(m_alpha/(1.-x));
  }
  switch(type) {
  case 1:
    ax*=(1.+x*x)/(1.-x);
    return -(1.+x)*(log(1.-x)-m_loga)+at+ax;
  case 2:
    ax*=(1.+sqr(1.-x))/x;
    return m_CF/m_CA*((1.+sqr(1.-x))/x*(log(1.-x)-m_loga)+ax);
  case 3:
    ax*=(1.-2.*x*(1.-x));
    return m_TR/m_CF*((x*x+sqr(1.-x))*(log(1.-x)-m_loga)+ax);
  case 4:
    ax*=x/(1.-x)+(1.-x)/x+x*(1.-x);
    return 2.*((1.-x)/x-1.+x*(1.-x))*(log(1.-x)-m_loga)+at+2.*ax;
  }
  return 0.;
}

double Flavour_Kernels::Kt4(int type,double x)
{
  if (type==2||type==3) return 0.;
  if (x<1.-m_alpha) return 0.;
  double l1x=log(1.-x);
  switch(type) {
  case 1:
    return -sqr(l1x)+sqr(m_loga);
  case 4:
    return -sqr(l1x)+sqr(m_loga);
  }
  return 0.;
}

double Flavour_Kernels::P1(int type,double x)
{
  switch(type) {
  case 1:
    return (1.+x*x)/(1.-x);
  case 4:
    return 2./(1.-x);
  }
  return 0.;
}

double Flavour_Kernels::P2(int type)
{
  if (type==4) return m_g2/m_CA;
  return 0.;
}

double Flavour_Kernels::P3(int type,double x)
{
  switch(type) {
  case 2:
    return m_CF/m_CA*(1.+sqr(1.-x))/x;
  case 3:
    return m_TR/m_CF*(x*x+sqr(1.-x));
  case 4:
    return 2.*((1.-x)/x-1.+x*(1.-x));
  }
  return 0.;
}

double Flavour_Kernels::P4(int type,double x)
{
  switch(type) {
  case 1:
    return -x-0.5*x*x-2.*log(1.-x);
  case 4:
    return -2.*log(1.-x);
  }
  return 0.;
}


