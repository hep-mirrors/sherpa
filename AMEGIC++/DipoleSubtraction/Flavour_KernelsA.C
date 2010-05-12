#include "AMEGIC++/DipoleSubtraction/Flavour_KernelsA.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;


Flavour_KernelsA::Flavour_KernelsA()
{
  m_cpldef = 0.0;
  p_cpl = NULL;
  
  Flavour hfl(kf_quark);
  m_nf = hfl.Size()/2;

  CSC.Init();
  m_g1 = 1.5*CSC.CF;
  m_g2 = 11./6.*CSC.CA-2./3.*CSC.TR*m_nf;

  m_alpha = 1.;
  m_loga = 0.;
}

void Flavour_KernelsA::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  if (cpls->find("Alpha_QCD")!=cpls->end()) p_cpl=cpls->find("Alpha_QCD")->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<"DipoleSplitting_Base:: alpha = "<<*p_cpl<<endl;
  m_cpldef = p_cpl->Default()/(2.*M_PI);
}

void Flavour_KernelsA::SetAlpha(double a)
{ m_alpha=a; m_loga=std::log(a); }

double Flavour_KernelsA::Kb1(int type,double x)
{
  switch(type) {
  case 1:
    return 2.*CSC.CF/(1.-x)*log((1.-x)/x);
  case 4:
    return 2.*CSC.CA/(1.-x)*log((1.-x)/x);
  }
  return 0.;
}

double Flavour_KernelsA::Kb2(int type)
{
  if (type==2||type==3) return 0.;
  switch(type) {
  case 1:
    return CSC.CF*(sqr(M_PI)-5.+2.*sqr(m_loga))-m_g1*(m_alpha-1-m_loga);
  case 4:
    return CSC.CA*(sqr(M_PI)-50./9.+2.*sqr(m_loga))-m_g2*(m_alpha-1-m_loga)+16./9.*CSC.TR*m_nf;
  }
  return 0.;
}

double Flavour_KernelsA::Kb3(int type,double x)
{
  double at=0.;
  if (m_alpha<1.&& (type==1||type==4)) {
    if (x<1.-m_alpha) at=log((1.-x)/(2.-x));
    at+=log(m_alpha*(2.-x)/(1.+m_alpha-x));
    at/=(1.-x);
  }
  switch(type) {
  case 1:
    return CSC.CF*(1.-x-(1.+x)*log(m_alpha*(1.-x)/x)+2.*at);
  case 2:
    return CSC.CF*((1.+sqr(1.-x))/x*log(m_alpha*(1.-x)/x)+x);
  case 3:
    return CSC.TR*((x*x+sqr(1.-x))*log(m_alpha*(1.-x)/x)+2.*x*(1.-x));
  case 4:
    return 2.*CSC.CA*(((1.-x)/x-1.+x*(1.-x))*log(m_alpha*(1.-x)/x)+at);
  }
  return 0.;
}

double Flavour_KernelsA::hfnc1(double x)
{
  double l1x=log(1.-x);
  return -0.5*sqr(l1x)+l1x*log(x)+DiLog(x);
}

double Flavour_KernelsA::Kb4(int type,double x)
{
  switch(type) {
  case 1:
    return 2.*CSC.CF*hfnc1(x);
  case 4:
    return 2.*CSC.CA*hfnc1(x);
  }
  return 0.;
}


double Flavour_KernelsA::t1(double x)
{
  if (x<1.-m_alpha) return 0.;
  return 1./(1.-x);
}

double Flavour_KernelsA::t2()
{
  return m_alpha;
}

double Flavour_KernelsA::t4(double x)
{
  if (x<1.-m_alpha) return 0.;
  return -log(1.-x)+m_loga;
}

double Flavour_KernelsA::ft(int type)
{
  switch(type) {
  case 1:
    return m_g1/CSC.CF;
  case 2:
    return m_g2/CSC.CA;
  }
  return 0.;
}


double Flavour_KernelsA::Kt1(int type,double x)
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

double Flavour_KernelsA::Kt2(int type)
{
  if (type==1||type==4) return -sqr(M_PI)/3.;
  return 0.;
}

double Flavour_KernelsA::Kt3(int type,double x)
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
    return CSC.CF/CSC.CA*((1.+sqr(1.-x))/x*(log(1.-x)-m_loga)+ax);
  case 3:
    ax*=(1.-2.*x*(1.-x));
    return CSC.TR/CSC.CF*((x*x+sqr(1.-x))*(log(1.-x)-m_loga)+ax);
  case 4:
    ax*=x/(1.-x)+(1.-x)/x+x*(1.-x);
    return 2.*((1.-x)/x-1.+x*(1.-x))*(log(1.-x)-m_loga)+at+2.*ax;
  }
  return 0.;
}

double Flavour_KernelsA::Kt4(int type,double x)
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

double Flavour_KernelsA::P1(int type,double x)
{
  switch(type) {
  case 1:
    return (1.+x*x)/(1.-x);
  case 4:
    return 2./(1.-x);
  }
  return 0.;
}

double Flavour_KernelsA::P2(int type)
{
  if (type==4) return m_g2/CSC.CA;
  return 0.;
}

double Flavour_KernelsA::P3(int type,double x)
{
  switch(type) {
  case 2:
    return CSC.CF/CSC.CA*(1.+sqr(1.-x))/x;
  case 3:
    return CSC.TR/CSC.CF*(x*x+sqr(1.-x));
  case 4:
    return 2.*((1.-x)/x-1.+x*(1.-x));
  }
  return 0.;
}

double Flavour_KernelsA::P4(int type,double x)
{
  switch(type) {
  case 1:
    return -x-0.5*x*x-2.*log(1.-x);
  case 4:
    return -2.*log(1.-x);
  }
  return 0.;
}


