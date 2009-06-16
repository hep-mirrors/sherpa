#include "AMEGIC++/DipoleSubtraction/Flavour_Kernels.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

Flavour_Kernels::Flavour_Kernels(MODEL::Model_Base *model) 
{
  m_cpl   = model->ScalarFunction(std::string("alpha_S"),sqr(rpa.gen.Ecms()));
  m_cpl /= 2.*M_PI;
  
  Flavour hfl(kf_quark);
  m_nf = hfl.Size()/2;

  CSC.Init();
  m_g1 = 1.5*CSC.CF;
  m_g2 = 11./6.*CSC.CA-2./3.*CSC.TR*m_nf;
}

double Flavour_Kernels::Kb1(int type,double x)
{
  switch(type) {
  case 1:
    return 2.*CSC.CF/(1.-x)*log((1.-x)/x);
  case 4:
    return 2.*CSC.CA/(1.-x)*log((1.-x)/x);
  }
  return 0.;
}

double Flavour_Kernels::Kb2(int type)
{
  switch(type) {
  case 1:
    return CSC.CF*(sqr(M_PI)-5.);
  case 4:
    return CSC.CA*(sqr(M_PI)-50./9.)+16./9.*CSC.TR*m_nf;
  }
  return 0.;
}

double Flavour_Kernels::Kb3(int type,double x)
{
  switch(type) {
  case 1:
    return CSC.CF*(1.-x-(1.+x)*log((1.-x)/x));
  case 2:
    return CSC.CF*((1.+sqr(1.-x))/x*log((1.-x)/x)+x);
  case 3:
    return CSC.TR*((x*x+sqr(1.-x))*log((1.-x)/x)+2.*x*(1.-x));
  case 4:
    return 2.*CSC.CA*((1.-x)/x-1.+x*(1.-x))*log((1.-x)/x);
  }
  return 0.;
}

double Flavour_Kernels::Kb4(int type,double x)
{
  if (type==2||type==3) return 0.;
  double l1x=log(1.-x);
  switch(type) {
  case 1:
    return 2.*CSC.CF*(-0.5*sqr(l1x)+l1x*log(x)+DiLog(x));
  case 4:
    return 2.*CSC.CA*(-0.5*sqr(l1x)+l1x*log(x)+DiLog(x));
  }
  return 0.;
}


double Flavour_Kernels::t1(double x)
{
  return 1./(1.-x);
}

double Flavour_Kernels::t2()
{
  return 1.;
}

double Flavour_Kernels::t4(double x)
{
  return -log(1.-x);
}

double Flavour_Kernels::ft(int type)
{
  switch(type) {
  case 1:
    return m_g1/CSC.CF;
  case 2:
    return m_g2/CSC.CA;
  }
  return 0.;
}


double Flavour_Kernels::Kt1(int type,double x)
{
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
  switch(type) {
  case 1:
    return -(1.+x)*log(1.-x);
  case 2:
    return CSC.CF/CSC.CA*(1.+sqr(1.-x))/x*log(1.-x);
  case 3:
    return CSC.TR/CSC.CF*(x*x+sqr(1.-x))*log(1.-x);
  case 4:
    return 2.*((1.-x)/x-1.+x*(1.-x))*log(1.-x);
  }
  return 0.;
}

double Flavour_Kernels::Kt4(int type,double x)
{
  if (type==2||type==3) return 0.;
  double l1x=log(1.-x);
  switch(type) {
  case 1:
    return -sqr(l1x);
  case 4:
    return -sqr(l1x);
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
  if (type==4) return m_g2/CSC.CA;
  return 0.;
}

double Flavour_Kernels::P3(int type,double x)
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


