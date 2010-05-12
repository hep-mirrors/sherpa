#include "AMEGIC++/DipoleSubtraction/Massive_Kernels.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

//for equations see hep-ph/0201036

Massive_Kernels::Massive_Kernels() 
{
  m_cpldef = 0.0;
  p_cpl = NULL;
  
  Flavour hfl(kf_quark);
  m_nf = hfl.Size()/2; //number of massless flavours

  CSC.Init();
  m_g1 = 1.5*CSC.CF;
  m_g2 = 11./6.*CSC.CA-2./3.*CSC.TR*m_nf;
  m_K1 = (3.5-sqr(M_PI)/6.)*CSC.CF;
  m_K2 = (67./18.-sqr(M_PI)/6.)*CSC.CA-10./9.*CSC.TR*m_nf;

  m_alpha = 1.;
  m_loga = 0.;

  int helpi,nfgs=m_nf;
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa.GetPath());
  reader.SetInputFile(rpa.gen.Variable("ME_DATA_FILE"));
  if (reader.ReadFromFile(helpi,"DIPOLE_NF_GSPLIT")) {
    nfgs = helpi;
    msg_Tracking()<<"Set number of flavours from gluon splitting="<<nfgs<<"."<<std::endl;
  }
  for (int i=1;i<=nfgs;i++) {
    Flavour flav((kf_code)(i));
    if (flav.IsMassive()) m_massflav.push_back(flav.Mass());
  }
  m_nmf=m_massflav.size();

  double helpd;
  m_kappa=2./3.;
  if (reader.ReadFromFile(helpd,"DIPOLE_KAPPA")) {
    m_kappa = helpd;
    msg_Tracking()<<"Set massive dipole kappa="<<m_kappa<<"."<<std::endl;
  }
}

void Massive_Kernels::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  if (cpls->find("Alpha_QCD")!=cpls->end()) p_cpl=cpls->find("Alpha_QCD")->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<"DipoleSplitting_Base:: alpha = "<<*p_cpl<<endl;
  m_cpldef = p_cpl->Default()/(2.*M_PI);
}

void Massive_Kernels::SetAlpha(double a)
{ 
  if (a!=1.) msg_Error()<<"Warning: DIPOLE_ALPHA!=1 not implemented for the massive case!"<<endl; 
  //m_alpha=a; m_loga=std::log(a); 
}

double Massive_Kernels::Lambda(double x,double y, double z)
{
  return sqr(x)+sqr(y)+sqr(z)-2.*(x*y+x*z+y*z);
}

void Massive_Kernels::CalcVS(double s,double mj,double mk)
// V^(S) as in (6.20)
{
  p_VS[0]=p_VS[1]=p_VS[2]=0.;
  if (mj>0.&&mk>0.) {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    double vjk=sqrt(Lambda(Q2,mj2,mk2))/s;
    double lrhoj=log(sqrt(((1.-vjk)*s+2.*mj2)/((1.+vjk)*s+2.*mj2)));
    double lrhok=log(sqrt(((1.-vjk)*s+2.*mk2)/((1.+vjk)*s+2.*mk2)));
    p_VS[1]=(lrhoj+lrhok)/vjk;
    p_VS[0]=(-sqr(lrhoj)-sqr(lrhok)-sqr(M_PI)/6.+(lrhoj+lrhok)*log(Q2/s))/vjk;
    return;
  }
  double m=mj+mk;
  if (m>0.) {
    double m2=sqr(m);
    double Q2=s+m2;
    double lms=log(m2/s);
    p_VS[2]=.5;
    p_VS[1]=.5*lms;
    p_VS[0]=-.25*sqr(lms)-sqr(M_PI)/12.-.5*log(s/Q2)*(lms+log(m2/Q2));
    return;
  }
  p_VS[2]=1.;
}

void Massive_Kernels::CalcVNS(int t,double s,double mj,double mk,bool ini)
{
  p_VNS[0]=p_VNS[1]=p_VNS[2]=0.;
  if (t==1) CalcVNSq(s,mj,mk);
  if (t==2) CalcVNSg(s,mk,ini);
}

void Massive_Kernels::CalcVNSq(double s,double mj,double mk)
//V^(NS)_q as defined in (6.21)-(6.23)
{
  if (mj==0.&&mk==0.) return;
  if (mj==0.) {
    double Q2=s+sqr(mj)+sqr(mk);
    double Q=sqrt(Q2);
    p_VNS[0] = m_g1/CSC.CF*(log(s/Q2)-2.*log(1.-mk/Q)-2.*mk/(Q+mk))
      +sqr(M_PI)/6.-DiLog(s/Q2);
    return;
  }
  if (mk==0.) {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    p_VNS[0] = (m_g1/CSC.CF-2.)*log(s/Q2)+sqr(M_PI)/6.
      -DiLog(s/Q2)-mj2/s*log(mj2/Q2);
    return;
  }
  {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    double Q=sqrt(Q2);
    double vjk=sqrt(Lambda(Q2,mj2,mk2))/s;
    double rhoj2=((1.-vjk)*s+2.*mj2)/((1.+vjk)*s+2.*mj2);
    double rhok2=((1.-vjk)*s+2.*mk2)/((1.+vjk)*s+2.*mk2);
    double rho2=rhoj2*rhok2;
    p_VNS[0] = m_g1/CSC.CF*log(s/Q2)
      +(log(rho2)*log(1.+rho2)+2.*DiLog(rho2)-DiLog(1-rhoj2)-DiLog(1-rhok2)-sqr(M_PI)/6.)/vjk
      +log(1.-mk/Q)-2.*log((sqr(Q-mk)-mj2)/Q2)-2.*mj2/s*log(mj/(Q-mk))
      -mk/(Q-mk)+2*mk*(2.*mk-Q)/s+0.5*sqr(M_PI);
  }
}

void Massive_Kernels::CalcVNSg(double s,double mk,bool ini)
//V^(NS_q) as defined in Eqs.(6.24) and (6.26); Q_aux-terms canceled with Gamma_g
{
  size_t nfjk=0;
  if (!ini) 
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) nfjk++;
  if (mk==0.) {
    for (size_t i=0;i<nfjk;i++) {
      double rho1=sqrt(1.-4.*sqr(m_massflav[i])/s);
      p_VS[0]+=log(0.5+0.5*rho1)-rho1*(1.+sqr(rho1)/3.)-0.5*log(sqr(m_massflav[i])/s);
    }
    p_VNS[0]*=4./3.*CSC.TR/CSC.CA;
    return;
  }
  {
    bool simplev=ini||IsEqual(m_kappa,2./3.);
    double Q2=s+sqr(mk);
    double Q=sqrt(Q2);
    p_VS[0]=m_g2/CSC.CA*(log(s/Q2)-2.*log(1.-mk/Q)-2.*mk/(Q+mk))+sqr(M_PI)/6.-DiLog(s/Q2);
    if (!simplev) p_VS[0]+=(m_kappa-2./3.)*sqr(mk)/s*((2.*CSC.TR*m_nf/CSC.CA-1.)*log(2.*mk/(Q+mk)));
    double nfc=0.;
    for (size_t i=0;i<nfjk;i++) {
      double rho1=sqrt(1.-4.*sqr(m_massflav[i])/sqr(Q-mk));
      nfc+=4./3.*(log(1.-mk/Q)+mk*rho1*rho1*rho1/(Q+mk)+log(0.5+0.5*rho1)-rho1*(1.+sqr(rho1)/3.)
		  -0.5*log(sqr(m_massflav[i])/Q2));
      if (!simplev) {
	double rho2=sqrt(1.-4.*sqr(m_massflav[i])/s);	
	nfc+=(m_kappa-2./3.)*2.*(rho2*rho2*rho2*log((rho2-rho1)/(rho2+rho1))-log((1.-rho1)/(1.+rho1))
		 -8.*rho1*sqr(m_massflav[i])/s);
      }
    }
    p_VNS[0]+=CSC.TR/CSC.CA*nfc;
  }
}

void Massive_Kernels::CalcGamma(int t,double mu2,double m)
{
  p_Gamma[2]=0.;
  if (t==2) {
    p_Gamma[0]=0.;
    p_Gamma[1]=m_g2;
    return;
  }
  if (t==1) {
    if (m==0.) {
      p_Gamma[0]=0.;
      p_Gamma[1]=m_g1;
      return;
    }
    p_Gamma[1]=CSC.CF;
    p_Gamma[0]=CSC.CF*(0.5*log(sqr(m)/mu2)-2.);
  }
}
    
void Massive_Kernels::Calculate(int t,double mu2,double s,double mj,double mk, bool ini)
{
  CalcVS(s,mj,mk);
  CalcVNS(t,s,mj,mk,ini);
  CalcGamma(t,mu2,mj);
  double lmus=log(mu2/s);
  p_Gamma[0]-=lmus*p_Gamma[1];
  if (t==1) {
    p_Gamma[0]+=m_g1*(1.+lmus)+m_K1;
    p_Gamma[0]/=CSC.CF;
    p_Gamma[1]/=CSC.CF;
  }
  if (t==2) {
    p_Gamma[0]+=m_g2*(1.+lmus)+m_K2;
    p_Gamma[0]/=CSC.CA;
    p_Gamma[1]/=CSC.CA;
  }
}

double Massive_Kernels::I_Fin()
{
  return p_VS[0]+p_VNS[0]+p_Gamma[0]-sqr(M_PI)/3.;
}

double Massive_Kernels::I_E1()
{
  return p_VS[1]+p_Gamma[1];
}

double Massive_Kernels::I_E2()
{
  return p_VS[2];
}

//muq2=m_j^2/s_ja

double Massive_Kernels::t1(int type,int spin,double muq2,double x)
// g(x)
{
  if (type==2||type==3) return 0.;
  x=1.-x;
  switch (spin) {
  case 1: 
    return 2./x*(1.+log(x+muq2)-log(x))-0.5*x/sqr(x+muq2); 
  case 2:
    return m_g2/CSC.CA/x; 
  }
  return 0.;
}

double Massive_Kernels::t2(int type,int spin,double muq2)
// h; in case of gluon muq2 must be s_ja!!
{
  if (type==2||type==3) return 0.;
  switch (spin) {
  case 1: {
    double mx=muq2/(1.+muq2);
    return m_g1/CSC.CF*(1.-2.*(log(sqrt(1.+muq2)-sqrt(muq2))+1./(sqrt(1./mx)+1.)))
      -muq2*log(mx)-0.5*mx;
  }
  case 2: {
    double mgs=0.;
    for (size_t i=0;i<m_nmf;i++) {
      double xp=1.-sqr(2.*m_massflav[i])/muq2;
      if (xp>0.) mgs+=pow(xp,1.5);
    }
    return (m_g2-CSC.TR*2./3.*mgs)/CSC.CA;
  }
  }
  return 0.;
}

double Massive_Kernels::t3(int type,int spin,double muq2,double x)
// k(x)
{
  if (spin==2) return 0.;
  double mx=log((1.-x)/(1.-x+muq2));
  switch (type) {
  case 1:
    return (1.+x)*mx;
  case 2:
    return -CSC.CF/CSC.CA*((1.+sqr(1.-x))*mx-2.*muq2*log((1.-x)/muq2+1.))/x;
  case 3:
    return -CSC.TR/CSC.CF*(x*x+sqr(1-x))*mx;
  case 4:
    return -2.*((1./x-2.+x*(1.-x))*mx-muq2/x*log((1.-x)/muq2+1.));
  }
  return 0.;
}

double Massive_Kernels::t4(int type,int spin,double muq2,double x)
// G(x)
{
  if (type==2||type==3) return 0.;
  double y=1.-x;
  double lny=log(x);
  switch (spin) {
  case 1: 
    return sqr(lny)+2.*(DiLog(-y/muq2)-DiLog(-1./muq2)-log(muq2)*lny)
      +0.5*(muq2*x/((1.+muq2)*(y+muq2))-log((1.+muq2)/(y+muq2)))-2.*lny; 
  case 2:
    return -m_g2/CSC.CA*lny; 
  }
  return 0.;
}

double Massive_Kernels::t5(int type,double x,double xp)
// g^{xp}(x)
{
  if (type==2||type==3) return 0.;
  if (x>xp) return 0.;
  x=1.-x;
  xp=1.-xp;
  return -2./3.*CSC.TR/CSC.CA*(x+0.5*xp)/sqr(x)*sqrt(1.-xp/x);
}

double Massive_Kernels::t6(int type,double xp)
// h^{xp}
{
  if (type==2||type==3) return 0.;
  double sxp=sqrt(xp);
  return -2./3.*CSC.TR/CSC.CA*(log((1.-sxp)/(1.+sxp))+sxp/3.*(6.-xp));
}

double Massive_Kernels::t7(int type,double x,double xp)
// G^{xp}(x)
{
  if (type==2||type==3) return 0.;
  if (x>xp) x=xp;
  return -2./3.*CSC.TR/CSC.CA*((sqrt(1.-(1.-xp)/(1.-x))*(5.+(1.-xp)/(1.-x))-sqrt(xp)*(6.-xp))/3.
			       -log((1.+xp)/2.-x+sqrt((1.-x)*(xp-x)))+log((1.+xp)/2.+sqrt(xp)));
}
