#include "AMEGIC++/DipoleSubtraction/Massive_KernelsA.H"
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

Massive_KernelsA::Massive_KernelsA() 
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
  m_aterm = 0.;

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

void Massive_KernelsA::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  if (cpls->find("Alpha_QCD")!=cpls->end()) p_cpl=cpls->find("Alpha_QCD")->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<"DipoleSplitting_Base:: alpha = "<<*p_cpl<<endl;
  m_cpldef = p_cpl->Default()/(2.*M_PI);
}

void Massive_KernelsA::SetAlpha(double a)
{ 
  m_alpha=a; m_loga=std::log(a); 
}

double Massive_KernelsA::Lambda(double x,double y, double z)
{
  return sqr(x)+sqr(y)+sqr(z)-2.*(x*y+x*z+y*z);
}

void Massive_KernelsA::CalcVS(double s,double mj,double mk)
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

void Massive_KernelsA::CalcVNS(int t,double s,double mj,double mk,bool ini)
{
  p_VNS[0]=p_VNS[1]=p_VNS[2]=0.;
  if (t==1) CalcVNSq(s,mj,mk);
  if (t==2) CalcVNSg(s,mk,ini);
}

void Massive_KernelsA::CalcVNSq(double s,double mj,double mk)
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

void Massive_KernelsA::CalcVNSg(double s,double mk,bool ini)
//V^(NS_q) as defined in Eqs.(6.24) and (6.26); Q_aux-terms canceled with Gamma_g
{
  size_t nfjk=0;
  if (!ini) 
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) nfjk++;
  if (mk==0.) {
    for (size_t i=0;i<nfjk;i++) {
      double rho1=sqrt(1.-4.*sqr(m_massflav[i])/s);
      p_VNS[0]+=log(0.5+0.5*rho1)-rho1*(1.+sqr(rho1)/3.)-0.5*log(sqr(m_massflav[i])/s);
    }
    p_VNS[0]*=4./3.*CSC.TR/CSC.CA;
    return;
  }
  else {
    bool simplev=ini||IsEqual(m_kappa,2./3.);
    double Q2=s+sqr(mk);
    double Q=sqrt(Q2);
    p_VNS[0]=m_g2/CSC.CA*(log(s/Q2)-2.*log(1.-mk/Q)-2.*mk/(Q+mk))+sqr(M_PI)/6.-DiLog(s/Q2);
    if (!simplev) p_VNS[0]+=(m_kappa-2./3.)*sqr(mk)/s*((2.*CSC.TR*m_nf/CSC.CA-1.)*log(2.*mk/(Q+mk)));
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

void Massive_KernelsA::CalcGamma(int t,double mu2,double m)
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

void Massive_KernelsA::Calculate(int t,double mu2,double s,double mj,double mk, bool ini, bool ini2)
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
  m_aterm=0.;
  if (!ini && !ini2 && m_alpha!=1.) CalcAterms( t,mu2,s,mj,mk);
}

double Massive_KernelsA::I_Fin()
{
  return p_VS[0]+p_VNS[0]+p_Gamma[0] - sqr(M_PI)/3. +m_aterm;
}

double Massive_KernelsA::I_E1()
{
  return p_VS[1]+p_Gamma[1];
}

double Massive_KernelsA::I_E2()
{
  return p_VS[2];
}

//muq2=m_j^2/s_ja

double Massive_KernelsA::t1(int type,int spin,double muq2,double x)
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
double Massive_KernelsA::t1p(int type,int spin,double muq2,double x)
// g(x) alpha terms
{
  if (type==2||type==3) return 0.;
  double aterm(0.);
  if (m_alpha<1.) aterm = -at1( type, spin, muq2, x);
  return aterm;
}

double Massive_KernelsA::t2(int type,int spin,double muq2)
// h; in case of gluon muq2 must be s_ja!!
{
  if (type==2||type==3) return 0.;
  double aterm(0.);
  if (m_alpha<1.) aterm = -at2(type, spin, muq2);
  switch (spin) {
  case 1: {
    double mx=muq2/(1.+muq2);
    if (IsZero(muq2)) return m_g1/CSC.CF+ aterm;
    return m_g1/CSC.CF*(1.-2.*(log(sqrt(1.+muq2)-sqrt(muq2))+1./(sqrt(1./mx)+1.)))
      -muq2*log(mx)-0.5*mx + aterm;
  }
  case 2: {
    double mgs=0.;
    for (size_t i=0;i<m_nmf;i++) {
      double xp=1.-sqr(2.*m_massflav[i])/muq2;
      if (xp>0.) mgs+=pow(xp,1.5);
    }
    return (m_g2-CSC.TR*2./3.*mgs)/CSC.CA + aterm;
  }
  }
  return aterm;
}

double Massive_KernelsA::t3(int type,int spin,double muq2,double x)
// k(x)
{
  double aterm(0.);
  if (m_alpha<1.) aterm = -at3( type, spin, muq2, x);
  if (IsZero(muq2)) return aterm;
  if (spin==2) return aterm;
  double mx=log((1.-x)/(1.-x+muq2));
  switch (type) {
  case 1:
    return (1.+x)*mx + aterm;
  case 2:
    return -CSC.CF/CSC.CA*((1.+sqr(1.-x))*mx-2.*muq2*log((1.-x)/muq2+1.))/x + aterm;
  case 3:
    return -CSC.TR/CSC.CF*(x*x+sqr(1.-x))*mx+ aterm;
  case 4:
    return -2.*((1./x-2.+x*(1.-x))*mx-muq2/x*log((1.-x)/muq2+1.)) + aterm;
  }
  return aterm;
}

double Massive_KernelsA::t4(int type,int spin,double muq2,double x)
// G(x)
{
  if (type==2||type==3) return 0.;
  double aterm(0.);
  if (m_alpha<1.) aterm = -at4( type, spin, muq2, x);
  double y=1.-x;
  double lny=log(y);
  if (IsZero(muq2)) {
  switch (spin) {
  case 1: 
    return -m_g1/CSC.CF*log(1.-x) + aterm; 
  case 2:
    return -m_g2/CSC.CA*log(1.-x) + aterm; 
  }
  }
  switch (spin) {
  case 1: 
    return sqr(lny)+2.*(DiLog(-y/muq2)-DiLog(-1./muq2)-log(muq2)*lny)
      +0.5*(muq2*x/((1.+muq2)*(y+muq2))-log((1.+muq2)/(y+muq2)))-2.*lny + aterm; 
  case 2:
    return -m_g2/CSC.CA*lny + aterm; 
  }
  return aterm;
}

double Massive_KernelsA::t5(int type,double x,double xp)
// g^{xp}(x)
{
  if (type==2||type==3) return 0.;
  if (x>xp) return 0.;
  x=1.-x;
  xp=1.-xp;
  return -2./3.*CSC.TR/CSC.CA*(x+0.5*xp)/sqr(x)*sqrt(1.-xp/x);
}

double Massive_KernelsA::t6(int type,double xp)
// h^{xp}
{
  if (type==2||type==3) return 0.;
  double sxp=sqrt(xp);
  return -2./3.*CSC.TR/CSC.CA*(log((1.-sxp)/(1.+sxp))+sxp/3.*(6.-xp));
}

double Massive_KernelsA::t7(int type,double x,double xp)
// G^{xp}(x)
{
  if (type==2||type==3) return 0.;
  if (x>xp) x=xp;
  return -2./3.*CSC.TR/CSC.CA*((sqrt(1.-(1.-xp)/(1.-x))*(5.+(1.-xp)/(1.-x))-sqrt(xp)*(6.-xp))/3.
			       -log((1.+xp)/2.-x+sqrt((1.-x)*(xp-x)))+log((1.+xp)/2.+sqrt(xp)));
}


/// alpha terms.

void Massive_KernelsA::CalcAterms(int t,double mu2,double s,double mj,double mk)
{
  m_aterm =0.;
  double Q2 = s+sqr(mj)+sqr(mk);
  double Q = sqrt(Q2);
  double muj2 = mj*mj/Q2;
  double muk2 = mk*mk/Q2;
  double muk = sqrt(muk2);
  double loga = log(m_alpha);
  if (t==1) {
    if (mj==0.) {
      if (mk == 0.) {
        m_aterm+= -sqr(log(m_alpha)) -3./2.*(log(m_alpha)+1.-m_alpha);
        return;
      }
      else{
        double yp = (1.-muk)/(1.+muk);
        double xp = yp*(1.-m_alpha) + sqrt((1.-m_alpha)*(1.-m_alpha*yp*yp));
        m_aterm += sqr(log((1.-yp*yp+2.*xp*yp)/(1.+yp-xp)/(1.-yp+xp))) - 2.*sqr(log((1.+yp-xp)/(1.+yp)))
          +4.*(log((1.+yp)/2.)*log((1.-yp+xp)/(1.-yp)) + log((1.+yp)/(2.*yp))*log((1.-yp*yp+2.*xp*yp)/(1.-yp*yp))
          +DiLog((1.-yp)/(1.+yp)) - DiLog((1.-yp*yp+2.*xp*yp)/sqr(1.+yp)) +DiLog((1.-yp+xp)/2.)
          - DiLog((1.-yp)/2.));
        m_aterm += -1.5*(log(m_alpha) +yp*(1.-m_alpha));
        return;
      }
    }
    else{
      if (mk == 0.) {
        double muq2 = muj2;
        double lmq2 = log(muj2);
        m_aterm+= 2.*(-log(m_alpha)*lmq2 - DiLog((muq2 -1.)/muq2) + DiLog((m_alpha*(muq2 -1.))/muq2));
        m_aterm+= 3./2.*(m_alpha-1.) - (3.-muq2)/(2.-2.*muq2)*log(m_alpha + (1.-m_alpha)*muq2)
                   - 0.5*m_alpha/(m_alpha + (1.-m_alpha)*muq2) +0.5 - 2.*log(m_alpha) 
                   + 2./(1.-muq2)*log(m_alpha+(1.-m_alpha)*muq2);
        return;
      }
      else{
        double mjmk1 = 1.-muj2 - muk2;
        double yp = 1. - 2.*muk*(1.-muk)/mjmk1;
        double ap = m_alpha*yp;
        double a = 2.*muk/(1.- muj2-muk2);
        double b = 2.*(1.-muk)/(1.- muj2-muk2);
        double c = 2.*muk*(1.-muk)/(1.- muj2-muk2);
        double d = 0.5*(1.- muj2-muk2);
        double x = yp - ap +sqrt((yp-ap)*(1./yp - ap +4.*muj2*muk2/((muj2 - sqr(1.-muk))*(1.- muj2-muk2))));
        double xp = ((1.-muk)*(1.-muk) - muj2 +sqrt(Lambda(1.,muj2,muk2)))/(1.- muj2-muk2);
        double xm = ((1.-muk)*(1.-muk) - muj2 -sqrt(Lambda(1.,muj2,muk2)))/(1.- muj2-muk2);
        double vjk = sqrt(Lambda(1.,muj2,muk2))/(1.-muj2-muk2);
        m_aterm +=  1.5*(1.+ap) +1./(1.-muk) - 2.*(2. - 2.*muj2 - muk)/mjmk1 +(1.-ap)*muj2/(2.*(muj2+ap*mjmk1))
        - 2.*log(ap*mjmk1/(sqr(1.-muk) - muj2)) + (1.+muj2-muk2)/(2.*mjmk1)*log((muj2+ap*mjmk1)/sqr(1.-muk)); ///eq A20
        m_aterm += 2.*(-DiLog((a+x)/(a+xp)) +DiLog(a/(a+xp)) +DiLog((xp-x)/(xp-b)) - DiLog(xp/(xp-b))
          +DiLog((c+x)/(c+xp)) - DiLog(c/(c+xp)) + DiLog((xm-x)/(xm+a)) - DiLog(xm/(xm+a))
          -DiLog((b-x)/(b-xm)) + DiLog(b/(b-xm)) - DiLog((xm-x)/(xm+c)) + DiLog(xm/(xm+c))
          +DiLog((b-x)/(b+a)) - DiLog(b/(b+a)) - DiLog((c+x)/(c-a)) + DiLog(c/(c-a))
          +log(c+x)*log((a-c)*(xp-x)/((a+x)*(c+xp))) - log(c)*log((a-c)*xp/(a*(c+xp)))
          +log(b-x)*log((a+x)*(xm-b)/((a+b)*(xm-x))) - log(b)*log(a*(xm-b)/(xm*(a+b)))
          -log((a+x)*(b-xp))*log(xp-x) + log(a*(b-xp))*log(xp)
          +log(d)*log((a+x)*xp*xm/(a*(xp-x)*(xm-x))) + log((xm-x)/xm)*log((c+xm)/(a+xm))
          +0.5*log((a+x)/a)*log(a*(a+x)*(a+xp)*(a+xp)))/vjk;
        return;
      }
    }
  }
  else{
    if (IsZero(mk)) {
      for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) {
     double muj2 = m_massflav[i]*m_massflav[i]/Q2;
     double muj4 = 4.*muj2*muj2;
     double a = sqrt(sqr(m_alpha*(1.-2.*muj2)) -muj4);
     double b = sqrt(1.-4.*muj2);
     m_aterm -= CSC.TR*2./3.*(2.*a/(2.*(m_alpha-1.)-m_alpha) +a + (2.*muj2 -1.)*
        (-log(-2.*(a+m_alpha*(2.*muj2-1.))) + 2.*atan(2.*muj2/a)
        +log(-2.*(2.*muj2 +b -1.)) - 2.*atan(2.*muj2/b)) + b); ///ref 0 eq A9
      }
      m_aterm+=-sqr(loga)+11./6.*(m_alpha-1.-loga) -m_nf*CSC.TR/CSC.CA*2./3.*(m_alpha-1.-loga);
    }
    else{
      double yp = 1. - 2.*muk*(1.-muk)/(1.-muk2);
      double ap = m_alpha*yp;
      double kappa = m_kappa;
      double yl = 1. + m_alpha*sqr(1. -muk) - muk2 -
                       sqrt(sqr(1.-muk)*(sqr(m_alpha*(1.-muk)) + sqr(1.+muk) -2.*(m_alpha*(1.+muk2))));
      m_aterm += -((11.*sqr(-2.+2.*muk+yl))/
                 ((-1.+muk2)*(-2.+yl))
                 -44.*log(2.-2.*muk) - 22.*log(muk) + 24.*log(2./(1. + muk))*
                 (log(2./(1. + muk)) + 2.*log(1. + muk)) +
                 (2.*((-11. + 15.*muk2)*log(2.*muk) + 
                 4.*muk2*(-log(-8.*(-1.+muk)*muk2)+
                 log(sqr(-2. + yl)+4.*muk2*(-1. + yl))) + 
                 (11. - 15.*muk2)*log(2. - yl)))/(-1. + muk2) + 
                 22.*log(2. - 2.*muk2 - yl) + 
                 22.*log(yl)-12.*(4.*log(1.-yl/2.)*log(-(yl/(-1.+muk2)))-
                 sqr(log(-(yl/(-1.+muk2)))) + 
                 sqr(log((-2.*(-2. + 2.*muk2 + yl))/
                 ((-1. + muk2)*(-2. + yl)))) + 
                 2.*log(-(yl/(-1.+muk2)))*(
                 log((-2.*(-2.+2.*muk2+yl))/((-1.+muk2)*(-2.+yl))) - 
                 2.*log(1. + yl/(-2. + 2.*muk2)))) + 48.*DiLog(1. - muk) -
                 48.*DiLog(1./(1. + muk)) - 
                 48.*DiLog(yl/2.) + 48.*DiLog(yl/(2. - 2.*muk2)))/12.;
      m_aterm += CSC.TR/CSC.CA*m_nf*(2./3.*((1.-muk-ap*(1.+muk))/(1.+muk) +log(ap*(1.+muk)/(1.-muk)))
          + (kappa - 2./3.)*2.*muk2/(1.-muk2)*log((1.-ap)*(1.+muk)/(2.*muk)));/// ref 42 eq 21. 
      for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) {
        double muj2 = m_massflav[i]*m_massflav[i]/Q2;
        double muj4 = 4.*muj2*muj2;
        double a = sqrt(1.-muk2);
        double c = -1. + 2.*muj2 +muk2;
        double c2 = c*c;
        double b = sqrt(c2*yp*yp - muj4);
        double d = sqrt(sqr(m_alpha*c*yp) - muj4);
        double e = sqrt(c2-muj4);
        yp = 1. - 2.*muk*(1.-muk)/(1.-2.*muj2 - muk2);
        m_aterm -= CSC.TR*(b*d*((-8.*muj2*muj2 + 2.*c2 +2.*c+4.*muj2)*log((m_alpha*c2*yp-d*e -muj4)/(-b*e+c2*yp-muj4))
              +2.*(c2+c-muj4+2.*muj2)*log((1.-yp)/(1.-m_alpha*yp)) + (-3.*c2 +4.*c*muj2 - 2.*c)*log((m_alpha*c*yp+d)/(b+c*yp))
              +2.*(c2-2.*(c+1.)*muj2+muj4)*(atan(2.*muj2/d) - atan(2.*muj2/b))) 
              +c*(c2*yp*(m_alpha*m_alpha*b*yp - 2.*m_alpha*b - d*(yp-2.)) + 4.*c*muj2*(b*(m_alpha*yp -1.) - d*yp +d) +muj4*(b-d)))
              /(3.*c*a*a*b*d);///ref 0 eq A6 reformulated. 
      }
    }
  }
}


double Massive_KernelsA::at1(int type,int spin,double muq2,double x)
// g(x)
{
///when doing x-dependent part of plus-func, muq2 is sja as defined in CS. When doing endpoints, muq2= (muq2)tilde
  if (type==2||type==3) return 0.;
  double res(0.);
  if (spin == 1) {  ///final state is quark
    if (IsZero(muq2)) {
      if (x<1.-m_alpha) res = 2.*log(1.-x)/(1.-x) + 1.5/(1.-x);
    }
    else {
      if (x<1.-m_alpha) res -= 2.*(log((1.+muq2)/muq2) - 1.)/(1.-x);
    }
  }
  else {
    /// final state is gluon - sum over all possible splittings
    if (x<1.-m_alpha) res -= CSC.TR/CSC.CA*m_nf*(2./3./(1.-x));
   if (x<1.-m_alpha) res -=( -2./(1.-x)*log(1.-x) - 11./6./(1.-x)); 
    size_t nfjk=0;
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i])<muq2) nfjk++;
    for(size_t i=0; i<nfjk;i++) {
      double muQ2 = (m_massflav[i]*m_massflav[i])/muq2;
      if (x<1.-m_alpha) res +=2./3.*((1.-x+2.*muQ2)/sqr(1.-x))*sqrt(1.-4.*muQ2/(1.-x));
    }
  }
  return res;
}

double Massive_KernelsA::at2(int type,int spin,double muq2)
// h; in case of gluon muq2 must be s_ja!!
{
  if (type==2||type==3) return 0.;
  double res(0.);
  double loga = log(m_alpha);
  if (spin == 1) {  ///final state is quark
    if (IsZero(muq2)) {
        res += (-1.5*loga - loga*loga);
    }
    else {
      res += 2.*log(m_alpha)*(log((1.+muq2)/muq2)-1.);
    }
  }
  else {
    /// final state is gluon - sum over all possible splittings
    res += CSC.TR/CSC.CA*m_nf*(loga*2./3.);
 res -= /*2.**/(loga*loga + 11./6.*loga); 
    double muQ2, a, b, c;
    c=sqrt(m_alpha);
    size_t nfjk=0;
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i])<muq2) nfjk++;
    for(size_t i=0; i<nfjk;i++) {
      muQ2 = (m_massflav[i]*m_massflav[i])/muq2;
      a = sqrt(1.-4.*muQ2);
      b = sqrt(m_alpha-4.*muQ2);
      res +=2./9.*(-4.*muQ2*(b/m_alpha/c +4./a) - 5.*b/c - sqr(4.*muQ2)/a +5./a +6.*log(b+c) - 6.*log(a+1.));
    }
  }
  return res;
}


double Massive_KernelsA::at3(int type,int spin,double muq2,double x)
// k(x)
{
  if (spin==2) muq2 = muq2*x;
  else muq2 = muq2/x; ///changing mu-notilde back to mutilde.
  double res(0.);
  if (type!=2 && type!=3){
  if (spin == 1) {  ///final state is quark
    if (IsZero(muq2)) {
      if (x<1.-m_alpha) {
        res += -2.*log(2.-x)/(1.-x);
      }
      
    }
    else {
      if (x<1.-m_alpha) {
        res -= (1.-x)/(2.*sqr(1.-x+muq2*x)) +2./(1.-x)*log((2.-x+muq2*x)*muq2/(1.+muq2)/(1.-x+muq2*x));
      }
    }
  }
  else {
    /// final state is gluon - sum over all possible splittings (only one non-zero here is g->gg)
  if (x<1.-m_alpha) res -= (2.*log(2.-x)/(1.-x));
  }
  }
  double zp=(1.-x)/(1.-x+muq2*x);
  if (spin==2) {
    zp = 1.;
  }
  switch (type) {
  case 1:
    if (zp>m_alpha) res-=2./(1.-x)*log(zp*(1.-x+m_alpha)/m_alpha/(1.-x+zp)) - (1.+x)*log(zp/m_alpha);
    break;
  case 2:
    if (zp>m_alpha) {
      if (zp!=1.) res+= -CSC.CF/CSC.CA*((1.+sqr(1.-x))/x*(log(zp/m_alpha)) +2.*muq2*log((1.-zp)/(1.-m_alpha)));
      else res += -CSC.CF/CSC.CA*(2.-2.*x+x*x)/x*log(zp/m_alpha);
    }
    break;
  case 3:
    if (zp>m_alpha) res += -CSC.TR/CSC.CF*(x*x+sqr(1.-x))*log(zp/m_alpha);
    break;
  case 4:
    if (zp>m_alpha) {
      if (zp!=1.) res += -2.*((1./x-2.+x*(1.-x))*log(zp/m_alpha) -muq2*log((1.-zp)/(1.-m_alpha)) - log(m_alpha*(1.-x+zp)/(zp*(1.-x+m_alpha)))/(1.-x));
      else res+= -2.*((1./x-2.+x*(1.-x))*log(zp/m_alpha) - log(m_alpha*(1.-x+zp)/(zp*(1.-x+m_alpha)))/(1.-x));
    }
    break;
  }
  return res;
}

double Massive_KernelsA::at4(int type,int spin,double muq2,double x)
// G(x)
{
  if (type==2||type==3) return 0.;
  double res(0.);
  if (spin == 1) {  ///final state is quark
    if (IsZero(muq2)) {
      if (x>1.-m_alpha) res -= sqr(m_loga) + 1.5*m_loga;
      if (x<1.-m_alpha) res -= sqr(log(1.-x)) + 1.5*log(1.-x);
    }
    else {
      if (x>1.-m_alpha) res -= - 2.*(log((1.+muq2)/muq2) - 1.)*m_loga;
      if (x<1.-m_alpha) res -= - 2.*(log((1.+muq2)/muq2) - 1.)*log(1.-x);
    }
  }
  else {
    /// final state is gluon - sum over all possible splittings
 if (x>1.-m_alpha) res -=((-CSC.TR/CSC.CA*m_nf*2./3.+11./6.)*m_loga + sqr(m_loga));
 if (x<1.-m_alpha) res-=( (-CSC.TR/CSC.CA*m_nf*2./3.+11./6.)*log(1.-x) + sqr(log(1.-x)));
    size_t nfjk=0;
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i])<muq2) nfjk++;
    for(size_t i=0; i<nfjk;i++) {
      double muQ2 = (m_massflav[i]*m_massflav[i])/muq2;
      double rt = sqrt(1.-4.*muQ2);
      double rta = sqrt(1.-4.*muQ2/m_alpha);
      if (x>1.-m_alpha) res +=2./3.*log(2.*m_alpha*(rta +1.) - 4.*muQ2) - 2./9./m_alpha*rta*(4.*muQ2 +5.*m_alpha)
                                          +2./9.*(4.*rt*muQ2 + 5.*rt - 3.*log(-2.*muQ2+rt+1.) - log(8.));
      rta = sqrt(1.-4.*muQ2/(1.-x));
      if (x<1.-m_alpha) res +=2./3.*log(2.*(1.-x)*(rta +1.) - 4.*muQ2) - 2./9./(1.-x)*rta*(4.*muQ2 +5.*(1.-x))
                                          +2./9.*(4.*rt*muQ2 + 5.*rt - 3.*log(-2.*muQ2+rt+1.) - log(8.));
    }
  }
  return res;
}



double Massive_KernelsA::Kt1(int type,double x)
{
  switch(type) {
  case 1:
    return 2./(1.-x)*log(1.-x);
  case 4:
    return 2./(1.-x)*log(1.-x);
  }
  return 0.;
}

double Massive_KernelsA::Kt2(int type)
{
  if (type==1||type==4) return -sqr(M_PI)/3.;
  return 0.;
}

double Massive_KernelsA::Kt3(int type,double x)
{
  double at=0.,ax=0.;
  if (m_alpha<1.) {
    if (m_alpha<(1.-x)) ax=log(m_alpha/(1.-x));
  }
  switch(type) {
  case 1:
    ax*=(1.+x*x)/(1.-x);
    return -(1.+x)*(log(1.-x))+at+ax;
  case 2:
    ax*=(1.+sqr(1.-x))/x;
    return CSC.CF/CSC.CA*((1.+sqr(1.-x))/x*(log(1.-x))+ax);
  case 3:
    ax*=(1.-2.*x*(1.-x));
    return CSC.TR/CSC.CF*((x*x+sqr(1.-x))*(log(1.-x))+ax);
  case 4:
    ax*=x/(1.-x)+(1.-x)/x+x*(1.-x);
    return 2.*((1.-x)/x-1.+x*(1.-x))*(log(1.-x))+at+2.*ax;
  }
  return 0.;
}

double Massive_KernelsA::Kt4(int type,double x)
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

