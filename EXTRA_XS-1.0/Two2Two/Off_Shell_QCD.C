#include "Off_Shell_QCD.H"

#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"
#include "Data_Reader.H"
#include "Random.H"
#include "Flow.H"

#define NC 3

using namespace EXTRAXS;

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_gg_qqb>(const size_t nin,const size_t nout,
						   const ATOOLS::Flavour *flavours,
						   const size_t nqed, const size_t nqcd)
{
  if (flavours[2].IsQuark() && flavours[3]==flavours[2].Bar() && 
      flavours[0].IsGluon() && flavours[1].IsGluon()){ 
    if (nqcd==2 && nqed==0) {
      return new Off_Shell_gg_qqb(nin,nout,flavours); 
    }
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
  p_addflavours = new ATOOLS::Flavour[2];
  p_addflavours[0]=ATOOLS::kf::gluon;
  p_addflavours[1]=ATOOLS::kf::gluon;
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  m_nstrong=4;
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
						  const ATOOLS::Flavour *flavours,
					      const size_t nqed, const size_t nqcd)
{
  if (flavours[2].IsGluon() && flavours[3].IsGluon() && 
      flavours[0].IsGluon() && flavours[1].IsGluon()){ 
    if (nqcd==2 && nqed==0) {
      return new Off_Shell_gg_gg(nin,nout,flavours); 
    }
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
  p_addflavours = new ATOOLS::Flavour[2];
  p_addflavours[0]=ATOOLS::kf::gluon;
  p_addflavours[1]=ATOOLS::kf::gluon;
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  m_zkey[0].Assign("z_1",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_zkey[1].Assign("z_2",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  if (!reader->ReadFromFile(m_jets,"LIPATOV_JETS","")) m_jets=2;
  delete reader;
  m_nstrong=4;
}

#include "MyComplex.H"

Complex App(const ATOOLS::Vec4D *p)
{
  double kp[4];
  Complex kt[4], qt[3];
  for (short unsigned int i=2;i<4;++i) {
    kp[i]=p[i].PPlus();
    double spt=p[i].PPerp();
    kt[i]=Complex(spt*p[i].CosPhi(),spt*p[i].SinPhi());
    ATOOLS::Vec4D pi=p[i-2];
    if (i==3) pi=-1.0*pi;
    spt=pi.PPerp();
    qt[i-1]=Complex(spt*pi.CosPhi(),spt*pi.SinPhi());
  }
  Complex ab23=sqrt(kp[3]/kp[2])*kt[2]-sqrt(kp[2]/kp[3])*kt[3];
  return 2.0*std::conj(qt[1])*qt[2]/kt[2]*sqrt(kp[2]/kp[3])/ab23;
}

Complex Apm(const ATOOLS::Vec4D *p)
{
  double kp[4], km[4];
  Complex kt[4], qt[3];
  for (short unsigned int i=2;i<4;++i) {
    kp[i]=p[i].PPlus();
    km[i]=p[i].PMinus();
    double spt=p[i].PPerp();
    kt[i]=Complex(spt*p[i].CosPhi(),spt*p[i].SinPhi());
    ATOOLS::Vec4D pi=p[i-2];
    if (i==3) pi=-1.0*pi;
    spt=pi.PPerp();
    qt[i-1]=Complex(spt*pi.CosPhi(),spt*pi.SinPhi());
  }
  double qt12=p[0].PPerp2(), qt22=p[1].PPerp2();
  double s23=(p[2]+p[3]).Abs2();
  double s3bb=(p[3]-p[1]).Abs2();
  return -2.0*std::conj(kt[2])/kt[2]*
    (-1.0/s23*(kt[3]*kt[3]*qt12/((km[2]+km[3])*kp[3])+
	       kt[2]*kt[2]*qt22/((kp[2]+kp[3])*km[2])+
	       s3bb*kt[2]*kt[3]/(km[2]*kp[3]))
     +(qt[2]+kt[3])*(qt[2]+kt[3])/s3bb
     -(qt[2]+kt[3])/s23*((km[2]+km[3])/km[2]*kt[2]-
			 (kp[2]+kp[3])/kp[3]*kt[3]));
}

Complex Lpp(const ATOOLS::Vec4D *p)
{
  double kp[4], km[4];
  Complex kt[4], qt[3];
  for (short unsigned int i=2;i<4;++i) {
    kp[i]=p[i].PPlus();
    km[i]=p[i].PMinus();
    double spt=p[i].PPerp();
    kt[i]=Complex(spt*p[i].CosPhi(),spt*p[i].SinPhi());
    ATOOLS::Vec4D pi=p[i-2];
    if (i==3) pi=-1.0*pi;
    spt=pi.PPerp();
    qt[i-1]=Complex(spt*pi.CosPhi(),spt*pi.SinPhi());
  }
  return 2.0*std::conj(qt[1])*(p[0]-p[2]).PPerp2()/(p[0]-p[2]).Abs2()*qt[2]
    /kt[2]/kt[3];
}

// #define USING__NLO_LEV

double Off_Shell_gg_gg::operator()(double s,double t,double u) 
{
  ATOOLS::Vec4D *const p=p_momenta;
#ifdef USING__NLO_LEV
  // nlo lipatov vertex
  Complex app34=App(p);
  Complex apm34=Apm(p);
  std::swap<ATOOLS::Vec4D>(p[2],p[3]);
  Complex app43=App(p);
  Complex amp43=Apm(p);
  std::swap<ATOOLS::Vec4D>(p[2],p[3]);
  double Mpp=std::abs(app34*std::conj(app34)+app43*std::conj(app43)
		      +0.5*(app34*std::conj(app43)+std::conj(app34)*app43));
  double Mpm=std::abs(apm34*std::conj(apm34)+amp43*std::conj(amp43));
  double M=2.0*(Mpp+Mpm);
#else
  Complex lpp34=Lpp(p);
  double M=2.0*(2.0*std::abs(lpp34*std::conj(lpp34)));
#endif
  return ATOOLS::sqr(4.0*M_PI*m_alphas)*(NC*NC)/(NC*NC-1)/(2.0*4.0)*M;
}

bool Off_Shell_gg_gg::SetColours(double s,double t,double u) 
{ 
  return true; 
}

bool Off_Shell_gg_gg::Trigger(const ATOOLS::Vec4D *const momenta) 
{
  bool result=true;
  if (m_jets<4) result=Integrable_Base::Trigger(momenta);
  else {
    ATOOLS::Vec4D temp[6];
    for (int i=2;i<4;++i) {
      temp[i-2]=momenta[i+2];
      temp[i+2]=momenta[i];
      temp[i]=p_addmomenta[i-2];
    }
    p_selector->SetNOut(m_nvector-m_nin);
    p_selector->SetNTot(m_nvector);
    result=PHASIC::Integrable_Base::Trigger(temp);
    p_selector->SetNTot(m_nin+m_nout);
    p_selector->SetNOut(m_nout);
  }
  return result;
}
