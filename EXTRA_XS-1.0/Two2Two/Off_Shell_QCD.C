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
  m_zkey[0].Assign("z_1",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_zkey[1].Assign("z_2",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  if (!reader->ReadFromFile(m_jets,"LIPATOV_JETS","")) m_jets=4;
  delete reader;
}

double Off_Shell_gg_gg::operator()(double s,double t,double u) 
{
  const ATOOLS::Vec4D *p=p_momenta;
  const ATOOLS::Vec4D &p1=p[4], &p2=p[5];
  const ATOOLS::Vec4D &k1=p[0], k2=-1.*p[1], &p3=p[2], &p4=p[3];
  ATOOLS::Vec4D k=k1-p3, l=k1-p4;;
  double S=2.*p1*p2, S2=2./S;
  double a=p2*k*S2, b=p1*k*S2;
  double a1=p2*k1*S2, b2=p1*k2*S2;
  const ATOOLS::Vec4D epslp(0.0,1.0,0.0,0.0), epsrp(0.0,0.0,1.0,0.0);
  double kp12=k1.PPerp2(), kp22=k2.PPerp2(), kp2=k.PPerp2();
  ATOOLS::Vec4D L1=(a1+2.0*kp12/(b*S))*p1+(b+2.0*kp2/(a1*S))*p2-(k1+k).Perp();
  ATOOLS::Vec4D L2=(a+2.0*kp2/(b2*S))*p1+(b2+2.0*kp22/(a*S))*p2-(k+k2).Perp();
  ATOOLS::Vec4D eps1=epsrp-(k1*epsrp)/(k1*p2)*p2;
  ATOOLS::Vec4D eps2=epslp-(k2*epslp)/(k2*p1)*p1;
  double A=ATOOLS::sqr((L1*eps1)*(L2*eps2)/kp2);
  double Xi1=p3.PMinus()/p3.PPlus(), Xi2=p4.PPlus()/p4.PMinus();
  const double m_cut=Xi1*Xi2;
  if (m_cut>=1.) return 0.;
  if (kp12/kp22<m_cut || kp22/kp12<m_cut) return 0.0;
  if (kp2/kp12<m_cut || kp12/kp2<m_cut) return 0.0;
  if (kp2/kp22<m_cut || kp22/kp2<m_cut) return 0.0;
  return ATOOLS::sqr(4.*M_PI*m_alphas)*(NC*NC)/(NC*NC-1)/(2.*4.)*A;
}

bool Off_Shell_gg_gg::SetColours(double s,double t,double u) 
{ 
  return true; 
}

bool Off_Shell_gg_gg::Trigger(const ATOOLS::Vec4D *const momenta) 
{
  if (m_jets<4) return Integrable_Base::Trigger(momenta);
  ATOOLS::Vec4D *temp = new ATOOLS::Vec4D[6];
  for (int i=4;i<6;++i) temp[i-4]=momenta[i];
  for (int i=2;i<4;++i) temp[i+2]=momenta[i];
  double x1=momenta[0]*momenta[5]/(momenta[4]*momenta[5]);
  double x2=momenta[1]*momenta[4]/(momenta[4]*momenta[5]);
  double b1=m_zkey[0][2]/(1.-m_zkey[0][2])*
    momenta[0].PPerp2()/x1/p_isrhandler->Pole();
  double b2=m_zkey[1][2]/(1.-m_zkey[1][2])*
    momenta[1].PPerp2()/x2/p_isrhandler->Pole();
  temp[2]=x1*(1./m_zkey[0][2]-1.)*momenta[4]+b1*momenta[5]-momenta[0].Perp();
  temp[3]=x2*(1./m_zkey[1][2]-1.)*momenta[5]+b2*momenta[4]-momenta[1].Perp();
  p_selector->SetNOut(m_nvector-m_nin);
  p_selector->SetNTot(m_nvector);
  bool result=p_selector->Trigger(temp);
  p_selector->SetNTot(m_nin+m_nout);
  p_selector->SetNOut(m_nout);
  delete temp;
  return result;
}
