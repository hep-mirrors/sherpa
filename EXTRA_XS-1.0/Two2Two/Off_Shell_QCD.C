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
  p_addflavours = new ATOOLS::Flavour[2];
  p_addflavours[0]=ATOOLS::kf::gluon;
  p_addflavours[1]=ATOOLS::kf::gluon;
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
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
  p_addflavours = new ATOOLS::Flavour[2];
  p_addflavours[0]=ATOOLS::kf::gluon;
  p_addflavours[1]=ATOOLS::kf::gluon;
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  m_zkey[0].Assign("z_1",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_zkey[1].Assign("z_2",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  if (!reader->ReadFromFile(m_jets,"LIPATOV_JETS","")) m_jets=4;
  delete reader;
}

#include "MyComplex.H"

#define USING__NLO_Lipatov_Vertex

Complex App(double k1p,double k1m,Complex kp1,
	    double k2p,double k2m,Complex kp2,
	    double p3p,double p3m,Complex pp3,
	    double p4p,double p4m,Complex pp4,
	    double s,double s3bb)
{
  return 2.0*std::conj(kp1)*kp2/pp3*sqrt(p3p/p4p/s);
}

Complex Apm(double k1p,double k1m,Complex kp1,
	    double k2p,double k2m,Complex kp2,
	    double p3p,double p3m,Complex pp3,
	    double p4p,double p4m,Complex pp4,
	    double s,double s3bb)
{
  double kp12=ATOOLS::sqr(std::abs(k1p));
  double kp22=ATOOLS::sqr(std::abs(k2p));
  return -2.0*std::conj(kp2)/kp2*
    (-1.0/s*(pp4*pp4*kp12/((p3m+p4m)*p4p)+
	     pp3*pp3*kp22/((p3p+p4p)*p3m)+
	     s3bb*pp3*pp4/(p3m*p4p))
     +(kp2+pp4)*(kp2+pp4)/s3bb
     -(kp2+pp4)/s*((p3m+p4m)/p3m*pp3-(p3p+p4p)/p4p*pp4));
}

#define USING__NLO_LEV
// #define USING__No_Turn

double Off_Shell_gg_gg::operator()(double s,double t,double u) 
{
#ifdef USING__NLO_LEV
  // use nlo lipatov vertex
  ATOOLS::Vec4D *const p=p_momenta;
#ifndef USING__No_Turn
  double c1=(p_addmomenta[0]+p[2]+p[3]).PPlus()/p[4].PPlus();
  double c2=(p_addmomenta[1]+p[2]+p[3]).PMinus()/p[5].PMinus();
  bool turn=false;
  if (c2>c1) {
    ATOOLS::Vec4D plus(1.0,0.0,0.0,1.0);
    ATOOLS::Vec4D minus(1.0,0.0,0.0,-1.0);
    ATOOLS::Poincare rot(plus,minus);
    for (short unsigned int i=0;i<6;++i) rot.Rotate(p[i]);
    for (short unsigned int i=0;i<2;++i) rot.Rotate(p_addmomenta[i]);
    std::swap<ATOOLS::Vec4D>(p[0],p[1]);
    std::swap<ATOOLS::Vec4D>(p[2],p[3]);
    std::swap<ATOOLS::Vec4D>(p[4],p[5]);
    std::swap<ATOOLS::Vec4D>(p_addmomenta[0],p_addmomenta[1]);
    turn=true;
  }
#endif
  double pp[6], pm[6];
  Complex pt[6];
  for (short unsigned int i=0;i<6;++i) {
    ATOOLS::Vec4D pi=p[i];
    if (i==1) pi=-1.0*pi;
    pp[i]=pi.PPlus();
    pm[i]=pi.PMinus();
    double spt=pi.PPerp();
    pt[i]=Complex(spt*pi.CosPhi(),spt*pi.SinPhi());
    //    std::cout<<i<<" "<<pi<<" "<<pp[i]<<" "<<pm[i]<<" "<<pt[i]<<std::endl;
  }
  double s42=(p[3]+p_addmomenta[1]+p[2]).Abs2();
  Complex app34=App(pp[0],pm[0],pt[0],pp[1],pm[1],pt[1],
		    pp[2],pm[2],pt[2],pp[3],pm[3],pt[3],s,s42);
  Complex apm34=Apm(pp[0],pm[0],pt[0],pp[1],pm[1],pt[1],
		    pp[2],pm[2],pt[2],pp[3],pm[3],pt[3],s,s42);
  Complex app43=App(pp[0],pm[0],pt[0],pp[1],pm[1],pt[1],
		    pp[3],pm[3],pt[3],pp[2],pm[2],pt[2],s,s42);
  Complex apm43=Apm(pp[0],pm[0],pt[0],pp[1],pm[1],pt[1],
		    pp[3],pm[3],pt[3],pp[2],pm[2],pt[2],s,s42);
  double Mpp=ATOOLS::sqr(std::abs(app34))+ATOOLS::sqr(std::abs(app43))
    +0.5*std::abs(app34*std::conj(app43)+std::conj(app34)*app43);
  double Mpm=ATOOLS::sqr(std::abs(apm34))+ATOOLS::sqr(std::abs(apm43))
    +0.5*std::abs(apm34*std::conj(apm43)+std::conj(apm34)*apm43);
  double M=(Mpp+Mpm);
  // check angular ordering of emitted gluons
  const ATOOLS::Vec4D &g1=p_addmomenta[0];
  const ATOOLS::Vec4D &g2=p_addmomenta[1];
  double Xi1=pm[2]/pp[2], Xi2=pp[3]/pm[3];
  double xi1=g1.PMinus()/g1.PPlus(), xi2=g2.PPlus()/g2.PMinus();
  if (xi1>=Xi1 || xi2>=Xi2) return 0.;
#ifndef USING__No_Turn
  if (turn) {
    ATOOLS::Vec4D plus(1.0,0.0,0.0,1.0);
    ATOOLS::Vec4D minus(1.0,0.0,0.0,-1.0);
    ATOOLS::Poincare rot(plus,minus);
    for (short unsigned int i=0;i<6;++i) rot.Rotate(p[i]);
    for (short unsigned int i=0;i<2;++i) rot.Rotate(p_addmomenta[i]);
    std::swap<ATOOLS::Vec4D>(p[0],p[1]);
    std::swap<ATOOLS::Vec4D>(p[2],p[3]);
    std::swap<ATOOLS::Vec4D>(p[4],p[5]);
    std::swap<ATOOLS::Vec4D>(p_addmomenta[0],p_addmomenta[1]);
  }
#endif
  return ATOOLS::sqr(4.0*M_PI*m_alphas)
    *(NC*NC)/(NC*NC-1)/(2.0*4.0)*M/(2.0*p[4]*p[5]);
#else
  const ATOOLS::Vec4D *p=p_momenta;
  const ATOOLS::Vec4D &p1=p[4], &p2=p[5];
  const ATOOLS::Vec4D &k1=p[0], k2=-1.*p[1], &p3=p[2], &p4=p[3];
  ATOOLS::Vec4D k=k1-p3, l=k1-p4, m=k1+k2;
  double S=2.*p1*p2, S2=2./S;
  double a=p2*k*S2, b=p1*k*S2;
  double a1=p2*k1*S2, b2=p1*k2*S2;
  double kp12=k1.PPerp2(), kp22=k2.PPerp2(), kp2=k.PPerp2();
#ifndef USING__LEV
  // use equivalence of lipatov vertex and triple gluon vertex
  ATOOLS::Vec4D A1=
    (k1.Perp()*k.Perp())*(k1+k)
    +(k1.Perp()*(p3-k))*k.Perp()
    -(k.Perp()*(k1+p3))*k1.Perp();
  ATOOLS::Vec4D A2=
    (k.Perp()*k2.Perp())*(k+k2)
    +(k.Perp()*(p4-k2))*k2.Perp()
    -(k2.Perp()*(k+p4))*k.Perp();
  double M1=2*(A1*p3)*(A1*p4)/(p3*p4)-A1*A1;
  double M2=2*(A2*p4)*(A2*p3)/(p4*p3)-A2*A2;
  double M=(M1*M2)*ATOOLS::sqr(S2*S2/(kp2*a1*a*b*b2));
  double Xi1=p3.PMinus()/p3.PPlus(), Xi2=p4.PPlus()/p4.PMinus();
  double m_cut=Xi1*Xi2;
  if (m_cut>=1.) return 0.0;
#else
  // use lipatov vertex 
  ATOOLS::Vec4D L1=
    (a1+2.0*kp12/(b*S))*p1
    +(b+2.0*kp2/(a1*S))*p2
    -(k1+k).Perp();
  ATOOLS::Vec4D L2=
    (a+2.0*kp2/(b2*S))*p1
    +(b2+2.0*kp22/(a*S))*p2
    -(k+k2).Perp();
  double M1=2*L1*p3*L1*p4/(p4*p3)-L1*L1;
  double M2=2*L2*p4*L2*p3/(p3*p4)-L2*L2;
  double M=(M1*M2)/ATOOLS::sqr(kp2);
  double Xi1=p3.PMinus()/p3.PPlus(), Xi2=p4.PPlus()/p4.PMinus();
  double m_cut=Xi1*Xi2;
  if (m_cut>=1.) return 0.0;
  // check kperp similarity
  if (kp12/kp22<m_cut || kp22/kp12<m_cut) return 0.0;
  if (kp2/kp12<m_cut || kp12/kp2<m_cut) return 0.0;
  if (kp2/kp22<m_cut || kp22/kp2<m_cut) return 0.0;
#endif
#endif
  return ATOOLS::sqr(4.0*M_PI*m_alphas)*(NC*NC)/(NC*NC-1)/(2.0*4.0)*M;
}

bool Off_Shell_gg_gg::SetColours(double s,double t,double u) 
{ 
  return true; 
}

bool Off_Shell_gg_gg::Trigger(const ATOOLS::Vec4D *const momenta) 
{
//   ATOOLS::Vec4D *temp = new ATOOLS::Vec4D[6];
//   for (int i=4;i<6;++i) temp[i-4]=momenta[i];
//   for (int i=2;i<4;++i) temp[i+2]=momenta[i];
//   double x1=momenta[0]*momenta[5]/(momenta[4]*momenta[5]);
//   double x2=momenta[1]*momenta[4]/(momenta[4]*momenta[5]);
//   double b1=m_zkey[0][2]/(1.-m_zkey[0][2])*
//     momenta[0].PPerp2()/x1/p_isrhandler->Pole();
//   double b2=m_zkey[1][2]/(1.-m_zkey[1][2])*
//     momenta[1].PPerp2()/x2/p_isrhandler->Pole();
//   temp[2]=x1*(1./m_zkey[0][2]-1.)*momenta[4]+b1*momenta[5]-momenta[0].Perp();
//   temp[3]=x2*(1./m_zkey[1][2]-1.)*momenta[5]+b2*momenta[4]-momenta[1].Perp();
//   if (!(temp[2]==p_addmomenta[0])) std::cout<<temp[2]<<p_addmomenta[0]<<std::endl;
//   if (!(temp[3]==p_addmomenta[1])) std::cout<<temp[3]<<p_addmomenta[1]<<std::endl;
  bool result=true;
  if (m_jets<4) result=Integrable_Base::Trigger(momenta);
  else {
    ATOOLS::Vec4D temp[6];
    for (int i=4;i<6;++i) temp[i-4]=momenta[i];
    for (int i=2;i<4;++i) temp[i+2]=p_addmomenta[i-2];
    p_selector->SetNOut(m_nvector-m_nin);
    p_selector->SetNTot(m_nvector);
    result=p_selector->Trigger(temp);
    p_selector->SetNTot(m_nin+m_nout);
    p_selector->SetNOut(m_nout);
  }
//   p_addmomenta[0]=temp[2];
//   p_addmomenta[1]=temp[3];
// //   for (size_t i=0;i<6;++i) std::cout<<momenta[i]<<std::endl;
// //   for (size_t i=0;i<6;++i) std::cout<<temp[i]<<std::endl;
//   delete [] temp;
  return result;
}
