#include "Off_Shell_QCD.H"

#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"
#include "Data_Reader.H"
#include "Random.H"
#include "Running_AlphaS.H"
#include "Flow.H"
#include "MathTools.H"
#include "Combined_Selector.H"
#include "Standard_Selector.H"

#define NC 3.0

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;

namespace EXTRAXS {

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
  double S2=p[4]*p[5], M2=p[2].Abs2();
  double a3=p[5]*p[2], a4=p[5]*p[3];
  double b3=p[4]*p[2], b4=p[4]*p[3];
  double k12=-p[0].Abs2(), k22=-p[1].Abs2();
  std::swap<double>(t,u);
  double delta=
    2.*S2*(2.*(a4*b3-a3*b4)/S2
	  -k12*b4/(p[4]*p[1])+k22*a4/(p[5]*p[0])+p[3]*(p[0]-p[1]));
  double abelian=
    4.*(S2*S2/((t-M2)*(u-M2))
	-1./(k12*k22)*ATOOLS::sqr(S2+2.*(a3*b4/(t-M2)+b3*a4/(u-M2))));
  double nonabelian=
    4.*S2*S2*(S2/(s*(p[4]*p[1])*(p[5]*p[0]))-1./((t-M2)*(u-M2))
	    -(1/(t-M2)-1/(u-M2))/s*(b3/(p[4]*p[1])-a3/(p[5]*p[0])))
    +2./(k12*k22)*(S2+4.*a3*b4/(t-M2)-delta/s)*(S2+4.*a4*b3/(u-M2)+delta/s);
  std::swap<double>(t,u);
  return ATOOLS::sqr(4.*M_PI*m_alphas*(p[0]*p[5])*(p[1]*p[4])/(S2*S2))*
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

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_gg_gg>(const size_t nin,const size_t nout,
						  const ATOOLS::Flavour *flavours,
					      const size_t nqed, const size_t nqcd)
{
  for (size_t i=0;i<nin+nout;++i) if (!flavours[i].IsGluon()) return NULL;
  if (nqcd==nin+nout-2 && nqed==0) 
    return new Off_Shell_gg_gg(nin,nout,flavours); 
  return NULL;
}

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
  m_nstrong=m_nin+m_nout;
}

double Off_Shell_gg_gg::operator()(double s,double t,double u) 
{
  ATOOLS::Vec4D *const p=p_momenta;
  double colfac=NC/(NC*NC-1.0);
  if (m_nout>0) for (size_t i=0;i<m_nout-1;++i) colfac*=NC;
  double scale=sqrt((p[m_nin+m_nout]-p[0]).PPerp2()*
		    (p[m_nin+m_nout+1]-p[1]).PPerp2());
  double asmean((*MODEL::as)(scale));
  if (p[0].PPerp()<m_qtcut || p[1].PPerp()<m_qtcut) return 0.0;
  double lasty((p[m_nin+m_nout]-p[0]).Y());
  Vec4D prop(p[0]);
  if (p_selector->Name()=="Combined_Selector") {
    Selector_Base *ipt = 
      ((Combined_Selector*)p_selector)->GetSelector("BFKL_PT_Selector");
    if (ipt==NULL) THROW(fatal_error,"BFKL ME needs BFKL_PT selector");
    if (!ipt->GetValue("qcut",m_qtcut)) 
      THROW(fatal_error,"Selector corrupted");
  }
  double M=4.0*sqr((p[0]+p[1]).Abs2());
  for (size_t i=0;i<m_nout;++i) {
    double y(p[m_nin+i].Y());
    if (y>=lasty) return 0.0;
    M*=2/p[m_nin+i].PPerp2();
    if (i>0) {
      double sud(exp(-asmean*(lasty-y)*log(prop.PPerp2()/m_qtcut)));
      if (sud>1.0) 
	msg.Error()<<"Off_Shell_gg_gg::operator()(..):"
		   <<"Sudakov "<<i<<" exceeds unity, \\Delta = "
		   <<sud<<" = exp(-"<<asmean<<"*("<<lasty<<"-"<<y<<")*log("
		   <<prop.PPerp2()<<"/"<<m_qtcut<<"))"<<std::endl;
      M*=sud;
    }
    prop-=p[m_nin+i];
    lasty=y;
  }
  double y((p[m_nin+m_nout+1]-p[1]).Y());
  if (y>=lasty) return 0.0;
  return pow(4.0*M_PI*m_alphas,m_nout)*colfac*M/4.0;
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
