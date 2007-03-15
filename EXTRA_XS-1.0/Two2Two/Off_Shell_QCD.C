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
#include "Doubly_Unintegrated_PDF.H"

#define NC 3.0

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PDF;

namespace EXTRAXS {

  template <> Single_XS *
  Single_XS::GetProcess<Off_Shell_gg_qqb>(const size_t nin,const size_t nout,
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
  msg_Info()<<METHOD<<"(): Init matrix element according to "
	    <<"Nucl.Phys.B 366(1991)135"<<std::endl;
  m_alphas=(*MODEL::as)(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms()));
  m_nvector=m_nvector+2;
  CreateMomenta(m_nvector);
  p_addflavours = new ATOOLS::Flavour[4];
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  m_nstrong=4;
  p_pdfs[1]=p_pdfs[0]=NULL;
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

bool Off_Shell_gg_qqb::SetPDFFlavours()
{
  if (p_pdfs[0]==NULL || p_pdfs[1]==NULL) {
    for (short unsigned int i(0);i<m_nin;++i) {
      p_pdfs[i]=dynamic_cast<Doubly_Unintegrated_PDF*>(p_isrhandler->PDF(i));
      if (p_pdfs[i]==NULL || 
	  p_pdfs[i]->Type().find("DUPDF")==std::string::npos)
	THROW(fatal_error,"Low-x ME needs UPDF.");
    }
  }
  // dice pdf jet flavours
  if (!p_pdfs[0]->SelectJetFlavour
      (p_addflavours[2],p_addflavours[0],ran.Get())) return false;
  if (!p_pdfs[1]->SelectJetFlavour
      (p_addflavours[3],p_addflavours[1],ran.Get())) return false;
  return true;
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

