#include "XS_QCD.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Flow.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;

namespace EXTRAXS {

  template <> Single_XS *Single_XS::GetProcess<XS_pp_q1qbar1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model,
   const size_t nqed,const size_t nqcd)
  {
    if (!model->IncludesModel("QED")) return NULL;
    if (flavours[0].IsPhoton() && flavours[1].IsPhoton() && 
	(flavours[2].IsQuark() || flavours[2].IsLepton()) && 
	flavours[3]==flavours[2].Bar()) { 
      if (nqcd==0 && nqed==2) {
	return new XS_pp_q1qbar1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
  
}

XS_pp_q1qbar1::XS_pp_q1qbar1(const size_t nin,const size_t nout, 
			     const ATOOLS::Flavour *fl,
			     XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  m_r=fl[2].IsAnti();
  m_g=std::abs(p_model->Constant("g_1"));
  m_eq=p_flavours[2].Charge();
  m_m2=sqr(p_flavours[2].Mass());
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_nstrong=2;
}

double XS_pp_q1qbar1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double tp(t-m_m2), up(u-m_m2);
  double mtt(2.0*(tp*(up-2.0*m_m2)-4.0*m_m2*m_m2)/(tp*tp));
  double muu(2.0*(up*(tp-2.0*m_m2)-4.0*m_m2*m_m2)/(up*up));
  double mtu(2.0*m_m2*(s-4.0*m_m2)/(tp*up));
  return sqr(csqr(m_g).real()*sqr(m_eq))*
    (p_flavours[2].Strong()?3.0:1.0)*(mtt+muu+2.0*mtu); 
}


bool XS_pp_q1qbar1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = Min(t,u);
  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = Min(t,u);
  p_colours[2+m_r][0]=p_colours[3-m_r][1]=Flow::Counter();
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_q1q2_q1q2>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsQuark() && flavours[1].IsQuark() && 
	flavours[0]!=flavours[1] &&
	((flavours[2]==flavours[0] && flavours[3]==flavours[1]) ||
	 (flavours[3]==flavours[0] && flavours[2]==flavours[1]))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1q2_q1q2(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_q1q2_q1q2::XS_q1q2_q1q2(const size_t nin,const size_t nout, 
			   const ATOOLS::Flavour *fl,
			   XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=fl[1].IsAnti();
  m_g=std::abs(p_model->Constant("g_3"));
  m_m12=sqr(p_flavours[0].Mass());
  m_m22=sqr(p_flavours[1].Mass());
  m_nstrong=4;
}

double XS_q1q2_q1q2::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double Mt(sqr(u-m_m12-m_m22)+sqr(s-m_m12-m_m22)+2.0*t*(m_m12+m_m22));
  return sqr(m_g*m_g)*4.0/9.0*Mt/(t*t);
}

bool XS_q1q2_q1q2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  int r = !(p_flavours[0] == p_flavours[2]);
  if (m_a==m_p) {
    /*
      0-----\   /-----2, if fl[0]==fl[2]
             \ /
	      X  t
             / \
      1-----/   \-----3, if fl[1]==fl[3]

      shower scale is u
    */
    p_colours[0][m_a] = p_colours[3-r][m_a] = Flow::Counter();
    p_colours[1][m_a] = p_colours[2+r][m_a] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
      pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq'->qq', set scale u "<<u<<"\n";
  }
  else {
    /*
      0-----+ +-----2
            | |
	    | |  t
            | |
      1-----+ +-----3

      shower scale is s
    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();
    p_colours[2+r][m_a] = p_colours[3-r][m_p] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
      pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb'->qqb', set scale s "<<s<<"\n";
  }
  if (swap) SwapInOrder();
  return 1; 
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_q1qbar1_q2qbar2>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model,
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
	flavours[2].IsQuark() && flavours[3]==flavours[2].Bar() &&
	flavours[0]!=flavours[2]) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar1_q2qbar2(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_q1qbar1_q2qbar2::XS_q1qbar1_q2qbar2(const size_t nin,const size_t nout, 
				       const ATOOLS::Flavour *fl,
				       XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_g=std::abs(p_model->Constant("g_3"));
  m_m12=sqr(p_flavours[0].Mass());
  m_m32=sqr(p_flavours[2].Mass());
  m_nstrong=4;
}

double XS_q1qbar1_q2qbar2::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double Ms(sqr(u-m_m12-m_m32)+sqr(t-m_m12-m_m32)+2.0*s*(m_m12+m_m32));
  return sqr(m_g*m_g)*4.0/9.0*Ms/(s*s); 
}

bool XS_q1qbar1_q2qbar2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  int r = !(p_flavours[0].IsAnti() == p_flavours[2].IsAnti());
  /*
    0\         /2, if fl[0].IsAnti()==fl[2].IsAnti()
      \   s   /
       ======= 
      /       \
    1/         \3, if fl[0].IsAnti()==fl[2].IsAnti()

    shower scale is t
  */
  p_colours[0][m_a] = p_colours[2+r][m_a] = Flow::Counter();
  p_colours[1][m_p] = p_colours[3-r][m_p] = Flow::Counter();

  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
  msg_Debugging()<<"xs: qqb->q'qb', set scale t "<<t<<"\n";
  if (swap) SwapInOrder();
  return 1; 
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_q1q1_q1q1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsQuark() && flavours[1]==flavours[0] &&
	flavours[2]==flavours[0] && flavours[3]==flavours[0]) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1q1_q1q1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_q1q1_q1q1::XS_q1q1_q1q1(const size_t nin,const size_t nout, 
			   const ATOOLS::Flavour *fl,
			   XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_g=std::abs(p_model->Constant("g_3"));
  m_m12=sqr(p_flavours[0].Mass());
  m_nstrong=4;
}

double XS_q1q1_q1q1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double Mt(sqr(u-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*t*m_m12); 
  double Mu(sqr(t-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*u*m_m12); 
  double Mtu(-2.0/3.0*(s*s-8.0*(s-2.0*m_m12)*m_m12-4.0*m_m12*m_m12));
  return sqr(m_g*m_g)*4.0/9.0*(Mt/(t*t)+Mu/(u*u)+Mtu/(t*u))/2.0;
}

bool XS_q1q1_q1q1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double Mt(sqr(u-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*t*m_m12); 
  double Mu(sqr(t-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*u*m_m12); 
  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->qq, set scale u "<<u<<"\n";
    /*
      0----\   /----2
            \ /
             X  t
            / \
      1----/   \----3

      shower scale is u
    */
    p_colours[3][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[2][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->qq, set scale t "<<t<<"\n";
    /*
      0----\   /----2
            \ /
             =  u
            / \
      1----/   \----3

      shower scale is t
    */
    p_colours[2][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[3][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_q1qbar1_q1qbar1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
	((flavours[2]==flavours[0] && flavours[3]==flavours[1]) ||
	 (flavours[3]==flavours[0] && flavours[2]==flavours[1]))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar1_q1qbar1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_q1qbar1_q1qbar1::XS_q1qbar1_q1qbar1(const size_t nin,const size_t nout, 
				       const ATOOLS::Flavour *fl,
				       XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_r=fl[0]!=fl[2];
  m_g=std::abs(p_model->Constant("g_3"));
  m_m12=sqr(p_flavours[0].Mass());
  m_nstrong=4;
}

double XS_q1qbar1_q1qbar1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double Mt(sqr(s-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*t*m_m12); 
  double Ms(sqr(t-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*s*m_m12); 
  double Mts(-2.0/3.0*(u*u-8.0*(u-2.0*m_m12)*m_m12-4.0*m_m12*m_m12));
  return sqr(m_g*m_g)*4.0/9.0*(Mt/(t*t)+Ms/(s*s)+Mts/(t*s));
}


bool XS_q1qbar1_q1qbar1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double Mt(sqr(s-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*t*m_m12); 
  double Ms(sqr(t-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*s*m_m12); 
  if (Ms >  (Mt+Ms) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->qqb, set scale t "<<t<<"\n";
    /*
      0\         /2, if fl[0]==fl[2]
        \   s   /
         =======
        /       \
      1/         \3, if fl[0]==fl[2]

      shower scale is t
    */
    p_colours[0][m_a] = p_colours[2+m_r][m_a] = Flow::Counter();	
    p_colours[1][m_p] = p_colours[3-m_r][m_p] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->qqb, set scale s "<<s<<"\n";
    /*
      0----+ +----2
           | |
           | | t
           | |
      1----+ +----3

      shower scale is s
    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();	
    p_colours[2+m_r][m_a] = p_colours[3-m_r][m_p] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_q1qbar1_gg>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
	flavours[2].IsGluon() && flavours[3].IsGluon()) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar1_gg(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_q1qbar1_gg::XS_q1qbar1_gg(const size_t nin,const size_t nout, 
			     const ATOOLS::Flavour *fl,
			     XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_m12=sqr(p_flavours[0].Mass());
  m_g=std::abs(p_model->Constant("g_3"));
  m_nstrong=4;
}

double XS_q1qbar1_gg::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double tp(t-m_m12), up(u-m_m12);
  double Mt(32.0/27.0*(tp*up-m_m12*(4.0*(m_m12+tp)+m_m12*tp/s))/(tp*tp));
  double Mu(32.0/27.0*(up*tp-m_m12*(4.0*(m_m12+up)+m_m12*up/s))/(up*up));
  double Ms(-8.0/3.0*(tp*tp+up*up+4.0*m_m12*s)/(s*s));
  return sqr(m_g*m_g)*(Mt+Mu+Ms)/2.0; 
}


bool XS_q1qbar1_gg::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double tp(t-m_m12), up(u-m_m12);
  double Mt(32.0/27.0*(tp*up-m_m12*(4.0*(m_m12+tp)+m_m12*tp/s))/(tp*tp));
  double Mu(32.0/27.0*(up*tp-m_m12*(4.0*(m_m12+up)+m_m12*up/s))/(up*up));
  p_colours[0][m_a] = Flow::Counter();
  p_colours[1][m_p] = Flow::Counter();
  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->gg, set scale s/t "<<s<<"/"<<t<<"\n";
    /*
      0------+====2
             |
             | t
             |
      1------+====3

      shower scale is s / t
    */
    p_colours[2][m_a] = p_colours[0][m_a];
    p_colours[3][m_p] = p_colours[1][m_p];
    p_colours[2][m_p] = p_colours[3][m_a] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->gg, set scale s/u "<<s<<"/"<<u<<"\n";
    /*
      0----\ +-==2
            \|/
             | u
            /|\
      1----/ +-==3

      shower scale is s / u
    */
    p_colours[3][m_a] = p_colours[0][m_a];
    p_colours[2][m_p] = p_colours[1][m_p];
    p_colours[3][m_p] = p_colours[2][m_a] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_gg_q1qbar1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsGluon() && flavours[1].IsGluon() && 
	flavours[2].IsQuark() && flavours[3]==flavours[2].Bar()) { 
      if (nqcd==2 && nqed==0) {
	return new XS_gg_q1qbar1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_gg_q1qbar1::XS_gg_q1qbar1(const size_t nin,const size_t nout, 
			     const ATOOLS::Flavour *fl,
			     XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_r=fl[2].IsAnti();
  m_g=std::abs(p_model->Constant("g_3"));
  m_m32=sqr(p_flavours[2].Mass());
  m_nstrong=4;
}

double XS_gg_q1qbar1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double tp(t-m_m32), up(u-m_m32);
  double Mt(1.0/6.0*(tp*up-m_m32*(4.0*(m_m32+tp)+m_m32*tp/s))/(tp*tp));
  double Mu(1.0/6.0*(up*tp-m_m32*(4.0*(m_m32+up)+m_m32*up/s))/(up*up));
  double Ms(-3.0/8.0*(tp*tp+up*up+4.0*m_m32*s)/(s*s));
  return sqr(m_g*m_g)*(Mt+Mu+Ms); 
}


bool XS_gg_q1qbar1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double tp(t-m_m32), up(u-m_m32);
  double Mt(1.0/6.0*(tp*up-m_m32*(4.0*(m_m32+tp)+m_m32*tp/s))/(tp*tp));
  double Mu(1.0/6.0*(up*tp-m_m32*(4.0*(m_m32+up)+m_m32*up/s))/(up*up));
  p_colours[0][0] = Flow::Counter();
  p_colours[0][1] = Flow::Counter();
  if (Mt*(1-m_r) +Mu*m_r > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: gg->qqb, set scale t/s "<<t<<"/"<<s<<"\n";
    /*
      0====+------2
           |
           |  t
           |
      1====+------3

      shower scale is s / t
    */
    p_colours[2+m_r][0] = p_colours[0][0];
    p_colours[3-m_r][1] = p_colours[1][1] = Flow::Counter();
    p_colours[1][0] = p_colours[0][1];
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: gg->qqb, set scale u/s "<<u<<"/"<<s<<"\n";
    /*
      0==-+ /----2
         \|/
          |  u
         /|\
      1==-+ \----3

      shower scale is u / s
    */
    p_colours[2+m_r][0] = p_colours[1][0] = Flow::Counter();
    p_colours[3-m_r][1] = p_colours[0][1];
    p_colours[1][1] = p_colours[0][0];
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_q1g_q1g>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (((flavours[0].IsQuark() && flavours[1].IsGluon()) && 
	 ((flavours[2]==flavours[0] && flavours[3].IsGluon()) || 
	  (flavours[3]==flavours[0] && flavours[2].IsGluon()))) ||
	((flavours[1].IsQuark() && flavours[0].IsGluon()) && 
	 ((flavours[2]==flavours[1] && flavours[3].IsGluon()) || 
	  (flavours[3]==flavours[1] && flavours[2].IsGluon()))))  { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1g_q1g(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_q1g_q1g::XS_q1g_q1g(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl,
		       XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_iniq=0;
  m_swaput=0;
  if (fl[1].IsQuark()){
    m_iniq=1;
    m_swaput=1;
  }
  m_finq=2;
  if (fl[3].IsQuark()) {
    m_finq=3;
    if (m_swaput) m_swaput=0;
    else m_swaput=1;
  }
  m_a=fl[m_iniq].IsAnti();
  m_p=1-m_a;
  m_mq2=sqr(p_flavours[m_iniq].Mass());
  m_g=std::abs(p_model->Constant("g_3"));
  m_nstrong=4;
}

double XS_q1g_q1g::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  if (m_swaput) std::swap<double>(t,u);
  double sp(s-m_mq2), up(u-m_mq2);
  double Ms(4.0/9.0*(sp*up-m_mq2*(4.0*(m_mq2+sp)+m_mq2*sp/t))/(sp*sp));
  double Mu(4.0/9.0*(up*sp-m_mq2*(4.0*(m_mq2+up)+m_mq2*up/t))/(up*up));
  double Mt(-(sp*sp+up*up+4.0*m_mq2*t)/(t*t));
  return -sqr(m_g*m_g)*(Ms+Mu+Mt); 
}

bool XS_q1g_q1g::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  if (m_swaput) std::swap<double>(t,u);
  double sp(s-m_mq2), up(u-m_mq2);
  double Ms(4.0/9.0*(sp*up-m_mq2*(4.0*(m_mq2+sp)+m_mq2*sp/t))/(sp*sp));
  double Mu(4.0/9.0*(up*sp-m_mq2*(4.0*(m_mq2+up)+m_mq2*up/t))/(up*up));
  p_colours[m_iniq][m_a] = Flow::Counter();
  p_colours[m_finq][m_a] = Flow::Counter();
  if (Mu > (Ms+Mu) * ran.Get()) {
    /*
      1====+----2, if fl[2].IsQuark() 
           |
           |  u
           |
      0----+====3, if fl[0].IsQuark()

      shower scale is t/u
    */
    p_colours[5-m_finq][m_a] = p_colours[m_iniq][m_a];
    p_colours[5-m_finq][m_p] = p_colours[1-m_iniq][m_p] = Flow::Counter();
    p_colours[1-m_iniq][m_a] = p_colours[m_finq][m_a];
    if (dabs(t)>dabs(u)) {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale t "<<t<<"\n";
    }
    else {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale u "<<u<<"\n";
    }
  }
  else {
    /*
      0\        /2, if fl[0].IsQuark && fl[2].IsQuark() 
        \      /
        | +--+ | 
        //    \\
      1//      \\3

      shower scale is s/t
    */
    p_colours[5-m_finq][m_p] = p_colours[m_finq][m_a];
    p_colours[1-m_iniq][m_a] = p_colours[5-m_finq][m_a] = Flow::Counter();
    p_colours[1-m_iniq][m_p] = p_colours[m_iniq][m_a];
    if (dabs(t)>s) {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale t "<<t<<"\n";
    }
    else {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale s "<<s<<"\n";
    }
  }
  if (swap) SwapInOrder();
  return true;
}
      
namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_gg_gg>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("QCD")) return NULL;
    if (flavours[0].IsGluon() && flavours[1].IsGluon() &&
	flavours[2].IsGluon() && flavours[3].IsGluon()) { 
      if (nqcd==2 && nqed==0) {
	return new XS_gg_gg(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }

}

XS_gg_gg::XS_gg_gg(const size_t nin,const size_t nout, 
		   const ATOOLS::Flavour *fl,
		   XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_g=std::abs(p_model->Constant("g_3"));
  m_nstrong=4;
}

double XS_gg_gg::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double Ms(1.0-t*u/(s*s));
  double Mt(1.0-s*u/(t*t));
  double Mu(1.0-s*t/(u*u));
  return sqr(m_g*m_g)*9.0/2.0*(Ms+Mt+Mu)/2.0;
}
  
bool XS_gg_gg::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
  msg_Debugging()<<"xs: gg->gg, set scale s "<<s<<"\n";
  p_colours[0][0] = Flow::Counter();
  p_colours[1][1] = Flow::Counter();
  double Mu(1.0+t*t/(u*s)-s*t/(u*u)-t*u/(s*s));
  double Ms(1.0+s*s/(t*u)-s*t/(u*u)-u*s/(t*t));
  double Mt(1.0+u*u/(s*t)-u*s/(t*t)-t*u/(s*s));
  double rr = ran.Get() * (Ms+Mt+Mu);
  if (rr-Mt < 0.) {
    /*
      0====++====2
           ||
           ||  t
           ||
      1====++====3

      shower scale is s/t/u
    */
    p_colours[2][0] = p_colours[0][0];
    p_colours[3][1] = p_colours[1][1];
    p_colours[0][1] = p_colours[1][0] = Flow::Counter();
    p_colours[2][1] = p_colours[3][0] = Flow::Counter();
  }
  else {
    if (rr-Mu-Mt < 0.) {
      /*
	0====+\---==3
             ||\ /
             || X u
             ||/ \
	1====+/---==2
	   
	shower scale is s/t/u
      */
      p_colours[3][0] = p_colours[0][0];
      p_colours[2][1] = p_colours[1][1];
      p_colours[0][1] = p_colours[1][0] = Flow::Counter();
      p_colours[3][1] = p_colours[2][0] = Flow::Counter();
    }
    else {
      /*
	0\\       //3
          \\  s  //
          | ===== |
          //     \\
	1//       \\2
	   
	shower scale is s/t/u
      */
      p_colours[2][0] = p_colours[0][0];
      p_colours[3][1] = p_colours[0][1] = Flow::Counter();
      p_colours[2][1] = p_colours[1][1];
      p_colours[3][0] = p_colours[1][0] = Flow::Counter();
    }
  }
  if (swap) SwapInOrder();
  return true;
}
    
