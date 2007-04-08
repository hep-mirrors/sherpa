#include "XS_QCD.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Flow.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;


/* 
   In all the differential cross sections the factor 1/16 Pi^2 is cancelled
   by the factor 4 Pi for each alpha
*/

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_pp_ffbar>(const size_t nin,const size_t nout,
					      const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (nqed!=2 || nqcd!=0) return NULL;
  if (flavours[2].IsFermion() && flavours[3]==flavours[2].Bar() &&
      flavours[0].IsPhoton()  && flavours[1]==flavours[0]) { 
    return new XS_pp_ffbar(nin,nout,flavours); 
  }
  return NULL;
}

}

XS_pp_ffbar::XS_pp_ffbar(const size_t nin,const size_t nout,
			 const ATOOLS::Flavour *fl) :
  Single_XS(nin,nout,fl)  
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (p_flavours[2].Strong()) {
    int a  = p_flavours[2].IsAnti();
    p_colours[2][a] = 500;
    p_colours[3][1-a] = 500;
  }
}
double XS_pp_ffbar::operator()(double s,double t,double u) 
{ 
  return 0.0; 
}

bool XS_pp_ffbar::SetColours(double s,double t,double u) 
{ 
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = s;
  return 1; 
}

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_q1q2_q1q2>(const size_t nin,const size_t nout,
					       const ATOOLS::Flavour *flavours, 
					       const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsQuark() && flavours[1].IsQuark() && flavours[0]!=flavours[1] &&
      ((flavours[2]==flavours[0] && flavours[3]==flavours[1]) ||
       (flavours[3]==flavours[0] && flavours[2]==flavours[1]))) { 
    if (nqcd==2 && nqed==0) {
      return new XS_q1q2_q1q2(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_q1q2_q1q2::XS_q1q2_q1q2(const size_t nin,const size_t nout, const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  a  = fl[0].IsAnti();
  p = fl[1].IsAnti();
  aS = (*as)(sqr(rpa.gen.Ecms()));

  m_nstrong=4;
}

double XS_q1q2_q1q2::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  return sqr(4.*M_PI*aS)* 4. * (s*s + u*u) / ( 9. * t*t);
}

bool XS_q1q2_q1q2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 2.0*(s*t*u)/(s*s+t*t+u*u);
  int r = !(p_flavours[0] == p_flavours[2]);
  if (a==p) {
    /*
      0-----\   /-----2, if fl[0]==fl[2]
             \ /
	      X  t
             / \
      1-----/   \-----3, if fl[1]==fl[3]

      shower scale is u
    */
    p_colours[0][a] = p_colours[3-r][a] = Flow::Counter();
    p_colours[1][a] = p_colours[2+r][a] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
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
    p_colours[0][a]   = p_colours[1][p]   = Flow::Counter();
    p_colours[2+r][a] = p_colours[3-r][p] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb'->qqb', set scale s "<<s<<"\n";
  }
  if (swap) SwapInOrder();
  return 1; 
}
bool XS_q1q2_q1q2::SetColours()                           { return 1; }

//----------------------------------------------------------------------

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_q1qbar1_q2qbar2>(const size_t nin,const size_t nout,
						     const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
      flavours[2].IsQuark() && flavours[3]==flavours[2].Bar() &&
      flavours[0]!=flavours[2]) { 
    if (nqcd==2 && nqed==0) {
      return new XS_q1qbar1_q2qbar2(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_q1qbar1_q2qbar2::XS_q1qbar1_q2qbar2(const size_t nin,const size_t nout, 
				       const ATOOLS::Flavour *fl)  : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a = fl[0].IsAnti();
  p = 1-a;

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_q1qbar1_q2qbar2::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  return sqr(4.*M_PI*aS)* 4. * (t*t + u*u) / ( 9. * s*s); 
}

bool XS_q1qbar1_q2qbar2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 2.0*(s*t*u)/(s*s+t*t+u*u);
  int r = !(p_flavours[0].IsAnti() == p_flavours[2].IsAnti());
  /*
    0\         /2, if fl[0].IsAnti()==fl[2].IsAnti()
      \   s   /
       ======= 
      /       \
    1/         \3, if fl[0].IsAnti()==fl[2].IsAnti()

    shower scale is t
  */
  p_colours[0][a] = p_colours[2+r][a] = Flow::Counter();
  p_colours[1][p] = p_colours[3-r][p] = Flow::Counter();

  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->q'qb', set scale t "<<t<<"\n";
  if (swap) SwapInOrder();
  return 1; 
}

bool XS_q1qbar1_q2qbar2::SetColours()                           { return 1; }

//----------------------------------------------------------------------
// Note : Combinatorical factor of 2 for identical outgoing particles explicitly added

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_q1q1_q1q1>(const size_t nin,const size_t nout,
					       const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsQuark() && flavours[1]==flavours[0] &&
      flavours[2]==flavours[0] && flavours[3]==flavours[0]) { 
    if (nqcd==2 && nqed==0) {
      return new XS_q1q1_q1q1(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_q1q1_q1q1::XS_q1q1_q1q1(const size_t nin,const size_t nout, const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a  = fl[0].IsAnti();

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_q1q1_q1q1::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  Mt    = (s*s + u*u) / (t*t);
  Mu    = (s*s + t*t) / (u*u);
  return sqr(4.*M_PI*aS) * 4./9.*(Mt + Mu - 2./3. * (s*s) / (u*t)) /2.;
}

bool XS_q1q1_q1q1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 2.0*(s*t*u)/(s*s+t*t+u*u);
  
  Mt      = 1. - 2.*(u*s) / (t*t);
  Mu      = 1. - 2.*(s*t) / (u*u);
  
  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->qq, set scale u "<<u<<"\n";
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->qq, set scale t "<<t<<"\n";
  }
  bool result=SetColours();
  if (swap) SwapInOrder();
  return result;
}


bool XS_q1q1_q1q1::SetColours() 
{
  if (Mt > (Mt+Mu) * ran.Get()) {
    /*
      0----\   /----2
            \ /
             X  t
            / \
      1----/   \----3

      shower scale is u
    */
    p_colours[3][a] = p_colours[0][a] = Flow::Counter();
    p_colours[2][a] = p_colours[1][a] = Flow::Counter();
  }
  else {
    /*
      0----\   /----2
            \ /
             =  u
            / \
      1----/   \----3

      shower scale is t
    */
    p_colours[2][a] = p_colours[0][a] = Flow::Counter();
    p_colours[3][a] = p_colours[1][a] = Flow::Counter();
  }
  return 1;
}

//----------------------------------------------------------------------

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_q1qbar1_q1qbar1>(const size_t nin,const size_t nout,
						     const ATOOLS::Flavour *flavours, 
						     const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
      ((flavours[2]==flavours[0] && flavours[3]==flavours[1]) ||
       (flavours[3]==flavours[0] && flavours[2]==flavours[1]))) { 
    if (nqcd==2 && nqed==0) {
      return new XS_q1qbar1_q1qbar1(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_q1qbar1_q1qbar1::XS_q1qbar1_q1qbar1(const size_t nin,const size_t nout, 
				       const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a  = fl[0].IsAnti();
  p  = 1-a;
  r = !(fl[0] == fl[2]);

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_q1qbar1_q1qbar1::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  Mt = 1. - 2.*(u*s)/(t*t); 
  Ms = 1. - 2.*(t*u)/(s*s); 
  return sqr(4.*M_PI*aS)*4./9.*(Mt + Ms - 2./3. * (u*u) / (s*t));
}


bool XS_q1qbar1_q1qbar1::SetColours(double s, double t, double u) {
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 2.0*(s*t*u)/(s*s+t*t+u*u);

  Mt = 1. - 2.*(u*s)/(t*t); 
  Ms = 1. - 2.*(t*u)/(s*s); 

  if (Ms >  (Mt+Ms) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->qqb, set scale t "<<t<<"\n";
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->qqb, set scale s "<<s<<"\n";
  }
  bool result=SetColours();
  if (swap) SwapInOrder();
  return result;
}

bool XS_q1qbar1_q1qbar1::SetColours() 
{
  if (Ms >  (Mt+Ms) * ran.Get()) {
    /*
      0\         /2, if fl[0]==fl[2]
        \   s   /
         =======
        /       \
      1/         \3, if fl[0]==fl[2]

      shower scale is t
    */
    p_colours[0][a] = p_colours[2+r][a] = Flow::Counter();	
    p_colours[1][p] = p_colours[3-r][p] = Flow::Counter();
  }
  else {
    /*
      0----+ +----2
           | |
           | | t
           | |
      1----+ +----3

      shower scale is s
    */
    p_colours[0][a]   = p_colours[1][p]   = Flow::Counter();	
    p_colours[2+r][a] = p_colours[3-r][p] = Flow::Counter();
  }
  return 1;
}

//----------------------------------------------------------------------
// Note : Combinatorical factor of 2 for identical outgoing particles explicitly added

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_q1qbar1_gg>(const size_t nin,const size_t nout,
						const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
      flavours[2].IsGluon() && flavours[3].IsGluon()) { 
    if (nqcd==2 && nqed==0) {
      return new XS_q1qbar1_gg(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_q1qbar1_gg::XS_q1qbar1_gg(const size_t nin,const size_t nout, 
			     const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a = fl[0].IsAnti();
  p = 1-a;

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_q1qbar1_gg::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  Mt = u/t;
  Mu = t/u;
  return sqr(4.*M_PI*aS)* (32./27.*(Mt+Mu) - 8./3.*(t*t +u*u)/ (s*s)) /2.;
}


bool XS_q1qbar1_gg::SetColours(double s, double t, double u) {
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 2.0*(s*t*u)/(s*s+t*t+u*u);

  Mt    = u/t;
  Mu    = t/u;

  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->gg, set scale s/t "<<s<<"/"<<t<<"\n";
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->gg, set scale s/u "<<s<<"/"<<u<<"\n";
  }
  bool result=SetColours();
  if (swap) SwapInOrder();
  return result;
}

bool XS_q1qbar1_gg::SetColours() 
{
  p_colours[0][a] = Flow::Counter();
  p_colours[1][p] = Flow::Counter();

  if (Mt > (Mt+Mu) * ran.Get()) {
    /*
      0------+====2
             |
             | t
             |
      1------+====3

      shower scale is s / t
    */
    p_colours[2][a] = p_colours[0][a];
    p_colours[3][p] = p_colours[1][p];
    p_colours[2][p] = p_colours[3][a] = Flow::Counter();
  }
  else {
    /*
      0----\ +-==2
            \|/
             | u
            /|\
      1----/ +-==3

      shower scale is s / u
    */
    p_colours[3][a] = p_colours[0][a];
    p_colours[2][p] = p_colours[1][p];
    p_colours[3][p] = p_colours[2][a] = Flow::Counter();
  }
  return 1;
}

//----------------------------------------------------------------------

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_gg_q1qbar1>(const size_t nin,const size_t nout,
						const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsGluon() && flavours[1].IsGluon() && 
      flavours[2].IsQuark() && flavours[3]==flavours[2].Bar()) { 
    if (nqcd==2 && nqed==0) {
      return new XS_gg_q1qbar1(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_gg_q1qbar1::XS_gg_q1qbar1(const size_t nin,const size_t nout, const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  r = fl[2].IsAnti();

  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_gg_q1qbar1::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  Mt = u/t;
  Mu = t/u;
  return sqr(4.*M_PI*aS)*(1./6.* ( Mt + Mu )   - 3./8. *(t*t +u*u)/ (s*s)); 
}


bool XS_gg_q1qbar1::SetColours(double s, double t, double u) {
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = (2.*s*t*u)/(s*s+t*t+u*u);

  Mt      = u/t;
  Mu      = t/u;

  if (Mt*(1-r) +Mu*r > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: gg->qqb, set scale t/s "<<t<<"/"<<s<<"\n";
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: gg->qqb, set scale u/s "<<u<<"/"<<s<<"\n";
  }
  bool result=SetColours();
  if (swap) SwapInOrder();
  return result;
}

bool XS_gg_q1qbar1::SetColours() 
{
  p_colours[0][0] = Flow::Counter();
  p_colours[0][1] = Flow::Counter();

  if (Mt*(1-r) +Mu*r > (Mt+Mu) * ran.Get()) {
    /*
      0====+------2
           |
           |  t
           |
      1====+------3

      shower scale is s / t
    */
    p_colours[2+r][0] = p_colours[0][0];
    p_colours[3-r][1] = p_colours[1][1] = Flow::Counter();
    p_colours[1][0] = p_colours[0][1];
  }
  else {
    /*
      0==-+ /----2
         \|/
          |  u
         /|\
      1==-+ \----3

      shower scale is u / s
    */
    p_colours[2+r][0] = p_colours[1][0] = Flow::Counter();
    p_colours[3-r][1] = p_colours[0][1];
    p_colours[1][1] = p_colours[0][0];
  }
  return 1;
}

//----------------------------------------------------------------------


namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_q1g_q1g>(const size_t nin,const size_t nout,
					     const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (((flavours[0].IsQuark() && flavours[1].IsGluon()) && 
       ((flavours[2]==flavours[0] && flavours[3].IsGluon()) || 
	(flavours[3]==flavours[0] && flavours[2].IsGluon()))) ||
      ((flavours[1].IsQuark() && flavours[0].IsGluon()) && 
       ((flavours[2]==flavours[1] && flavours[3].IsGluon()) || 
	(flavours[3]==flavours[1] && flavours[2].IsGluon()))))  { 
    if (nqcd==2 && nqed==0) {
      return new XS_q1g_q1g(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_q1g_q1g::XS_q1g_q1g(const size_t nin,const size_t nout, const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  ini_q=0;
  swap_ut=0;
  if (fl[1].IsQuark()){
    ini_q=1;
    swap_ut=1;
  }

  fin_q=2;
  if (fl[3].IsQuark()) {
    fin_q=3;
    if (swap_ut) swap_ut=0;
    else swap_ut=1;
  }

  a = fl[ini_q].IsAnti();
  p = 1-a;

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_q1g_q1g::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  if (swap_ut) {
    Ms = t/s;
    Mu = s/t;
    return  sqr(4.*M_PI*aS)*(-4./9. * (Ms + Mu) +  (s*s + t*t)/(u*u));
  }
  Ms = u/s;
  Mu = s/u;
  return  sqr(4.*M_PI*aS)*(-4./9. * (Ms + Mu) +  (s*s + u*u)/(t*t));
}

bool XS_q1g_q1g::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = (2.*s*t*u)/(s*s+t*t+u*u);
  if (swap_ut) {
    Ms      = t/s;
    Mu      = s/t;
  }
  else {
    Ms      = u/s;
    Mu      = s/u;
  }
  if (Mu > (Ms+Mu) * ran.Get()) {
    if (dabs(t)>dabs(u)) {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale t "<<t<<"\n";
    }
    else {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale u "<<u<<"\n";
    }
  }
  else {
    if (dabs(t)>s) {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale t "<<t<<"\n";
    }
    else {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->qg, set scale s "<<s<<"\n";
    }
  }
  bool result=SetColours();
  if (swap) SwapInOrder();
  return result;
}
      
bool XS_q1g_q1g::SetColours() 
{
  p_colours[ini_q][a] = Flow::Counter();
  p_colours[fin_q][a] = Flow::Counter();

  if (Mu > (Ms+Mu) * ran.Get()) {
    /*
      1====+----2, if fl[2].IsQuark() 
           |
           |  u
           |
      0----+====3, if fl[0].IsQuark()

      shower scale is t/u
    */
    p_colours[5-fin_q][a] = p_colours[ini_q][a];
    p_colours[5-fin_q][p] = p_colours[1-ini_q][p] = Flow::Counter();
    p_colours[1-ini_q][a] = p_colours[fin_q][a];
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
    p_colours[5-fin_q][p] = p_colours[fin_q][a];
    p_colours[1-ini_q][a] = p_colours[5-fin_q][a] = Flow::Counter();
    p_colours[1-ini_q][p] = p_colours[ini_q][a];
  }

  return 1;
}

//----------------------------------------------------------------------
// Note : Combinatorical factor of 2 for identical outgoing particles explicitly added

namespace EXTRAXS {

template <> 
Single_XS *Single_XS::GetProcess<XS_gg_gg>(const size_t nin,const size_t nout,
					   const ATOOLS::Flavour *flavours, 
					      const size_t nqed, const size_t nqcd)
{
  if (flavours[0].IsGluon() && flavours[1].IsGluon() &&
      flavours[2].IsGluon() && flavours[3].IsGluon()) { 
    if (nqcd==2 && nqed==0) {
      return new XS_gg_gg(nin,nout,flavours); 
    }
  }
  return NULL;
}

}

XS_gg_gg::XS_gg_gg(const size_t nin,const size_t nout, const ATOOLS::Flavour *fl) : 
  Single_XS(nin,nout,fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  aS = (*as)(sqr(rpa.gen.Ecms()));
  m_nstrong=4;
}

double XS_gg_gg::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  Ms = 1 - t*u/(s*s);
  Mt = 1 - s*u/(t*t);
  Mu = 1 - s*t/(u*u);
  return sqr(4.*M_PI*aS)*9./2. * ( Ms + Mt + Mu )/2.;
}
  
bool XS_gg_gg::SetColours(double s, double t, double u) {
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = (2.*s*t*u)/(s*s+t*t+u*u);
  
  Mu      = 1 + t*t/(u*s) - s*t/(u*u) - t*u/(s*s);
  Ms      = 1 + s*s/(t*u) - s*t/(u*u) - u*s/(t*t);
  Mt      = 1 + u*u/(s*t) - u*s/(t*t) - t*u/(s*s);

  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
  msg_Debugging()<<"xs: gg->gg, set scale s "<<s<<"\n";

  bool result=SetColours();
  if (swap) SwapInOrder();
  return result;
}
    
bool XS_gg_gg::SetColours() 
{
  p_colours[0][0] = Flow::Counter();
  p_colours[1][1] = Flow::Counter();

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
  return 1;
}





