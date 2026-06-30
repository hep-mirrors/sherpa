#include "AMISIC++/Perturbative/QCD_Processes.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/LDME.H"
#include "EXTRA_XS/Two2Two/XS_Quarkonia.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

//////////////////////////////////////////////////////////////////////////
//
// All colour distributions below according to Webber & Marchesini
//
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// gg -> gg
//////////////////////////////////////////////////////////////////////////

gg_gg::gg_gg() : XS_Base() { m_name = string("gg->gg"); }

void gg_gg::Calc(const double & s,const double & t,const double & u) {
  m_Ms = 9./4. * (1. + s*s/(t*u) - s*t/(u*u) - u*s/(t*t));
  m_Mt = 9./4. * (1. + u*u/(s*t) - u*s/(t*t) - t*u/(s*s));
  m_Mu = 9./4. * (1. + t*t/(u*s) - s*t/(u*u) - t*u/(s*s));
  m_lastxs = (m_Mt+m_Mu+m_Ms)/2.;
}

bool gg_gg::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  m_colours[0][0] = Flow::Counter();
  m_colours[1][1] = Flow::Counter();
  double rr = ran->Get() * (m_Ms+m_Mt+m_Mu);
  if (rr-m_Mt<0.) {
    /*
      
      0====++====2
           ||
           ||  t
           ||
      1====++====3
      
    */
    m_colours[1][0] = Flow::Counter();
    m_colours[3][0] = Flow::Counter();
    m_colours[2][0] = m_colours[0][0];
    m_colours[3][1] = m_colours[1][1];
    m_colours[0][1] = m_colours[1][0];
    m_colours[2][1] = m_colours[3][0];
  }
  else if (rr-m_Mt-m_Mu<0.) {
    /*
      
	0====++====3
             ||
             ||
             ||
	1====++====2
	   
    */
    m_colours[1][0] = Flow::Counter();
    m_colours[2][0] = Flow::Counter();
    m_colours[3][0] = m_colours[0][0];
    m_colours[2][1] = m_colours[1][1];
    m_colours[0][1] = m_colours[1][0];
    m_colours[3][1] = m_colours[2][0];
  }
  else {
    /*
    
	0====++====2
             ||
             ||
             ||
	3====++====1
	   
    */
    m_colours[0][1] = Flow::Counter();
    m_colours[1][0] = Flow::Counter();
    m_colours[2][0] = m_colours[0][0];
    m_colours[3][1] = m_colours[0][1];
    m_colours[2][1] = m_colours[1][1];
    m_colours[3][0] = m_colours[1][0];
  }
  return true; 
}

DECLARE_XSBASE_GETTER(gg_gg,"XS_gg_gg")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,gg_gg>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsGluon() && flavs[1].IsGluon() &&
      flavs[2].IsGluon() && flavs[3].IsGluon()) {
    return new gg_gg();
  }
  return NULL;
}


//////////////////////////////////////////////////////////////////////////
// gg -> qqbar
//////////////////////////////////////////////////////////////////////////

gg_qqbar::gg_qqbar() : XS_Base() { m_name = string("gg->qqbar"); }

void gg_qqbar::Calc(const double & s,const double & t,const double & u) {
  double help = 3./16. * ((t*t+u*u)/(s*s) - 1./9.);
  // h^c(s,t,u) + h^c(s,u,t)
  m_Mt = u/t*help;
  m_Mu = t/u*help;
  m_Ms = 0.;
  m_lastxs = (m_Mt+m_Mu+m_Ms);
}

bool gg_qqbar::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  m_colours[0][0] = Flow::Counter();
  m_colours[0][1] = Flow::Counter();
  if ((m_Mt+m_Mu)*ran->Get()-m_Mt<0.) {
    /*
    
      0====+------2
           |
           |  t
           |
      1====+------3

    */
    m_colours[1][1] = Flow::Counter();
    m_colours[2][0] = m_colours[0][0];
    m_colours[3][1] = m_colours[1][1];
    m_colours[1][0] = m_colours[0][1];
  }
  else {
    /*

      0====+------3
           |
           |  u
           |
      1====+------2

    */
    m_colours[1][0] = Flow::Counter();
    m_colours[2][0] = m_colours[1][0];
    m_colours[3][1] = m_colours[0][1];
    m_colours[1][1] = m_colours[0][0];
  }
  m_colours[2][1] = m_colours[3][0] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(gg_qqbar,"XS_gg_qqbar")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,gg_qqbar>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsGluon() && flavs[1].IsGluon() &&
      flavs[2].IsQuark() && !flavs[2].IsAnti() && 
      flavs[3]==flavs[2].Bar()) {
    if (flavs[2].Mass()<=1.e-6) return new gg_qqbar();
    THROW(fatal_error,"Massive quarks not yet enabled in gg -> QQbar");
  }
  return NULL;
}



//////////////////////////////////////////////////////////////////////////
// qqbar -> gg
//////////////////////////////////////////////////////////////////////////
qqbar_gg::qqbar_gg() : XS_Base() { m_name = string("qqbar_gg"); }

void qqbar_gg::Calc(const double & s,const double & t,const double & u) {
  double help = 4./3. * ((t*t+u*u)/(s*s) - 1./9.);
  // h^c(s,t,u) + h^c(s,u,t)
  m_Mt = u/t*help;
  m_Mu = t/u*help;
  m_Ms = 0.;
  m_lastxs = (m_Mt+m_Mu+m_Ms)/2.;
}

bool qqbar_gg::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t anti = size_t(flavs[0].IsAnti()), part = 1-anti; 
  m_colours[0][anti] = Flow::Counter();
  m_colours[1][part] = Flow::Counter();
  m_colours[0][part] = m_colours[1][anti] = 0;
  if ((m_Mt+m_Mu) * ran->Get() < m_Mt) {
     /*
    
      0------+====2
             |
             | t
             |
      1------+====3

    */
    m_colours[3][anti] = Flow::Counter();
    m_colours[2][anti] = m_colours[0][anti];
    m_colours[3][part] = m_colours[1][part];
    m_colours[2][part] = m_colours[3][anti];
  }
  else {
     /*

      0------+====3
             |
             | u
             |
      1------+====2

    */
    m_colours[2][anti] = Flow::Counter();
    m_colours[3][anti] = m_colours[0][anti];
    m_colours[2][part] = m_colours[1][part];
    m_colours[3][part] = m_colours[2][anti];
  }
  return true;
}

DECLARE_XSBASE_GETTER(qqbar_gg,"XS_qqbar_gg")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,qqbar_gg>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[0]==flavs[1].Bar() &&
      flavs[2]==Flavour(kf_gluon) && flavs[3]==Flavour(kf_gluon)) {
    if (flavs[0].Mass()<1.e-6) return new qqbar_gg();
    THROW(fatal_error,"Massive quarks not yet enabled in QQbar->gg");
  }
  return NULL;
}



//////////////////////////////////////////////////////////////////////////
// qg -> qg
//////////////////////////////////////////////////////////////////////////

qg_qg::qg_qg() : XS_Base() { m_name = string("qg->qg"); }

void qg_qg::Calc(const double & s,const double & t,const double & u) {
  // h^c(t,s,u) + h^c(t,u,s)
  m_Mt = -4./3. * u/s * ((s*s+u*u)/(t*t) - 1./9.);
  m_Mu = -4./3. * s/u * ((s*s+u*u)/(t*t) - 1./9.);    
  m_Ms = 0.;
  m_lastxs = 3./8.*(m_Mt+m_Mu+m_Ms);
}

bool qg_qg::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t quark = flavs[0].IsQuark()?0:1, gluon = 1-quark;
  size_t trip = flavs[quark].IsAnti()?1:0, anti=1-trip;
  m_colours[quark][trip]   = Flow::Counter();
  m_colours[quark+2][trip] = Flow::Counter();
  m_colours[quark][anti]   = m_colours[quark+2][anti] = 0;
  if ((m_Mt+m_Mu) * ran->Get() < m_Mt) {
     /*
    
      0------++----1
             ||
             || t
             ||
      2======++====3

    */
    m_colours[gluon+2][trip] = Flow::Counter();
    m_colours[gluon][anti]   = m_colours[quark][trip];
    m_colours[gluon+2][anti] = m_colours[quark+2][trip];
    m_colours[gluon][trip]   = m_colours[gluon+2][trip];
  }
  else {
     /*

      0------+====3
             |
             | u
             |
      1======+----2

    */
    m_colours[gluon+2][anti] = Flow::Counter();
    m_colours[gluon+2][trip] = m_colours[quark][trip];
    m_colours[gluon][trip]   = m_colours[quark+2][trip];
    m_colours[gluon][anti]   = m_colours[gluon+2][anti];
  }
  return true;
}

DECLARE_XSBASE_GETTER(qg_qg,"XS_qg_qg")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,qg_qg>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[0]==flavs[2] &&
      flavs[1]==Flavour(kf_gluon) && flavs[3]==Flavour(kf_gluon)) {
    if (flavs[0].Mass()<1.e-6) return new qg_qg();
    THROW(fatal_error,"Massive quarks not yet enabled in Qg->Qg");
  }
  else if (flavs[1].IsQuark() && flavs[1]==flavs[3] &&
	   flavs[0]==Flavour(kf_gluon) && flavs[2]==Flavour(kf_gluon)) {
    if (flavs[1].Mass()<1.e-6) return new qg_qg();
    THROW(fatal_error,"Massive quarks not yet enabled in gQ->gQ");
  }
  return NULL;
}



//////////////////////////////////////////////////////////////////////////
// qq -> qq
//////////////////////////////////////////////////////////////////////////

qq_qq::qq_qq() : XS_Base() { m_name = string("qq->qq"); }

void qq_qq::Calc(const double & s,const double & t,const double & u) {
  // h^b(s,t,u) + h^b(s,u,t)
  m_Mt     = 4./9. * (s*s+t*t)/(u*u) + 8./27. * s/u;
  m_Mu     = 4./9. * (s*s+u*u)/(t*t) + 8./27. * s/t;
  m_Ms     = 0.;
  m_lastxs = (m_Ms + m_Mt + m_Mu)/2.;
}

bool qq_qq::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t anti(flavs[0].IsAnti());
  for (size_t i=0;i<4;i++) m_colours[i][1-anti] = 0;
  m_colours[3][anti] = Flow::Counter();
  m_colours[2][anti] = Flow::Counter();
  if ((m_Mt+m_Mu) * ran->Get() < m_Mt) {
     /*
    
      0------++-----2
             ||
             || t
             ||
      1------++-----3

    */
    m_colours[0][anti] = m_colours[3][anti];
    m_colours[1][anti] = m_colours[2][anti];
  }
  else {
     /*

      0------++-----3
             ||
             || u
             ||
      1------++-----2

    */
    m_colours[0][anti] = m_colours[2][anti];
    m_colours[1][anti] = m_colours[3][anti];
  }
  return true;
}

DECLARE_XSBASE_GETTER(qq_qq,"XS_qq_qq")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,qq_qq>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0] &&
      flavs[2]==flavs[0] && flavs[3]==flavs[0]) {
    if (flavs[0].Mass()<1.e-6) return new qq_qq();
    THROW(fatal_error,"Massive quarks not yet enabled in QQ->QQ");
  }
  return NULL;
}



//////////////////////////////////////////////////////////////////////////
// qqbar -> qqbar
//////////////////////////////////////////////////////////////////////////

qqbar_qqbar::qqbar_qqbar() : XS_Base() { m_name = string("qqbar->qqbar"); }

void qqbar_qqbar::Calc(const double & s,const double & t,const double & u) {
  m_Ms     = 4./9. * (t*t+u*u)/(s*s) - 8./27. * u/s;
  m_Mt     = 4./9. * (s*s+u*u)/(t*t) - 8./27. * u/t;
  m_Mu     = 0.;
  m_lastxs = m_Ms + m_Mt + m_Mu;
}

bool qqbar_qqbar::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t anti0 = size_t(flavs[0].IsAnti()), anti1 = 1-anti0;
  if (m_Ms >  (m_Mt+m_Ms) * ran->Get()) {
    /*
    
      0\         /2, if fl[0]==fl[2]
        \   s   /
         =======
        /       \
      1/         \3, if fl[0]==fl[2]

    */
    m_colours[2][anti0] = Flow::Counter();
    m_colours[3][anti1] = Flow::Counter();
    m_colours[0][anti0] = m_colours[2][anti0];
    m_colours[1][anti1] = m_colours[3][anti1];
  }
  else {
    /*

      0----+ +----2
           | |
           | | t
           | |
      1----+ +----3

    */
    m_colours[1][anti1] = Flow::Counter();
    m_colours[3][anti1] = Flow::Counter();
    m_colours[0][anti0] = m_colours[1][anti1];
    m_colours[2][anti0] = m_colours[3][anti1];
  }
  m_colours[0][1-anti0] = m_colours[2][1-anti0] = 0;
  m_colours[1][1-anti1] = m_colours[3][1-anti1] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(qqbar_qqbar,"XS_qqbar_qqbar")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,qqbar_qqbar>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0].Bar() &&
      flavs[2]==flavs[0] && flavs[3]==flavs[2].Bar()) {
    if (flavs[0].Mass()<1.e-6) return new qqbar_qqbar();
    THROW(fatal_error,"Massive quarks not yet enabled in QQbar->QQbar");
  }
  return NULL;
}



//////////////////////////////////////////////////////////////////////////
// q1q2 -> q1q2
//////////////////////////////////////////////////////////////////////////

q1q2_q1q2::q1q2_q1q2() : XS_Base() { m_name = string("q1q2_q1q2"); }

void q1q2_q1q2::Calc(const double & s,const double & t,const double & u) {
  m_lastxs = 4./9. * (s*s+u*u)/(t*t);
}

bool q1q2_q1q2::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t anti0(flavs[0].IsAnti()), anti1(flavs[1].IsAnti());
  if (anti0==anti1) {
    // particle--particle or anti-particle--anti-particle: u channel colours
    m_colours[3][anti0] = Flow::Counter();
    m_colours[2][anti1] = Flow::Counter();
    m_colours[0][anti0] = m_colours[3][anti0];
    m_colours[1][anti1] = m_colours[2][anti1];
    m_colours[0][1-anti0] = m_colours[3][1-anti0] = 0;
    m_colours[1][1-anti1] = m_colours[2][1-anti1] = 0;
  }
  else {
    // particle--anti-particle: t-channel colours
    m_colours[1][anti1] = Flow::Counter();
    m_colours[3][anti1] = Flow::Counter();
    m_colours[0][anti0] = m_colours[1][anti1];
    m_colours[2][anti0] = m_colours[3][anti1];
    m_colours[0][1-anti0] = m_colours[1][1-anti1] = 0;
    m_colours[2][1-anti0] = m_colours[3][1-anti1] = 0;
  }
  return true; 
}

DECLARE_XSBASE_GETTER(q1q2_q1q2,"XS_q1q2_q1q2")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,q1q2_q1q2>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1].IsQuark() &&
      flavs[0]!=flavs[1] && flavs[2]==flavs[0] && flavs[3]==flavs[1]) {
    if (flavs[0].Mass()<1.e-6 && flavs[1].Mass()<1.e-6) return new q1q2_q1q2();
    THROW(fatal_error,"Massive quarks not yet enabled in QQ'->QQ'");
  }
  return NULL;
}



//////////////////////////////////////////////////////////////////////////
// q1q1bar -> q2q2bar
//////////////////////////////////////////////////////////////////////////

q1q1bar_q2q2bar::q1q1bar_q2q2bar() : XS_Base() {
  m_name = string("q1q1bar->q2q2bar");
}

void q1q1bar_q2q2bar::Calc(const double & s,const double & t,const double & u) {
  m_lastxs = 4./9. * (t*t+u*u)/(s*s);
}

bool q1q1bar_q2q2bar::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t anti(flavs[0].IsAnti()), part=1-anti;
  m_colours[2][anti] = Flow::Counter();
  m_colours[3][part] = Flow::Counter();
  m_colours[0][anti] = m_colours[2][anti];
  m_colours[1][part] = m_colours[3][part];
  m_colours[0][1-anti] = m_colours[2][1-anti] = 0;
  m_colours[1][1-part] = m_colours[3][1-part] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(q1q1bar_q2q2bar,"XS_q1q1bar_q2q2bar")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,q1q1bar_q2q2bar>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0].Bar() &&
      flavs[2].IsQuark() && flavs[2]!=flavs[0] && flavs[3]==flavs[2].Bar()) {
    if (flavs[0].Mass()<1.e-6 && flavs[1].Mass()<1.e-6) return new q1q1bar_q2q2bar();
    THROW(fatal_error,"Massive quarks not yet enabled in QQbar'->QQbar'");
  }
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////
//      QUARKONIA
//
// singlet m_pref 4.*sqr(4.*M_PI)/sqr(4.*M_PI)*R_0^2
// octet m_pref = pow(4.*M_PI,3)/sqr(4.*M_PI)*LDME. The division by 16 pi is because 
// of a conversion factor later (assumed only for 2nd order processes) 
// from alpha_s to g^2 so we account for it here.
// Otherwise, same logic as in EXTRA_XS
/////////////////////////////////////////////////////////////////////////////
XS_gg_g3S1::XS_gg_g3S1(const ATOOLS::Flavour_Vector& flavs): XS_Base(), Q(flavs)
{
  m_name = string("gg->g3S1");
}

void XS_gg_g3S1::Calc(const double & s,const double & t,const double & u) 
{
  //there is a 1/(16*M_PI) difference to EXTRA_XS because of the coupling treatment
  m_lastxs = Q.R02*4*Quarkonium::Calc::gg_g3S1(s,t,u,Q.mass2);
}

bool XS_gg_g3S1::SetColours(const ATOOLS::Flavour_Vector& flavs) 
{
  Quarkonium::ColourFlowSetter::gg_gSinglet(Q.index,m_colours);
  return true;
}


// OCTETS
// 1S0
DECLARE_XSBASE_GETTER(XS_qg_q1S0_oct,"XS_qg_q1S0_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qg_q1S0_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qg_qQuarkonium(flavs, {kf_1S0_c})) return new XS_qg_q1S0_oct(flavs);
  return NULL;
}

XS_qg_q1S0_oct::XS_qg_q1S0_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsGluon()) m_g = i;
    if (i<2 && flavs[i].IsQuark()) m_q = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_a      = flavs[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qg_q1S0_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qg_q1S0_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qg_q1S0_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qg_qOctet(m_a,m_g,Q.index, m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_qqbar_g1S0_oct,"XS_qqbar_g1S0_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qqbar_g1S0_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qqbar_gQuarkonium(flavs, {kf_1S0_c})) return new XS_qqbar_g1S0_oct(flavs);
  return NULL;
}

XS_qqbar_g1S0_oct::XS_qqbar_g1S0_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsAnti())  m_a = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qqbar_g1S0_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qqbar_g1S0_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qqbar_g1S0_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qqbar_gOctet(Q.index,m_a,m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_gg_g1S0_oct,"XS_gg_g1S0_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g1S0_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::gg_gQuarkonium(flavs, {kf_1S0_c})) return new XS_gg_g1S0_oct(flavs);
  return NULL;
}

XS_gg_g1S0_oct::XS_gg_g1S0_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs){
  m_name = string("gg->g1S0_oct");
  const ATOOLS::Flavour_Vector& fl = flavs;
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_gg_g1S0_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::gg_g1S0_oct(s,t,u,Q.mass2);
  m_lastxs =  all*m_pref;
}

bool XS_gg_g1S0_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::gg_gOctet(Q.index,m_colours);
  return true;
}

///// 3S1

DECLARE_XSBASE_GETTER(XS_qg_q3S1_oct,"XS_qg_q3S1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qg_q3S1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qg_qQuarkonium(flavs, {kf_3S1_c})) return new XS_qg_q3S1_oct(flavs);
  return NULL;
}

XS_qg_q3S1_oct::XS_qg_q3S1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsGluon()) m_g = i;
    if (i<2 && flavs[i].IsQuark()) m_q = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_a      = flavs[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qg_q3S1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qg_q3S1_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qg_q3S1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qg_qOctet(m_a,m_g,Q.index, m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_qqbar_g3S1_oct,"XS_qqbar_g3S1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qqbar_g3S1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qqbar_gQuarkonium(flavs, {kf_3S1_c})) return new XS_qqbar_g3S1_oct(flavs);
  return NULL;
}

XS_qqbar_g3S1_oct::XS_qqbar_g3S1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsAnti())  m_a = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qqbar_g3S1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qqbar_g3S1_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qqbar_g3S1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qqbar_gOctet(Q.index,m_a,m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_gg_g3S1_oct,"XS_gg_g3S1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g3S1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::gg_gQuarkonium(flavs, {kf_3S1_c})) return new XS_gg_g3S1_oct(flavs);
  return NULL;
}

XS_gg_g3S1_oct::XS_gg_g3S1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs){
  m_name = string("gg->g3S1_oct");
  const ATOOLS::Flavour_Vector& fl = flavs;
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_gg_g3S1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::gg_g3S1_oct(s,t,u,Q.mass2);
  m_lastxs =  all*m_pref;
}

bool XS_gg_g3S1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::gg_gOctet(Q.index,m_colours);
  return true;
}
// 3P0

DECLARE_XSBASE_GETTER(XS_qg_q3P0_oct,"XS_qg_q3P0_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qg_q3P0_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qg_qQuarkonium(flavs, {kf_3P0_c})) return new XS_qg_q3P0_oct(flavs);
  return NULL;
}

XS_qg_q3P0_oct::XS_qg_q3P0_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsGluon()) m_g = i;
    if (i<2 && flavs[i].IsQuark()) m_q = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_a      = flavs[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qg_q3P0_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qg_q3P0_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qg_q3P0_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qg_qOctet(m_a,m_g,Q.index, m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_qqbar_g3P0_oct,"XS_qqbar_g3P0_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qqbar_g3P0_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qqbar_gQuarkonium(flavs, {kf_3P0_c})) return new XS_qqbar_g3P0_oct(flavs);
  return NULL;
}

XS_qqbar_g3P0_oct::XS_qqbar_g3P0_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsAnti())  m_a = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qqbar_g3P0_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qqbar_g3P0_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qqbar_g3P0_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qqbar_gOctet(Q.index,m_a,m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_gg_g3P0_oct,"XS_gg_g3P0_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g3P0_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::gg_gQuarkonium(flavs, {kf_3P0_c})) return new XS_gg_g3P0_oct(flavs);
  return NULL;
}

XS_gg_g3P0_oct::XS_gg_g3P0_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs){
  m_name = string("gg->g3P0_oct");
  const ATOOLS::Flavour_Vector& fl = flavs;
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_gg_g3P0_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::gg_g3P0_oct(s,t,u,Q.mass2);
  m_lastxs =  all*m_pref;
}

bool XS_gg_g3P0_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::gg_gOctet(Q.index,m_colours);
  return true;
}
// 3P1

DECLARE_XSBASE_GETTER(XS_qg_q3P1_oct,"XS_qg_q3P1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qg_q3P1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qg_qQuarkonium(flavs, {kf_3P1_c})) return new XS_qg_q3P1_oct(flavs);
  return NULL;
}

XS_qg_q3P1_oct::XS_qg_q3P1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsGluon()) m_g = i;
    if (i<2 && flavs[i].IsQuark()) m_q = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_a      = flavs[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qg_q3P1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qg_q3P1_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qg_q3P1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qg_qOctet(m_a,m_g,Q.index, m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_qqbar_g3P1_oct,"XS_qqbar_g3P1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qqbar_g3P1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qqbar_gQuarkonium(flavs, {kf_3P1_c})) return new XS_qqbar_g3P1_oct(flavs);
  return NULL;
}

XS_qqbar_g3P1_oct::XS_qqbar_g3P1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsAnti())  m_a = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qqbar_g3P1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qqbar_g3P1_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qqbar_g3P1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qqbar_gOctet(Q.index,m_a,m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_gg_g3P1_oct,"XS_gg_g3P1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g3P1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::gg_gQuarkonium(flavs, {kf_3P1_c})) return new XS_gg_g3P1_oct(flavs);
  return NULL;
}

XS_gg_g3P1_oct::XS_gg_g3P1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs){
  m_name = string("gg->g3P1_oct");
  const ATOOLS::Flavour_Vector& fl = flavs;
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_gg_g3P1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::gg_g3P1_oct(s,t,u,Q.mass2);
  m_lastxs =  all*m_pref;
}

bool XS_gg_g3P1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::gg_gOctet(Q.index,m_colours);
  return true;
}

// 3P2

DECLARE_XSBASE_GETTER(XS_qg_q3P2_oct,"XS_qg_q3P2_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qg_q3P2_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qg_qQuarkonium(flavs, {kf_3P2_c})) return new XS_qg_q3P2_oct(flavs);
  return NULL;
}

XS_qg_q3P2_oct::XS_qg_q3P2_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsGluon()) m_g = i;
    if (i<2 && flavs[i].IsQuark()) m_q = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_a      = flavs[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qg_q3P2_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qg_q3P2_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qg_q3P2_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qg_qOctet(m_a,m_g,Q.index, m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_qqbar_g3P2_oct,"XS_qqbar_g3P2_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_qqbar_g3P2_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::qqbar_gQuarkonium(flavs, {kf_3P2_c})) return new XS_qqbar_g3P2_oct(flavs);
  return NULL;
}

XS_qqbar_g3P2_oct::XS_qqbar_g3P2_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs)
{
  for (short int i=0;i<4;i++) {
    if (i<2 && flavs[i].IsAnti())  m_a = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_qqbar_g3P2_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::qqbar_g3P2_oct(s,t,u,Q.mass2);
  m_lastxs = m_pref*all;
}

bool XS_qqbar_g3P2_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::qqbar_gOctet(Q.index,m_a,m_colours);
  return true;
}

DECLARE_XSBASE_GETTER(XS_gg_g3P2_oct,"XS_gg_g3P2_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g3P2_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  if (flavs.size()!=4) return NULL;
  if(Quarkonium::IsCandidate::gg_gQuarkonium(flavs, {kf_3P2_c})) return new XS_gg_g3P2_oct(flavs);
  return NULL;
}

XS_gg_g3P2_oct::XS_gg_g3P2_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(), Q(flavs){
  m_name = string("gg->g3P2_oct");
  const ATOOLS::Flavour_Vector& fl = flavs;
  m_pref = 4.*M_PI*Q.ldme;
}

void XS_gg_g3P2_oct::Calc(const double & s,const double & t,const double & u) 
{
  double all = Quarkonium::Calc::gg_g3P2_oct(s,t,u,Q.mass2);
  m_lastxs =  all*m_pref;
}

bool XS_gg_g3P2_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  Quarkonium::ColourFlowSetter::gg_gOctet(Q.index,m_colours);
  return true;
}

/////////////////////////////////////////////////////////////
DECLARE_XSBASE_GETTER(XS_gg_g3S1,"XS_gg_g3S1")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g3S1>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  const Flavour_Vector fl = flavs;
  if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
  fl[3].IsMeson() ) {
    kf_code kfc = fl[3].Kfcode();
    if (kfc==kf_J_psi_1S || kfc == kf_psi_2S || kfc==kf_Upsilon_1S || kfc == kf_Upsilon_2S || kfc == kf_Upsilon_3S){
      
      if(!IsZero(GetLDME(kfc))) return new XS_gg_g3S1(flavs);
      WarnZeroLDME(kfc);
    }    
  }
  if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
  fl[2].IsMeson() ) {
    kf_code kfc = fl[2].Kfcode();
    if (kfc==kf_J_psi_1S || kfc == kf_psi_2S || kfc==kf_Upsilon_1S || kfc == kf_Upsilon_2S || kfc == kf_Upsilon_3S){
      
      if(!IsZero(GetLDME(kfc))) return new XS_gg_g3S1(flavs);
      WarnZeroLDME(kfc);
    }      
  }
return NULL;
}
