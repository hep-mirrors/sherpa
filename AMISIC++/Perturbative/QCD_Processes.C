#include "AMISIC++/Perturbative/QCD_Processes.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/LDME.H"

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
/////////////////////////////////////////////////////////////////////////////

DECLARE_XSBASE_GETTER(XS_gg_g3S1_oct,"XS_gg_g3S1_oct")
XS_Base * ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,XS_gg_g3S1_oct>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  
  // const Flavour_Vector fl = args.Flavours();
  if (flavs.size()!=4) return NULL;
  if ( flavs[0].IsGluon() && flavs[1].IsGluon() && flavs[2].IsGluon()  &&
  flavs[3].IsOctetMeson() ) {
    kf_code kfc = flavs[3].Kfcode();
    if (kfc==kf_3S1_c) 
    {
      std::cout<<"FLAVOURS: "<<flavs[0].Kfcode()<<" "<<flavs[1].Kfcode()<<" "<<flavs[2].Kfcode()<<" "<<flavs[3].Kfcode()<<"\n";
      if(!IsZero(GetTotalLDME(kfc))) return new XS_gg_g3S1_oct(flavs);
    }
  }
  if ( flavs[0].IsGluon() && flavs[1].IsGluon() && flavs[3].IsGluon()  &&
  flavs[2].IsOctetMeson() ) {
    kf_code kfc = flavs[2].Kfcode();
    if (kfc==kf_3S1_c) 
    {
      std::cout<<"FLAVOURS: "<<flavs[0].Kfcode()<<" "<<flavs[1].Kfcode()<<" "<<flavs[2].Kfcode()<<" "<<flavs[3].Kfcode()<<"\n";
      if(!IsZero(GetTotalLDME(kfc))) return new XS_gg_g3S1_oct(flavs);
    }
  }
  
  return NULL;
}

XS_gg_g3S1_oct::XS_gg_g3S1_oct(const ATOOLS::Flavour_Vector & flavs): XS_Base(){
  m_name = string("gg->g3S1_oct");
  const ATOOLS::Flavour_Vector& fl = flavs;
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  const kf_code kfc = fl[m_S].Kfcode();
  m_mass = ATOOLS::Flavour((kfc / 100) % 10).Mass(true) +
  ATOOLS::Flavour((kfc / 10) % 10).Mass(true);
  m_mass2 = sqr(m_mass);
  LDME = GetTotalLDME(fl[m_S].Kfcode());
  m_pref = pow(4.*M_PI,3)*LDME;
}

void XS_gg_g3S1_oct::Calc(const double & s,const double & t,const double & u) 
{
  double sM2 = sqr(s-m_mass2);
  double tM2 = sqr(t-m_mass2);
  double uM2 = sqr(u-m_mass2);

  double heq0 = 2.*s*m_mass2*(sqr(t)+sqr(u))*t*u;
  double heq1 = sqr(s)*(sqr(sM2)+pow(t,4)+pow(u,4)+2.*sqr(m_mass2)*sqr(t*u/s));
  double nom = 27.*(s*t+t*u+u*s)-19.*sqr(m_mass2);
  double dnom = sM2*sM2*tM2*uM2;
  m_lastxs =  -1./(144.*pow(m_mass,3))*(heq0+heq1)*m_pref*nom/dnom;
}

bool XS_gg_g3S1_oct::SetColours(const ATOOLS::Flavour_Vector & flavs) 
{
  for (size_t i = 0; i<4;i++){
    if (flavs[i].IsOctetMeson()) m_S = i;
  }
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  size_t cross = ran->Get()<0.5 ? 0 : 1;
  if (cross==0) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[1][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[m_S][1-bit] = m_colours[5-m_S][bit] = Flow::Counter();
  }
  if (cross==1) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[1][bit] = m_colours[5-m_S][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[m_S][1-bit] = Flow::Counter();
  }
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


XS_gg_g3S1::XS_gg_g3S1(const ATOOLS::Flavour_Vector& flavs): XS_Base()
{
  m_name = string("gg->g3S1");
  const ATOOLS::Flavour_Vector& fl = flavs;
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  const kf_code kfc = fl[m_S].Kfcode();
  m_mass = ATOOLS::Flavour((kfc / 100) % 10).Mass(true) +
  ATOOLS::Flavour((kfc / 10) % 10).Mass(true);
  m_mass2 = sqr(m_mass);
  
  m_R02 =  GetLDME(fl[m_S].Kfcode())*2.*M_PI/(3.*3.);
}

void XS_gg_g3S1::Calc(const double & s,const double & t,const double & u) 
{
  double M2 = s+t+u;
  double sM2 = sqr(s-M2);
  double tM2 = sqr(t-M2);
  double uM2 = sqr(u-M2);
  double all = sqr(s)/(uM2*tM2)+sqr(t)/(uM2*sM2)+sqr(u)/(sM2*tM2);
  m_pref   = (5./9.)*sqr(4.*M_PI)*sqrt(M2)*m_R02;
  m_lastxs = m_pref*all;
}

bool XS_gg_g3S1::SetColours(const ATOOLS::Flavour_Vector& flavs) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}
