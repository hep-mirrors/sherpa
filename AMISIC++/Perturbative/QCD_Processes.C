#include "AMISIC++/Perturbative/QCD_Processes.H"
#include "ATOOLS/Math/Random.H"

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
  double rr = ran->Get() * (m_Ms+m_Mt+m_Mu);
  if (rr-m_Mt<0.) {
    /*
      
      0====++====2
           ||
           ||  t
           ||
      1====++====3
      
    */
    m_colours[2][0] = m_colours[0][0];
    m_colours[3][1] = m_colours[1][1];
    m_colours[0][1] = m_colours[1][0] = Flow::Counter();
    m_colours[2][1] = m_colours[3][0] = Flow::Counter();
  }
  else if (rr-m_Mt-m_Mu<0.) {
    /*
      
	0====++====3
             ||
             ||
             ||
	1====++====2
	   
    */
    m_colours[3][0] = m_colours[0][0];
    m_colours[2][1] = m_colours[1][1];
    m_colours[0][1] = m_colours[1][0] = Flow::Counter();
    m_colours[3][1] = m_colours[2][0] = Flow::Counter();
  }
  else {
    /*
    
	0====++====2
             ||
             ||
             ||
	3====++====1
	   
    */
    m_colours[2][0] = m_colours[0][0];
    m_colours[3][1] = m_colours[0][1] = Flow::Counter();
    m_colours[2][1] = m_colours[1][1];
    m_colours[3][0] = m_colours[1][0] = Flow::Counter();
  }
  return true; 
}

DECLARE_XSBASE_GETTER(gg_gg,"XS_gg_gg")
XS_Base * Getter<XS_Base,Flavour_Vector,gg_gg>::
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
    m_colours[2][0] = m_colours[0][0];
    m_colours[3][1] = m_colours[1][1] = Flow::Counter();
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
    m_colours[2][0] = m_colours[1][0] = Flow::Counter();
    m_colours[3][1] = m_colours[0][1];
    m_colours[1][1] = m_colours[0][0];
  }
  m_colours[2][1] = m_colours[3][0] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(gg_qqbar,"XS_gg_qqbar")
XS_Base * Getter<XS_Base,Flavour_Vector,gg_qqbar>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsGluon() && flavs[1].IsGluon() &&
      flavs[2].IsQuark() && !flavs[2].IsAnti() && 
      flavs[3]==flavs[2].Bar()) {
    if (flavs[2].Mass()<=1.e-6) return new gg_qqbar();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for gg -> QQbar "
		 <<"with Q = "<<flavs[2]<<".\n";
      exit(1);
    }
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
    m_colours[2][anti] = m_colours[0][anti];
    m_colours[3][part] = m_colours[1][part];
    m_colours[2][part] = m_colours[3][anti] = Flow::Counter();
  }
  else {
     /*

      0------+====3
             |
             | u
             |
      1------+====2

    */
    m_colours[3][anti] = m_colours[0][anti];
    m_colours[2][part] = m_colours[1][part];
    m_colours[3][part] = m_colours[2][anti] = Flow::Counter();
  }
  return true;
}

DECLARE_XSBASE_GETTER(qqbar_gg,"XS_qqbar_gg")
XS_Base * Getter<XS_Base,Flavour_Vector,qqbar_gg>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[0]==flavs[1].Bar() &&
      flavs[2]==Flavour(kf_gluon) && flavs[3]==Flavour(kf_gluon)) {
    if (flavs[0].Mass()<1.e-6) return new qqbar_gg();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for QQbar->gg "
		 <<"with Q = "<<flavs[0]<<".\n";
      exit(1);
    }
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
    m_colours[gluon][anti]   = m_colours[quark][trip];
    m_colours[gluon+2][anti] = m_colours[quark+2][trip];
    m_colours[gluon][trip]   = m_colours[gluon+2][trip] = Flow::Counter();
  }
  else {
     /*

      0------+====3
             |
             | u
             |
      1======+----2

    */
    m_colours[gluon+2][trip] = m_colours[quark][trip];
    m_colours[gluon][trip]   = m_colours[quark+2][trip];
    m_colours[gluon][anti]   = m_colours[gluon+2][anti] = Flow::Counter();
  }
  return true;
}

DECLARE_XSBASE_GETTER(qg_qg,"XS_qg_qg")
XS_Base * Getter<XS_Base,Flavour_Vector,qg_qg>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[0]==flavs[2] &&
      flavs[1]==Flavour(kf_gluon) && flavs[3]==Flavour(kf_gluon)) {
    if (flavs[0].Mass()<1.e-6) {
      return new qg_qg();
    }
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for Qg->Qg "
		 <<"with Q = "<<flavs[0]<<".\n";
      exit(1);
    }
  }
  else if (flavs[1].IsQuark() && flavs[1]==flavs[3] &&
	   flavs[0]==Flavour(kf_gluon) && flavs[2]==Flavour(kf_gluon)) {
    if (flavs[1].Mass()<1.e-6) return new qg_qg();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for gQ->gQ "
		 <<"with Q = "<<flavs[1]<<".\n";
      exit(1);
    }
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
  if ((m_Mt+m_Mu) * ran->Get() < m_Mt) {
     /*
    
      0------++-----2
             ||
             || t
             ||
      1------++-----3

    */
    m_colours[0][anti] = m_colours[3][anti] = Flow::Counter();
    m_colours[1][anti] = m_colours[2][anti] = Flow::Counter();
  }
  else {
     /*

      0------++-----3
             ||
             || u
             ||
      1------++-----2

    */
    m_colours[0][anti] = m_colours[2][anti] = Flow::Counter();
    m_colours[1][anti] = m_colours[3][anti] = Flow::Counter();
  }
  return true;
}

DECLARE_XSBASE_GETTER(qq_qq,"XS_qq_qq")
XS_Base * Getter<XS_Base,Flavour_Vector,qq_qq>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0] &&
      flavs[2]==flavs[0] && flavs[3]==flavs[0]) {
    if (flavs[0].Mass()<1.e-6) return new qq_qq();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for QQ->QQ "
		 <<"with Q = "<<flavs[0]<<".\n";
      exit(1);
      vector<double> masses;
      masses.resize(4);
      for (size_t i=0;i<4;i++) masses[i] = flavs[i].Mass();
    }
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
    m_colours[0][anti0] = m_colours[2][anti0] = Flow::Counter();	
    m_colours[1][anti1] = m_colours[3][anti1] = Flow::Counter();
  }
  else {
    /*

      0----+ +----2
           | |
           | | t
           | |
      1----+ +----3

    */
    m_colours[0][anti0] = m_colours[1][anti1] = Flow::Counter();	
    m_colours[2][anti0] = m_colours[3][anti1] = Flow::Counter();
  }
  m_colours[0][1-anti0] = m_colours[2][1-anti0] = 0;
  m_colours[1][1-anti1] = m_colours[3][1-anti1] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(qqbar_qqbar,"XS_qqbar_qqbar")
XS_Base * Getter<XS_Base,Flavour_Vector,qqbar_qqbar>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0].Bar() &&
      flavs[2]==flavs[0] && flavs[3]==flavs[2].Bar()) {
    if (flavs[0].Mass()<1.e-6) return new qqbar_qqbar();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for QQBAR->QQBAR "
		 <<"with Q = "<<flavs[0]<<".\n";
      exit(1);
      vector<double> masses;
      masses.resize(4);
      for (size_t i=0;i<4;i++) masses[i] = flavs[i].Mass();
    }
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
    m_colours[0][anti0] = m_colours[3][anti0] = Flow::Counter();
    m_colours[1][anti1] = m_colours[2][anti1] = Flow::Counter();
    m_colours[0][1-anti0] = m_colours[3][1-anti0] = 0;
    m_colours[1][1-anti1] = m_colours[2][1-anti1] = 0;
  }
  else {
    // particle--anti-particle: t-channel colours
    m_colours[0][anti0] = m_colours[1][anti1] = Flow::Counter();
    m_colours[2][anti0] = m_colours[3][anti1] = Flow::Counter();
    m_colours[0][1-anti0] = m_colours[1][1-anti1] = 0;
    m_colours[2][1-anti0] = m_colours[3][1-anti1] = 0;
  }
  return true; 
}

DECLARE_XSBASE_GETTER(q1q2_q1q2,"XS_q1q2_q1q2")
XS_Base * Getter<XS_Base,Flavour_Vector,q1q2_q1q2>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1].IsQuark() &&
      flavs[0]!=flavs[1] && flavs[2]==flavs[0] && flavs[3]==flavs[1]) {
    if (flavs[0].Mass()<1.e-6 && flavs[1].Mass()<1.e-6) return new q1q2_q1q2();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for qq'->qq' "
		 <<"with q, q' = "<<flavs[0]<<", "<<flavs[1]<<".\n";
      exit(1);
      vector<double> masses;
      masses.resize(4);
      for (size_t i=0;i<4;i++) masses[i] = flavs[i].Mass();
      //return new q1q2_q1q2(masses);
    }
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
  m_colours[0][anti] = m_colours[2][anti] = Flow::Counter();
  m_colours[1][part] = m_colours[3][part] = Flow::Counter();
  m_colours[0][1-anti] = m_colours[2][1-anti] = 0;
  m_colours[1][1-part] = m_colours[3][1-part] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(q1q1bar_q2q2bar,"XS_q1q1bar_q2q2bar")
XS_Base * Getter<XS_Base,Flavour_Vector,q1q1bar_q2q2bar>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0].Bar() &&
      flavs[2].IsQuark() && flavs[2]!=flavs[0] && flavs[3]==flavs[2].Bar()) {
    if (flavs[0].Mass()<1.e-6 && flavs[1].Mass()<1.e-6)
      return new q1q1bar_q2q2bar();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for qqbar->q'q'bar"
		 <<"with q, q' = "<<flavs[0]<<", "<<flavs[2]<<".\n";
      exit(1);
      vector<double> masses;
      masses.resize(4);
      for (size_t i=0;i<4;i++) masses[i] = flavs[i].Mass();
      //return new q1q1bar_q2q2bar(masses);
    }
  }
  return NULL;
}

