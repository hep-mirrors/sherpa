#include "XS_QCD.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

/* 
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha
*/

XS_q1q2_q1q2::XS_q1q2_q1q2(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  a  = _fl[0].IsAnti();
  p_colours[0][a] = p_colours[2][a] = 500;
  p_colours[1][a] = p_colours[3][a] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_q1q2_q1q2::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  return sqr(4.*M_PI*aS)* 4. * (s*s + u*u) / ( 9. * t*t);
}

bool XS_q1q2_q1q2::SetColours(double s,double t,double u) { 
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}
bool XS_q1q2_q1q2::SetColours()                           { return 1; }

//----------------------------------------------------------------------

XS_q1qbar1_q2qbar2::XS_q1qbar1_q2qbar2(int _nin,int _nout, 
				       Flavour * _fl)  : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a = _fl[0].IsAnti();
  p = 1-a;
  p_colours[0][a] = p_colours[1][p] = 500;
  p_colours[2][0] = p_colours[3][1] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_q1qbar1_q2qbar2::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  return sqr(4.*M_PI*aS)* 4. * (t*t + u*u) / ( 9. * s*s); 
}

bool XS_q1qbar1_q2qbar2::SetColours(double s,double t,double u) { 
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}
bool XS_q1qbar1_q2qbar2::SetColours()                           { return 1; }

//----------------------------------------------------------------------
// Note : Combinatorical factor of 2 for identical outgoing particles explicitly added

XS_q1q1_q1q1::XS_q1q1_q1q1(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a  = _fl[0].IsAnti();

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_q1q1_q1q1::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  Mt    = (s*s + u*u) / (t*t);
  Mu    = (s*s + t*t) / (u*u);
  return sqr(4.*M_PI*aS) * 4./9.*(Mt + Mu - 2./3. * (s*s) / (u*t)) /2.;
}

bool XS_q1q1_q1q1::SetColours(double s, double t, double u) 
{
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  Mt      = (s*s + u*u) / (t*t);
  Mu      = (s*s + t*t) / (u*u);
  return SetColours();
}


bool XS_q1q1_q1q1::SetColours() {
  if (Mt > (Mt+Mu) * ran.Get()) {
    p_colours[3][a] = p_colours[0][a] = 500;
    p_colours[2][a] = p_colours[1][a] = 501;
  }
  else {
    p_colours[2][a] = p_colours[0][a] = 500;
    p_colours[3][a] = p_colours[1][a] = 501;
  }
  return 1;
}

//----------------------------------------------------------------------

XS_q1qbar1_q1qbar1::XS_q1qbar1_q1qbar1(int _nin,int _nout, 
				       Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a  = _fl[0].IsAnti();
  p  = 1-a;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_q1qbar1_q1qbar1::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  Mt = (s*s + u*u) / (t*t);
  Ms = (t*t + u*u) / (s*s);
  return sqr(4.*M_PI*aS)*4./9.*(Mt + Ms - 2./3. * (u*u) / (s*t));
}


bool XS_q1qbar1_q1qbar1::SetColours(double s, double t, double u) {
  Mt    = (s*s + u*u) / (t*t);
  Ms    = (t*t + u*u) / (s*s);
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_q1qbar1_q1qbar1::SetColours() {
  if (Mt >  (Mt+Ms) * ran.Get()) {
    p_colours[0][a] = p_colours[2][a] = 500;	
    p_colours[1][p] = p_colours[3][p] = 501;
  }
  else {
    p_colours[0][a] = p_colours[3][p] = 500;	
    p_colours[1][p] = p_colours[2][a] = 501;
  }
  return 1;
}

//----------------------------------------------------------------------
// Note : Combinatorical factor of 2 for identical outgoing particles explicitly added

XS_q1qbar1_gg::XS_q1qbar1_gg(int _nin,int _nout, 
			     Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  a = _fl[0].IsAnti();
  p = 1-a;

  p_colours[0][a] = 500;
  p_colours[1][p] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_q1qbar1_gg::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  Mt = u/t;
  Mu = t/u;
  return sqr(4.*M_PI*aS)* (32./27.*(Mt+Mu) - 8./3.*(t*t +u*u)/ (s*s)) /2.;
}


bool XS_q1qbar1_gg::SetColours(double s, double t, double u) {
  Mt    = u/t;
  Mu    = t/u;
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_q1qbar1_gg::SetColours() {
  if (Mt > (Mt+Mu) * ran.Get()) {
    p_colours[2][a] = p_colours[0][a];
    p_colours[3][p] = p_colours[1][p];
    p_colours[2][p] = p_colours[3][a] = 502;
  }
  else {
    p_colours[3][a] = p_colours[0][a];
    p_colours[2][p] = p_colours[1][p];
    p_colours[3][p] = p_colours[2][a] = 502;
  }
  return 1;
}

//----------------------------------------------------------------------

XS_gg_q1qbar1::XS_gg_q1qbar1(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

  p_colours[0][0] = 500;
  p_colours[0][1] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_gg_q1qbar1::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  Mt = u/t;
  Mu = t/u;
  return sqr(4.*M_PI*aS)*(1./6.* ( Mt + Mu )   - 3./8. *(t*t +u*u)/ (s*s)); 
}


bool XS_gg_q1qbar1::SetColours(double s, double t, double u) {
  Mt      = u/t;
  Mu      = t/u;
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_gg_q1qbar1::SetColours() {
  if (Mt > (Mt+Mu) * ran.Get()) {
    p_colours[2][0] = p_colours[0][0];
    p_colours[3][1] = p_colours[1][1] = 502;
    p_colours[1][0] = p_colours[0][1];
  }
  else {
    p_colours[2][0] = p_colours[1][0] = 502;
    p_colours[3][1] = p_colours[0][1];
    p_colours[1][1] = p_colours[1][0];
  }
  return 1;
}

//----------------------------------------------------------------------

XS_q1g_q1g::XS_q1g_q1g(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  a = _fl[0].IsAnti();
  p = 1-a;

  p_colours[0][a] = 500;
  p_colours[2][a] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_q1g_q1g::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  Ms = u/s;
  Mu = s/u;
  return  sqr(4.*M_PI*aS)*(-4./9. * (Ms + Mu) +  (s*s + u*u)/(t*t));
}

bool XS_q1g_q1g::SetColours(double s, double t, double u) {
  Ms      = u/s;
  Mu      = s/u;
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}
      
bool XS_q1g_q1g::SetColours() {
  if (Ms > (Ms+Mu) * ran.Get()) {
    p_colours[3][a] = p_colours[0][a];
    p_colours[3][p] = p_colours[1][p] = 502;
    p_colours[1][a] = p_colours[2][a];
  }
  else {
    p_colours[3][p] = p_colours[2][a];
    p_colours[1][a] = p_colours[3][a] = 502;
    p_colours[1][p] = p_colours[0][a];
  }
  return 1;
}

//----------------------------------------------------------------------
// Note : Combinatorical factor of 2 for identical outgoing particles explicitly added

XS_gg_gg::XS_gg_gg(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  p_colours[0][0] = 500;
  p_colours[1][1] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
}

double XS_gg_gg::operator()(double s,double t,double u) {
  if (s<m_thres) return 0.;
  Ms = 1 - t*u/(s*s);
  Mt = 1 - s*u/(t*t);
  Mu = 1 - s*t/(u*u);
  return sqr(4.*M_PI*aS)*9./2. * ( Ms + Mt + Mu )/2.;
}
  
bool XS_gg_gg::SetColours(double s, double t, double u) {
  Ms      = 1 - t*u/(s*s);
  Mt      = 1 - s*u/(t*t);
  Mu      = 1 - s*t/(u*u);
  m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}
    
bool XS_gg_gg::SetColours() {
  double rr = ran.Get() * (Ms+Mt+Mu);
  if (rr-Mt < 0.) {
    p_colours[2][0] = p_colours[0][0];
    p_colours[3][1] = p_colours[1][1];
    p_colours[0][1] = p_colours[1][0] = 502;
    p_colours[2][1] = p_colours[3][0] = 503;
  }
  else {
    if (rr-Mu-Mt < 0.) {
      p_colours[3][0] = p_colours[0][0];
      p_colours[2][1] = p_colours[1][1];
      p_colours[0][1] = p_colours[1][0] = 502;
      p_colours[3][1] = p_colours[2][0] = 503;
    }
    else {
      p_colours[2][0] = p_colours[0][0];
      p_colours[3][1] = p_colours[0][1] = 502;
      p_colours[2][1] = p_colours[1][1];
      p_colours[3][0] = p_colours[1][0] = 503;
    }
  }
  return 1;
}




