#include "XS_QCD.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

/* 
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha
   
   AlphaS is set to a fixed value since we intend to calculate it really 'running' 
   which can only be guaranteed, if we obtain the value on each sprime and y
   independently. 
   This statement is not a paradox, you simply have to set 
   
    KFactorScheme = 1
    Schalescheme  = 1

   in 'Run.dat'
   This will perform appropriate changes for calculations with XS' as well as for those with 
   generated Amplitudes.
*/

XS_q1q2_q1q2::XS_q1q2_q1q2(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  m12 = sqr(_fl[0].Mass());
  m22 = sqr(_fl[1].Mass());
  if ((m12 != 0) && (m22 != 0)) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;
  a  = _fl[0].IsAnti();
  colours[0][a] = colours[2][a] = 500;
  colours[1][a] = colours[3][a] = 501;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_q1q2_q1q2::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) return sqr(4.*M_PI*aS)* (4./9.) * ( s*s + u*u ) / ( t*t );
  else return sqr(4.*M_PI*aS) * (4./9.) * ( sqr(s-m12-m22) + sqr(u-m12-m22) - 2.*t*(m12+m22) + 16.*m12*m22 ) / ( t*t );
};

bool XS_q1q2_q1q2::SetColours(double s,double t,double u) { 
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}
bool XS_q1q2_q1q2::SetColours()                           { return 1; }

//----------------------------------------------------------------------


XS_q1qbar1_q2qbar2::XS_q1qbar1_q2qbar2(int _nin,int _nout, 
				       Flavour * _fl)  : 
  Single_XS(_nin,_nout,_fl) 
{
  m12 = sqr(_fl[0].Mass());
  m32 = sqr(_fl[2].Mass());
  if ((m12 != 0) && (m32 != 0)) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a = _fl[0].IsAnti();
  p = 1-a;
  colours[0][a] = colours[1][p] = 500;
  colours[2][0] = colours[3][1] = 501;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_q1qbar1_q2qbar2::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) return sqr(4.*M_PI*aS) * (4./9.) * ( t*t + u*u ) / ( s*s ); 
  else return sqr(4.*M_PI*aS) * (4./9.) * ( sqr(t-m12-m32) + sqr(u-m12-m32) + 2.*s*(m12+m32) ) / ( s*s );
};

bool XS_q1qbar1_q2qbar2::SetColours(double s,double t,double u) { 
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}
bool XS_q1qbar1_q2qbar2::SetColours()                           { return 1; }

//----------------------------------------------------------------------

XS_q1q1_q1q1::XS_q1q1_q1q1(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  m2 = sqr(_fl[0].Mass());
  if (m2 != 0) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a  = _fl[0].IsAnti();

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_q1q1_q1q1::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) {
    Mt    = ( s*s + u*u ) / ( t*t );
    Mu    = ( s*s + t*t ) / ( u*u );
    Mtu   = (2./3.) * ( s*s ) /( u*t );
  }
  else {
    Mt    = ( sqr(s-2.*m2) + sqr(u-2.*m2) - 4.*t*m2 + 16.*m2*m2 ) / ( t*t );
    Mu    = ( sqr(s-2.*m2) + sqr(t-2.*m2) - 4.*u*m2 + 16.*m2*m2 ) / ( u*u );
    Mtu   = (2./3.) * ( sqr(s-2.*m2) + 8.*m2*m2 ) / ( t*u );
  }
  return sqr(4.*M_PI*aS) * (4./9.) * ( Mt + Mu - Mtu);
};

bool XS_q1q1_q1q1::SetColours(double s, double t, double u) 
{
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  Mt    = (s*s + u*u) / (t*t);
  Mu    = (s*s + t*t) / (u*u);
  return SetColours();
}


bool XS_q1q1_q1q1::SetColours() {
  if (Mt > (Mt+Mu) * ran.Get()) {
    colours[3][a] = colours[0][a] = 500;
    colours[2][a] = colours[1][a] = 501;
  }
  else {
    colours[2][a] = colours[0][a] = 500;
    colours[3][a] = colours[1][a] = 501;
  }
  return 1;
}

//----------------------------------------------------------------------

XS_q1qbar1_q1qbar1::XS_q1qbar1_q1qbar1(int _nin,int _nout, 
				       Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  m2 = sqr(_fl[0].Mass());
  if (m2 != 0) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a  = _fl[0].IsAnti();
  p  = 1-a;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_q1qbar1_q1qbar1::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) {
    Mt  = ( s*s + u*u ) / ( t*t );
    Ms  = ( t*t + u*u ) / ( s*s );
    Mts = (2./3.) * ( u*u ) / ( s*t );
  }
  else {
    Mt    = ( sqr(s-2.*m2) + sqr(u-2.*m2) - 4.*t*m2 + 16.*m2*m2 ) / ( t*t );
    Ms    = ( sqr(t-2.*m2) + sqr(u-2.*m2) + 4.*s*m2 ) / ( s*s );
    Mts   = (2./3.) * ( sqr(s-2.*m2) + 8.*m2*m2 ) / ( s*t );
  }
  return sqr(4.*M_PI*aS) * (4./9.) * ( Mt + Ms - Mts );
};


bool XS_q1qbar1_q1qbar1::SetColours(double s, double t, double u) {
  Mt    = (s*s + u*u) / (t*t);
  Ms    = (t*t + u*u) / (s*s);
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_q1qbar1_q1qbar1::SetColours() {
  if (Mt >  (Mt+Ms) * ran.Get()) {
    colours[0][a] = colours[2][a] = 500;	
    colours[1][p] = colours[3][p] = 501;
  }
  else {
    colours[0][a] = colours[3][p] = 500;	
    colours[1][p] = colours[2][a] = 501;
  }
  return 1;
}

//----------------------------------------------------------------------

XS_q1qbar1_gg::XS_q1qbar1_gg(int _nin,int _nout, 
			     Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  m2 = sqr(_fl[0].Mass());
  if (m2 != 0) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a = _fl[0].IsAnti();
  p = 1-a;

  colours[0][a] = 500;
  colours[1][p] = 501;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_q1qbar1_gg::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) {
    Mt  = u / t;
    Mu  = t / u;
    Mtu = 0;
    Ms  = ( t*t + u*u ) / ( s*s );
  }
  else {
    Mt  = ( (t-m2)*(u+m2) - 16.*m2*(u+m2) + 2.*m2*(s-m2) ) / sqr(u-m2);
    Mu  = ( (u-m2)*(t+m2) - 16.*m2*(t+m2) + 2.*m2*(s-m2) ) / sqr(t-m2);
    Mtu = ( -2.*m2*(s-2.*m2) + 4.*m2*m2 ) / ( (t-m2) * (u-m2) );
    Ms  = ( sqr(t-m2) + sqr(u-m2) + 2.*s*m2 ) / ( s*s );
  }
  return sqr(4.*M_PI*aS) * ( (32./27.) * ( Mt + Mu ) - (8./3.) * Ms );
};


bool XS_q1qbar1_gg::SetColours(double s, double t, double u) {
  Mt    = u/t;
  Mu    = t/u;
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_q1qbar1_gg::SetColours() {
  if (Mt > (Mt+Mu) * ran.Get()) {
    colours[2][a] = colours[0][a];
    colours[3][p] = colours[1][p];
    colours[2][p] = colours[3][a] = 502;
  }
  else {
    colours[3][a] = colours[0][a];
    colours[2][p] = colours[1][p];
    colours[3][p] = colours[2][a] = 502;
  }
  return 1;
}

//----------------------------------------------------------------------

XS_gg_q1qbar1::XS_gg_q1qbar1(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  m2 = sqr(_fl[2].Mass());
  if (m2 != 0) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  colours[0][0] = 500;
  colours[0][1] = 501;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_gg_q1qbar1::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) {
    Mt  = u / t;
    Mu  = t / u;
    Mtu = 0;
    Ms  = ( t*t + u*u ) / ( s*s );
  }
  else {
    Mt  = ( (u-m2)*(t+m2) - 16.*m2*(t+m2) + 2.*m2*(s-m2) ) / sqr(t-m2);
    Mu  = ( (t-m2)*(u+m2) - 16.*m2*(u+m2) + 2.*m2*(s-m2) ) / sqr(u-m2);
    Mtu = ( -2.*m2*(s-2.*m2) + 4.*m2*m2 ) / ( (t-m2) * (u-m2) );
    Ms  = ( sqr(t-m2) + sqr(u-m2) + 2.*s*m2 ) / ( s*s );
  }
  return sqr(4.*M_PI*aS) * ( (1./6.) * ( Mt + Mu + Mtu ) - (3./8.) * Ms );
};


bool XS_gg_q1qbar1::SetColours(double s, double t, double u) {
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_gg_q1qbar1::SetColours() {
  if (Mt > (Mt+Mu) * ran.Get()) {
    colours[2][0] = colours[0][0];
    colours[3][1] = colours[1][1] = 502;
    colours[1][0] = colours[0][1];
  }
  else {
    colours[2][0] = colours[1][0] = 502;
    colours[3][1] = colours[0][1];
    colours[1][1] = colours[1][0];
  }
  return 1;
}

//----------------------------------------------------------------------

XS_q1g_q1g::XS_q1g_q1g(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (int i=0; i<2; i++) if ((m2 = sqr(_fl[2].Mass())) != 0) break;
  if (m2 != 0) massive = 1;
  else massive = 0;

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;
  a = _fl[0].IsAnti();

  colours[0][a] = 500;
  colours[2][a] = 501;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_q1g_q1g::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  if (!massive) {
    Ms  = u / s;
    Mu  = s / u;
    Msu = 0;
    Mt  = ( s*s + u*u ) / ( t*t );
  }
  else {
    Ms  = ( (u-m2)*(s+m2) - 16.*m2*(s+m2) + 2.*m2*(t-m2) ) / sqr(s-m2);
    Mu  = ( (s-m2)*(u+m2) - 16.*m2*(u+m2) + 2.*m2*(t-m2) ) / sqr(u-m2);
    Msu = ( -2.*m2*(t-2.*m2) + 4.*m2*m2 ) / ( (s-m2) * (u-m2) );
    Mt  = ( sqr(s-m2) + sqr(u-m2) + 2.*t*m2 ) / ( t*t );
  }
  return  sqr(4.*M_PI*aS)*( (-4./9.) * ( Ms + Mu + Msu ) + Mt );
};

bool XS_q1g_q1g::SetColours(double s, double t, double u) {
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}
      
bool XS_q1g_q1g::SetColours() {
  if (Ms > (Ms+Mu) * ran.Get()) {
    colours[3][a] = colours[0][a];
    colours[3][p] = colours[1][p] = 502;
    colours[1][a] = colours[2][a];
  }
  else {
    colours[3][p] = colours[2][a];
    colours[1][a] = colours[3][a] = 502;
    colours[1][p] = colours[0][a];
  }
  return 1;
}

//----------------------------------------------------------------------

XS_gg_gg::XS_gg_gg(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;
  colours[0][0] = 500;
  colours[1][1] = 501;

  aS = as->AsFixed();

  SetISRTypes(_fl);
};

double XS_gg_gg::operator()(double s,double t,double u) {
  if (s<Thres()) return 0.;
  Ms = 1 - t*u/(s*s);
  Mt = 1 - s*u/(t*t);
  Mu = 1 - s*t/(u*u);  
  return sqr(4.*M_PI*aS)*9./2. * ( Ms + Mt + Mu );
};
  
bool XS_gg_gg::SetColours(double s, double t, double u) {
  Ms    = 1 - t*u/(s*s);
  Mt    = 1 - s*u/(t*t);
  Mu    = 1 - s*t/(u*u);
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}
    
bool XS_gg_gg::SetColours() {
  double rr = ran.Get() * (Ms+Mt+Mu);
  if (rr-Mt < 0.) {
    colours[2][0] = colours[0][0];
    colours[3][1] = colours[1][1];
    colours[0][1] = colours[1][0] = 502;
    colours[2][1] = colours[3][0] = 503;
  }
  else {
    if (rr-Mu-Mt < 0.) {
      colours[3][0] = colours[0][0];
      colours[2][1] = colours[1][1];
      colours[0][1] = colours[1][0] = 502;
      colours[3][1] = colours[2][0] = 503;
    }
    else {
      colours[2][0] = colours[0][0];
      colours[3][1] = colours[0][1] = 502;
      colours[2][1] = colours[1][1];
      colours[3][0] = colours[1][0] = 503;
    }
  }
  return 1;
}




