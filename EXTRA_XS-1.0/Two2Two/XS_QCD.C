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
*/

XS_q1q2_q1q2::XS_q1q2_q1q2(int _nin,int _nout, Flavour * _fl) : 
  Single_XS(_nin,_nout,_fl) 
{
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;
  a  = _fl[0].isanti();
  colours[0][a] = colours[2][a] = 500;
  colours[1][a] = colours[3][a] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_q1q2_q1q2::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  return sqr(4.*M_PI*aS)* 4. * (s*s + u*u) / ( 9. * t*t);
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
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a = _fl[0].isanti();
  p = 1-a;
  colours[0][a] = colours[1][p] = 500;
  colours[2][0] = colours[3][1] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_q1qbar1_q2qbar2::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  return sqr(4.*M_PI*aS)* 4. * (t*t + u*u) / ( 9. * s*s); 
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
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a  = _fl[0].isanti();

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_q1q1_q1q1::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  Mt    = (s*s + u*u) / (t*t);
  Mu    = (s*s + t*t) / (u*u);
  return sqr(4.*M_PI*aS) * 4./9.*(Mt + Mu - 2./3. * (s*s) / (u*t));
};

bool XS_q1q1_q1q1::SetColours(double s, double t, double u) 
{
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  Mt    = (s*s + u*u) / (t*t);
  Mu    = (s*s + t*t) / (u*u);
  return SetColours();
}


bool XS_q1q1_q1q1::SetColours() {
  if (Mt > (Mt+Mu) * Ran.get()) {
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
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a  = _fl[0].isanti();
  p  = 1-a;

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_q1qbar1_q1qbar1::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  Mt = (s*s + u*u) / (t*t);
  Ms = (t*t + u*u) / (s*s);
  return sqr(4.*M_PI*aS)*4./9.*(Mt + Ms - 2./3. * (u*u) / (s*t));
};


bool XS_q1qbar1_q1qbar1::SetColours(double s, double t, double u) {
  Mt    = (s*s + u*u) / (t*t);
  Ms    = (t*t + u*u) / (s*s);
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_q1qbar1_q1qbar1::SetColours() {
  if (Mt >  (Mt+Ms) * Ran.get()) {
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
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  a = _fl[0].isanti();
  p = 1-a;

  colours[0][a] = 500;
  colours[1][p] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_q1qbar1_gg::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  Mt = u/t;
  Mu = t/u;
  return sqr(4.*M_PI*aS)* (32./27.* ( Mt + Mu )   - 8./3. *(t*t +u*u)/ (s*s));
};


bool XS_q1qbar1_gg::SetColours(double s, double t, double u) {
  Mt    = u/t;
  Mu    = t/u;
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_q1qbar1_gg::SetColours() {
  if (Mt > (Mt+Mu) * Ran.get()) {
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
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;

  colours[0][0] = 500;
  colours[0][1] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_gg_q1qbar1::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  Mt = u/t;
  Mu = t/u;
  return sqr(4.*M_PI*aS)*(1./6.* ( Mt + Mu )   - 3./8. *(t*t +u*u)/ (s*s)); 
};


bool XS_gg_q1qbar1::SetColours(double s, double t, double u) {
  Mt    = u/t;
  Mu    = t/u;
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}

bool XS_gg_q1qbar1::SetColours() {
  if (Mt > (Mt+Mu) * Ran.get()) {
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
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;
  a = _fl[0].isanti();

  colours[0][a] = 500;
  colours[2][a] = 501;

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_q1g_q1g::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  Ms = u/s;
  Mu = s/u;
  return  sqr(4.*M_PI*aS)*(-4./9. * (Ms + Mu) +  (s*s + u*u)/(t*t));
};

bool XS_q1g_q1g::SetColours(double s, double t, double u) {
  Ms    = u/s;
  Mu    = s/u;
  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return SetColours();
}
      
bool XS_q1g_q1g::SetColours() {
  if (Ms > (Ms+Mu) * Ran.get()) {
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

  aS = (*as)(sqr(rpa.gen.Ecms()));
};

double XS_gg_gg::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
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
  double rr = Ran.get() * (Ms+Mt+Mu);
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




