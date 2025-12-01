#include "HADRON_RESCATTERING/XSecs/PiPi.H"
#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

PiPi::PiPi() :
  m_mpi(Flavour(kf_pi_plus).HadMass()), m_mpi2(sqr(m_mpi)),
  m_mK(Flavour(kf_K_plus).HadMass()),   m_mK2(sqr(m_mK))
{
  Test();
  exit(1);
}

const double PiPi::
XStot(const double & s,const int & I1,const int & I2) const {
  if (s<sqr(1.42)) {
    return ( 16.*M_PI/sigma(s)/sqrt(sqr(s-2.*m_mpi2)-4.*m_mpi2*m_mpi2) *
	     ( hatXSel(s,I1,I2) + hatXSinel(s,I1,I2) ) );
  }
  return ( 4.*sqr(M_PI)/sqrt(sqr(s-2.*m_mpi2)-4.*m_mpi2*m_mpi2) *
	   hatXStot(s,I1,I2) );
}

const double PiPi::
hatXStot(const double & s,const int & I1,const int & I2) const {
  Complex ampl = Complex(0.,0.);
  int I = I1+I2;
  if (I==0 && I1==0) {
    ampl = ( 1./3.  * sigma_I0(s,0.) +
	     1./4.  * sigma_I1(s,0.) +
	     2./3.  * sigma_I2(s,0.) );
  }
  else if (I==0 && (I1==1 || I1==-1)) {
    ampl = ( 1./3.  * sigma_I0(s,0.) +
	              sigma_I1(s,0.) +
	     1./6.  * sigma_I2(s,0.) );
  }
  else if (I==1 || I==-1) {
    ampl = ( 1./3. * sigma_I0(s,0.) +
	     1./6. * sigma_I1(s,0.) -
	     1./3. * sigma_I2(s,0.) );
  }
  else if (I==2 || I==-2) {
    ampl = ( 1./3.  * sigma_I0(s,0.) -
	     1./2.  * sigma_I1(s,0.) +
	     1./6.  * sigma_I2(s,0.) );
  }
  return ampl.imag();
}

const double PiPi::
hatXSel(const double & s,const int & I1,const int & I2) const {
  double xsel = 0.;
  int I = I1+I2;
  if (I==0 && I1==0) {
    xsel =  ( 1./3. * (xsel_partial(0,0,s)+
		       xsel_partial(0,2,s))  +
	      2./3. * (xsel_partial(2,0,s) +
		       xsel_partial(2,2,s)) );
  }
  else if (I==0 && (I1==1 || I1==-1)) {
    // pi^+ pi^- -> pi^+ pi^-
    xsel =  ( 1./3. * (xsel_partial(0,0,s)+
		       xsel_partial(0,2,s)) +
	      1./2. *  xsel_partial(1,1,s)  +
	      1./6. * (xsel_partial(2,0,s)+
		       xsel_partial(2,2,s)) ); 
  }
  else if (I==1 || I==-1) {
    xsel = ( 1./2. *  xsel_partial(1,1,s)  +
	     1./2. * (xsel_partial(2,0,s)+
		      xsel_partial(2,2,s)) );
  }
  else if (I==2 || I==-2) {
    xsel = (xsel_partial(2,0,s)+
	    xsel_partial(2,2,s) );
  }
  return 2.*xsel;
}

const double PiPi::
hatXSinel(const double & s,const int & I1,const int & I2) const {
  double xsinel = 0.;
  int I = I1+I2;
  if (I==0 && I1==0) {
    xsinel = ( 1./3. * (xsinel_partial(0,0,s)+
			xsinel_partial(0,2,s)) +
	       2./3. * (xsinel_partial(2,0,s)+
			xsinel_partial(2,2,s)) );
  }
  else if (I==0 && (I1==1 || I1==-1)) {
    xsinel = ( 1./3. * (xsinel_partial(0,0,s)+
			xsinel_partial(0,2,s))  +
	       1./2. *  xsinel_partial(1,1,s)   +
	       1./6. * (xsinel_partial(2,0,s)+
			xsinel_partial(2,2,s)) ); 
  }
  else if (I==1 || I==-1) {
    xsinel = ( 1./2. *  xsinel_partial(1,1,s)   +
	       1./2. * (xsinel_partial(2,0,s)+
			xsinel_partial(2,2,s)) );
  }
  else if (I==2 || I==-2) {
    xsinel = (xsinel_partial(2,0,s)+
	      xsinel_partial(2,2,s));
  }
  return 2.*xsinel;
}


const double PiPi::eta(const size_t & I,const size_t & l,
		       const double & s) const {
  if (s < 4.*m_mK2) return 1.;
  /////////////////////////////////////////////////////////////////////////
  // I = 0, S-wave (A.4 & A.5)
  /////////////////////////////////////////////////////////////////////////
  if (I==0 && l==0) {
    double k2 = 1./4.-m_mK2/s, k = sqrt(k2);
    return 1.-(6.4*k - 16.8*k2) * (sqr(1.500)-s)/s;
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 2, S-wave (A.8c)
  /////////////////////////////////////////////////////////////////////////
  else if (I==2 && l==0 && s>sqr(1.05)) {
    return 1.-0.17*pow(1.-sqr(1.05)/s,1.5);
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 1, P-wave (A.10)
  /////////////////////////////////////////////////////////////////////////
  else if (I==1 && l==1) {
    return 1.-0.30*(sqrt(s)-1.);
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 0, D-wave  (A.11c)
  /////////////////////////////////////////////////////////////////////////
  else if (I==0 && l==2) {
    double k2 = s/4.-m_mK2;
    return ( 1. - 0.262 * sqrt(k2/(sqr(1.275)/4.-m_mK2)) );
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 2, D-wave (text in Sec A.8)
  /////////////////////////////////////////////////////////////////////////
  else if (I==2 && l==2) {
    return 1.;
  }
  return 1.;
}

const double PiPi::delta(const size_t & I,const size_t & l,
			 const double & s) const {
  double delta = 0.;
  double sqrt_s = sqrt(s);
  /////////////////////////////////////////////////////////////////////////
  // I = 0, S-wave (A.1) & (A.3, A.5)
  // Add pi to the result for delta at higher masses to end up on the right
  // sheet that reproduces data. 
  /////////////////////////////////////////////////////////////////////////
  if (I==0 && l==0) {
    if (s < 4.*m_mK2) {
      double sqrt_s0s = sqrt(4.*m_mK2-s);
      double cotdelta = ( 1./sigma(s) * m_mpi2/(s-sqr(0.195)/2.) *
			  (1.-s/sqr(0.790)) *
			  (17.4 + 4.3*(sqrt_s-sqrt_s0s)/(sqrt_s+sqrt_s0s)) );
      delta = M_PI/2.-atan(cotdelta);
    }
    else {
      double k2 = s/4.-m_mK2/s, k = sqrt(k2);
      double cotdelta = ( (1.3 * (s-sqr(0.920)) * (sqr(1.320)-s) * k)/
			  (sqr(1.320) * sqrt(s) * k2) );
      delta = 3.*M_PI/2.-atan(cotdelta);
    }
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 2, S-wave (A.6, A.7) & 
  /////////////////////////////////////////////////////////////////////////
  else if (I==2 && l==0) {
    if (s < 4.*m_mK2) {
      double sqrt_s0s = sqrt(sqr(1.05)-s);
      double cotdelta = ( 1./sigma(s) * m_mpi2/(s-2.*sqr(0.147)) *
			  (-80.8 - 77.0*(sqrt_s-sqrt_s0s)/(sqrt_s+sqrt_s0s)) );
      delta = atan(1./cotdelta);
    }
    else {
      double sqrt_s0s = sqrt(sqr(1.45)-s);
      double cotdelta = ( 1./sigma(s) * m_mpi2/(s-2.*m_mpi2) *
			  (-125.0 - 119.0*(sqrt_s-sqrt_s0s)/(sqrt_s+sqrt_s0s)) );
      delta = atan(1./cotdelta);
    }
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 1, P-wave (A.9a, A.9b) & (A.10)
  /////////////////////////////////////////////////////////////////////////
  else if (I==1 && l==1) {
    if (s < 4.*m_mK2) {
      double sqrt_s0s = sqrt(sqr(1.05)-s);
      double cotdelta = ( 1./(sigma(s)*k(s)) * (sqr(0.776)-s) *
			  (1.064 - 0.170*(sqrt_s-sqrt_s0s)/
			   (sqrt_s+sqrt_s0s)) );
      delta = M_PI/2.-atan(cotdelta);
    }
    else delta = 2.69 + 1.1 * (sqrt(s)-1.);
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 0, D-wave (A.11a, A.11b)
  /////////////////////////////////////////////////////////////////////////
  else if (I==0 && l==2) {
    double sqrt_s0s = sqrt(sqr(1.45)-s);
    double cotdelta = ( 1./(sigma(s) * sqr(k(s))) *
			m_mpi2*(sqr(1.275)-s)*
			(23.5 + 24.8*(sqrt_s-sqrt_s0s)/(sqrt_s+sqrt_s0s)) );
    delta = M_PI/2.-atan(cotdelta);
  }
  /////////////////////////////////////////////////////////////////////////
  // I = 2, D-wave (A.12a, A.12b)
  /////////////////////////////////////////////////////////////////////////
  else if (I==2 && l==2) {
    double sqrt_s0s = sqrt(sqr(1.45)-s);
    double w_s      = (sqrt_s-sqrt_s0s)/(sqrt_s+sqrt_s0s);
    double cotdelta = ( 1./(sigma(s) * sqr(k(s))) * sqr(m_mpi2) *
			s/(4.*(m_mpi2+sqr(0.212))-s) *
			(2900. + 7300.*w_s + 25400.*sqr(w_s)) );
    delta = atan(1./cotdelta);
  }
  return delta;
}

const Complex PiPi::sigma_I0(const double & s,const double & t) const {
  return ( hrpars->Xi_pom() * sqr(hrpars->beta_pi_pom(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_pom(t)) +
	   hrpars->Xi_f2()  * sqr(hrpars->beta_pi_f2(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_f2(t)) );
}

const Complex PiPi::sigma_I1(const double & s,const double & t) const {
  return ( hrpars->Xi_rho() * sqr(hrpars->beta_pi_rho(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_rho(t)) );
}

const Complex PiPi::sigma_I2(const double & s,const double & t) const {
  return ( hrpars->Xi_I2() * sqr(hrpars->beta_pi_I2(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_I2(t)) );
}



void PiPi::Test() {
  size_t Nsteps = 100;
  double step = (1.45-2.*m_mpi)/double(Nsteps);

  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";
  msg_Out()<<"PIPI test called";

  for (size_t i=1;i<Nsteps;i++) {
    double  s       = sqr(2.*m_mpi+i*step);
    double  eta00   = eta(0,0,s);
    double  delta00 = delta(0,0,s);
    double  eta20   = eta(2,0,s);
    double  delta20 = delta(2,0,s);
    double  eta11   = eta(1,1,s);
    double  delta11 = delta(1,1,s);
    double  eta02   = eta(0,2,s);
    double  delta02 = delta(0,2,s);
    double  eta22   = eta(2,2,s);
    double  delta22 = delta(2,2,s);
    msg_Out()<<setw(8)<<setprecision(4)<<sqrt(s)
      // <<setw(12)<<setprecision(4)<<(XStot(s)*rpa->Picobarn()*1.e-12)
	     <<setw(12)<<setprecision(4)<<eta00
	     <<setw(12)<<setprecision(4)<<(delta00*180./M_PI)
	     <<setw(12)<<setprecision(4)<<eta20
	     <<setw(12)<<setprecision(4)<<(delta20*180./M_PI)
	     <<setw(12)<<setprecision(4)<<eta11
	     <<setw(12)<<setprecision(4)<<(delta11*180./M_PI)
	     <<setw(12)<<setprecision(4)<<eta02
	     <<setw(12)<<setprecision(4)<<(delta02*180./M_PI)
	     <<setw(12)<<setprecision(4)<<eta22
	     <<setw(12)<<setprecision(4)<<(delta22*180./M_PI)
	     <<"\n";
  }
  msg_Out()<<"----------------------------------------------------------------\n"
	   <<setw(8)<<setprecision(4)<<"E"
	   <<setw(12)<<setprecision(4)<<"00"
	   <<setw(12)<<setprecision(4)<<"+-"
	   <<setw(12)<<setprecision(4)<<"+0"
	   <<setw(12)<<setprecision(4)<<"++"
	   <<"\n";
  for (size_t i=1;i<Nsteps;i++) {
    double  s       = sqr(2.*m_mpi+i*step);
    msg_Out()<<setw(8)<<setprecision(4)<<sqrt(s)
	     <<setw(12)<<setprecision(4)<<(XStot(s,0,0)*rpa->Picobarn()*1.e-9)
	     <<setw(12)<<setprecision(4)<<(XStot(s,1,-1)*rpa->Picobarn()*1.e-9)
	     <<setw(12)<<setprecision(4)<<(XStot(s,1,0)*rpa->Picobarn()*1.e-9)
	     <<setw(12)<<setprecision(4)<<(XStot(s,1,1)*rpa->Picobarn()*1.e-9)
	     <<"\n";
  }
  msg_Out()<<"----------------------------------------------------------------\n"
	   <<setw(8)<<setprecision(4)<<"E"
	   <<setw(12)<<setprecision(4)<<"00"
	   <<setw(12)<<setprecision(4)<<"+-"
	   <<setw(12)<<setprecision(4)<<"+0"
	   <<setw(12)<<setprecision(4)<<"++"
	   <<"\n";
  double start = 0.5, inc = 0.1;
  while (start<21.) {
    double  s       = sqr(start);
    start          += inc;
    if (start>1.5) { inc = .5; }
    if (start>10.) { inc = 1.; }
    msg_Out()<<setw(8)<<setprecision(4)<<sqrt(s)
	     <<setw(12)<<setprecision(4)<<(XStot(s,0,0)*rpa->Picobarn()*1.e-9)
	     <<setw(12)<<setprecision(4)<<(XStot(s,1,-1)*rpa->Picobarn()*1.e-9)
	     <<setw(12)<<setprecision(4)<<(XStot(s,1,0)*rpa->Picobarn()*1.e-9)
	     <<setw(12)<<setprecision(4)<<(XStot(s,1,1)*rpa->Picobarn()*1.e-9)
	     <<"\n";
  }
  msg_Out()<<"----------------------------------------------------------------\n"
	   <<setw(8)<<setprecision(4)<<"E"
	   <<setw(12)<<setprecision(4)<<"00"
	   <<setw(12)<<setprecision(4)<<"+-"
	   <<setw(12)<<setprecision(4)<<"+0"
	   <<setw(12)<<setprecision(4)<<"++"
	   <<"\n";
  msg_Out()<<setw(8)<<setprecision(4)<<sqrt(1.69)
	   <<setw(12)<<setprecision(4)<<(XStot(1.69,0,0)*rpa->Picobarn()*1.e-9)
	   <<setw(12)<<setprecision(4)<<(XStot(1.69,1,-1)*rpa->Picobarn()*1.e-9)
	   <<setw(12)<<setprecision(4)<<(XStot(1.69,1,0)*rpa->Picobarn()*1.e-9)
	   <<setw(12)<<setprecision(4)<<(XStot(1.69,1,1)*rpa->Picobarn()*1.e-9)
	   <<"\n";
  msg_Out()<<setw(8)<<setprecision(4)<<sqrt(2.25)
	   <<setw(12)<<setprecision(4)<<(XStot(2.25,0,0)*rpa->Picobarn()*1.e-9)
	   <<setw(12)<<setprecision(4)<<(XStot(2.25,1,-1)*rpa->Picobarn()*1.e-9)
	   <<setw(12)<<setprecision(4)<<(XStot(2.25,1,0)*rpa->Picobarn()*1.e-9)
	   <<setw(12)<<setprecision(4)<<(XStot(2.25,1,1)*rpa->Picobarn()*1.e-9)
	   <<"\n";
}


