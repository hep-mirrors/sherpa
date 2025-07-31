#include "HADRON_RESCATTERING/XSecs/PiK.H"
#include "HADRON_RESCATTERING/XSecs/PiPi.H"
#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

PiK::PiK() :
  m_mpi(Flavour(kf_pi_plus).HadMass()), m_mpi2(sqr(m_mpi)),
  m_mK(Flavour(kf_K_plus).HadMass()),   m_mK2(sqr(m_mK)),
  m_meta(Flavour(kf_eta).HadMass()),    m_meta2(sqr(m_meta)),
  m_sum2(m_mK2+m_mpi2),                 m_Delta2(m_mK2-m_mpi2),
  m_sKeta(sqr(m_meta+m_mK))
{
  Test();
}

const double PiK::
XStot(const double & s,const int & I1,const int & I2) const {
  if (s<sqr(1.8)) {
    return ( 16.*M_PI/sigma(s)/sqrt(sqr(s-m_mK2-m_mpi2)-4.*m_mK2*m_mpi2) *
	     ( hatXSel(s,I1,I2) + hatXSinel(s,I1,I2) ) );
  }
  return ( 4.*sqr(M_PI)/sqrt(sqr(s-m_mK2-m_mpi2)-4.*m_mK2*m_mpi2) *
	   hatXStot(s,I1,I2) );
}

const double PiK::
hatXStot(const double & s,const int & I1,const int & I2) const {
  Complex ampl = Complex(0.,0.);
  ampl = ( 1./3.  * sigma_I0(s,0.) +
	   1./3.  * sigma_I1(s,0.) );
  return ampl.imag();
}

const double PiK::
hatXSel(const double & s,const int & I1,const int & I2) const {
  double xsel = 0.;
  int I = I1+I2;
  return 2.*xsel;
}

const double PiK::
hatXSinel(const double & s,const int & I1,const int & I2) const {
  double xsinel = 0.;
  int I = I1+I2;
  return 2.*xsinel;
}


const double PiK::eta(const size_t & I,const size_t & l,
		       const double & s) const {
  return 1.;
}

const Complex PiK::ampl(const double & s,const int & I, const int & l) const {
  if (I==3 || (I==1 && s<m_sKeta)) {
    double delta_l = delta(3,l,s), cos_d = cos(delta_l), sin_d = sin(delta_l);
    return Complex(cos_d,sin_d)*sin_d/sigma(s);
  }
  else if (I==1 && s>=m_sKeta) {
    if (l==0) return ampl_10(s);
  }
  return Complex(0.,0.);
}

const Complex PiK::ampl_10(const double & s) const {    
  double sr1     = sqr(1.399), e1 = 1.,    G1 = 0.499;
  double sr2     = sqr(1.815), e2 = 0.184, G2 = 0.29;
  double qpiK    = sqrt(k2(s));
  double qetaK   = sqrt(keta2(s));
  double qpiKhat = sqrt(k2(m_sKeta));
  
  double qpiKr1  = sqrt(k2(sr1));
  double p1_quot = p1(qpiK)/p1(qpiKr1);
  double diffs1  = (qpiK-qpiKhat)/(qpiKr1-qpiKhat);
  double thres1  = sqrt((sqr(s-m_meta2-m_mK2)-4.*m_meta2*m_mK2)/
			(sqr(sr1-m_meta2-m_mK2)-4.*m_meta2*m_mK2));
  double qpiKr2  = sqrt(k2(sr2));
  double p2_quot = p2(qpiK)/p2(qpiKr2);
  double diffs2  = (qpiK-qpiKhat)/(qpiKr2-qpiKhat);
  double thres2  = sqrt((sqr(s-m_meta2-m_mK2)-4.*m_meta2*m_mK2)/
		      (sqr(sr2-m_meta2-m_mK2)-4.*m_meta2*m_mK2));

  double P1    = m_beta*(sr1-s) + e1 * G1 * p1_quot * diffs1;
  double Q1    =             (1.-e1) * G1 * p1_quot * thres1;
  double P2    =                  e2 * G2 * p2_quot * diffs2; 
  double Q2    =             (1.-e2) * G2 * p2_quot * thres2;

  double S0arg = 2.*qetaK*(-0.2+4.76*sqr(qetaK)); 
  Complex S0b  = Complex(cos(S0arg),sin(S0arg));
  return ( Complex(0.,-1./(2.*sigma(s)) )*
	   ( S0b*Sn(s,sr1,P1,Q1)*Sn(s,sr2,P2,Q2)- Complex(1.,0.)) );
}



const double PiK::delta(const size_t & I,const size_t & l,
			const double & s) const {
  double delta = 0.;
  if (I==3) {
    double s0      = sqr(1.84);
    double y0      = sqr((s0-m_Delta2)/(s0+m_Delta2));
    double omega_s = omega(s,y0,1.45), cotdelta = 0.;
    if (l==0) {
      /////////////////////////////////////////////////////////////////////////
      // I = 3/2, S-wave
      // Eq. (11) with parameters from UFD fit, Table 1, both in 1602.08404
      /////////////////////////////////////////////////////////////////////////
      double s_adler  = m_sum2;
      cotdelta = (1./(sigma(s) * (s-s_adler)) *
		  (2.25 + 4.21*omega_s + 2.45*sqr(omega_s)) );
    }
    else if (l==1) {
      /////////////////////////////////////////////////////////////////////////
      // I = 3/2, P-wave
      // Eq. (26) with parameters from UFD fit, Table 5, both in 1602.08404
      /////////////////////////////////////////////////////////////////////////
      cotdelta = 1./(sigma(s) * k2(s)) * (-14.8 + 2.7*omega_s);
    }
    else if (l==2) {
      /////////////////////////////////////////////////////////////////////////
      // I = 3/2, D-wave
      // Eq. (35) with parameters from UFD fit, Table 8, both in 1602.08404
      /////////////////////////////////////////////////////////////////////////
      cotdelta = ( 1./(sigma(s) * sqr(k2(s))) *
		   (-1.70 - 6.5*omega_s - 36.*sqr(omega_s)) );
    }
    delta = M_PI/2.-atan(cotdelta);
  }
  else if (I==1) {
    if (s<m_sKeta) {
      /////////////////////////////////////////////////////////////////////////
      // I = 1/2, S-wave, elastic regime
      // Eq. (13) with parameters from UFD fit, Table 2, both in 1602.08404
      /////////////////////////////////////////////////////////////////////////
      double s_adler  = (m_sum2+2.*sqrt(sqr(m_Delta2)+m_mK2*m_mpi2))/5.;
      double s0       = sqr(1.1), y0 = sqr((s0-m_Delta2)/(s0+m_Delta2));
      double omega_s  = omega(s,y0,1.15);
      double cotdelta = (1./sigma(s) * 1./(s-s_adler) *
			 (0.411 + 0.181*omega_s) );
      delta = M_PI/2.-atan(cotdelta);
    }
    else {
    }
  }
  return delta;
}

const Complex PiK::sigma_I0(const double & s,const double & t) const {
  return ( hrpars->Xi_pom() * sqr(hrpars->beta_K_pom(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_pom(t)) +
	   hrpars->Xi_f2()  * sqr(hrpars->beta_K_f2(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_f2(t)) );
}

const Complex PiK::sigma_I1(const double & s,const double & t) const {
  return ( hrpars->Xi_rho() * sqr(hrpars->beta_K_rho(t)) *
	   pow(s/hrpars->s0(),hrpars->alpha_rho(t)) );
}



void PiK::Test() {
  size_t Nsteps = 100;
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


