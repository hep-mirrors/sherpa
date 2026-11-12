#include "HADRON_RESCATTERING/XSecs/PiPi.H"
#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

PiPi::PiPi() :
  m_mpi(Flavour(kf_pi_plus).HadMass()), m_mpi2(sqr(m_mpi)),
  m_mK(Flavour(kf_K_plus).HadMass()),   m_mK2(sqr(m_mK)), 
  m_alpha(1.), m_s0(4.*m_mK2), m_sp(Complex(0.996, -0.025))
{
  Test();
  exit(1);
}

const Complex PiPi::t00(const double & s) const {
  return t00_conf(s);
  Complex conf = t00_conf(s), f0 = t0_f0(s);
  return conf + f0 + 2.*Complex(0.,1.)*sigma(s)*conf*f0;
}

const Complex PiPi::t00_conf(const double & s) const {
  return 1./sigma(s) * 1./(Phi00(s)-Complex(0.,1.));
}

const Complex PiPi::t0_f0(const double & s) const {
  Complex f_sp     = f(m_sp);
  Complex sig_sp   = sqrt(1.-4.*m_mpi2/m_sp);
  Complex JBar_pi  = JBar(m_sp,m_mpi2);
  Complex JBar_K   = JBar(m_sp,m_mK2);
  double  f_spR    = f_sp.real(),    f_spI    = f_sp.imag();   
  double  sig_spR  = sig_sp.real(),  sig_spI  = sig_sp.imag();   
  double  JBar_piR = JBar_pi.real(), JBar_piI = JBar_pi.imag();   
  double  JBar_KR  = JBar_K.real(),  JBar_KI  = JBar_K.imag();
  double  d = ( JBar_piI * m_sp.real() + JBar_piR * m_sp.imag() +
		2. * ( sig_spI * m_sp.imag() -
		       sig_spR * m_sp.real() ) );
  double  M = ( ( ( f_spI*JBar_KR + f_spR*JBar_KI ) *
		  ( m_sp.imag() * (JBar_piI - 2.*sig_spR) -
		    m_sp.real() * (JBar_piR - 2.*sig_spI) ) +
		  ( JBar_piI - 2.*sig_spR) *
		  ( sqr(m_sp.imag()) + sqr(m_sp.real())) )/d -
		( f_spI * JBar_KI - f_spR * JBar_KR) );
  double  G = -( f_spI * JBar_KR - f_spR * JBar_KI + m_sp.imag())/d;

  return s*G/(M-s-JBar(s,m_mpi2)*s*G-JBar(s,m_mK2)*m_mK2*f(s));
}

const Complex PiPi::Phi00(const double & s) const {
  double prefactor = m_mpi2/(sigma(s)*(s-s_z02_S/2.));
  Complex sum = Complex(s_z02_S/(m_mpi*sqrt(s)), 0.);
  for (size_t n=0;n<6;n++) {
    sum += s_B_S[n] * pow(omega(s),n);
  }
  return prefactor * sum;
}

const Complex PiPi::f(const Complex & s) const {
  Complex term2 = 2.*s*s-Complex(1.,0.);
  Complex term3 = 2.*s*term2 - s;
  return ( s_K_S[0] * Complex(1.,0.) + s_K_S[1] * s +
	   s_K_S[2] * term2 + s_K_S[3] * term3 );
}

const Complex PiPi::JBar(const Complex & s,const double & mi2) const {
  Complex sig = sqrt(1.-4.*mi2/s);
  return (2.+sig*log((sig-1.)/(sig+1.)))/M_PI;
}

double PiPi::XStot(const double & s) {
  if (s<sqr(1.4)) {
    return 4.*M_PI/q(s) * t00(s).imag();
  }
  return 0.;    
}

double PiPi::XSel(const double & s) {
  return 0.;
}

void PiPi::Test() {
  size_t Nsteps = 100;
  double step = (1.4-2.*m_mpi)/double(Nsteps);
  for (size_t i=1;i<Nsteps;i++) {
    double  s       = sqr(2.*m_mpi+i*step);
    Complex term    = Complex(1.,0.)+Complex(0.,2.*sigma(s))*t00(s);
    double  eta00   = abs(term);
    double  delta00 = 1./2.*arg(term);
    Complex delta   = cos(Phi00(s))/sin(Phi00(s));
    msg_Out()<<setw(8)<<setprecision(4)<<sqrt(s)
	     <<setw(24)<<setprecision(4)<<term
	     <<setw(12)<<setprecision(4)<<(XStot(s)*rpa->Picobarn()*1.e-12)
	     <<setw(12)<<setprecision(4)<<eta00
	     <<setw(12)<<setprecision(4)<<(delta00*180./M_PI)
	     <<setw(12)<<setprecision(4)<<(delta*180.*M_PI)
	     <<"\n";
  }
}


double PiPi::s_z02_S    = sqr(0.137);
double PiPi::s_B_S[6]   = { 12.2, -0.9, 15.9, -5.7, -22.5, 6.9 };
double PiPi::s_K_S[4]   = { 5.25, -4.40, 0.175, -0.28 };
double PiPi::s_d_S[3]   = { -5.4, 0., 0. };
double PiPi::s_eps_S[3] = { 10.3, 0., 0. };

