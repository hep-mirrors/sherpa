#include "NLL_Branching_Probabilities.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Exception.H"

namespace SHERPA {
  const double   NC    = 3.;
  const double   CA    = NC;
  const double   CF    = (NC*NC-1.)/(2.*NC);
  const int      NF    = 5;
  const double   TR    =  1./2.;
  const double   BETA0 = (11.*CA-2.*NF)/3.;
  const double   BETA1 = (17.*CA*CA- 3.*CF*NF-5.*CA*NF)/3.;
  const double   KAPPA = CA*(67.0/18.0-ATOOLS::sqr(M_PI)/6.0)-NF*TR*10.0/9.0;
}

using namespace SHERPA;
using namespace ATOOLS;
using namespace MODEL;

Gamma_Lambda_Base::
Gamma_Lambda_Base(bpt::code type,bpm::code mode,double lambda, 
		  MODEL::Running_AlphaS * runas,double qmass,double asfac): 
  p_runas(runas), m_type(type), m_mode(mode), 
  m_colfac(0.), m_dlog(0.), m_slog(0.), m_qmass(qmass),
  m_lambda(lambda), m_as_factor(asfac), m_kfac(0.), m_pca(0)
{ 
  m_power[4]=m_power[3]=m_power[2]=m_power[1]=m_power[0]=0.0;
}

double Gamma_Lambda_Base::AlphaS(double t)                
{ 
  if (p_runas) return (*p_runas)(dabs(t)*m_as_factor);
  return 4.*M_PI/(BETA0*log(dabs(t)/sqr(m_lambda)));
}

void Gamma_Lambda_Base::SetZBounds
(const double &q, const double &Q,
 const double &m,const double &m1,const double &m2)
{
  m_z[1]=m_z[0]=q/Q;
  if (!(m_mode&bpm::z_shift)) return;
  m_z[0]=Max(m1/Q,m_z[0]);
  m_z[1]=Max(m2/Q,m_z[1]);
}

double Gamma_Lambda_Base::Gamma(double q, double Q) 
{
  if (Q<m_qmass) return 0.0;
  m_lastas=AlphaS(sqr(q));
  double val(2.*m_colfac*m_lastas/M_PI/q);
  if (m_mode & bpm::power_corrs) {
    double lf(m_pca==0?log((1.0-m_z[0])/m_z[1]):log(1.0/m_z[1]));
    val*=m_dlog*(1.+m_kfac*m_lastas/(2.*M_PI))*lf+m_slog
      +m_z[1]*(m_power[0]+m_z[1]*(m_power[1]+m_z[1]*m_power[2]));
  }
  else {
    val*=m_dlog*(1.+m_kfac*m_lastas/(2.*M_PI))*log(1.0/m_z[1])+m_slog;
  }
  return val;
}

double Gamma_Lambda_Base::IntGamma(double Q0, double Q) 
{
  if (Q0<m_lambda) {
    msg_Error()<<"Gamma_Lambda_Base::IntGamma("<<Q0<<","<<Q<<"): \n"
	       <<"   Lower bound below lambda_QCD = "<<m_lambda
	       <<" GeV.\n   Return 1e6."<<std::endl;
    return 1.e6;
  }
  double xi0(log(Q0/m_lambda)), xi1(log(Q/m_lambda));
  double xic(Max(xi1+m_slog/m_dlog,xi0));
  double result((m_dlog*xi1+m_slog)*log(dabs(xic/xi0))+m_dlog*(xi0-xic));
  if (m_power[0]>0.) result += m_power[0]*m_lambda/Q * 
    (ReIncompleteGamma0(-xi0)-ReIncompleteGamma0(-xi1));  
  if (m_qmass!=0.0) THROW(fatal_error,"No massive sudakov in analytic mode.");
  return 4.*m_colfac/BETA0*result;
}


GammaQ_QG_Lambda::GammaQ_QG_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double qmass, double asfac) : 
  Gamma_Lambda_Base(bpt::gamma_q2qg,mode,lambda,runas,qmass,asfac) 
{
  m_colfac = CF;
  m_dlog = 1.0;
  m_slog = -3./4.;
  if (m_mode & bpm::power_corrs) {
    if (m_mode & bpm::fs) {
      m_power[0] = 1.5;
    }
    else {
      m_power[0] = 1.0;
      m_power[1] = -0.25;
      m_pca=1;
    }
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaQ_QG_Lambda::Gamma(double q, double Q) 
{
  SetZBounds(q,Q,m_qmass,m_qmass,0.0);
  double val(Gamma_Lambda_Base::Gamma(q,Q));
  if (m_qmass==0.0 || 
      !(m_mode&bpm::massive || m_mode&bpm::dead_cone)) return Max(val,0.0);
  if (m_mode&bpm::dead_cone) {
    if (q<m_qmass) val+=2.0*m_colfac*m_lastas/M_PI/q*log(q/m_qmass);
  }
  else {
    if (!(m_mode & bpm::power_corrs)) {
      val+=m_colfac*m_lastas/M_PI/q* 
	(0.5-q/m_qmass*atan(m_qmass/q)-
	 (1.0-0.5*sqr(q/m_qmass))*log(1.0+sqr(m_qmass/q)));
    }
    else {
      double l(m_z[1]), m(m_z[0]), ll(l+m-1.0), a(q/m_qmass);
      double lt(log((a*a+sqr(1.0-m))/(a*a+l*l)));
      double at(atan((1.0-m)/a)-atan(l/a));
      double cmp(0.5*(3.0*(1.0-2.0*l)-4.0*ll-(2.0-a*a)*lt)-a*at);
      double cmm(ll/(a*a+l*l)/(a*a+sqr(1.0-m))*
		 (a*a*(3.0*(a*a+1.0)+2.0*(m*m+l*l)+m*(l-5.0))
		  +2.0*l*l*sqr(1.0-m)));
      val+=m_colfac*m_lastas/M_PI/q*(cmp+cmm);
    }
  }
  return Max(val,0.0);
}

GammasQ_sQG_Lambda::GammasQ_sQG_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double qmass, double asfac) : 
  Gamma_Lambda_Base(bpt::gamma_q2qg,mode,lambda,runas,qmass,asfac) 
{
  m_colfac = CF;
  m_dlog = 1.0;
  m_slog = -1.;
}

double GammasQ_sQG_Lambda::Gamma(double q, double Q) 
{
  //no power corrections supported
  SetZBounds(q,Q,m_qmass,m_qmass,0.0);
  double val(Gamma_Lambda_Base::Gamma(q,Q));
  if (m_qmass==0.0 || 
      !(m_mode&bpm::massive || m_mode&bpm::dead_cone)) return Max(val,0.0);
  if (m_mode&bpm::dead_cone) {
    if (q<m_qmass) val+=2.0*m_colfac*m_lastas/M_PI/q*log(q/m_qmass);
  }
  else {
    val+=m_colfac*m_lastas/M_PI/q* 
	(0.5-q/m_qmass*atan(m_qmass/q)-
	 (1.0-0.5*sqr(q/m_qmass))*log(1.0+sqr(m_qmass/q)));
  }
  return Max(val,0.0);
}

GammaQ_GQ_Lambda::GammaQ_GQ_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double qmass,double asfac) : 
  Gamma_Lambda_Base(bpt::gamma_q2gq,mode,lambda,runas,qmass,asfac) 
{
  m_colfac = CF;
  m_dlog = 1.0;
  m_slog = -3./4.;
  if (m_mode & bpm::power_corrs) {
    if (m_mode & bpm::fs) {
      m_power[0] = 1.5;
    }
    else {
      m_power[0] = 1.0;
      m_power[1] = -0.25;
      m_pca=1;
    }
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaQ_GQ_Lambda::Gamma(double q, double Q) 
{
  SetZBounds(q,Q,m_qmass,0.0,m_qmass);
  double val(Gamma_Lambda_Base::Gamma(q,Q));
  if (m_qmass==0.0 || 
      !(m_mode&bpm::massive || m_mode&bpm::dead_cone)) return Max(val,0.0);
  if (m_mode&bpm::dead_cone) {
    if (q<m_qmass) val+=2.0*m_colfac*m_lastas/M_PI/q*log(q/m_qmass);
  }
  else {
    if (!(m_mode & bpm::power_corrs)) {
      val+=m_colfac*m_lastas/M_PI/q* 
	(0.5-q/m_qmass*atan(m_qmass/q)-
	 (1.0-0.5*sqr(q/m_qmass))*log(1.0+sqr(m_qmass/q)));
    }
    else {
      double l(m_z[1]), m(m_z[0]), ll(l+m-1.0), a(q/m_qmass);
      double lt(log((a*a+sqr(1.0-m))/(a*a+l*l)));
      double at(atan((1.0-m)/a)-atan(l/a));
      double cmp(0.5*(3.0*(1.0-2.0*l)-4.0*ll-(2.0-a*a)*lt)-a*at);
      double cmm(ll/(a*a+l*l)/(a*a+sqr(1.0-m))*
		 (a*a*(3.0*(a*a+1.0)+2.0*(m*m+l*l)+m*(l-5.0))
		  +2.0*l*l*sqr(1.0-m)));
      val+=m_colfac*m_lastas/M_PI/q*(cmp+cmm);
    }
  }
  return Max(val,0.0);
}

GammaG_GG_Lambda::GammaG_GG_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double asfac):   
  Gamma_Lambda_Base(bpt::gamma_g2gg,mode,lambda,runas,0.0,asfac) 
{
  m_colfac = CA;
  m_dlog = 1.0;
  m_slog = -11./12.;
  if (m_mode & bpm::power_corrs) {
    m_power[0] = 2.0;
    m_power[1] = -0.5;
    m_power[2] = 1.0/3.0;
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaG_GG_Lambda::Gamma(double q, double Q) 
{
  SetZBounds(q,Q,0.0,0.0,0.0);
  return Max(Gamma_Lambda_Base::Gamma(q,Q),0.0);
}

GammaG_QQ_Lambda::GammaG_QQ_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double qmass,double asfac) : 
  Gamma_Lambda_Base(bpt::gamma_g2qq,
		    bpm::code((mode&(bpm::is|bpm::fs))|bpm::power_corrs),
		    lambda,runas,qmass,asfac) 
{
  m_colfac = TR;
  m_dlog   = 0.;
  m_slog = 1./3.;
  if (m_mode & bpm::power_corrs) {
    if (m_mode & bpm::fs) {
      m_power[0] = -1.0;
      m_power[1] = 1.0;
      m_power[2] = -2.0/3.0;
    }
    else {
      m_power[0] = -0.5;
      m_power[1] = 0.5;
      m_power[2] = -1.0/3.0;
      m_pca=1;
    }
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaG_QQ_Lambda::Gamma(double q, double Q) 
{
  SetZBounds(q,Q,0.0,m_qmass,m_qmass);
  if (m_qmass==0.0 || !(m_mode&bpm::massive) ||
      m_mode&bpm::dead_cone) return Max(Gamma_Lambda_Base::Gamma(q,Q),0.0);
  m_lastas=AlphaS(sqr(q));
  double val(m_colfac*m_lastas/M_PI*q/(q*q+sqr(m_qmass)));
  if (!(m_mode & bpm::power_corrs)) {
    val*=1.0-1.0/3.0/(1.0+sqr(m_qmass/q));
  }
  else {
    double l(m_z[0]), a(q/m_qmass);
    val*=2.0*((1.0-2.0*l*l*l)/3.0-l*(l-1.0))+
      ((1.0+4.0*l*l*l)/3.0-2.0*l*l)/(a*a+1.0);
  }
  return Max(val,0.0);
}

