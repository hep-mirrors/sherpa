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
  m_type(type), m_mode(mode), 
  m_colfac(0.), m_dlog(0.), m_slog(0.), m_qmass(qmass),
  m_lambda(lambda), m_as_factor(asfac), p_runas(runas), 
  m_kfac(0.)
{ 
  m_power[2]=m_power[1]=m_power[0]=0.0;
}

double Gamma_Lambda_Base::AlphaS(double t)                
{ 
  if (p_runas) return (*p_runas)(dabs(t)/m_as_factor);
  return 4.*M_PI/(BETA0*log(dabs(t)/sqr(m_lambda)));
}

double Gamma_Lambda_Base::Gamma(double q, double Q) 
{
  if (Q<m_qmass) return 0.0;
  m_lastas=AlphaS(sqr(q));
  double val(2.*m_colfac*m_lastas/M_PI/q);
  if (m_mode & bpm::power_corrs) {
    double qr(q/Q);
    val*=m_dlog*(1.+m_kfac*m_lastas/(2.*M_PI))*log(1.0/qr-1.0)+m_slog
      +qr*(m_power[0]+qr*(m_power[1]+qr*m_power[2]));
  }
  else {
    val*=m_dlog*(1.+m_kfac*m_lastas/(2.*M_PI))*log(Q/q)+m_slog;
  }
  return Max(val,0.0);
}

double Gamma_Lambda_Base::IntGamma(double Q0, double Q) 
{
  if (Q0<m_lambda) {
    msg.Error()<<"Gamma_Lambda_Base::IntGamma("<<Q0<<","<<Q<<"): \n"
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
  m_dlog   = 1.;
  if (m_mode & (bpm::linear_term | bpm::power_corrs)) {
    m_slog   = -3./4.;
    if (m_mode & bpm::power_corrs) m_power[0] = 1.5;
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaQ_QG_Lambda::Gamma(double q, double Q) 
{
  double val(Gamma_Lambda_Base::Gamma(q,Q));
  if (m_qmass==0.0 || !(m_mode&bpm::massive)) return val;
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
  m_dlog   = 1.;
  if (m_mode & (bpm::linear_term | bpm::power_corrs)) {
    m_slog = -3./4.;
    if (m_mode & bpm::power_corrs) m_power[0] = 1.5;
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaQ_GQ_Lambda::Gamma(double q, double Q) 
{
  double val(Gamma_Lambda_Base::Gamma(q,Q));
  if (m_qmass==0.0 || !(m_mode&bpm::massive)) return val;
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

GammaG_GG_Lambda::GammaG_GG_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double asfac):   
  Gamma_Lambda_Base(bpt::gamma_g2gg,mode,lambda,runas,0.0,asfac) 
{
  m_colfac = CA;
  m_dlog   = 1.;
  if (m_mode & (bpm::linear_term | bpm::power_corrs)) {
    m_slog = -11./12.;
    if (m_mode & bpm::power_corrs) {
      m_power[0] = 2.0;
      m_power[1] = -0.5;
      m_power[2] = 1.0/3.0;
    }
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

GammaG_QQ_Lambda::GammaG_QQ_Lambda(bpm::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,
				   double qmass,double asfac) : 
  Gamma_Lambda_Base(bpt::gamma_g2qq,mode,lambda,runas,qmass,asfac) 
{
  m_colfac = TR;
  m_dlog   = 0.;
  if (m_mode & (bpm::linear_term | bpm::power_corrs)) {
    m_slog = 1./3.;
    if (m_mode & bpm::power_corrs) {
      m_power[0] = -1.0;
      m_power[1] = 1.0;
      m_power[2] = -2.0/3.0;
    }
  }
  if (m_mode & bpm::soft_kfac) m_kfac=KAPPA;
}

double GammaG_QQ_Lambda::Gamma(double q, double Q) 
{
  double val(Gamma_Lambda_Base::Gamma(q,Q));
  if (m_qmass==0.0 || !(m_mode&bpm::massive) ||
      m_mode&bpm::dead_cone) return val;
  val+=m_colfac*m_lastas/M_PI/q/(1.0+sqr(m_qmass/q))* 
    (1.0-1.0/3.0/(1.0+sqr(m_qmass/q)));
  return Max(val,0.0);
}

