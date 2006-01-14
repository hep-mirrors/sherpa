#include "NLL_Branching_Probabilities.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"

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
Gamma_Lambda_Base(BPType::code type,BPMode::code mode,double lambda, 
		  MODEL::Running_AlphaS * runas, int nf,double asfac): 
  m_type(type), m_mode(mode), 
  m_colfac(0.), m_dlog(0.), m_slog(0.), m_power(0.),
  m_lambda(lambda), m_as_factor(asfac), p_runas(runas), 
  m_kfac(0.)
{ }

double Gamma_Lambda_Base::AlphaS(double t)                
{ 
  if (p_runas) return (*p_runas)(dabs(t)/m_as_factor);
  return 4.*M_PI/(BETA0*log(dabs(t)/sqr(m_lambda)));
}

double Gamma_Lambda_Base::Gamma(double q, double Q) 
{
  double as_q(AlphaS(sqr(q)));
  double val = 2.*m_colfac* as_q/M_PI/q * 
    (m_dlog * (1.+m_kfac*as_q/(2.*M_PI)) * log(Q/q) + 
     m_slog + m_power*q/Q);
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
  if (m_power>0.) result += m_power*m_lambda/Q * 
    (ReIncompleteGamma0(-xi0)-ReIncompleteGamma0(-xi1));  
  return 4.*m_colfac/BETA0*result;
}


GammaQ_QG_Lambda::GammaQ_QG_Lambda(BPMode::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,double asfac) : 
  Gamma_Lambda_Base(BPType::gamma_q2qg,mode,lambda,runas,-1,asfac) 
{
  m_colfac = CF;
  m_dlog   = 1.;
  if (m_mode & (BPMode::linear_term | BPMode::power_corrs)) {
    m_slog   = -3./4.;
    if (m_mode & BPMode::power_corrs) m_power  = 1.;
  }
  if (m_mode & BPMode::soft_kfac) m_kfac=KAPPA;
}

GammaQ_GQ_Lambda::GammaQ_GQ_Lambda(BPMode::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,double asfac) : 
  Gamma_Lambda_Base(BPType::gamma_q2gq,mode,lambda,runas,-1,asfac) 
{
  m_colfac = CF;
  m_dlog   = 1.;
  if (m_mode & (BPMode::linear_term | BPMode::power_corrs)) {
    m_slog = -3./4.;
    if (m_mode & BPMode::power_corrs) m_power  = 1.;
  }
  if (m_mode & BPMode::soft_kfac) m_kfac=KAPPA;
}

GammaG_GG_Lambda::GammaG_GG_Lambda(BPMode::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,double asfac) : 
  Gamma_Lambda_Base(BPType::gamma_g2gg,mode,lambda,runas,-1,asfac) 
{
  m_colfac = CA;
  m_dlog   = 1.;
  if (m_mode & (BPMode::linear_term | BPMode::power_corrs)) {
    m_slog = -11./12.;
    if (m_mode & BPMode::power_corrs) m_power  = 1.;
  }
  if (m_mode & BPMode::soft_kfac) m_kfac=KAPPA;
}

GammaG_QQ_Lambda::GammaG_QQ_Lambda(BPMode::code mode, double lambda, 
				   MODEL::Running_AlphaS * runas,double asfac) : 
  Gamma_Lambda_Base(BPType::gamma_g2qq,mode,lambda,runas,-1,asfac) 
{
  m_colfac = TR;
  m_dlog   = 0.;
  if (m_mode & (BPMode::linear_term | BPMode::power_corrs)) {
    m_slog = 1./3.;
    if (m_mode & BPMode::power_corrs) m_power  = 0.;
  }
  if (m_mode & BPMode::soft_kfac) m_kfac=KAPPA;
}

