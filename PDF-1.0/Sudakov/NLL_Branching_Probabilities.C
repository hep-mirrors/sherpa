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
}

using namespace SHERPA;
using namespace ATOOLS;
using namespace MODEL;

Gamma_Lambda_Base::Gamma_Lambda_Base(BP::code mode,double lambda, 
				     MODEL::Running_AlphaS * runas, int nf): 
  m_mode(mode), m_lambda(lambda), p_runas(runas) 
{
  m_powercorr=0;
  m_f2=0.;
  m_f3=0.;
  m_f4=0.;
  m_qlimit=sqrt(1./3.) * rpa.gen.Ecms();
  
  if (nf==-1) nf=NF;

  m_nlo    = 0;
  m_kfac = 0;

  if (m_mode & BP::gamma_kfac) {
    m_kfac=  CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*NF;
    m_nlo =1;
  }
//   std::cout<<" Kfactor = "<<m_kfac<<std::endl;
//   std::cout<<" m_mode="<<m_mode<<std::endl;

  switch (m_mode & 11) {
  case BP::gammaq_powercorr : 
    m_powercorr=1;
    m_f2=1.;
  case BP::gammaq : 
    m_cc=CF;
    m_f1=-3./4.;
    break;
  case BP::gammag_powercorr : 
    m_powercorr=1;
    m_f2=1.;
  case BP::gammag : 
    m_cc=CA;
    m_f1=-11./12.;
    break;
  case BP::gammaf_powercorr : 
    m_powercorr=1;
  case BP::gammaf : 
    m_cc=nf/6.;
    m_f1=0.;
    break;
  default:
    m_cc=CF;
    m_f1=-3./4.;
    break;
  }
  if ((m_mode & BP::gamma_corr) && ((m_mode & 3)==BP::gammaq)) {
//     std::cout<<" ME corrections "<<std::endl;
//     std::cout<<" m_mode="<<m_mode<<std::endl;
    m_f2=0.973527; // c1
    m_f3=0.306434; // c2
    m_f4=-2.;      // c4/2
  }
}

double Gamma_Lambda_Base::AlphaS(double t)                
{ 
  if (t<0.) t = -t;

  if (p_runas) {
    return (*p_runas)(t);
  }

  return 4.*M_PI/(BETA0*log(t/sqr(m_lambda)));

  // 2 M_PI /(BETA0* log(q/m_lambda));
}

double Gamma_Lambda_Base::GammaF(double Q0, double Q) {
  if ((m_mode & BP::gamma_kinlim) && (Q0>m_qlimit)) return 0.;
  //  double val =(4.*m_cc)/(BETA0*Q0*log(Q0/m_lambda));
  double a    = AlphaS(sqr(Q0));
  double val = 2.*m_cc/Q0 * a/M_PI;

  if ((m_mode & BP::gamma_cut) && (val<0.)) return 0.;
  return val;
}

double Gamma_Lambda_Base::IntGammaF(double Q0, double Q) {
  double xi0 = log(Q0/m_lambda);
  double xi1 = log(Q/m_lambda);
  return 4.*m_cc/BETA0*log(xi1/xi0);
}

double Gamma_Lambda_Base::Gamma(double q, double Q) {
//   if (q==20. && Q==2000.)  std::cout<<"Gamma("<<q<<","<<Q<<",  m="<<m_mode<<")"<<std::endl;
  if ((m_mode & BP::gamma_kinlim) && (q>m_qlimit)) return 0.;
//   double val=(4.*m_cc*(m_f1 + m_f2*q/Q  + log(Q/q)))/
//     (BETA0*q*log(q/m_lambda));

  double a   = AlphaS(sqr(q));
  //  double val = 2.*m_cc/q * a/M_PI * (m_f1 + m_f2*q/Q  + (1.+a/(2.*M_PI)*m_kfac)*log(Q/q));
  double val = 2.*m_cc* a/M_PI *(1./q * (m_f1 + m_f2*q/Q  + m_f3*sqr(q/Q)  + (1.+a/(2.*M_PI)*m_kfac)*log(Q/q))
				 + m_f4 * q/sqr(Q) * log(Q/q));
  if ((m_mode & BP::gamma_cut) && (val<0.)) return 0.;
  return val;
}

double Gamma_Lambda_Base::IntGamma(double Q0, double Q) {
//   std::cout<<" IntGamma "<<m_mode<<" called "<<std::endl;

  double xi0 = log(Q0/m_lambda);
  double xi1 = log(Q/m_lambda);
  double fac = 4.*m_cc/BETA0;
  
  double part1= fac* (log(Q0/Q) + (m_f1 + xi1)*log(dabs(xi1/xi0)));
  double part2= 0;
  if (m_powercorr) 
    part2=fac* m_lambda/Q *(ReIncompleteGamma0(-xi0)-ReIncompleteGamma0(-xi1));
  
  if (xi0<=0.) {
    part1=1.e3;
    std::cout<<" WARNING: Q0 smaller than Lambda! "<<std::endl;
  }  
  return part1+part2;
}


GammaQ_Lambda::GammaQ_Lambda(BP::code mode, double lambda, MODEL::Running_AlphaS * runas)
  : Gamma_Lambda_Base(BP::gammaq|mode,lambda,runas) {}

GammaG_Lambda::GammaG_Lambda(BP::code mode, double lambda, MODEL::Running_AlphaS * runas)
  : Gamma_Lambda_Base(BP::gammag|mode,lambda,runas) {}

GammaF_Lambda::GammaF_Lambda(BP::code mode, double lambda, MODEL::Running_AlphaS * runas, int nf)
  : Gamma_Lambda_Base(BP::gammaf|mode,lambda,runas,nf) {}


double GammaF_Lambda::Gamma(double q, double Q) {
  return GammaF(q,Q);
}
double GammaF_Lambda::IntGamma(double q, double Q) {
  return IntGammaF(q,Q);
}


// ==================================================

Gamma_AlphaS_Base::Gamma_AlphaS_Base(BP::code mode, double asmu, double mu2, int nf) :
  m_mode(mode), m_asmu(asmu), m_mu2(mu2)
{
  if (nf==-1) nf=NF;

  m_nlo    = 0;

  m_kfac = 0;
  if (m_nlo) m_kfac=  CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*NF;

  m_f2=0.;

  switch (m_mode & 11) {
  case BP::gammaq_powercorr : 
  case BP::gammaq : 
    m_cc=CF;
    m_f1=-3./4.;
    break;
  case BP::gammag_powercorr : 
  case BP::gammag : 
    m_cc=CA;
    m_f1=-11./12.;
    break;
  case BP::gammaf_powercorr : 
  case BP::gammaf : 
    m_cc=nf/6.;
    m_f1=0.;
    break;
  default:
    m_cc=CF;
    m_f1=-3./4.;
    break;
  }
}


double Gamma_AlphaS_Base::AlphaS(double t)                
{ 
  if (t<0.) t = -t;
  // - lambda - parametrisation:
  //  return 4.*M_PI/(beta0*log(t/lambda2));
  // - first order:
  double   w = 1.-BETA0*m_asmu/(4.*M_PI)*log(m_mu2/t);
  double   a = m_asmu/w;
  if (m_nlo) a *= (1. - (BETA1*m_asmu)/(BETA0*2.*M_PI*w)*log(w));
  return a;
}


double Gamma_AlphaS_Base::Gamma(double q, double Q)
{
  double a = AlphaS(q*q);
  return 2.*m_cc/M_PI * a/q * ((1.+a/(2.*M_PI)*m_kfac)*log(Q/q)+m_f1); 
}

double Gamma_AlphaS_Base::IntGamma(double q, double Q)
{
  if ((m_nlo==0) && (m_kfac==0.)) {
    double balpi = BETA0*m_asmu/(4.*M_PI) ;
    double eta0  = balpi * log(sqr(q)/m_mu2);
    double eta1  = balpi * log(sqr(Q)/m_mu2);
    return 8.*M_PI * m_cc/(sqr(BETA0)* m_asmu) *
      ( balpi * log(sqr(q/Q)) + (1.+ balpi*(log(sqr(Q)/m_mu2) + 2.*m_f1))*log((1+eta1)/(1+eta0)));
  } 
  else {
    msg.Error()<<"Error in Gamma_AlphaS_Base::IntGamma"<<std::endl
	       <<"    Analytic version for higher orders not implemented yet."<<std::endl;
    return -1.;
  }
}

double Gamma_AlphaS_Base::GammaF(double q, double Q)
{
  return 2.*m_cc/M_PI * AlphaS(q*q)/q;
}

double Gamma_AlphaS_Base::IntGammaF(double q, double Q)
{
  if ((m_nlo==0) && (m_kfac==0.)) {
    double eta0 = BETA0*m_asmu/(4.*M_PI) * log(sqr(q)/m_mu2);
    double eta1 = BETA0*m_asmu/(4.*M_PI) * log(sqr(Q)/m_mu2);
    return 2.*NF/(3.*BETA0) * log((1+eta1)/(1+eta0));
  }
  else {
    msg.Error()<<"ERROR in NLL_Sudakov::IntGammaQ."<<std::endl
	       <<"    Analytic version for higher orders not implemented yet."<<std::endl;
    return -1.;
  }
}


GammaQ_AlphaS::GammaQ_AlphaS(BP::code mode, double asmu,double mu2):
  Gamma_AlphaS_Base(BP::gammaq|mode,asmu,mu2) {}

GammaG_AlphaS::GammaG_AlphaS(BP::code mode, double asmu,double mu2):
  Gamma_AlphaS_Base(BP::gammag|mode,asmu,mu2) {}

GammaF_AlphaS::GammaF_AlphaS(BP::code mode, double asmu,double mu2, int nf):
  Gamma_AlphaS_Base(BP::gammaf|mode,asmu,mu2,nf) {}


double GammaF_AlphaS::Gamma(double q, double Q) {
  return GammaF(q,Q);
}
double GammaF_AlphaS::IntGamma(double q, double Q) {
  return IntGammaF(q,Q);
}


// ==================================================

Gamma_Lambda_Massive::Gamma_Lambda_Massive(BP::code mode, double lambda, 
					   MODEL::Running_AlphaS * runas, ATOOLS::Flavour fl) :
  Gamma_Lambda_Base(mode,lambda,runas), m_mass(fl.PSMass())
{

}

double Gamma_Lambda_Massive::GammaQ(double q, double Q) 
{
  double massless=Gamma_Lambda_Base::Gamma(q,Q);
  if (m_mass==0.) return massless;

  double massive =0;

  //  double pref = 2./(BETA0*q*log(q/m_lambda)) * m_cc;

  double a    = AlphaS(sqr(q));
  double pref = m_cc/q * a/M_PI;
  double rest=0;

  if (m_mass/q<5.e-2) {
    double x = m_mass/q;
    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    double x8 = x4*x4;
    rest = (-11.*x2)/12. + (7.*x4)/15. - (53.*x6)/168. + (43.*x8)/180.;
    massive = pref * rest;
  }
  else {
    double m= m_mass;
    double m2=m*m;
    double q2=q*q;
    rest = (0.5 - q/m * atan(m/q) - (2. *m2   - q2)/(2.*m2) *log((m2+q2)/q2));
    massive = pref * rest;
  }

  return massless + massive;
}

double Gamma_Lambda_Massive::IntGammaQ(double q, double Q) 
{
  if (m_mass==0.) return Gamma_Lambda_Base::IntGamma(q,Q);
  return -1.;
}


double Gamma_Lambda_Massive::GammaF(double q, double Q) 
{
  if (m_mass==0.) return Gamma_Lambda_Base::GammaF(q,Q);
  double mfac = 1./(1.+sqr(m_mass/q));
  double rest = mfac * (1. - 1./3.*mfac);
  //  double pref = 2./(BETA0*q*log(q/m_lambda));
  double a    = AlphaS(sqr(q));
  double pref = m_cc/q * a/M_PI;
  return pref*rest;
}


double Gamma_Lambda_Massive::IntGammaF(double q, double Q) 
{
  if (m_mass==0.) return Gamma_Lambda_Base::IntGammaF(q,Q) ;
  return -1.;
}


GammaG_Lambda_Massive::GammaG_Lambda_Massive(BP::code mode, double lambda, 
					     MODEL::Running_AlphaS * runas):
  Gamma_Lambda_Massive(BP::gammag|mode,lambda,runas,Flavour(kf::gluon)) {}

double GammaG_Lambda_Massive::Gamma(double q, double Q) 
{
  return Gamma_Lambda_Base::Gamma(q,Q); 
}

double GammaG_Lambda_Massive::IntGamma(double q, double Q) 
{
  return Gamma_Lambda_Base::IntGamma(q,Q);
}


GammaQ_Lambda_Massive::GammaQ_Lambda_Massive(BP::code mode, double lambda, MODEL::Running_AlphaS * runas, Flavour fl):
  Gamma_Lambda_Massive(BP::gammaq|mode,lambda,runas,fl) {}

double GammaQ_Lambda_Massive::Gamma(double q, double Q) 
{
  return Gamma_Lambda_Massive::GammaQ(q,Q);
}

double GammaQ_Lambda_Massive::IntGamma(double q, double Q) 
{
  return Gamma_Lambda_Massive::IntGammaQ(q,Q);
}


GammaF_Lambda_Massive::GammaF_Lambda_Massive(BP::code mode, double lambda, 
					     MODEL::Running_AlphaS * runas, Flavour fl):
  Gamma_Lambda_Massive(BP::gammaf|mode,lambda,runas,fl) 
{
  m_cc=TR;
}

double GammaF_Lambda_Massive::Gamma(double q, double Q) 
{
  return Gamma_Lambda_Massive::GammaF(q,Q);
}

double GammaF_Lambda_Massive::IntGamma(double q, double Q) 
{
  return Gamma_Lambda_Massive::IntGammaF(q,Q);
}

